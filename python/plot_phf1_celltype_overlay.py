#!/usr/bin/env python
# ---------------------------------------------------------------------------
# plot_phf1_celltype_overlay.py
#
# Figure panels: 1E
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
"""
plot_phf1_celltype_overlay.py
=============================

Paper figures: a two-channel composite of the post-stain images (PHF1 ptau + DAPI
nuclei) with CosMx cells overlaid as their actual segmentation polygons coloured by
celltype. Per sample and per view (whole section + N zoomed single FOVs) it writes:
  <sample>_stain[_zoom_fovK].png      composite only (PHF1 + DAPI) + scale bar
  <sample>_grey[_zoom_fovK].png       composite + cells in uniform grey
  <sample>_celltype[_zoom_fovK].png   composite + cells coloured by celltype (+ legend)

Channels are additively blended: DAPI -> --dapi-color (default blue), PHF1 -> --phf1-color
(default red); overlap reads magenta. Cells/polygons are placed with the validated
transform the extractor uses (FOV offsets + mm_per_px from stitched_slides/*/.zattrs;
M from PHF1/Transformed; DAPI shares PHF1 pixel coords). A scale bar (image px size from
M) is on every panel. Zoom FOVs default to the N FOVs with the most manual PHF1+ cells
(the `PHF1` column of the cells CSV); falls back to ptau-densest by stain, or use --zoom-fov.

Inputs:
  --cells-csv   cell_id, sample_id, celltype, x_FOV_px, y_FOV_px [, PHF1]  (default PHF1/seu_coords.csv)
  --palette-csv celltype,hex
Run under spatial_env (numpy, tifffile, matplotlib).

Usage:
  python plot_phf1_celltype_overlay.py --polygons --zoom --n-zoom 3
"""

import argparse
import csv
import gzip
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import extract_phf1_intensity as E

try:
    import tifffile
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.lines import Line2D
    from matplotlib.collections import PolyCollection
except ImportError as e:
    sys.exit(f"Missing dependency ({e}); run under spatial_env.")

POLY_CSV = {
    2: "flatFiles/flatFiles/IGFQ002102_matthews_1312026_CosMx_RNA_1b/"
       "IGFQ002102_matthews_1312026_CosMx_RNA_1b-polygons.csv.gz",
    1: "flatFiles/flatFiles/IGFQ002102_matthews_1312026_CosMx_RNA_2a/"
       "IGFQ002102_matthews_1312026_CosMx_RNA_2a-polygons.csv.gz",
}


def load_palette(path):
    return {r["celltype"]: r["hex"] for r in csv.DictReader(open(path, newline=""))}


def load_cells(path):
    by_sample = {}
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh)
        for col in ("cell_id", "sample_id", "celltype", "x_FOV_px", "y_FOV_px"):
            if col not in rdr.fieldnames:
                sys.exit(f"{path} lacks required column '{col}'.")
        has_phf1 = "PHF1" in rdr.fieldnames
        for r in rdr:
            d = by_sample.setdefault(r["sample_id"],
                                     {"cell_id": [], "celltype": [], "xl": [], "yl": [], "phf1": []})
            d["cell_id"].append(r["cell_id"]); d["celltype"].append(r["celltype"])
            d["xl"].append(float(r["x_FOV_px"])); d["yl"].append(float(r["y_FOV_px"]))
            d["phf1"].append(str(r["PHF1"]).strip().upper() in ("TRUE", "1") if has_phf1 else False)
    for d in by_sample.values():
        d["cell_id"] = np.array(d["cell_id"], dtype=object)
        d["celltype"] = np.array(d["celltype"], dtype=object)
        d["xl"] = np.array(d["xl"]); d["yl"] = np.array(d["yl"])
        d["phf1"] = np.array(d["phf1"], dtype=bool)
    return by_sample, has_phf1


# --------------------------------------------------------------------------- #
# transforms + geometry
# --------------------------------------------------------------------------- #
def _local_to_phf1px(xl, yl, mm_per_px, X_mm, Y_mm, Minv):
    w0 = Y_mm + mm_per_px * yl
    w1 = -X_mm + mm_per_px * xl
    dst = Minv @ np.vstack([w0, w1, np.ones(len(xl))])
    return dst[0], dst[1]


def cells_to_phf1_px(base, transform_dir, sample, cell_ids, xl, yl):
    geom = {}
    Minv = np.linalg.inv(E.load_napari_affine(os.path.join(transform_dir, f"{sample}.txt")))
    rows = np.empty(len(xl)); cols = np.empty(len(xl)); by_fov = {}
    for i, cid in enumerate(cell_ids):
        pfx, fov, _ = E.parse_cell_id(cid)
        by_fov.setdefault((pfx, fov), []).append(i)
    for (pfx, fov), idx in by_fov.items():
        slide = E.PREFIX_TO_SLIDE[pfx]
        if slide not in geom:
            geom[slide] = E.load_slide_geometry(base, slide)
        mm_per_px, offsets = geom[slide]
        X_mm, Y_mm = offsets[fov]
        idx = np.array(idx)
        r, c = _local_to_phf1px(xl[idx], yl[idx], mm_per_px, X_mm, Y_mm, Minv)
        rows[idx] = r; cols[idx] = c
    return rows, cols


def load_polygons(base, prefix, cell_ids_set):
    path = os.path.join(base, POLY_CSV[prefix])
    if not os.path.exists(path):
        return {}
    acc = {}
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh):
            cid = r["cell"]
            if cid not in cell_ids_set:
                continue
            d = acc.setdefault(cid, {"x": [], "y": [], "fov": int(r["fov"])})
            d["x"].append(float(r["x_local_px"])); d["y"].append(float(r["y_local_px"]))
    return {cid: (np.column_stack([d["x"], d["y"]]), d["fov"]) for cid, d in acc.items()}


def polygons_to_phf1px(base, transform_dir, sample, polys):
    geom = {}
    Minv = np.linalg.inv(E.load_napari_affine(os.path.join(transform_dir, f"{sample}.txt")))
    out = {}
    for cid, (verts, fov) in polys.items():
        pfx, _, _ = E.parse_cell_id(cid)
        slide = E.PREFIX_TO_SLIDE[pfx]
        if slide not in geom:
            geom[slide] = E.load_slide_geometry(base, slide)
        mm_per_px, offsets = geom[slide]
        X_mm, Y_mm = offsets[fov]
        r, c = _local_to_phf1px(verts[:, 0], verts[:, 1], mm_per_px, X_mm, Y_mm, Minv)
        out[cid] = np.column_stack([r, c])
    return out


def px_size_um_from_M(M):
    return float(np.sqrt(abs(np.linalg.det(M[:2, :2])))) * 1000.0


def fov_size(base, slide):
    cm = json.load(open(os.path.join(base, E.ZATTRS.format(slide=slide))))["CosMx"]
    return int(cm.get("fov_width", 4256)), int(cm.get("fov_height", 4256))


def fov_bbox_phf1(base, transform_dir, sample, prefix, fov):
    slide = E.PREFIX_TO_SLIDE[prefix]
    mm_per_px, offsets = E.load_slide_geometry(base, slide)
    X_mm, Y_mm = offsets[fov]
    Minv = np.linalg.inv(E.load_napari_affine(os.path.join(transform_dir, f"{sample}.txt")))
    fw, fh = fov_size(base, slide)
    cx = np.array([0, fw - 1, 0, fw - 1], float); cy = np.array([0, 0, fh - 1, fh - 1], float)
    r, c = _local_to_phf1px(cx, cy, mm_per_px, X_mm, Y_mm, Minv)
    return int(r.min()), int(r.max()), int(c.min()), int(c.max())


def rank_fovs_by_phf1(fovs, phf1_bool, n):
    cnt = {}
    for f, p in zip(fovs, phf1_bool):
        if p:
            cnt[int(f)] = cnt.get(int(f), 0) + 1
    ranked = sorted(cnt.items(), key=lambda x: -x[1])[:n]
    return [f for f, _ in ranked], cnt


def rank_fovs_by_stain(img, base, transform_dir, sample, prefix, fovs_present, n):
    H, W = img.shape; scored = []
    for f in fovs_present:
        try:
            r0, r1, c0, c1 = fov_bbox_phf1(base, transform_dir, sample, prefix, int(f))
        except Exception:
            continue
        r0, r1 = max(0, r0), min(H, r1); c0, c1 = max(0, c0), min(W, c1)
        if r1 <= r0 or c1 <= c0:
            continue
        scored.append((float(np.asarray(img[r0:r1, c0:c1]).mean()), int(f)))
    scored.sort(reverse=True)
    return [f for _, f in scored[:n]]


# --------------------------------------------------------------------------- #
# image processing + drawing
# --------------------------------------------------------------------------- #
def subtract_background(img, p):
    if not p or p <= 0:
        return img, 0.0
    nz = img[img > 0]
    bg = float(np.percentile(nz, p)) if nz.size else 0.0
    out = img.astype(np.float32) - bg
    out[out < 0] = 0
    return out, bg


def composite(phf1_sub, dapi_sub, pcol, dcol):
    """Additive RGB composite of two single-channel images, each contrast-stretched."""
    pv = max(1.0, np.percentile(phf1_sub, 99.5))
    dv = max(1.0, np.percentile(dapi_sub, 99.5))
    p = np.clip(phf1_sub / pv, 0, 1); dd = np.clip(dapi_sub / dv, 0, 1)
    pr, pg, pb = mcolors.to_rgb(pcol); dr, dg, db = mcolors.to_rgb(dcol)
    rgb = np.stack([p * pr + dd * dr, p * pg + dd * dg, p * pb + dd * db], axis=-1)
    return np.clip(rgb, 0, 1)


def add_scalebar(ax, Ws, Hs, px_size_um, down, fixed_um=None, lw=2.5, show_label=True):
    field_um = Ws * down * px_size_um
    if field_um <= 0:
        return
    if fixed_um:
        bar_um = fixed_um
    else:
        target = field_um * 0.2
        k = 10 ** int(np.floor(np.log10(target)))
        bar_um = next((m * k for m in (1, 2, 5, 10) if m * k >= target), 10 * k)
    bar_disp = bar_um / px_size_um / down
    x0, y = 0.05 * Ws, 0.95 * Hs
    ax.plot([x0, x0 + bar_disp], [y, y], color="white", lw=lw, solid_capstyle="butt")
    if show_label:
        lbl = f"{bar_um:g} µm" if bar_um < 1000 else f"{bar_um/1000:g} mm"
        ax.text(x0 + bar_disp / 2, y - 0.015 * Hs, lbl, color="white", ha="center", va="bottom", fontsize=7)


def render(sample, phf1, dapi, outdir, down, pt, suffix, palette,
           rows, cols, celltypes, cell_polys, poly_lw, poly_alpha, px_size_um, pcol, dcol,
           fig_width_cm=None, scalebar_um=None, legend_cts=None,
           scalebar_lw=2.5, scalebar_label=True, dpi=300):
    rgb = composite(np.asarray(phf1[::down, ::down]), np.asarray(dapi[::down, ::down]), pcol, dcol)
    Hs, Ws = rgb.shape[:2]
    use_poly = cell_polys is not None and any(p is not None for p in cell_polys)
    LEG_FRAC = 0.66   # image share of the figure width when a legend is present (fixed-width mode)

    def base_ax(with_legend=False):
        if fig_width_cm:                       # exact figure width (image + legend) in cm
            w_in = fig_width_cm / 2.54
            f = LEG_FRAC if with_legend else 1.0
            fig = plt.figure(figsize=(w_in, w_in * f * Hs / Ws))
            ax = fig.add_axes([0, 0, f, 1])
        else:
            fig, ax = plt.subplots(figsize=(Ws / 300.0, Hs / 300.0))
        ax.imshow(rgb, interpolation="nearest")
        ax.set_xlim(0, Ws); ax.set_ylim(Hs, 0); ax.axis("off")
        add_scalebar(ax, Ws, Hs, px_size_um, down, scalebar_um, lw=scalebar_lw, show_label=scalebar_label)
        return fig, ax

    def save(fig, path):
        if fig_width_cm:
            fig.savefig(path, dpi=dpi)                                   # keep exact figsize
        else:
            fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.02)
        plt.close(fig)

    def cell_layer(ax, mode):
        if use_poly:
            verts, faces, present = [], [], []
            for i, p in enumerate(cell_polys):
                if p is None:
                    continue
                verts.append(np.column_stack([p[:, 1] / down, p[:, 0] / down]))
                ct = celltypes[i]
                faces.append("#BBBBBB" if mode == "grey" else palette.get(ct, "#777777"))
                present.append(ct)
            ax.add_collection(PolyCollection(verts, facecolors=faces, edgecolors=faces,
                              linewidths=poly_lw, alpha=(min(0.4, poly_alpha) if mode == "grey" else poly_alpha)))
            return present
        rs, cs = rows / down, cols / down
        inb = (rs >= 0) & (rs < Hs) & (cs >= 0) & (cs < Ws)
        if mode == "grey":
            ax.scatter(cs[inb], rs[inb], s=pt, c="#BBBBBB", linewidths=0, alpha=0.6)
            return list(celltypes[inb])
        for c in [x for x in palette if x in set(celltypes[inb].tolist())]:
            m = celltypes[inb] == c
            ax.scatter(cs[inb][m], rs[inb][m], s=pt, c=palette[c], linewidths=0, alpha=0.9)
        return list(celltypes[inb])

    fig, ax = base_ax(False); save(fig, os.path.join(outdir, f"{sample}_stain{suffix}.png"))
    fig, ax = base_ax(False); cell_layer(ax, "grey"); save(fig, os.path.join(outdir, f"{sample}_grey{suffix}.png"))

    fig, ax = base_ax(True); present = cell_layer(ax, "celltype")
    pres = legend_cts if legend_cts is not None else [c for c in palette if c in set(present)]
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=3,
                      markerfacecolor=palette[c], markeredgewidth=0, label=c) for c in pres]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=4, frameon=False, labelspacing=0.25, handletextpad=0.2)
    save(fig, os.path.join(outdir, f"{sample}_celltype{suffix}.png"))


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-dir", default="<PROJECT_ROOT>/phf1_v2")
    ap.add_argument("--cells-csv", default=None)
    ap.add_argument("--palette-csv", default=None)
    ap.add_argument("--phf1-dir", default=None)
    ap.add_argument("--dapi-dir", default=None)
    ap.add_argument("--transform-dir", default=None)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--samples", nargs="*", default=None)
    ap.add_argument("--phf1-color", default="red")
    ap.add_argument("--dapi-color", default="blue")
    ap.add_argument("--downsample", type=int, default=6)
    ap.add_argument("--point-size", type=float, default=1.0)
    ap.add_argument("--polygons", action="store_true")
    ap.add_argument("--poly-lw", type=float, default=0.2)
    ap.add_argument("--poly-alpha", type=float, default=0.5)
    ap.add_argument("--bg-percentile", type=float, default=50.0, help="PHF1 background subtraction percentile")
    ap.add_argument("--dapi-bg-percentile", type=float, default=10.0, help="DAPI background subtraction percentile")
    ap.add_argument("--zoom", action="store_true")
    ap.add_argument("--n-zoom", type=int, default=3, help="number of zoom FOVs (most manual PHF1+ cells)")
    ap.add_argument("--zoom-fov", default=None, help="comma list of explicit FOV numbers to zoom")
    ap.add_argument("--zoom-downsample", type=int, default=1)
    ap.add_argument("--zoom-point-size", type=float, default=8.0)
    ap.add_argument("--zoom-width-cm", type=float, default=9.38, help="zoom figure width (image+legend) in cm")
    ap.add_argument("--zoom-scalebar-um", type=float, default=100.0, help="zoom scale-bar length in microns")
    ap.add_argument("--skip-whole", action="store_true", help="only (re)render the zoom FOVs")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--scalebar-lw", type=float, default=2.5)
    ap.add_argument("--hide-scalebar-text", action="store_true")
    args = ap.parse_args()
    base = args.base_dir
    args.cells_csv     = args.cells_csv     or os.path.join(base, "PHF1/seu_coords.csv")
    args.palette_csv   = args.palette_csv   or os.path.join(base, "PHF1/celltype_palette.csv")
    args.phf1_dir      = args.phf1_dir      or os.path.join(base, "PHF1/PHF1")
    args.dapi_dir      = args.dapi_dir      or os.path.join(base, "PHF1/DAPI")
    args.transform_dir = args.transform_dir or os.path.join(base, "PHF1/Transformed")
    args.out_dir       = args.out_dir       or os.path.join(base, "plots/phf1_celltype_overlay")
    os.makedirs(args.out_dir, exist_ok=True)

    palette = load_palette(args.palette_csv)
    cells, has_phf1 = load_cells(args.cells_csv)
    samples = args.samples or sorted(cells.keys())

    # fixed, complete legend = every celltype present anywhere in the dataset (palette order),
    # so every panel shows the same full legend regardless of which cells are in view.
    all_cts = set()
    for s in samples:
        if s in cells:
            all_cts.update(cells[s]["celltype"].tolist())
    legend_cts = [c for c in palette if c in all_cts]

    poly_cache = {}
    if args.polygons:
        for pfx in {E.parse_cell_id(cells[s]["cell_id"][0])[0] for s in samples if s in cells}:
            ids = set()
            for s in samples:
                if s in cells and E.parse_cell_id(cells[s]["cell_id"][0])[0] == pfx:
                    ids.update(cells[s]["cell_id"].tolist())
            poly_cache[pfx] = load_polygons(base, pfx, ids)
            print(f"[overlay] polygons slide prefix {pfx}: {len(poly_cache[pfx])} cells")

    for s in samples:
        if s not in cells:
            print(f"[overlay] {s}: not in cells CSV, skipping"); continue
        d = cells[s]
        pfx = E.parse_cell_id(d["cell_id"][0])[0]
        Mp = os.path.join(args.transform_dir, f"{s}.txt")
        pp = os.path.join(args.phf1_dir, f"{s}.tif"); dp = os.path.join(args.dapi_dir, f"{s}.tif")
        if not (os.path.exists(Mp) and os.path.exists(pp) and os.path.exists(dp)):
            print(f"[overlay] {s}: missing transform/PHF1/DAPI, skipping"); continue
        px_um = px_size_um_from_M(E.load_napari_affine(Mp))
        rows, cols = cells_to_phf1_px(base, args.transform_dir, s, d["cell_id"], d["xl"], d["yl"])

        cell_polys = None
        if args.polygons:
            raw = {cid: poly_cache.get(pfx, {})[cid] for cid in d["cell_id"] if cid in poly_cache.get(pfx, {})}
            ph = polygons_to_phf1px(base, args.transform_dir, s, raw)
            cell_polys = [ph.get(cid) for cid in d["cell_id"]]

        phf1 = tifffile.imread(pp); phf1 = phf1 if phf1.ndim == 2 else phf1[..., 0]
        dapi = tifffile.imread(dp); dapi = dapi if dapi.ndim == 2 else dapi[..., 0]
        phf1, pbg = subtract_background(phf1, args.bg_percentile)
        dapi, dbg = subtract_background(dapi, args.dapi_bg_percentile)
        H, W = phf1.shape
        fovs = np.array([E.parse_cell_id(cid)[1] for cid in d["cell_id"]])

        zoom_fovs = []
        if args.zoom:
            if args.zoom_fov:
                zoom_fovs = [int(v) for v in args.zoom_fov.split(",")]
            elif has_phf1 and d["phf1"].any():
                zoom_fovs, cnt = rank_fovs_by_phf1(fovs, d["phf1"], args.n_zoom)
                print(f"[overlay] {s}: zoom FOVs (manual PHF1+ counts) = "
                      f"{[(f, cnt[f]) for f in zoom_fovs]}")
            else:
                zoom_fovs = rank_fovs_by_stain(phf1, base, args.transform_dir, s, pfx, np.unique(fovs), args.n_zoom)
                print(f"[overlay] {s}: zoom FOVs (ptau-densest stain) = {zoom_fovs}")

        if not args.skip_whole:
            render(s, phf1, dapi, args.out_dir, args.downsample, args.point_size, "", palette,
                   rows, cols, d["celltype"], cell_polys, args.poly_lw, args.poly_alpha, px_um,
                   args.phf1_color, args.dapi_color, legend_cts=legend_cts,
                   scalebar_lw=args.scalebar_lw, scalebar_label=not args.hide_scalebar_text, dpi=args.dpi)

        for zf in zoom_fovs:
            r0, r1, c0, c1 = fov_bbox_phf1(base, args.transform_dir, s, pfx, zf)
            r0, r1 = max(0, r0), min(H, r1); c0, c1 = max(0, c0), min(W, c1)
            sel = fovs == zf; idx = np.nonzero(sel)[0]
            cpolys = None
            if cell_polys is not None:
                cpolys = [None if cell_polys[i] is None else cell_polys[i] - np.array([r0, c0]) for i in idx]
            render(s, phf1[r0:r1, c0:c1], dapi[r0:r1, c0:c1], args.out_dir, args.zoom_downsample,
                   args.zoom_point_size, f"_zoom_fov{zf}", palette, rows[sel] - r0, cols[sel] - c0,
                   d["celltype"][sel], cpolys, args.poly_lw, args.poly_alpha, px_um,
                   args.phf1_color, args.dapi_color,
                   fig_width_cm=args.zoom_width_cm, scalebar_um=args.zoom_scalebar_um, legend_cts=legend_cts,
                   scalebar_lw=args.scalebar_lw, scalebar_label=not args.hide_scalebar_text, dpi=args.dpi)
        print(f"[overlay] {s}: whole + {len(zoom_fovs)} zoom FOV(s) (PHF1 bg {pbg:.0f}, DAPI bg {dbg:.0f}, {px_um:.3f} um/px)")
        del phf1, dapi


if __name__ == "__main__":
    main()
