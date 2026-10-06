#!/usr/bin/env python
# ---------------------------------------------------------------------------
# plot_tangle_cell_crops.py
#
# Figure panels: 2C, 2D
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
"""
plot_tangle_cell_crops.py
=========================

Single-cell crops centred on manually-annotated tangle-bearing cells (PHF1+), for a
gallery figure. Per sample it renders up to N crops, each: PHF1 (red) + DAPI (blue)
composite, CosMx segmentation OUTLINES (no fill), transcripts as small dots, and a
scale bar. Saved size_cm x size_cm, no legend.

  <sample>_tangle<K>_<cell_id>.png

Options:
  --celltype       restrict tangle cells to a celltype (e.g. Exc-IT-L2-3-CBLN2-HOPX)
  --celltypes      comma/semicolon list of celltypes (e.g. all neuron labels)
  --rank-col       per-cell column to rank on (default phf1_intensity_p95; use
                   phf1_intensity_p95_z for the per-sample robust z). Cells with no
                   value in that column are dropped rather than ranked last.
  --palette-csv    colour the segmentation outlines by celltype (else white)
  --no-transcripts skip the ~1.1 GB/slide tx flat file; needed for large galleries
  --index-tsv      one row per rendered crop (rank, cell, intensity, file)
  --pad-rank       zero-pad the rank in the filename so the listing sorts by brightness
  --highlight-genes  comma/semicolon list; these transcripts are coloured (highlight),
                     all others drawn faint grey (context). Without it, all RNA one colour.

Tangle cells are the qualifying manual PHF1+ cells with the highest phf1_intensity_p95.
Ranking on any other column, or on manual PHF1-negative cells, is what --rank-col and
--phf1-status neg are for; --only-cells renders a given list of cell_ids instead.
Polygons / transcripts come from the flatFiles, transformed to PHF1-image px with the
validated map (reuses plot_phf1_celltype_overlay + extract_phf1_intensity).

Run under spatial_env (numpy, pandas, tifffile, matplotlib).

Usage:
  python plot_tangle_cell_crops.py --celltype Exc-IT-L2-3-CBLN2-HOPX --n-cells 10 \
     --palette-csv PHF1/celltype_palette.csv --highlight-genes "NDUFB8;COX4I1;..." \
     --out-dir plots/phf1_tangle_exc_oxphos
"""

import argparse
import csv
import gzip
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import extract_phf1_intensity as E
import plot_phf1_celltype_overlay as O

try:
    import tifffile
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
except ImportError as e:
    sys.exit(f"Missing dependency ({e}); run under spatial_env.")

TX_CSV = {
    2: "flatFiles/flatFiles/IGFQ002102_matthews_1312026_CosMx_RNA_1b/"
       "IGFQ002102_matthews_1312026_CosMx_RNA_1b_tx_file.csv.gz",
    1: "flatFiles/flatFiles/IGFQ002102_matthews_1312026_CosMx_RNA_2a/"
       "IGFQ002102_matthews_1312026_CosMx_RNA_2a_tx_file.csv.gz",
}
_prefix_cache = {}
_GEOM = {}
_MINV = {}


#: numeric per-cell columns read from the cells CSV when present; any of these
#: can be used as --rank-col. Absent columns become all-NaN so a cells CSV without
#: them (e.g. a slim PHF1/seu_coords.csv) still loads.
NUMERIC_COLS = ("phf1_intensity_p95", "phf1_intensity_p95_z", "phf1_intensity_mean",
                "phf1_intensity_p95_pct", "phf1_frac_in_bounds")


def load_cells(path):
    by_sample = {}
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh)
        num = [c for c in NUMERIC_COLS if c in rdr.fieldnames]
        has_ct = "celltype" in rdr.fieldnames
        for r in rdr:
            d = by_sample.get(r["sample_id"])
            if d is None:
                d = {"cell_id": [], "xl": [], "yl": [], "phf1": [], "celltype": []}
                d.update({c: [] for c in NUMERIC_COLS})
                by_sample[r["sample_id"]] = d
            d["cell_id"].append(r["cell_id"])
            d["xl"].append(float(r["x_FOV_px"])); d["yl"].append(float(r["y_FOV_px"]))
            d["phf1"].append(str(r["PHF1"]).strip().upper() in ("TRUE", "1"))
            d["celltype"].append(r["celltype"] if has_ct else "")
            for c in NUMERIC_COLS:
                v = r.get(c, "") if c in num else ""
                d[c].append(float(v) if v not in ("", "NA", "NaN") else np.nan)
    for d in by_sample.values():
        for k in ("xl", "yl") + NUMERIC_COLS:
            d[k] = np.array(d[k], float)
        d["p95"] = d["phf1_intensity_p95"]              # back-compat alias
        d["cell_id"] = np.array(d["cell_id"], dtype=object)
        d["celltype"] = np.array(d["celltype"], dtype=object)
        d["phf1"] = np.array(d["phf1"], dtype=bool)
    return by_sample


def pick_tangles(d, n, celltype=None, status="pos", select="top", only=None,
                 rank_col="phf1_intensity_p95", celltypes=None):
    if only:
        return [i for i, c in enumerate(d["cell_id"]) if c in only]
    if status == "pos":
        idx = np.where(d["phf1"])[0]
    elif status == "neg":
        idx = np.where(~d["phf1"])[0]
    else:
        idx = np.arange(len(d["phf1"]))
    if celltype:
        idx = idx[d["celltype"][idx] == celltype]
    if celltypes:
        idx = idx[np.isin(d["celltype"][idx], list(celltypes))]
    if idx.size == 0:
        return []
    if rank_col not in d:
        sys.exit(f"--rank-col {rank_col} is not a column of the cells CSV "
                 f"(known numeric columns: {', '.join(NUMERIC_COLS)})")
    key = d[rank_col][idx]
    if not np.any(np.isfinite(key)):
        # Column absent from the CSV entirely (e.g. a slim PHF1/seu_coords.csv):
        # return cells in file order.
        print(f"[tangle] WARNING: no finite {rank_col} in this sample; taking cells in file order")
        return list(idx[:n])
    idx = idx[np.isfinite(key)]                         # unmeasured cells are dropped, not ranked last
    if idx.size == 0:
        return []
    order = idx[np.argsort(-d[rank_col][idx])]          # brightest first
    if select == "spread" and len(order) > n:           # representative sample across the range
        return list(order[np.linspace(0, len(order) - 1, n).round().astype(int)])
    return list(order[:n])


def load_polys_for_fovs(base, prefix, fovs):
    path = os.path.join(base, O.POLY_CSV[prefix])
    acc = {}
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh):
            if int(r["fov"]) not in fovs:
                continue
            cid = r["cell"]
            dd = acc.setdefault(cid, {"x": [], "y": [], "fov": int(r["fov"])})
            dd["x"].append(float(r["x_local_px"])); dd["y"].append(float(r["y_local_px"]))
    return {cid: (np.column_stack([v["x"], v["y"]]), v["fov"]) for cid, v in acc.items()}


def load_tx_for_fovs(base, prefix, fovs):
    path = os.path.join(base, TX_CSV[prefix])
    keep = []
    for chunk in pd.read_csv(path, usecols=["fov", "x_local_px", "y_local_px", "target"],
                             chunksize=3_000_000, compression="gzip"):
        keep.append(chunk[chunk["fov"].isin(fovs)])
    return pd.concat(keep, ignore_index=True) if keep else pd.DataFrame(columns=["fov", "x_local_px", "y_local_px", "target"])


def _geom_for(base, sample, fov):
    slide = E.PREFIX_TO_SLIDE[_prefix_cache[sample]]
    if slide not in _GEOM:
        _GEOM[slide] = E.load_slide_geometry(base, slide)
    mm_per_px, offsets = _GEOM[slide]
    X_mm, Y_mm = offsets[int(fov)]
    return mm_per_px, X_mm, Y_mm


def _minv(base, transform_dir, sample):
    if sample not in _MINV:
        _MINV[sample] = np.linalg.inv(E.load_napari_affine(os.path.join(transform_dir, f"{sample}.txt")))
    return _MINV[sample]


def to_phf1px_byfov(base, transform_dir, sample, fov_arr, xl, yl):
    Minv = _minv(base, transform_dir, sample)
    fov_arr = np.asarray(fov_arr); xl = np.asarray(xl); yl = np.asarray(yl)
    rows = np.empty(len(xl)); cols = np.empty(len(xl))
    for fov in np.unique(fov_arr):
        m = fov_arr == fov
        try:
            mm_per_px, X_mm, Y_mm = _geom_for(base, sample, fov)
        except KeyError:
            rows[m] = np.nan; cols[m] = np.nan; continue
        r, c = O._local_to_phf1px(xl[m], yl[m], mm_per_px, X_mm, Y_mm, Minv)
        rows[m] = r; cols[m] = c
    return rows, cols


def render_crop(path, phf1_crop, dapi_crop, polys_crop, tx_hi, tx_other, px_um, scalebar_um,
                pcol, dcol, outline_lw, hi_color, hi_size, hi_alpha,
                other_color, other_size, other_alpha, show_other, size_cm, dpi=300,
                scalebar_lw=2.5, scalebar_label=True, target_poly=None,
                target_color=None, target_lw=1.2):
    rgb = O.composite(phf1_crop, dapi_crop, pcol, dcol)
    Hs, Ws = rgb.shape[:2]
    fig = plt.figure(figsize=(size_cm / 2.54, size_cm / 2.54))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(rgb, interpolation="bilinear")
    ax.set_xlim(0, Ws); ax.set_ylim(Hs, 0); ax.axis("off")
    if show_other and len(tx_other):
        ax.scatter(tx_other[:, 1], tx_other[:, 0], s=other_size, c=other_color, linewidths=0, alpha=other_alpha)
    if polys_crop:
        verts = [np.column_stack([p[:, 1], p[:, 0]]) for p, _ in polys_crop]
        cols_ = [c for _, c in polys_crop]
        ax.add_collection(PolyCollection(verts, facecolors="none", edgecolors=cols_, linewidths=outline_lw))
    if target_color and target_poly is not None:
        # the ranked cell, so a gallery reader knows which of the outlined cells
        # carries the intensity value it was selected on
        ax.add_collection(PolyCollection([np.column_stack([target_poly[:, 1], target_poly[:, 0]])],
                                         facecolors="none", edgecolors=target_color,
                                         linewidths=target_lw))
    if len(tx_hi):
        ax.scatter(tx_hi[:, 1], tx_hi[:, 0], s=hi_size, c=hi_color, linewidths=0, alpha=hi_alpha)
    O.add_scalebar(ax, Ws, Hs, px_um, 1, scalebar_um, lw=scalebar_lw, show_label=scalebar_label)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-dir", default="<PROJECT_ROOT>/phf1_v2")
    ap.add_argument("--cells-csv", default=None)
    ap.add_argument("--palette-csv", default=None, help="colour outlines by celltype (else white)")
    ap.add_argument("--phf1-dir", default=None)
    ap.add_argument("--dapi-dir", default=None)
    ap.add_argument("--transform-dir", default=None)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--samples", nargs="*", default=None)
    ap.add_argument("--celltype", default=None, help="restrict tangle cells to this celltype")
    ap.add_argument("--celltypes", default=None,
                    help="comma/semicolon list of celltypes to keep (e.g. the neuron labels)")
    ap.add_argument("--rank-col", default="phf1_intensity_p95",
                    help=f"per-cell column to rank on; one of {', '.join(NUMERIC_COLS)}. "
                         "Cells with no value are dropped.")
    ap.add_argument("--n-cells", type=int, default=3)
    ap.add_argument("--phf1-status", choices=["pos", "neg", "any"], default="pos",
                    help="centre on PHF1+ (pos), PHF1- (neg) or any cells")
    ap.add_argument("--select", choices=["top", "spread"], default="top",
                    help="top = brightest p95; spread = representative sample across p95 range")
    ap.add_argument("--only-cells", default=None, help="comma list of cell_ids to render (overrides selection)")
    ap.add_argument("--crop-label", default="tangle", help="filename stem: <sample>_<label><k>_<cell_id>")
    ap.add_argument("--file-suffix", default="", help="appended to each filename (e.g. _highdef)")
    ap.add_argument("--pad-rank", type=int, default=0,
                    help="zero-pad the rank in the filename to this width (e.g. 3 -> negneuron001)")
    ap.add_argument("--index-tsv", default=None,
                    help="write one row per rendered crop (rank, cell, intensity, file) here")
    ap.add_argument("--no-transcripts", action="store_true",
                    help="skip loading the tx flat file entirely (~1.1 GB gz per slide); "
                         "necessary for large galleries spanning many FOVs")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--window-um", type=float, default=60.0)
    ap.add_argument("--size-cm", type=float, default=4.0)
    ap.add_argument("--scalebar-um", type=float, default=10.0)
    ap.add_argument("--scalebar-lw", type=float, default=2.5)
    ap.add_argument("--hide-scalebar-text", action="store_true")
    ap.add_argument("--phf1-color", default="red")
    ap.add_argument("--dapi-color", default="blue")
    ap.add_argument("--outline-lw", type=float, default=0.6)
    ap.add_argument("--outline-default", default="white", help="outline colour for cells with no celltype")
    ap.add_argument("--target-color", default=None,
                    help="redraw the ranked cell's outline in this colour so it is identifiable "
                         "among its neighbours (default: off)")
    ap.add_argument("--target-lw", type=float, default=1.2)
    ap.add_argument("--highlight-genes", default=None, help="comma/semicolon gene list to colour")
    ap.add_argument("--highlight-color", default="#FFE800")
    ap.add_argument("--highlight-size", type=float, default=2.0)
    ap.add_argument("--highlight-alpha", type=float, default=0.95)
    ap.add_argument("--other-rna-color", default="#8A8A8A")
    ap.add_argument("--other-rna-size", type=float, default=0.5)
    ap.add_argument("--other-rna-alpha", type=float, default=0.5)
    ap.add_argument("--hide-other-rna", action="store_true")
    ap.add_argument("--bg-percentile", type=float, default=50.0)
    ap.add_argument("--dapi-bg-percentile", type=float, default=10.0)
    args = ap.parse_args()
    base = args.base_dir
    args.cells_csv     = args.cells_csv     or os.path.join(base, "PHF1/seu_coords.csv")
    args.phf1_dir      = args.phf1_dir      or os.path.join(base, "PHF1/PHF1")
    args.dapi_dir      = args.dapi_dir      or os.path.join(base, "PHF1/DAPI")
    args.transform_dir = args.transform_dir or os.path.join(base, "PHF1/Transformed")
    args.out_dir       = args.out_dir       or os.path.join(base, "plots/phf1_tangle_crops")
    os.makedirs(args.out_dir, exist_ok=True)

    palette = O.load_palette(args.palette_csv) if args.palette_csv else {}
    ct_set = None
    if args.celltypes:
        ct_set = {c.strip() for c in args.celltypes.replace(";", ",").split(",") if c.strip()}
        print(f"[tangle] restricting to {len(ct_set)} celltypes")
    gene_set = set()
    if args.highlight_genes:
        gene_set = {g.strip() for g in args.highlight_genes.replace(";", ",").split(",") if g.strip()}
        print(f"[tangle] highlighting {len(gene_set)} genes")
    cells = load_cells(args.cells_csv)
    samples = args.samples or sorted(cells.keys())
    cid2ct = {}                       # cell_id -> celltype (for outline colour)
    for s in samples:
        if s in cells:
            for cid, ct in zip(cells[s]["cell_id"], cells[s]["celltype"]):
                cid2ct[cid] = ct

    chosen = {}
    fovs_by_prefix = {}
    for s in samples:
        if s not in cells:
            continue
        d = cells[s]
        _prefix_cache[s] = E.parse_cell_id(d["cell_id"][0])[0]
        picks = []
        only = {c.strip() for c in args.only_cells.split(",")} if args.only_cells else None
        for i in pick_tangles(d, args.n_cells, args.celltype, args.phf1_status, args.select, only,
                              args.rank_col, ct_set):
            cid = d["cell_id"][i]; fov = E.parse_cell_id(cid)[1]
            picks.append(dict(cell_id=cid, fov=fov, xl=d["xl"][i], yl=d["yl"][i],
                              celltype=d["celltype"][i],
                              **{c: d[c][i] for c in NUMERIC_COLS}))
            fovs_by_prefix.setdefault(_prefix_cache[s], set()).add(fov)
        chosen[s] = picks
        note = "" if len(picks) == args.n_cells else f"  (only {len(picks)} qualify)"
        print(f"[tangle] {s}: {len(picks)} cells{note} -> FOVs {sorted({p['fov'] for p in picks})}")

    poly_by_prefix, tx_by_prefix = {}, {}
    for pfx, fovs in fovs_by_prefix.items():
        poly_by_prefix[pfx] = load_polys_for_fovs(base, pfx, fovs)
        tx_by_prefix[pfx] = None if args.no_transcripts else load_tx_for_fovs(base, pfx, fovs)
        n_tx = "skipped" if args.no_transcripts else len(tx_by_prefix[pfx])
        print(f"[tangle] slide prefix {pfx}: {len(poly_by_prefix[pfx])} polygons, "
              f"{n_tx} transcripts in {len(fovs)} FOVs")

    index_rows = []

    for s in samples:
        if s not in chosen or not chosen[s]:
            continue
        pfx = _prefix_cache[s]
        px_um = O.px_size_um_from_M(E.load_napari_affine(os.path.join(args.transform_dir, f"{s}.txt")))
        win_px = int(round(args.window_um / px_um))
        phf1 = tifffile.imread(os.path.join(args.phf1_dir, f"{s}.tif"))
        phf1 = phf1 if phf1.ndim == 2 else phf1[..., 0]
        dapi = tifffile.imread(os.path.join(args.dapi_dir, f"{s}.tif"))
        dapi = dapi if dapi.ndim == 2 else dapi[..., 0]
        phf1, _ = O.subtract_background(phf1, args.bg_percentile)
        dapi, _ = O.subtract_background(dapi, args.dapi_bg_percentile)
        H, W = phf1.shape

        polys = poly_by_prefix.get(pfx, {})
        poly_px = {cid: np.column_stack(O._local_to_phf1px(v[:, 0], v[:, 1], *_geom_for(base, s, fov),
                                                           Minv=_minv(base, args.transform_dir, s)))
                   for cid, (v, fov) in polys.items()}
        # Centroid index: the in-window test below runs once per crop, so with a
        # few hundred crops x ~10^5 slide polygons it has to be vectorised.
        poly_ids = list(poly_px.keys())
        poly_arrs = [poly_px[cid] for cid in poly_ids]
        if poly_ids:
            cent_r = np.array([p[:, 0].mean() for p in poly_arrs])
            cent_c = np.array([p[:, 1].mean() for p in poly_arrs])
            poly_cols = np.array([palette.get(cid2ct.get(cid, ""), args.outline_default)
                                  if palette else args.outline_default for cid in poly_ids], dtype=object)
        else:
            cent_r = cent_c = np.empty(0)
            poly_cols = np.empty(0, dtype=object)
        txdf = tx_by_prefix.get(pfx)
        if txdf is not None and len(txdf):
            trow, tcol = to_phf1px_byfov(base, args.transform_dir, s, txdf["fov"].values,
                                         txdf["x_local_px"].values, txdf["y_local_px"].values)
            t_is_hi = txdf["target"].isin(gene_set).values if gene_set else np.zeros(len(txdf), bool)
        else:
            trow = tcol = np.array([]); t_is_hi = np.array([], bool)

        for k, p in enumerate(chosen[s], 1):
            cr, cc = O._local_to_phf1px(np.array([p["xl"]]), np.array([p["yl"]]),
                                        *_geom_for(base, s, p["fov"]), Minv=_minv(base, args.transform_dir, s))
            cr, cc = float(cr[0]), float(cc[0])
            r0 = int(max(0, cr - win_px / 2)); r1 = int(min(H, cr + win_px / 2))
            c0 = int(max(0, cc - win_px / 2)); c1 = int(min(W, cc + win_px / 2))
            if len(cent_r):
                inbox = np.where((cent_r >= r0) & (cent_r < r1) & (cent_c >= c0) & (cent_c < c1))[0]
            else:
                inbox = np.empty(0, int)
            polys_crop = [(poly_arrs[j] - np.array([r0, c0]), poly_cols[j]) for j in inbox]
            tgt = poly_px.get(p["cell_id"])
            target_poly = None if tgt is None else tgt - np.array([r0, c0])
            if len(trow):
                inwin = (trow >= r0) & (trow < r1) & (tcol >= c0) & (tcol < c1)
                hi = inwin & t_is_hi
                oth = inwin & ~t_is_hi if gene_set else inwin   # no highlight list -> all "other"
                tx_hi = np.column_stack([trow[hi] - r0, tcol[hi] - c0])
                tx_other = np.column_stack([trow[oth] - r0, tcol[oth] - c0])
            else:
                tx_hi = np.empty((0, 2)); tx_other = np.empty((0, 2))
            rank = f"{k:0{args.pad_rank}d}" if args.pad_rank else str(k)
            out = os.path.join(args.out_dir, f"{s}_{args.crop_label}{rank}_{p['cell_id']}{args.file_suffix}.png")
            index_rows.append(dict(sample_id=s, rank=k, cell_id=p["cell_id"],
                                   celltype=p.get("celltype", ""), fov=p["fov"],
                                   x_FOV_px=p["xl"], y_FOV_px=p["yl"],
                                   crop_row=round(cr, 1), crop_col=round(cc, 1),
                                   window_um=args.window_um, px_um=round(px_um, 5),
                                   n_cells_in_crop=len(polys_crop),
                                   target_outline=target_poly is not None,
                                   file=os.path.basename(out),
                                   **{c: p.get(c, np.nan) for c in NUMERIC_COLS}))
            render_crop(out, phf1[r0:r1, c0:c1], dapi[r0:r1, c0:c1], polys_crop, tx_hi, tx_other,
                        px_um, args.scalebar_um, args.phf1_color, args.dapi_color, args.outline_lw,
                        args.highlight_color, args.highlight_size, args.highlight_alpha,
                        args.other_rna_color, args.other_rna_size, args.other_rna_alpha,
                        not args.hide_other_rna, args.size_cm, args.dpi,
                        args.scalebar_lw, not args.hide_scalebar_text,
                        target_poly, args.target_color, args.target_lw)
            print(f"[tangle] {s} tangle{k} {p['cell_id']}: {len(polys_crop)} cells, "
                  f"{len(tx_hi)} highlight + {len(tx_other)} other tx -> {out}")
        del phf1, dapi

    if args.index_tsv and index_rows:
        cols = ["sample_id", "rank", "cell_id", "celltype", "fov"] + list(NUMERIC_COLS) + \
               ["x_FOV_px", "y_FOV_px", "crop_row", "crop_col", "window_um", "px_um",
                "n_cells_in_crop", "target_outline", "file"]
        os.makedirs(os.path.dirname(os.path.abspath(args.index_tsv)), exist_ok=True)
        pd.DataFrame(index_rows)[cols].to_csv(args.index_tsv, sep="\t", index=False, na_rep="NA")
        print(f"[tangle] wrote {len(index_rows)} rows -> {args.index_tsv}")


if __name__ == "__main__":
    main()
