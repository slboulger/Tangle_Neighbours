#!/usr/bin/env python
# ---------------------------------------------------------------------------
# extract_phf1_intensity.py
#
# Upstream pipeline - builds the objects every panel reads
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
"""
extract_phf1_intensity.py
=========================

Steps 2-3 of the PHF1 immunofluorescence intensity pipeline.

For every CosMx cell, aggregate the PHF1 (phospho-tau) post-stain IF signal over
that cell's *segmentation mask* (CellLabels_F*.tif, uint16, pixel value = local cell
label) and write a per-cell table that R/add_phf1_intensity.R joins onto seu_PHF1.rds.

Coordinate transform (EXACT, deterministic -- no fitting, no image registration)
--------------------------------------------------------------------------------
The PHF1 images were aligned in napari with the napari_CosMx plugin, which loaded
the per-slide stitched folders in `stitched_slides/{1b,2a}/`. That plugin places all
CosMx layers in a millimetre "world" frame as

    world_layer = scale * mosaic_px + translate ,   scale = mm_per_px,
                  translate = (min Y_mm, -max X_mm)           (non-DASH instruments)

Working that through the plugin's FOV stitching collapses, for a cell pixel at
FOV-local (col,row) in FOV `f`, to:

    world0 =  Y_mm[f] + mm_per_px * row
    world1 = -X_mm[f] + mm_per_px * col

with `X_mm,Y_mm` the FOV offsets and `mm_per_px = scale_um/1000`, BOTH read from the
exact values the plugin used (`stitched_slides/{slide}/images/.zattrs`). The PHF1
transform `M` in `PHF1/Transformed/<sample>.txt` maps PHF1 px (row,col) -> world, so

    PHF1px(row,col) = M^-1 . [world0, world1, 1]

i.e. per FOV the full map is  T = M^-1 . A_fov , with

    A_fov = [[0,         mm_per_px,  Y_mm[f]],
             [mm_per_px, 0,         -X_mm[f]],
             [0,         0,          1      ]]   acting on [col, row, 1].

Validated against the manual PHF1+ calls (e.g. sample BBN_24895: AUC 0.79, median masked
intensity 52 in PHF1+ vs 12 in PHF1-). The CellLabels pixel frame equals the metadata
(CenterX_local_px, CenterY_local_px) frame (verified, no flip), and the napari layer
`rotate` is 0 by default. The per-sample QC reports AUC vs the manual calls so any
sample where these assumptions fail is flagged, not silently used.

Run under the `spatial_env` conda env (HPC). Needs numpy + tifffile (+ json/csv from
stdlib; matplotlib optional for diagnostic PNGs).

Usage:
    python extract_phf1_intensity.py --mode extract        # full mask aggregation + QC
    python extract_phf1_intensity.py --mode validate       # fast: per-sample centroid AUC + overlay PNGs
    python extract_phf1_intensity.py --mode both
    # optional: --samples BBN_24895 BBN00638047   --max-fovs 5   (quick tests)
"""

import argparse
import csv
import json
import os
import sys

import numpy as np

try:
    import tifffile
except ImportError:
    sys.exit("tifffile not available - run under the spatial_env conda environment.")


# ----------------------------------------------------------------------------- #
# Slide configuration. cell_id is c_<prefix>_<fov>_<cellID>; the prefix is INVERTED
# vs SlideLabel (slide 1b -> prefix 2, slide 2a -> prefix 1).
# ----------------------------------------------------------------------------- #
PREFIX_TO_SLIDE = {2: "1b", 1: "2a"}
CELLSTATS = {
    "1b": "DecodedFiles/DecodedFiles/IGFQ002102_matthews_1312026_CosMx_RNA_1b/"
          "20260114_151706_S1/CellStatsDir",
    "2a": "DecodedFiles/DecodedFiles/IGFQ002102_matthews_1312026_CosMx_RNA_2a/"
          "20260114_151706_S2/CellStatsDir",
}
ZATTRS = "stitched_slides/{slide}/images/.zattrs"   # exact fov_offsets + scale_um the plugin used

SAT_VALUE = 255          # 8-bit saturation
DAPI_ZERO_THRESH = 5     # cells with dapi mean below this count as "near-zero" (QC only)


# ----------------------------------------------------------------------------- #
# Helpers
# ----------------------------------------------------------------------------- #
def parse_cell_id(cid):
    """'c_2_1_345' -> (prefix=2, fov=1, cellID=345)."""
    p = cid.split("_")
    return int(p[1]), int(p[2]), int(p[3])


def load_napari_affine(path):
    """3x3 matrix mapping PHF1 px (row,col,1) -> napari world (w0,w1,1)."""
    M = np.loadtxt(path, delimiter=",")
    if M.shape != (3, 3):
        raise ValueError(f"{path}: expected 3x3 affine, got {M.shape}")
    return M


def load_slide_geometry(base_dir, slide):
    """Read the exact mm_per_px and per-FOV (X_mm, Y_mm) the napari_CosMx plugin used."""
    z = json.load(open(os.path.join(base_dir, ZATTRS.format(slide=slide))))
    cm = z["CosMx"]
    if "scale_um" in cm:
        mm_per_px = cm["scale_um"] / 1000.0
    elif "scale_mm" in cm:
        mm_per_px = cm["scale_mm"]
    else:
        raise ValueError(f"no scale_um/scale_mm in {slide} .zattrs")
    fo = cm["fov_offsets"]
    # fov_offsets is a dict of columns -> {row_index: value}; index FOV->(X_mm,Y_mm)
    fovs = fo["FOV"]; xs = fo["X_mm"]; ys = fo["Y_mm"]
    offsets = {int(fovs[k]): (float(xs[k]), float(ys[k])) for k in fovs}
    return mm_per_px, offsets


def fov_affine(mm_per_px, X_mm, Y_mm):
    """3x3 mapping FOV-local mask px [col,row,1] -> napari world [w0,w1,1]."""
    return np.array([[0.0,       mm_per_px,  Y_mm],
                     [mm_per_px, 0.0,       -X_mm],
                     [0.0,       0.0,        1.0]])


def auc_mann_whitney(scores, pos_mask):
    """ROC-AUC via the rank-sum statistic (tie-aware). No scipy/sklearn dependency."""
    scores = np.asarray(scores, dtype=float)
    pos_mask = np.asarray(pos_mask, dtype=bool)
    n_pos = int(pos_mask.sum()); n_neg = int((~pos_mask).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=float)
    ranks[order] = np.arange(1, len(scores) + 1)
    s_sorted = scores[order]
    i = 0
    while i < len(s_sorted):
        j = i
        while j + 1 < len(s_sorted) and s_sorted[j + 1] == s_sorted[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + 1 + j + 1) / 2.0
        i = j + 1
    return (ranks[pos_mask].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def mask_label_centroids(lab):
    """Returns (centroids dict label->(col_mean,row_mean), ys, xs, vals, n_tot, K)."""
    ys, xs = np.nonzero(lab)
    vals = lab[ys, xs].astype(np.int64)
    K = int(vals.max()) + 1 if vals.size else 1
    n = np.bincount(vals, minlength=K)
    sx = np.bincount(vals, weights=xs.astype(np.float64), minlength=K)
    sy = np.bincount(vals, weights=ys.astype(np.float64), minlength=K)
    cent = {k: (sx[k] / n[k], sy[k] / n[k]) for k in range(1, K) if n[k] > 0}
    return cent, ys, xs, vals, n, K


def _read_image(path):
    img = tifffile.imread(path)
    return img if img.ndim == 2 else img[..., 0]


def _write_qc(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


# ----------------------------------------------------------------------------- #
# seu_coords -> cohort maps + manual PHF1 labels
# ----------------------------------------------------------------------------- #
def load_seu_coords(path):
    """Returns:
       sample_of_fov[(prefix,fov)] -> sample_id
       fovs_of_sample[sample]      -> sorted [(prefix,fov)]
       phf1_of_cell[cell_id]       -> bool (manual PHF1 call)
       centroids_local[sample]     -> dict (prefix,fov,cellID) -> (x_FOV_px, y_FOV_px)  [validate]
    """
    sample_of_fov, phf1_of_cell, cent_local = {}, {}, {}
    have_local = None
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh)
        have_local = "x_FOV_px" in rdr.fieldnames and "y_FOV_px" in rdr.fieldnames
        for r in rdr:
            cid = r["cell_id"]; pfx, fov, cellID = parse_cell_id(cid)
            sample_of_fov[(pfx, fov)] = r["sample_id"]
            phf1_of_cell[cid] = str(r["PHF1"]).strip().upper() in ("TRUE", "1")
            if have_local:
                cent_local.setdefault(r["sample_id"], {})[(pfx, fov, cellID)] = (
                    float(r["x_FOV_px"]), float(r["y_FOV_px"]))
    fovs_of_sample = {}
    for key, samp in sample_of_fov.items():
        fovs_of_sample.setdefault(samp, []).append(key)
    for s in fovs_of_sample:
        fovs_of_sample[s].sort()
    return sample_of_fov, fovs_of_sample, phf1_of_cell, cent_local


# ----------------------------------------------------------------------------- #
# Validate (fast): centroid AUC + overlay PNG per sample, using the deterministic map
# ----------------------------------------------------------------------------- #
def validate(args, fovs_of_sample, phf1_of_cell, cent_local):
    geom = {}
    qc_rows = []
    os.makedirs(args.diagnostics_dir, exist_ok=True)
    samples = args.samples or sorted(fovs_of_sample.keys())
    for samp in samples:
        if samp not in cent_local:
            print(f"[validate] {samp}: no local centroids in seu_coords (re-export with x_FOV_px), skipping")
            continue
        phf1_path = os.path.join(args.phf1_dir, f"{samp}.tif")
        tf_path = os.path.join(args.transform_dir, f"{samp}.txt")
        if not (os.path.exists(phf1_path) and os.path.exists(tf_path)):
            print(f"[validate] {samp}: missing image/transform, skipping")
            continue
        Minv = np.linalg.inv(load_napari_affine(tf_path))
        img = _read_image(phf1_path); H, W = img.shape
        rows_px, cols_px, labels = [], [], []
        for (pfx, fov, cellID), (xl, yl) in cent_local[samp].items():
            slide = PREFIX_TO_SLIDE[pfx]
            if slide not in geom:
                geom[slide] = load_slide_geometry(args.base_dir, slide)
            mm_per_px, offsets = geom[slide]
            if fov not in offsets:
                continue
            X_mm, Y_mm = offsets[fov]
            w0 = Y_mm + mm_per_px * yl; w1 = -X_mm + mm_per_px * xl
            px = Minv @ np.array([w0, w1, 1.0])
            rows_px.append(px[0]); cols_px.append(px[1])
            labels.append(phf1_of_cell.get(f"c_{pfx}_{fov}_{cellID}", False))
        rp = np.array(rows_px); cp = np.array(cols_px); lab = np.array(labels, dtype=bool)
        inb = (rp >= 0) & (rp < H) & (cp >= 0) & (cp < W)
        ri = np.clip(np.rint(rp[inb]).astype(int), 0, H - 1)
        ci = np.clip(np.rint(cp[inb]).astype(int), 0, W - 1)
        auc = auc_mann_whitney(img[ri, ci].astype(float), lab[inb])
        qc_rows.append(dict(sample_id=samp, n_cells=len(lab), n_phf1_pos=int(lab.sum()),
                            frac_in_bounds=round(float(inb.mean()), 4),
                            centroid_auc=round(auc, 4) if not np.isnan(auc) else "NA"))
        flag = " [LOW AUC]" if (not np.isnan(auc) and auc < 0.65) else ""
        print(f"[validate] {samp}: frac_in={inb.mean():.3f} centroid_AUC={auc:.3f} nPHF1+={int(lab.sum())}{flag}")
        _overlay_png(args, samp, img, rp, cp, lab, inb, auc)
    _write_qc(os.path.join(args.out_dir, "phf1_validate_qc.csv"), qc_rows)
    print(f"[validate] wrote phf1_validate_qc.csv + diagnostics/*.png")


def _overlay_png(args, samp, img, rp, cp, lab, inb, auc):
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except Exception:
        return
    H, W = img.shape; step = max(1, max(H, W) // 1500); sub = img[::step, ::step]
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(sub, cmap="gray", vmax=max(1, np.percentile(sub, 99.5)))
    neg = np.where(inb & (~lab))[0]
    if len(neg) > 2000:
        neg = neg[:: max(1, len(neg) // 2000)]
    pos = inb & lab
    ax.scatter(cp[neg] / step, rp[neg] / step, s=1, c="cyan", alpha=0.3, label="PHF1- cells")
    ax.scatter(cp[pos] / step, rp[pos] / step, s=6, c="red", label="PHF1+ (manual)")
    ax.set_title(f"{samp}  centroid AUC={auc:.3f}"); ax.legend(loc="upper right", fontsize=7); ax.axis("off")
    fig.tight_layout(); fig.savefig(os.path.join(args.diagnostics_dir, f"{samp}_overlay.png"), dpi=120)
    plt.close(fig)


# ----------------------------------------------------------------------------- #
# Extract (full mask aggregation)
# ----------------------------------------------------------------------------- #
def extract(args, fovs_of_sample, phf1_of_cell):
    geom = {}
    rows, qc_rows = [], []
    samples = args.samples or sorted(fovs_of_sample.keys())

    for samp in samples:
        phf1_path = os.path.join(args.phf1_dir, f"{samp}.tif")
        dapi_path = os.path.join(args.dapi_dir, f"{samp}.tif")
        tf_path = os.path.join(args.transform_dir, f"{samp}.txt")
        if not (os.path.exists(phf1_path) and os.path.exists(tf_path)):
            print(f"[extract] {samp}: missing image/transform, skipping")
            continue
        Minv = np.linalg.inv(load_napari_affine(tf_path))
        phf1 = _read_image(phf1_path); H, W = phf1.shape
        dapi = _read_image(dapi_path) if os.path.exists(dapi_path) else None

        n_cells_samp, oob_cells = 0, 0
        val_mean, val_p95, val_lab = [], [], []   # per-sample acceptance AUC vs manual PHF1
        fovs = fovs_of_sample[samp]
        if args.max_fovs:
            fovs = fovs[: args.max_fovs]

        for (pfx, fov) in fovs:
            slide = PREFIX_TO_SLIDE[pfx]
            if slide not in geom:
                geom[slide] = load_slide_geometry(args.base_dir, slide)
            mm_per_px, offsets = geom[slide]
            if fov not in offsets:
                print(f"[extract] {samp} c_{pfx}_{fov}: FOV not in offsets, skipping")
                continue
            X_mm, Y_mm = offsets[fov]
            mask_path = os.path.join(args.base_dir, CELLSTATS[slide],
                                     f"FOV{fov:05d}", f"CellLabels_F{fov:05d}.tif")
            if not os.path.exists(mask_path):
                print(f"[extract] {samp} c_{pfx}_{fov}: mask not found, skipping FOV")
                continue
            lab = tifffile.imread(mask_path).astype(np.int64)
            _, ys, xs, vals, n_tot, K = mask_label_centroids(lab)
            if K <= 1:
                continue

            T = Minv @ fov_affine(mm_per_px, X_mm, Y_mm)     # [col,row,1] -> PHF1 px [row,col,1]
            src = np.vstack([xs.astype(np.float64), ys.astype(np.float64), np.ones(len(xs))])
            dst = T @ src
            rp = np.rint(dst[0]).astype(np.int64); cp = np.rint(dst[1]).astype(np.int64)
            inb = (rp >= 0) & (rp < H) & (cp >= 0) & (cp < W)

            pv = np.zeros(len(vals)); pv[inb] = phf1[rp[inb], cp[inb]].astype(np.float64)
            dv = np.zeros(len(vals))
            if dapi is not None:
                dv[inb] = dapi[rp[inb], cp[inb]].astype(np.float64)

            n_all = n_tot.astype(np.float64)                              # all mask px (incl. OOB)
            n_in = np.bincount(vals[inb], minlength=K).astype(np.float64)  # in-bounds px
            sum_p = np.bincount(vals, weights=pv, minlength=K)            # integrated DN (OOB=0)
            sumsq = np.bincount(vals, weights=pv * pv, minlength=K)
            sum_d = np.bincount(vals, weights=dv, minlength=K)

            # one per-label 8-bit histogram over in-bounds px -> percentiles / max / bright-fractions
            hist = np.zeros((K, 256))
            if inb.any():
                comb = vals[inb] * 256 + pv[inb].astype(np.int64)
                hist = np.bincount(comb, minlength=K * 256).reshape(K, 256).astype(np.float64)
            csum = np.cumsum(hist, axis=1)
            tot = csum[:, -1]                                             # == n_in

            def q_bin(qq):                                               # exact 8-bit percentile
                return np.where(tot > 0, np.argmax(csum >= (qq * tot)[:, None], axis=1).astype(float), np.nan)
            median_p = q_bin(0.50); p90 = q_bin(0.90); p95 = q_bin(0.95); p99 = q_bin(0.99)
            max_p = np.where(tot > 0, (255 - np.argmax(hist[:, ::-1] > 0, axis=1)).astype(float), np.nan)
            with np.errstate(invalid="ignore", divide="ignore"):
                mean_p = np.where(n_in > 0, sum_p / np.maximum(n_in, 1), np.nan)
                sd_p = np.sqrt(np.clip(np.where(n_in > 0, sumsq / np.maximum(n_in, 1) - mean_p ** 2, np.nan), 0, None))
                mean_d = np.where(n_in > 0, sum_d / np.maximum(n_in, 1), np.nan)
                frac_in_b = np.where(n_all > 0, n_in / np.maximum(n_all, 1), np.nan)
                frac_sat = np.where(tot > 0, hist[:, 255] / np.maximum(tot, 1), np.nan)
                frac_gt64 = np.where(tot > 0, (tot - csum[:, 64]) / np.maximum(tot, 1), np.nan)
                frac_gt128 = np.where(tot > 0, (tot - csum[:, 128]) / np.maximum(tot, 1), np.nan)

            def fmt(a, k, nd=4):
                return round(float(a[k]), nd) if not np.isnan(a[k]) else ""

            def int_or_blank(a, k):
                return "" if np.isnan(a[k]) else int(a[k])

            for k in range(1, K):
                if n_all[k] <= 0:
                    continue
                cid = f"c_{pfx}_{fov}_{k}"
                if n_in[k] == 0:
                    oob_cells += 1
                rows.append((
                    cid, samp, pfx, fov, k,
                    fmt(mean_p, k), int_or_blank(median_p, k), fmt(sd_p, k),
                    int_or_blank(p90, k), int_or_blank(p95, k), int_or_blank(p99, k), int_or_blank(max_p, k),
                    round(float(sum_p[k]), 1),
                    int(n_all[k]), fmt(frac_in_b, k), fmt(frac_sat, k),
                    fmt(frac_gt64, k), fmt(frac_gt128, k),
                    fmt(mean_d, k),
                ))
                n_cells_samp += 1
                if (cid in phf1_of_cell) and frac_in_b[k] >= 0.5 and not np.isnan(mean_p[k]):
                    val_mean.append(mean_p[k]); val_p95.append(p95[k]); val_lab.append(phf1_of_cell[cid])

        auc_mean = auc_mann_whitney(val_mean, val_lab) if val_lab else float("nan")
        auc_p95 = auc_mann_whitney(val_p95, val_lab) if val_lab else float("nan")
        qc_rows.append(dict(sample_id=samp, n_cells=n_cells_samp, n_fovs=len(fovs),
                            n_cells_fully_oob=oob_cells,
                            n_phf1_pos=int(np.sum(val_lab)) if val_lab else 0,
                            mask_mean_auc=round(auc_mean, 4) if not np.isnan(auc_mean) else "NA",
                            mask_p95_auc=round(auc_p95, 4) if not np.isnan(auc_p95) else "NA"))
        print(f"[extract] {samp}: {n_cells_samp} cells / {len(fovs)} FOVs (oob={oob_cells})  "
              f"AUC vs manual PHF1: mean={auc_mean:.3f} p95={auc_p95:.3f}")
        del phf1
        if dapi is not None:
            del dapi

    header = ["cell_id", "sample_id", "slide_id", "fov", "cell_ID",
              "phf1_intensity_mean", "phf1_intensity_median", "phf1_intensity_sd",
              "phf1_intensity_p90", "phf1_intensity_p95", "phf1_intensity_p99", "phf1_intensity_max",
              "phf1_intensity_sum", "phf1_n_mask_px", "phf1_frac_in_bounds", "phf1_frac_saturated",
              "phf1_frac_gt64", "phf1_frac_gt128", "dapi_intensity_mean"]
    with open(args.out_csv, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(header); w.writerows(rows)
    _write_qc(os.path.join(args.out_dir, "phf1_extraction_qc.csv"), qc_rows)
    print(f"[extract] wrote {args.out_csv} ({len(rows)} cells) and phf1_extraction_qc.csv")


# ----------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-dir", default="<PROJECT_ROOT>/phf1_v2")
    ap.add_argument("--mode", choices=["validate", "extract", "both"], default="both")
    ap.add_argument("--phf1-dir", default=None)
    ap.add_argument("--dapi-dir", default=None)
    ap.add_argument("--transform-dir", default=None)
    ap.add_argument("--seu-coords", default=None)
    ap.add_argument("--out-csv", default=None)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--diagnostics-dir", default=None)
    ap.add_argument("--samples", nargs="*", default=None)
    ap.add_argument("--max-fovs", type=int, default=0)
    args = ap.parse_args()

    base = args.base_dir
    args.phf1_dir        = args.phf1_dir        or os.path.join(base, "PHF1/PHF1")
    args.dapi_dir        = args.dapi_dir        or os.path.join(base, "PHF1/DAPI")
    args.transform_dir   = args.transform_dir   or os.path.join(base, "PHF1/Transformed")
    args.seu_coords      = args.seu_coords      or os.path.join(base, "PHF1/seu_coords.csv")
    args.out_dir         = args.out_dir         or os.path.join(base, "PHF1")
    args.out_csv         = args.out_csv         or os.path.join(base, "PHF1/phf1_intensity_per_cell.csv")
    args.diagnostics_dir = args.diagnostics_dir or os.path.join(base, "PHF1/diagnostics")
    os.makedirs(args.out_dir, exist_ok=True)

    sample_of_fov, fovs_of_sample, phf1_of_cell, cent_local = load_seu_coords(args.seu_coords)

    if args.mode in ("validate", "both"):
        validate(args, fovs_of_sample, phf1_of_cell, cent_local)
    if args.mode in ("extract", "both"):
        extract(args, fovs_of_sample, phf1_of_cell)


if __name__ == "__main__":
    main()
