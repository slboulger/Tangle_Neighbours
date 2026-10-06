#!/usr/bin/env python
# ---------------------------------------------------------------------------
# imc_phf1_nucleus_crops_singles.py
#
# Figure panels: 5D
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# imc_phf1_nucleus_crops_singles.py
#
# Plain single-channel exports of selected cells, for assembling figures by hand. Cells are
# selected from the candidate pairs in plots/imc_phf1_nucleus_crops/pairs/pairs_*.pdf (labels
# are ROI:object) and passed with --cells. One image per cell x channel, no text, 10 um scale
# bar bottom-left on every image:
#   DNA    greyscale
#   mask   the target nucleus' DNA Otsu mask, white fill on black (other nuclei not drawn)
#   PHF1   red
#   merge  PHF1 (red) + mask OUTLINE (white, drawn just inside the mask's pixel edges, so it
#          never covers signal outside the nucleus; an outline rather than a fill, so PHF1
#          over the nucleus stays visible).
#
# Display processing and crop windows are identical to imc_phf1_nucleus_crops.py and
# imc_phf1_nucleus_crops_pairs.py (shared code in imc_crop_utils.py): Supp1 recipe per ROI,
# 32 x 32 um window centred on the regionprops centroid, zero-padded at the image edge.
# Files: <roi>_<obj>_<cluster>_<PHF1pos|PHF1neg>_<channel>.png / .tif, nearest-neighbour
# upscaled x UP (TIFF resolution tag = UP px/um, so the physical scale is preserved).
#
#   python python/imc_phf1_nucleus_crops_singles.py \
#       --cells BBN110.26073_1:47,BBN110.26073_1:1345

import argparse, os, sys
import numpy as np, pandas as pd, tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import imc_crop_utils as U

ap = argparse.ArgumentParser()
ap.add_argument("--cells", default="BBN110.26073_1:47,BBN110.26073_1:1345")
ap.add_argument("--half", type=int, default=16)
ap.add_argument("--scalebar_um", type=float, default=10)
ap.add_argument("--up", type=int, default=20, help="nearest-neighbour upscale factor")
ap.add_argument("--outline_px", type=int, default=3, help="outline width in upscaled px")
ap.add_argument("--out", default=U.OUTROOT + "/singles")
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)
N = 2 * args.half

cells = [(c.split(":")[0], int(c.split(":")[1])) for c in args.cells.split(",")]
rois = sorted({r for r, _ in cells})
lab = U.load_labels(rois, f"{U.OUTROOT}/cell_labels_pairs.tsv")

R, log = {}, []
for roi, obj in cells:
    if roi not in R:
        R[roi] = U.load_roi(roi, lab)
        R[roi]["pad"] = dict(disp={k: np.pad(v, args.half) for k, v in R[roi]["disp"].items()},
                             mask=np.pad(R[roi]["mask"], args.half))
    row = lab[(lab.sample_id == roi) & (lab.ObjectNumber == obj)]
    if len(row) != 1:
        raise SystemExit(f"{roi}:{obj} not found in spe")
    row = row.iloc[0]
    y0, x0, inb, cy, cx = U.crop_window(R[roi], obj, args.half)
    P = R[roi]["pad"]
    t, m = U.tiles(dict(disp=P["disp"], mask=P["mask"]), obj, y0 + args.half, x0 + args.half, N)
    assert m.sum() == row.area, "mask area != measured area"

    clus = row.celltype_clusters.replace("Excitatory neuron cluster ", "Exc").split(" (")[0]
    marker = row.celltype_clusters.split("(")[-1].rstrip(")").replace(", ", "-").replace(" ", "")
    ph = "PHF1pos" if row.PHF1_Otsu == "PHF1_pos" else "PHF1neg"
    stem = f"{args.out}/{roi}_{obj}_{clus}_{marker}_{ph}"
    imgs = {
        "DNA":   U.upscale_with_outline(t["DNA"], None, args.up, args.outline_px),
        "mask":  U.upscale_with_outline(t["mask"], None, args.up, args.outline_px),
        "PHF1":  U.upscale_with_outline(t["PHF1"], None, args.up, args.outline_px),
        "merge": U.upscale_with_outline(t["PHF1"], m, args.up, args.outline_px),
    }
    for ch, big in imgs.items():
        big = U.burn_scalebar(big, args.scalebar_um, args.up)
        plt.imsave(f"{stem}_{ch}.png", big)
        tifffile.imwrite(f"{stem}_{ch}.tif", big, resolution=(args.up / U.PX_UM, args.up / U.PX_UM),
                         metadata={"unit": "um"})
    log.append(dict(roi=roi, ObjectNumber=obj, cell_id=row.cell_id,
                    celltype_clusters=row.celltype_clusters, PHF1_Otsu=row.PHF1_Otsu,
                    area_um2=row.area, centroid_row=cy, centroid_col=cx, crop_row0=y0,
                    crop_col0=x0, crop_um=N, edge_padded=not inb, scalebar_um=args.scalebar_um,
                    upscale=args.up, outline_px=args.outline_px))

src = f"{args.out}/source_data_singles.tsv"
new = pd.DataFrame(log)
if os.path.exists(src):  # accumulate across calls, newest export of a cell wins
    old = pd.read_csv(src, sep="\t")
    keep = ~old.set_index(["roi", "ObjectNumber"]).index.isin(new.set_index(["roi", "ObjectNumber"]).index)
    new = pd.concat([old[keep], new], ignore_index=True)
new.to_csv(src, sep="\t", index=False)
print(f"wrote {len(log)} cells x 4 channels -> {args.out}")
