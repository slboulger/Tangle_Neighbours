#!/usr/bin/env python
# ---------------------------------------------------------------------------
# imc_phf1_nucleus_crops_pairs.py
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
# imc_phf1_nucleus_crops_pairs.py
#
# Browsing PDF of PHF1- / PHF1+ nucleus pairs, matched on excitatory cluster, for choosing
# cells to export with imc_phf1_nucleus_crops_singles.py. Not a figure.
#
# PAIRING (deterministic): every PHF1_Otsu+ cell of the 4 excitatory clusters in the chosen
# ROIs is shown. Within each ROI x cluster, PHF1+ cells are sorted by nucleus area (largest
# first) and the i-th is paired with the i-th PHF1- cell of the SAME ROI and cluster ranked by
# closeness to the PHF1- group median area -- so each partner is a typical PHF1- nucleus, and
# both members of a pair share one per-ROI display range. PHF1- partners must have a full crop
# window; PHF1+ cells near the image edge are zero-padded so none is dropped.
#
# Each half-row: DNA (grey) | mask (white fill) | PHF1 (red) | merge (PHF1 red + white mask
# outline drawn just inside the mask's pixel edges). 10 um bar bottom-left on every tile.
# Labels give ROI:object (the --cells argument for the singles script) and area. Tiles are drawn
# at native 1 um/px with the outline and scale bar as VECTOR overlays, so they stay crisp at any
# zoom (a burned-in pixel outline aliases away when the tile is downsampled). Default ROIs: the three BBN110.26073 (Braak 3-4) ROIs; BBN_9389_2 is Braak 0-1
# with no PHF1+ cells.
#
# OUTPUTS plots/imc_phf1_nucleus_crops/pairs/
#   pairs_<rois>.pdf, source_data_pairs.tsv, stats_pairs_display_ranges.tsv
#
#   python python/imc_phf1_nucleus_crops_pairs.py

import argparse, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from skimage import measure

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import imc_crop_utils as U

ap = argparse.ArgumentParser()
ap.add_argument("--rois", default="BBN110.26073_1,BBN110.26073_3,BBN110.26073_4")
ap.add_argument("--half", type=int, default=16)
ap.add_argument("--scalebar_um", type=float, default=10)
ap.add_argument("--rows_per_page", type=int, default=8)
args = ap.parse_args()
ROIS = args.rois.split(",")
OUT = U.OUTROOT + "/pairs"
os.makedirs(OUT, exist_ok=True)
N = 2 * args.half

CLUSTERS = ["Excitatory neuron cluster 1 (MAP2)", "Excitatory neuron cluster 2 (CALB1)",
            "Excitatory neuron cluster 3 (RORB)", "Excitatory neuron cluster 4 (RELN)"]
short = lambda cl: cl.replace("Excitatory neuron cluster ", "Exc").replace(" (", " ").rstrip(")")

lab = U.load_labels(ROIS, f"{U.OUTROOT}/cell_labels_pairs.tsv")
R, ranges = {}, []
for r in ROIS:
    R[r] = U.load_roi(r, lab)
    ranges += R[r]["ranges"]
    # zero-pad so edge cells can still be centred; pad == half, coordinates shift by half
    R[r]["pad"] = dict(disp={k: np.pad(v, args.half) for k, v in R[r]["disp"].items()},
                       mask=np.pad(R[r]["mask"], args.half))

rows = []
for cl in CLUSTERS:
    for r in ROIS:
        sub = lab[(lab.sample_id == r) & (lab.celltype_clusters == cl) & lab.PHF1_Otsu.notna()]
        pos = sub[sub.PHF1_Otsu == "PHF1_pos"].sort_values("area", ascending=False)
        neg = sub[sub.PHF1_Otsu == "PHF1_neg"].copy()
        neg["inb"] = [U.crop_window(R[r], o, args.half)[2] for o in neg.ObjectNumber]
        neg = neg[neg.inb]
        med = neg.area.median()
        neg = neg.assign(dev=(neg.area - med).abs()).sort_values(["dev", "ObjectNumber"])
        if len(neg) < len(pos):
            sys.exit(f"{r} {cl}: fewer PHF1- partners than PHF1+ cells")
        for i, (p, q) in enumerate(zip(pos.itertuples(), neg.head(len(pos)).itertuples())):
            rec = dict(cluster=cl, roi=r, pair=i + 1, BraakGroup=p.BraakGroup,
                       neg_median_area=med, pos_median_area=pos.area.median(),
                       n_pos=len(pos), n_neg=len(sub) - len(pos))
            for tag, c in (("neg", q), ("pos", p)):
                y0, x0, inb, cy, cx = U.crop_window(R[r], c.ObjectNumber, args.half)
                rec.update({f"{tag}_obj": int(c.ObjectNumber), f"{tag}_cell_id": c.cell_id,
                            f"{tag}_area_um2": c.area, f"{tag}_eccentricity": c.eccentricity,
                            f"{tag}_crop_row0": y0, f"{tag}_crop_col0": x0,
                            f"{tag}_in_bounds": inb})
            rows.append(rec)
pairs = pd.DataFrame(rows)
pairs.insert(0, "pair_id", [f"{short(c).split()[0]}-{r}-{i:02d}"
                            for c, r, i in zip(pairs.cluster, pairs.roi, pairs.pair)])


def render(r, obj, y0, x0):
    P = R[r]["pad"]
    t, m = U.tiles(dict(disp=P["disp"], mask=P["mask"]), obj, y0 + args.half, x0 + args.half, N)
    t["merge"] = t["PHF1"]
    return t, m


def vector_outline(ax, m, lw=0.7):
    big = np.kron(np.pad(m, 1).astype(float), np.ones((10, 10)))
    for c in measure.find_contours(big, 0.5):
        ax.plot(c[:, 1] / 10 - 1.5, c[:, 0] / 10 - 1.5, color="white", lw=lw)


def vector_scalebar(ax):
    L = args.scalebar_um / U.PX_UM
    ax.add_patch(Rectangle((0.5, N - 2.5), L, 0.9, color="white", lw=0))


plt.rcParams.update({"pdf.fonttype": 42, "font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 6})
CHS = ["DNA", "mask", "PHF1", "merge"]
tw, gap, lab_h, W_in = 0.8, 0.25, 0.22, 8.27
pdf_path = f"{OUT}/pairs_{'_'.join(ROIS)}.pdf"
with PdfPages(pdf_path) as pdf:
    for cl in CLUSTERS:
        pc = pairs[pairs.cluster == cl]
        for start in range(0, len(pc), args.rows_per_page):
            chunk = pc.iloc[start:start + args.rows_per_page]
            H_in = 0.7 + len(chunk) * (tw + lab_h + 0.06)
            fig = plt.figure(figsize=(W_in, H_in))
            fig.text(0.5, 1 - 0.25 / H_in,
                     f"{cl}  --  PHF1- (left, nearest PHF1- median) vs PHF1+ (right), same ROI  "
                     f"[page {start // args.rows_per_page + 1}]", ha="center", va="center", fontsize=8)
            fig.text(0.5, 1 - 0.45 / H_in, "DNA | mask | PHF1 | merge (PHF1 + mask outline);  "
                     f"bar {args.scalebar_um:g} um;  label = ROI:object  nucleus area",
                     ha="center", va="center", fontsize=6, color="0.3")
            x_left = (W_in - (8 * tw + 3 * 0.04 * 2 + gap)) / 2
            for i, rr in enumerate(chunk.itertuples()):
                ytop = H_in - 0.6 - i * (tw + lab_h + 0.06)
                for side, tag in enumerate(("neg", "pos")):
                    obj = getattr(rr, f"{tag}_obj")
                    imgs, m = render(rr.roi, obj, getattr(rr, f"{tag}_crop_row0"),
                                     getattr(rr, f"{tag}_crop_col0"))
                    xs = x_left + side * (4 * tw + 3 * 0.04 + gap)
                    edge = "" if getattr(rr, f"{tag}_in_bounds") else "  [edge-padded]"
                    fig.text(xs / W_in, (ytop - lab_h * 0.45) / H_in,
                             f"{'PHF1+' if tag == 'pos' else 'PHF1-'}  {rr.roi}:{obj}   "
                             f"{getattr(rr, f'{tag}_area_um2'):.0f} um^2{edge}"
                             + (f"      {rr.pair_id}" if tag == "neg" else ""),
                             fontsize=6, va="center",
                             fontweight="bold" if tag == "pos" else "normal")
                    for j, ch in enumerate(CHS):
                        ax = fig.add_axes([(xs + j * (tw + 0.04)) / W_in, (ytop - lab_h - tw) / H_in,
                                           tw / W_in, tw / H_in])
                        ax.imshow(imgs[ch], interpolation="nearest")
                        if ch == "merge":
                            vector_outline(ax, m)
                        vector_scalebar(ax)
                        ax.set_xlim(-0.5, N - 0.5); ax.set_ylim(N - 0.5, -0.5)
                        ax.set_xticks([]); ax.set_yticks([])
                        for s in ax.spines.values():
                            s.set_visible(False)
            pdf.savefig(fig)
            plt.close(fig)

pairs.to_csv(f"{OUT}/source_data_pairs.tsv", sep="\t", index=False)
pd.DataFrame(ranges).to_csv(f"{OUT}/stats_pairs_display_ranges.tsv", sep="\t", index=False)
print(f"{len(pairs)} pairs -> {pdf_path}")
print(pairs.groupby(["cluster", "roi"]).size().to_string())
print("pairs with PHF1+ larger:", (pairs.pos_area_um2 > pairs.neg_area_um2).sum(), "/", len(pairs))
