#!/usr/bin/env python
# ---------------------------------------------------------------------------
# imc_phf1_nucleus_crops_figure.py
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
# imc_phf1_nucleus_crops_figure.py
#
# Paper figure panel (19 cm wide): IMC nuclei of PHF1- and PHF1+ excitatory neurons, matched for
# excitatory cluster. The displayed cells were chosen by eye from the candidate pairs sheet
# (pairs/pairs_BBN110.26073_1_BBN110.26073_3_BBN110.26073_4.pdf, from imc_phf1_nucleus_crops_pairs.py);
# each PHF1- cell has a PHF1+ partner of the SAME ROI and cluster:
#   BBN110.26073_1  523 (Exc3, PHF1-)  <->  653  (Exc3, PHF1+)
#   BBN110.26073_3  1060 (Exc4, PHF1-) <->  1098 (Exc4, PHF1+)
#   BBN110.26073_3  1079 (Exc2, PHF1-) <->  990  (Exc2, PHF1+)
# Columns left to right: 523, 1060, 1079 | 653, 1098, 990.
# The shared excitatory cluster of each pair is asserted; same-ROI is reported per pair in the
# stats log rather than asserted.
# Rows: DNA (grey) + nucleus mask outline; PHF1 (red) + nucleus mask outline.
# Under each column: measured nucleus area (um^2).
# 10 um scale bar bottom-left on every tile, no text.
#
# Display = Supp1 recipe per ROI (imc_crop_utils.py): partners share one display range; tiles
# from different ROIs do not. Outline = measured DNA Otsu mask (pixel count == spe area, asserted).
#
# Outputs plots/imc_phf1_nucleus_crops/figure/: plot_imc_phf1_nucleus_figure.pdf (+ .png preview),
# source_data_imc_phf1_nucleus_figure.tsv, stats_imc_phf1_nucleus_figure.txt
#
#   python python/imc_phf1_nucleus_crops_figure.py

import os, sys, datetime, platform
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from skimage import measure

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import imc_crop_utils as U

CELLS = [("BBN110.26073_1", 523), ("BBN110.26073_3", 1060), ("BBN110.26073_3", 1079),
         ("BBN110.26073_1", 653), ("BBN110.26073_3", 1098), ("BBN110.26073_3", 990)]
HALF = 16                    # 32 x 32 um crop
SCALEBAR_UM = 10
WIDTH_CM = 19.0
OUTLINE = "#56B4E9"          # Okabe-Ito sky blue: contrasts with both grey DNA and red PHF1
OUT = U.OUTROOT + "/figure"
os.makedirs(OUT, exist_ok=True)
N = 2 * HALF

rois = sorted({r for r, _ in CELLS})
lab = U.load_labels(rois, f"{U.OUTROOT}/cell_labels_pairs.tsv")

R, ranges, rows = {}, [], []
for r in rois:
    R[r] = U.load_roi(r, lab)
    ranges += R[r]["ranges"]
    R[r]["pad"] = dict(disp={k: np.pad(v, HALF) for k, v in R[r]["disp"].items()},
                       mask=np.pad(R[r]["mask"], HALF))

tiles = []
for r, obj in CELLS:
    row = lab[(lab.sample_id == r) & (lab.ObjectNumber == obj)]
    assert len(row) == 1, (r, obj)
    row = row.iloc[0]
    assert row.celltype_clusters.startswith("Excitatory"), (r, obj, row.celltype_clusters)
    y0, x0, inb, cy, cx = U.crop_window(R[r], obj, HALF)
    P = R[r]["pad"]
    t, m = U.tiles(dict(disp=P["disp"], mask=P["mask"]), obj, y0 + HALF, x0 + HALF, N)
    assert m.sum() == row.area, "mask area != measured area"
    short = row.celltype_clusters.replace("Excitatory neuron cluster ", "Exc").split(" (")[0]
    ph = "PHF1+" if row.PHF1_Otsu == "PHF1_pos" else "PHF1-"
    tiles.append(dict(t=t, m=m, short=short, ph=ph, area=row.area))
    rows.append(dict(column=len(rows) + 1, roi=r, patient_id=row.patient_id,
                     BraakGroup=row.BraakGroup, ObjectNumber=obj, cell_id=row.cell_id,
                     celltype_clusters=row.celltype_clusters, label=short, PHF1_Otsu=row.PHF1_Otsu,
                     nucleus_area_um2=row.area, eccentricity=row.eccentricity,
                     centroid_row=cy, centroid_col=cx, crop_row0=y0, crop_col0=x0, crop_um=N,
                     edge_padded=not inb))
src = pd.DataFrame(rows)
# the matching the figure claims: column i (PHF1-) and i+3 (PHF1+) share excitatory cluster;
# same ROI is reported, not required
for i in range(3):
    a, b = src.iloc[i], src.iloc[i + 3]
    assert a.PHF1_Otsu == "PHF1_neg" and b.PHF1_Otsu == "PHF1_pos"
    assert a.celltype_clusters == b.celltype_clusters, (a.ObjectNumber, b.ObjectNumber)


# ------------------------------------------------------------------ figure
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7,
                     "axes.linewidth": 0})
W = WIDTH_CM / 2.54
left, right = 0.30, 0.02           # inches: row labels on the left
gap, group_gap = 0.045, 0.16
tw = (W - left - right - 4 * gap - group_gap) / 6
top_hdr = 0.42                     # PHF1-/PHF1+ bracket + celltype labels
row_gap = 0.045
bot = 0.20                         # nucleus-area annotation under each column
H = top_hdr + 2 * tw + row_gap + bot

fig = plt.figure(figsize=(W, H))


def xpos(j):
    return left + j * (tw + gap) + (group_gap - gap if j >= 3 else 0)


def outline(ax, m):
    big = np.kron(np.pad(m, 1).astype(float), np.ones((10, 10)))
    for c in measure.find_contours(big, 0.5):
        ax.plot(c[:, 1] / 10 - 1.5, c[:, 0] / 10 - 1.5, color=OUTLINE, lw=0.8,
                solid_joinstyle="miter", solid_capstyle="butt")


def scalebar(ax):
    L = SCALEBAR_UM / U.PX_UM
    ax.add_patch(Rectangle((1.0, N - 2.6), L, 0.75, facecolor="white", edgecolor="none"))


for j, tl in enumerate(tiles):
    for i, ch in enumerate(["DNA", "PHF1"]):
        y = H - top_hdr - (i + 1) * tw - i * row_gap
        ax = fig.add_axes([xpos(j) / W, y / H, tw / W, tw / H])
        ax.imshow(tl["t"][ch], interpolation="nearest")
        outline(ax, tl["m"])
        scalebar(ax)
        ax.set_xlim(-0.5, N - 0.5); ax.set_ylim(N - 0.5, -0.5)
        ax.set_xticks([]); ax.set_yticks([])
        if j == 0:
            fig.text((left - 0.07) / W, (y + tw / 2) / H, ch, rotation=90, ha="center",
                     va="center", fontsize=7)
    fig.text((xpos(j) + tw / 2) / W, (H - top_hdr + 0.05) / H, tl["short"], ha="center",
             va="bottom", fontsize=7)
    # measured nucleus area (DNA Otsu mask pixel count, 1 um/px) under the column
    fig.text((xpos(j) + tw / 2) / W, (bot - 0.05) / H, f"{tl['area']:.0f} \u00b5m\u00b2",
             ha="center", va="top", fontsize=7)

# group brackets: tangle-free = PHF1_Otsu negative, tangle-bearing = PHF1_Otsu positive
for g, lab_txt in ((0, "Tangle-free"), (1, "Tangle-bearing")):   # PHF1- / PHF1+ (PHF1_Otsu)
    x0, x1 = xpos(3 * g), xpos(3 * g + 2) + tw
    yb = H - top_hdr + 0.24
    fig.add_artist(plt.Line2D([x0 / W, x1 / W], [yb / H, yb / H], color="black", lw=0.6))
    fig.text(((x0 + x1) / 2) / W, (yb + 0.035) / H, lab_txt, ha="center", va="bottom", fontsize=7)

base = f"{OUT}/plot_imc_phf1_nucleus_figure"
fig.savefig(base + ".pdf")
fig.savefig(base + ".png", dpi=600)
plt.close(fig)

# ------------------------------------------------------------------ source data + stats
src.to_csv(f"{OUT}/source_data_imc_phf1_nucleus_figure.tsv", sep="\t", index=False)
rg = pd.DataFrame(ranges)
with open(f"{OUT}/stats_imc_phf1_nucleus_figure.txt", "w") as f:
    P = lambda *a: print(*a, file=f)
    P(f"imc_phf1_nucleus_crops_figure.py  run {datetime.datetime.now():%Y-%m-%d %H:%M}")
    P(f"Figure width {WIDTH_CM} cm; crops {N} x {N} um (1 um/px); scale bar {SCALEBAR_UM} um on every tile.")
    P("\nIMAGE PANEL ONLY -- representative examples, no test.")
    P("\n=== Cells (left to right) ===")
    P(src[["column", "roi", "ObjectNumber", "label", "celltype_clusters", "PHF1_Otsu",
           "nucleus_area_um2", "eccentricity", "edge_padded"]].to_string(index=False))
    P("\nMatching: column i (PHF1-) and column i+3 (PHF1+) share excitatory cluster (asserted).")
    for i in range(3):
        a, b = src.iloc[i], src.iloc[i + 3]
        P(f"  pair {i + 1}: {a.roi}:{a.ObjectNumber} <-> {b.roi}:{b.ObjectNumber}  "
          + ("same ROI, shared display range" if a.roi == b.roi else
             "DIFFERENT ROIs -> different display ranges; cluster-matched only"))
    P("Chosen by eye from the candidate pairs sheet pairs/pairs_BBN110.26073_1_BBN110.26073_3_BBN110.26073_4.pdf.")
    P("\nWithin-pair nucleus area (PHF1+ minus PHF1-):")
    for i in range(3):
        a, b = src.iloc[i], src.iloc[i + 3]
        P(f"  {a.roi}:{a.ObjectNumber} -> {b.roi}:{b.ObjectNumber} {a.label}: {a.nucleus_area_um2:.0f} -> {b.nucleus_area_um2:.0f} um^2 "
          f"({b.nucleus_area_um2 - a.nucleus_area_um2:+.0f}, {100 * (b.nucleus_area_um2 / a.nucleus_area_um2 - 1):+.0f}%)")
    P("\n=== Display processing (Supp1/Save_imgs.ijm step 1) ===")
    P("Per channel over the whole ROI: 3x3 median (Despeckle) -> ImageJ Enhance Contrast saturated=0.35 -> 8-bit.")
    P("Partners share a display range; the three ROIs do not. DNA grey, PHF1 red; outline = DNA Otsu nucleus mask.")
    P(rg[rg.channel.isin(["DNA", "PHF1"])].to_string(index=False))
    P(f"\npython {platform.python_version()}  numpy {np.__version__}  pandas {pd.__version__}  "
      f"matplotlib {matplotlib.__version__}  platform {platform.platform()}")
print("wrote", base + ".pdf")
