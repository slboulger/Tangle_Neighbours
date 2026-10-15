#!/usr/bin/env python
# ---------------------------------------------------------------------------
# imc_l23_roi_phf1_p62_ubk48.py
#
# Figure panels: 2I
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# imc_l23_roi_phf1_p62_ubk48.py
#
# Paper panel: representative 400 x 400 um layer 2/3 subregion of IMC ROI BBN110.26073_1,
# PHF1 / p62 / ubK48 single channels (+ DNA in the merge).
# Processing, per channel: despeckle (3x3 median filter); ImageJ-style auto contrast with 0.35%
# saturated pixels computed on the whole ROI (256-bin histogram) to set the display maximum;
# display minimum = 1.2 x the crop median (at least 0.05 x, at most 0.8 x, the auto maximum);
# linear 8-bit mapping onto the channel LUT.
# Orientation: the full image is rotated 180 deg so layer 1 (pia) is at the
# top. Layer comes from the cell table (`layer`, a cut on depth from pia)
# joined to regionprops centroids, which are in the img/ pixel frame
# (spe Pos_X/Pos_Y are the 180-deg-rotated frame and are not used here).
# Crop: 400 x 400 um (1 um/px), lying wholly inside the L2/3 row band, placed
# on the window holding the most PHF1_Otsu+ cells.
# LUTs: PHF1 Red, p62 Cyan, ubK48 Yellow, DNA (191Ir) Navy (merge only).
# 100 um scale bar, 6 px thick, no text, on every panel.
#
# Inputs:
#   <USER_DATA_ROOT>/EC_IMC_Project/img/<ROI>.tiff          multichannel IMC image
#   <USER_DATA_ROOT>/EC_IMC_Project/regionprops/<ROI>.csv   steinbock regionprops (centroids)
#   <IMC_ROOT>/imc_cell_table.csv   CSV export of colData(spe), with at least the columns
#                                   sample_id, ObjectNumber, layer and PHF1_Otsu
# Outputs (<PROJECT_ROOT>/phf1_v2/plots/imc_l23_roi_phf1_p62_ubk48/):
#   01_PHF1.png, 02_p62.png, 03_ubK48.png, 04_merge_DNA.png (panels),
#   00_overview_ROI_box.png (orientation check, not a panel), export_log.txt
#
#   python python/imc_l23_roi_phf1_p62_ubk48.py
import os
import numpy as np, pandas as pd, tifffile
from scipy.ndimage import median_filter
from PIL import Image, ImageDraw

PROJ = "<USER_DATA_ROOT>/EC_IMC_Project/"
CELLS = os.path.expanduser("<IMC_ROOT>/imc_cell_table.csv")
SID, SAT, BAR, W = "BBN110.26073_1", 0.35, 100, 400
OUT = "<PROJECT_ROOT>/phf1_v2/plots/imc_l23_roi_phf1_p62_ubk48/"
CHANS = {"PHF1":  (5,  (1, 0, 0)),          # channel (1-based), LUT (Red)
         "p62":   (22, (0, 1, 1)),          # Cyan
         "ubK48": (1,  (1, 1, 0)),          # Yellow
         "DNA":   (35, (0, 0, 128 / 255))}  # 191Ir, Navy (#000080); merge only
SINGLES = ["PHF1", "p62", "ubK48"]


def ij_minmax(x, sat=SAT):
    """ImageJ Enhance Contrast min/max on a float image (256-bin histogram)."""
    lo, hi = float(x.min()), float(x.max())
    h, _ = np.histogram(x, bins=256, range=(lo, hi))
    t = x.size * sat / 200.0
    hmin = np.argmax(np.cumsum(h) > t)
    hmax = 255 - np.argmax(np.cumsum(h[::-1]) > t)
    bs = (hi - lo) / 256
    return lo + hmin * bs, lo + hmax * bs


def scalebar(rgb):
    out = rgb.copy()
    h, w = out.shape[:2]
    m = (w + h) // 100
    out[h - m - 6:h - m, m:m + BAR] = 255
    return out


stack = tifffile.TiffFile(PROJ + f"img/{SID}.tiff")
H, Wimg = stack.pages[0].shape

# Layer per cell in the rotated (L1-up) pixel frame
ct = pd.read_csv(CELLS, usecols=["sample_id", "ObjectNumber", "layer", "PHF1_Otsu"])
ct = ct[ct.sample_id == SID]
rp = pd.read_csv(PROJ + f"regionprops/{SID}.csv")
c = ct.merge(rp, left_on="ObjectNumber", right_on="Object", how="inner")
assert len(c) == len(ct), "cell table / regionprops join lost cells"
c["row"] = H - 1 - c["centroid-0"]          # 180 deg rotation
c["col"] = Wimg - 1 - c["centroid-1"]
assert c.loc[c.layer == "1", "row"].median() < c.loc[c.layer == "2/3", "row"].median(), "L1 not on top"
l23 = c[c.layer == "2/3"]
top, bot = int(np.ceil(c.loc[c.layer == "1", "row"].max())), int(np.floor(c.loc[c.layer == "4", "row"].min()))
assert bot - top >= W and Wimg >= W

# Window inside L2/3 with the most PHF1+ cells (ties -> nearest the band centre)
pos = l23[l23.PHF1_Otsu == "PHF1_pos"]
best = None
for y0 in range(top, bot - W + 1, 5):
    for x0 in range(0, Wimg - W + 1, 5):
        n = int(((pos.row >= y0) & (pos.row < y0 + W) & (pos.col >= x0) & (pos.col < x0 + W)).sum())
        d = abs(y0 + W / 2 - (top + bot) / 2) + abs(x0 + W / 2 - Wimg / 2)
        if best is None or (n, -d) > (best[0], -best[1]):
            best = (n, d, y0, x0)
n_pos, _, y0, x0 = best
os.makedirs(OUT, exist_ok=True)
log = [f"Source: {PROJ}img/{SID}.tiff ({Wimg} x {H} px, 1 um/px)",
       "Rotation: 180 deg (layer 1 at top of the full image)",
       f"L2/3 band (rotated frame, between last L1 and first L4 centroid): rows {top}:{bot}",
       f"Crop (rotated frame): y={y0}:{y0 + W} x={x0}:{x0 + W}; "
       f"unrotated img/ frame: y={H - y0 - W}:{H - y0} x={Wimg - x0 - W}:{Wimg - x0}",
       f"PHF1_Otsu+ cells in crop: {n_pos} (of {len(pos)} in L2/3, {int((c.PHF1_Otsu == 'PHF1_pos').sum())} in ROI)",
       "Despeckle; auto contrast saturated=0.35% on ROI; display min = max(1.2 x crop median, 0.05 x auto max)",
       "LUTs: PHF1 Red, p62 Cyan, ubK48 Yellow, DNA Navy #000080 (merge only)",
       f"Scale bar: {BAR} um, 6 px, no text (every panel)", ""]

rgb = {}
for name, (ch, lut) in CHANS.items():
    img = np.rot90(stack.pages[ch - 1].asarray(), 2)                 # L1 to top
    img = median_filter(img, size=3, mode="nearest")                 # Despeckle
    _, mx = ij_minmax(img)
    crop = img[y0:y0 + W, x0:x0 + W]
    mn = min(max(1.2 * float(np.median(crop)), 0.05 * mx), 0.8 * mx)
    v = np.clip(np.floor((crop - mn) * 256 / (mx - mn)), 0, 255)
    rgb[name] = np.stack([v * c_ for c_ in lut], -1)
    log.append(f"  {ch:2d} {name}: used {mn:.4f} - {mx:.4f}")

for j, name in enumerate(SINGLES, 1):
    Image.fromarray(scalebar(rgb[name].astype(np.uint8))).save(OUT + f"{j:02d}_{name}.png")
merge = np.clip(sum(rgb.values()), 0, 255).astype(np.uint8)
Image.fromarray(scalebar(merge)).save(OUT + f"{len(SINGLES) + 1:02d}_merge_DNA.png")

# Orientation check: full rotated DNA image with the crop boxed (not a panel)
dna = np.rot90(stack.pages[CHANS["DNA"][0] - 1].asarray(), 2)
lo, hi = ij_minmax(median_filter(dna, size=3, mode="nearest"))
ov = Image.fromarray(np.clip((dna - lo) * 255 / (hi - lo), 0, 255).astype(np.uint8)).convert("RGB")
dr = ImageDraw.Draw(ov)
dr.rectangle([x0, y0, x0 + W - 1, y0 + W - 1], outline=(255, 255, 0), width=4)
for r_ in (top, bot):
    dr.line([0, r_, Wimg, r_], fill=(255, 0, 255), width=3)
ov.save(OUT + "00_overview_ROI_box.png")
log.append("00_overview_ROI_box.png: full rotated DNA, crop in yellow, L2/3 band limits in magenta (check only)")

open(OUT + "export_log.txt", "w").write("\n".join(log) + "\n")
print("\n".join(log))
