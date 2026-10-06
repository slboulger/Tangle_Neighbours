# ---------------------------------------------------------------------------
# imc_crop_utils.py
#
# Shared utility - sourced by the scripts above, no panel of its own
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# imc_crop_utils.py
#
# Shared helpers for the IMC single-nucleus crop scripts (imc_phf1_nucleus_crops_pairs.py,
# imc_phf1_nucleus_crops_singles.py). See imc_phf1_nucleus_crops.py for provenance notes.
#
# - Display processing = Supp1/Save_imgs.ijm step 1: per channel over the WHOLE ROI, 3x3 median
#   (Despeckle) -> ImageJ Enhance Contrast saturated=0.35 -> linear 8-bit. Never per crop.
# - Masks = EC_IMC_Project/DNA_Otsu_Masks (gated: object pixel count == spe area).
# - Centroids from regionprops/<roi>.csv (spe spatialCoords are not in this pixel frame).

import os, subprocess
import numpy as np, pandas as pd, tifffile
from scipy import ndimage

PROJ = "<PROJECT_ROOT>/phf1_v2"
EC = "<USER_DATA_ROOT>/EC_IMC_Project/"
SPE = "<IMC_ROOT>/spe.rds"
OUTROOT = PROJ + "/plots/imc_phf1_nucleus_crops"

CH = {"PHF1": 4, "CALB1": 10, "Reelin": 30, "DNA": 34}   # 0-based slices of img/*.tiff
SATURATED = 0.35
PX_UM = 1.0


def imagej_enhance_contrast(a, saturated=SATURATED):
    """ImageJ ContrastEnhancer.getMinAndMax for a 32-bit image (256-bin histogram)."""
    a = a.ravel()
    lo, hi = float(a.min()), float(a.max())
    if hi <= lo:
        return lo, lo + 1
    nb = 256
    bs = (hi - lo) / nb
    hist, _ = np.histogram(a, bins=nb, range=(lo, hi))
    thr = int(a.size * saturated / 200.0)
    hmin = int(np.argmax(np.cumsum(hist) > thr))
    hmax = nb - 1 - int(np.argmax(np.cumsum(hist[::-1]) > thr))
    mn, mx = lo + hmin * bs, lo + hmax * bs
    return (mn, mx) if mx > mn else (mn, mn + 1)


def load_labels(rois, cache_path):
    """celltype_clusters + PHF1_Otsu per ObjectNumber from spe.rds (cached TSV)."""
    if os.path.exists(cache_path):
        lab = pd.read_csv(cache_path, sep="\t", dtype={"sample_id": str, "patient_id": str})
        if set(rois) <= set(lab.sample_id):
            return lab[lab.sample_id.isin(rois)]
    rcode = f'''
    suppressPackageStartupMessages(library(SingleCellExperiment))
    spe <- readRDS("{SPE}")
    k <- spe$sample_id %in% c({",".join(f'"{r}"' for r in rois)})
    cd <- as.data.frame(colData(spe)[k, c("sample_id","patient_id","BraakGroup","ObjectNumber",
                                          "area","eccentricity","celltype_clusters","PHF1_Otsu")])
    cd$cell_id <- colnames(spe)[k]
    write.table(cd, "{cache_path}", sep = "\\t", quote = FALSE, row.names = FALSE)
    '''
    env = dict(os.environ, PATH="/usr/local/bin:" + os.environ.get("PATH", ""))
    subprocess.run(["Rscript", "-e", rcode], check=True, env=env)
    return pd.read_csv(cache_path, sep="\t", dtype={"sample_id": str, "patient_id": str})


def load_roi(roi, lab=None, channels=("DNA", "PHF1")):
    """Processed 8-bit display channels, label mask, regionprops, display ranges for one ROI."""
    img = tifffile.imread(f"{EC}img/{roi}.tiff")
    assert img.shape[0] == 35, img.shape
    mask = np.squeeze(tifffile.imread(f"{EC}DNA_Otsu_Masks/{roi}.tiff")).astype(np.int64)
    assert mask.shape == img.shape[1:]
    disp, ranges = {}, []
    for nm in channels:
        i = CH[nm]
        d = ndimage.median_filter(img[i], size=3, mode="nearest")
        mn, mx = imagej_enhance_contrast(d)
        disp[nm] = np.clip((d - mn) * 256.0 / (mx - mn), 0, 255).astype(np.uint8)
        ranges.append(dict(roi=roi, channel=nm, slice=i + 1, display_min=mn, display_max=mx))
    rp = pd.read_csv(f"{EC}regionprops/{roi}.csv").set_index("Object")
    if lab is not None:
        ids, cnt = np.unique(mask[mask > 0], return_counts=True)
        lr = lab[lab.sample_id == roi]
        d = pd.Series(cnt, index=ids).reindex(lr.ObjectNumber).values - lr.area.values
        if np.isnan(d).any() or np.abs(d).max() != 0:
            raise SystemExit(f"{roi}: DNA_Otsu_Masks areas do not match spe area -- wrong mask set")
    return dict(disp=disp, mask=mask, rp=rp, shape=mask.shape, ranges=ranges)


def crop_window(R, obj, half):
    cy, cx = R["rp"].loc[obj, "centroid-0"], R["rp"].loc[obj, "centroid-1"]
    y0, x0 = int(round(cy)) - half, int(round(cx)) - half
    H, W = R["shape"]
    inb = y0 >= 0 and x0 >= 0 and y0 + 2 * half <= H and x0 + 2 * half <= W
    return y0, x0, inb, cy, cx


def tiles(R, obj, y0, x0, n):
    """8-bit RGB arrays for DNA (grey), mask (white fill), PHF1 (red), plus the boolean mask."""
    dna = R["disp"]["DNA"][y0:y0 + n, x0:x0 + n]
    phf = R["disp"]["PHF1"][y0:y0 + n, x0:x0 + n]
    m = R["mask"][y0:y0 + n, x0:x0 + n] == obj
    z = np.zeros_like(dna)
    return dict(DNA=np.stack([dna] * 3, -1),
                mask=np.stack([m.astype(np.uint8) * 255] * 3, -1),
                PHF1=np.stack([phf, z, z], -1)), m


def upscale_with_outline(rgb, m, up, thick):
    """Nearest-neighbour upscale; draw the mask outline in white, `thick` upscaled px wide,
    just INSIDE the pixel edges of the mask (so it never covers signal outside the nucleus)."""
    big = np.kron(rgb, np.ones((up, up, 1), dtype=np.uint8))
    if m is not None:
        bm = np.kron(m.astype(np.uint8), np.ones((up, up), dtype=np.uint8)).astype(bool)
        er = ndimage.binary_erosion(bm, iterations=thick, border_value=0)
        big[bm & ~er] = 255
    return big


def burn_scalebar(big, um, up, bar_h=None, margin_px=None):
    """White bar, bottom-left, length `um` (1 um/px), no text."""
    n = big.shape[0] // up
    L = int(round(um / PX_UM)) * up
    bh = bar_h or max(2, up // 2)
    mg = margin_px or up
    big[n * up - 2 * up - bh:n * up - 2 * up, mg:mg + L] = 255
    return big
