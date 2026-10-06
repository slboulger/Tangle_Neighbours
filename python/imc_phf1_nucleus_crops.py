#!/usr/bin/env python
# ---------------------------------------------------------------------------
# imc_phf1_nucleus_crops.py
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
# imc_phf1_nucleus_crops.py
#
# Representative IMC crops of single excitatory-neuron nuclei, PHF1+ vs PHF1-, matched for
# celltype, to illustrate the nucleus-area comparison. Image panel only -- no inference is made
# here.
#
# INPUTS (EC_IMC_Project, read-only)
#   img/<roi>.tiff              35-channel float32 steinbock stack, 1 um/px. This is the SAME file
#                               Supp1 was made from (Supp1/BBN110.26073_1.tiff is byte-identical in size).
#   DNA_Otsu_Masks/<roi>.tiff   nucleus label mask. Verified on every run: per-object pixel counts
#                               equal spe `area` exactly (max |diff| 0). SIMPLI_Masks/ does NOT
#                               match (off by up to 227 px) and is not used.
#   regionprops/<roi>.csv       centroids (centroid-0 = row, centroid-1 = col). spe spatialCoords
#                               are NOT in this pixel frame, so they are not used for cropping.
#   spe.rds (via Rscript)       celltype_clusters + PHF1_Otsu per ObjectNumber; cached as
#                               cell_labels.tsv in the output dir.
#
# DISPLAY PROCESSING = Supp1/Save_imgs.ijm, step 1, re-implemented:
#   per channel, on the WHOLE ROI: Despeckle (ImageJ 3x3 median) -> Enhance Contrast
#   saturated=0.35 (ImageJ ContrastEnhancer: 256-bin histogram, 0.175% clipped per tail) ->
#   linear map to 8-bit, Grays LUT. Ranges are computed per ROI (as Supp1), never per crop, so
#   every crop from one ROI -- PHF1+ and PHF1- alike -- shares the same display range. Crops
#   from different ROIs do not; that is stated in the stats log.
#   The nucleus outline is the measured DNA Otsu mask boundary, so the size comparison does not
#   depend on display contrast.
#
# SELECTION (deterministic, not hand-picked): within each ROI x cluster x PHF1 group, drop
#   cells whose crop window leaves the image, then take the k cells whose nucleus area is
#   closest to that group's median area. Override with --pick ROI:OBJ,... after browsing the
#   contact sheets (all candidates, ranked).
#
# NOTE BBN_9389_2 is Braak 0-1 and has NO PHF1_Otsu-positive cells, so every PHF1+ example comes
#   from BBN110.26073_1 (Braak 3-4). The headline pairing is within BBN110.26073_1; BBN_9389_2 supplies a
#   no-tangle PHF1- reference only.
#
# OUTPUTS plots/imc_phf1_nucleus_crops/
#   plot_imc_phf1_nucleus_crops_<cluster>.pdf (+ .png, 600 dpi)
#   contact_sheet_<cluster>.pdf          every candidate, ranked by |area - group median|
#   source_data_imc_phf1_nucleus_crops.tsv   one row per displayed cell
#   source_data_imc_phf1_nucleus_crops_candidates.tsv  all candidate cells
#   stats_imc_phf1_nucleus_crops.txt     selection rule, display ranges, descriptive effect sizes
#
# Run locally on the mount (base anaconda has numpy/scipy/skimage/tifffile/matplotlib):
#   python python/imc_phf1_nucleus_crops.py

import argparse, os, subprocess, sys, datetime, platform
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage, stats
from skimage import measure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

PROJ = "<PROJECT_ROOT>/phf1_v2"
EC = "<USER_DATA_ROOT>/EC_IMC_Project/"
SPE = "<IMC_ROOT>/spe.rds"
OUT = PROJ + "/plots/imc_phf1_nucleus_crops"

# 0-based slice indices in img/*.tiff (rownames(spe) order; Supp1 marker list)
CH = {"PHF1": 4, "CALB1": 10, "Reelin": 30, "DNA": 34}
CLUSTERS = {
    "Exc2_CALB1": ("Excitatory neuron cluster 2 (CALB1)", "CALB1"),
    "Exc4_RELN":  ("Excitatory neuron cluster 4 (RELN)", "Reelin"),
}
SATURATED = 0.35
PX_UM = 1.0  # IMC ablation 1 um / px

ap = argparse.ArgumentParser()
ap.add_argument("--rois", default="BBN_9389_2,BBN110.26073_1")
ap.add_argument("--k", type=int, default=3, help="cells per ROI x cluster x PHF1 group")
ap.add_argument("--half", type=int, default=16, help="crop half-width in px (= um)")
ap.add_argument("--scalebar_um", type=float, default=10)
ap.add_argument("--pick", default="", help="override: ROI:OBJ,ROI:OBJ,... (replaces auto picks for the groups those cells belong to)")
args = ap.parse_args()
ROIS = args.rois.split(",")
os.makedirs(OUT, exist_ok=True)


# ---------------------------------------------------------------- labels from spe
lab_path = f"{OUT}/cell_labels.tsv"
if not os.path.exists(lab_path):
    rcode = f'''
    suppressPackageStartupMessages(library(SingleCellExperiment))
    spe <- readRDS("{SPE}")
    k <- spe$sample_id %in% c({",".join(f'"{r}"' for r in ROIS)})
    cd <- as.data.frame(colData(spe)[k, c("sample_id","patient_id","BraakGroup","ObjectNumber",
                                          "area","eccentricity","celltype_clusters","PHF1_Otsu")])
    cd$cell_id <- colnames(spe)[k]
    write.table(cd, "{lab_path}", sep = "\\t", quote = FALSE, row.names = FALSE)
    '''
    subprocess.run(["Rscript", "-e", rcode], check=True)
lab = pd.read_csv(lab_path, sep="\t", dtype={"sample_id": str, "patient_id": str})
missing = set(ROIS) - set(lab.sample_id)
if missing:
    sys.exit(f"cell_labels.tsv lacks {missing}; delete it to regenerate")


# ---------------------------------------------------------------- Supp1 display processing
def imagej_enhance_contrast(a, saturated=SATURATED):
    """ImageJ ContrastEnhancer.getMinAndMax for a 32-bit image (256-bin histogram)."""
    a = a.ravel()
    lo, hi = float(a.min()), float(a.max())
    if hi <= lo:
        return lo, lo + 1
    nb = 256
    bin_size = (hi - lo) / nb
    hist, _ = np.histogram(a, bins=nb, range=(lo, hi))
    thr = int(a.size * saturated / 200.0)
    c = np.cumsum(hist)
    hmin = int(np.argmax(c > thr))
    cr = np.cumsum(hist[::-1])
    hmax = nb - 1 - int(np.argmax(cr > thr))
    mn, mx = lo + hmin * bin_size, lo + hmax * bin_size
    if mx <= mn:
        mx = mn + 1
    return mn, mx


def to8(a, mn, mx):
    return np.clip((a - mn) * 256.0 / (mx - mn), 0, 255).astype(np.uint8)


rois = {}
ranges = []
for r in ROIS:
    img = tifffile.imread(f"{EC}img/{r}.tiff")
    assert img.shape[0] == 35, img.shape
    mask = np.squeeze(tifffile.imread(f"{EC}DNA_Otsu_Masks/{r}.tiff")).astype(np.int64)
    assert mask.shape == img.shape[1:]
    disp = {}
    for nm, i in CH.items():
        d = ndimage.median_filter(img[i], size=3, mode="nearest")  # ImageJ Despeckle
        mn, mx = imagej_enhance_contrast(d)
        disp[nm] = to8(d, mn, mx)
        ranges.append(dict(roi=r, channel=nm, slice=i + 1, display_min=mn, display_max=mx))
    rp = pd.read_csv(f"{EC}regionprops/{r}.csv").set_index("Object")
    # gate: mask object areas must equal spe area (DNA Otsu masks are what was measured)
    ids, cnt = np.unique(mask[mask > 0], return_counts=True)
    area_mask = pd.Series(cnt, index=ids)
    lr = lab[lab.sample_id == r]
    d = (area_mask.reindex(lr.ObjectNumber).values - lr.area.values)
    if np.isnan(d).any() or np.abs(d).max() != 0:
        sys.exit(f"{r}: DNA_Otsu_Masks areas do not match spe area -- wrong mask set")
    rois[r] = dict(disp=disp, mask=mask, rp=rp, shape=mask.shape)
ranges = pd.DataFrame(ranges)


# ---------------------------------------------------------------- candidates + selection
cand = []
for key, (cl, marker) in CLUSTERS.items():
    for r in ROIS:
        H, W = rois[r]["shape"]
        rp = rois[r]["rp"]
        sub = lab[(lab.sample_id == r) & (lab.celltype_clusters == cl) & lab.PHF1_Otsu.notna()]
        for _, c in sub.iterrows():
            cy, cx = rp.loc[c.ObjectNumber, "centroid-0"], rp.loc[c.ObjectNumber, "centroid-1"]
            y0, x0 = int(round(cy)) - args.half, int(round(cx)) - args.half
            inb = y0 >= 0 and x0 >= 0 and y0 + 2 * args.half <= H and x0 + 2 * args.half <= W
            cand.append(dict(cluster_key=key, celltype_clusters=cl, marker=marker, roi=r,
                             patient_id=c.patient_id, BraakGroup=c.BraakGroup,
                             ObjectNumber=int(c.ObjectNumber), cell_id=c.cell_id,
                             PHF1=("PHF1+" if c.PHF1_Otsu == "PHF1_pos" else "PHF1-"),
                             area_um2=c.area * PX_UM ** 2, eccentricity=c.eccentricity,
                             centroid_row=cy, centroid_col=cx, crop_row0=y0, crop_col0=x0,
                             crop_in_bounds=inb))
cand = pd.DataFrame(cand)
cand["group"] = cand.roi + " " + cand.PHF1
g = cand.groupby(["cluster_key", "roi", "PHF1"])
cand["group_median_area"] = g.area_um2.transform("median")
cand["group_n"] = g.area_um2.transform("size")
cand["area_pctile_in_group"] = g.area_um2.rank(pct=True)
cand["dev_from_median"] = (cand.area_um2 - cand.group_median_area).abs()
cand = cand.sort_values(["cluster_key", "roi", "PHF1", "crop_in_bounds", "dev_from_median"],
                        ascending=[True, True, True, False, True]).reset_index(drop=True)
cand["rank"] = cand.groupby(["cluster_key", "roi", "PHF1"]).cumcount() + 1

ok = cand[cand.crop_in_bounds]
sel = ok.groupby(["cluster_key", "roi", "PHF1"]).head(args.k).copy()
sel["selection"] = "auto_median"
if args.pick:
    picks = [(p.split(":")[0], int(p.split(":")[1])) for p in args.pick.split(",")]
    pk = pd.concat([ok[(ok.roi == r) & (ok.ObjectNumber == o)] for r, o in picks])
    if len(pk) != len(picks):
        sys.exit("some --pick cells are not in-bounds candidates")
    for (ck, r, ph), _ in pk.groupby(["cluster_key", "roi", "PHF1"]):
        sel = sel[~((sel.cluster_key == ck) & (sel.roi == r) & (sel.PHF1 == ph))]
    pk = pk.assign(selection="manual_pick")
    sel = pd.concat([sel, pk])


# ---------------------------------------------------------------- drawing
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 6})
OUTLINE_GREY = "#56B4E9"   # Okabe-Ito sky blue on the grey DNA tile
UP = 10                    # outline drawn on a 10x nearest-upsampled mask -> follows pixel edges


def crop(row):
    R = rois[row.roi]
    y0, x0, s = row.crop_row0, row.crop_col0, 2 * args.half
    t = {k: v[y0:y0 + s, x0:x0 + s] for k, v in R["disp"].items()}
    t["mask"] = (R["mask"][y0:y0 + s, x0:x0 + s] == row.ObjectNumber)
    return t


def outline(ax, m, colour, lw=0.6):
    big = np.kron(np.pad(m, 1).astype(float), np.ones((UP, UP)))
    for c in measure.find_contours(big, 0.5):
        ax.plot(c[:, 1] / UP - 1 - 0.5, c[:, 0] / UP - 1 - 0.5, color=colour, lw=lw)


def merge_rgb(t, marker):
    p = t["PHF1"] / 255.0
    m = t[marker] / 255.0
    rgb = np.stack([p, m, p], axis=-1)  # PHF1 magenta, marker green (CB-safe pair)
    return np.clip(rgb, 0, 1)


def scalebar(ax, label):
    s = 2 * args.half
    L = args.scalebar_um / PX_UM
    ax.add_patch(Rectangle((s - 1.5 - L, s - 2.5), L, 1.0, color="white", lw=0))
    if label:
        ax.text(s - 1.5 - L / 2, s - 4.2, f"{args.scalebar_um:g} $\\mu$m", color="white",
                ha="center", va="bottom", fontsize=5)


def tile(ax, img, cmap=None):
    ax.imshow(img, cmap=cmap, vmin=0, vmax=255 if cmap else None, interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


group_order = [(r, "PHF1-") for r in ROIS] + [(r, "PHF1+") for r in ROIS]
for key, (cl, marker) in CLUSTERS.items():
    cols = []
    for r, ph in group_order:
        s = sel[(sel.cluster_key == key) & (sel.roi == r) & (sel.PHF1 == ph)]
        s = s.sort_values("area_um2")
        cols += [row for _, row in s.iterrows()]
    if not cols:
        continue
    rows = ["DNA", "PHF1", marker, "merge"]
    tw = 0.55
    fig, axes = plt.subplots(len(rows), len(cols), figsize=(tw * len(cols) + 0.45, tw * len(rows) + 0.5),
                             squeeze=False, gridspec_kw=dict(wspace=0.04, hspace=0.04))
    fig.subplots_adjust(left=0.45 / (tw * len(cols) + 0.45), right=0.995,
                        top=1 - 0.38 / (tw * len(rows) + 0.5), bottom=0.12 / (tw * len(rows) + 0.5))
    prev = None
    for j, row in enumerate(cols):
        t = crop(row)
        for i, ch in enumerate(rows):
            ax = axes[i, j]
            if ch == "merge":
                tile(ax, merge_rgb(t, marker))
                outline(ax, t["mask"], "white", 0.5)
                scalebar(ax, label=(j == 0))
            else:
                tile(ax, t[ch], cmap="gray")
                if ch == "DNA":
                    outline(ax, t["mask"], OUTLINE_GREY)
            if j == 0:
                ax.set_ylabel(ch, fontsize=6, rotation=0, ha="right", va="center", labelpad=3)
        axes[0, j].set_title(f"{row.area_um2:.0f} $\\mu$m$^2$", fontsize=5, pad=1.5)
        grp = (row.roi, row.PHF1)
        if grp != prev:
            n = sum(1 for c in cols if (c.roi, c.PHF1) == grp)
            x0 = axes[0, j].get_position().x0
            x1 = axes[0, j + n - 1].get_position().x1
            yt = axes[0, j].get_position().y1 + 0.13 / (tw * len(rows) + 0.5)
            braak = row.BraakGroup.replace("Braak_", "Braak ").replace("_", "-")
            fig.text((x0 + x1) / 2, yt, f"{row.PHF1}  {row.roi} ({braak})", ha="center",
                     va="bottom", fontsize=6, fontweight="bold" if row.PHF1 == "PHF1+" else "normal")
            fig.lines.append(plt.Line2D([x0, x1], [yt - 0.01, yt - 0.01], transform=fig.transFigure,
                                        color="black", lw=0.5))
            prev = grp
    base = f"{OUT}/plot_imc_phf1_nucleus_crops_{key}"
    fig.savefig(base + ".pdf")
    fig.savefig(base + ".png", dpi=600)
    plt.close(fig)

    # contact sheet: every candidate (merge + outline), ranked by closeness to group median
    cc = cand[cand.cluster_key == key]
    groups = [(r, ph) for r, ph in group_order if ((cc.roi == r) & (cc.PHF1 == ph)).any()]
    ncol = 12
    nrow_g = [int(np.ceil(((cc.roi == r) & (cc.PHF1 == ph)).sum() / ncol)) for r, ph in groups]
    fig, axes = plt.subplots(sum(nrow_g) + len(groups), ncol,
                             figsize=(ncol * 0.6, (sum(nrow_g) + len(groups)) * 0.62), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    rr = 0
    for (r, ph), nr in zip(groups, nrow_g):
        sub = cc[(cc.roi == r) & (cc.PHF1 == ph)]
        axes[rr, 0].text(0, 0.3, f"{cl} -- {ph} {r}  (n = {len(sub)}, median area "
                         f"{sub.area_um2.median():.0f} um^2; ranked by |area - median|)",
                         fontsize=7, transform=axes[rr, 0].transAxes)
        rr += 1
        for n, (_, row) in enumerate(sub.iterrows()):
            ax = axes[rr + n // ncol, n % ncol]
            ax.axis("on")
            if row.crop_in_bounds:
                t = crop(row)
                tile(ax, merge_rgb(t, marker))
                outline(ax, t["mask"], "white", 0.4)
            else:
                tile(ax, np.zeros((2 * args.half, 2 * args.half, 3)))
            ax.set_title(f"{r}:{row.ObjectNumber}  {row.area_um2:.0f}", fontsize=4.5, pad=1)
        rr += nr
    fig.tight_layout(pad=0.3)
    fig.savefig(f"{OUT}/contact_sheet_{key}.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- source data + stats
keep = ["cluster_key", "celltype_clusters", "marker", "roi", "patient_id", "BraakGroup",
        "ObjectNumber", "cell_id", "PHF1", "area_um2", "eccentricity", "group_n",
        "group_median_area", "area_pctile_in_group", "rank", "centroid_row", "centroid_col",
        "crop_row0", "crop_col0"]
sel[keep + ["selection"]].to_csv(f"{OUT}/source_data_imc_phf1_nucleus_crops.tsv", sep="\t", index=False)
cand[keep + ["crop_in_bounds"]].to_csv(f"{OUT}/source_data_imc_phf1_nucleus_crops_candidates.tsv",
                                       sep="\t", index=False)
ranges.to_csv(f"{OUT}/stats_imc_phf1_nucleus_crops_display_ranges.tsv", sep="\t", index=False)

with open(f"{OUT}/stats_imc_phf1_nucleus_crops.txt", "w") as f:
    P = lambda *a: print(*a, file=f)
    P(f"imc_phf1_nucleus_crops.py  run {datetime.datetime.now():%Y-%m-%d %H:%M}")
    P(f"ROIs: {ROIS}; crop {2*args.half} x {2*args.half} px = um; k = {args.k} per group; "
      f"scale bar {args.scalebar_um:g} um")
    P("\nIMAGE PANEL ONLY -- not an inferential test.")
    P("\n=== Inputs ===")
    P(f"image   {EC}img/<roi>.tiff  (same stack as Supp1)")
    P(f"mask    {EC}DNA_Otsu_Masks/<roi>.tiff  -- per-object pixel count == spe area for every cell (gate passed)")
    P(f"labels  {SPE}  (celltype_clusters, PHF1_Otsu)")
    P("\n=== Display processing (Supp1/Save_imgs.ijm step 1) ===")
    P(f"per channel, whole ROI: 3x3 median (Despeckle) -> ImageJ Enhance Contrast saturated={SATURATED}")
    P("(256-bin histogram, 0.175% clipped per tail) -> linear 8-bit. Range per ROI, shared by all crops of that ROI.")
    P("Crops from DIFFERENT ROIs have different display ranges (as in Supp1). Outline = measured mask boundary.")
    P("Merge: PHF1 magenta, celltype marker green, nucleus outline white.")
    P(ranges.to_string(index=False))
    P("\n=== Selection ===")
    P("Per ROI x cluster x PHF1 group: crop window inside image, then the k cells with area closest to the")
    P("group median (deterministic). selection == manual_pick marks --pick overrides.")
    if args.pick:
        P(f"--pick {args.pick}")
    P("\n=== Group sizes and nucleus area (um^2) ===")
    summ = cand.groupby(["cluster_key", "roi", "BraakGroup", "PHF1"]).area_um2.agg(
        n="size", mean="mean", median="median", q25=lambda x: x.quantile(.25),
        q75=lambda x: x.quantile(.75)).reset_index()
    P(summ.to_string(index=False))
    P("\nNOTE BBN_9389_2 (Braak 0-1) has no PHF1_Otsu-positive cells; all PHF1+ examples are BBN110.26073_1.")
    P("\n=== Effect sizes: PHF1+ vs PHF1- within BBN110.26073_1 (single ROI, descriptive) ===")
    es = []
    for key in CLUSTERS:
        for r in ROIS:
            a = cand[(cand.cluster_key == key) & (cand.roi == r)]
            pos, neg = a[a.PHF1 == "PHF1+"].area_um2.values, a[a.PHF1 == "PHF1-"].area_um2.values
            if len(pos) < 2 or len(neg) < 2:
                continue
            diff = pos.mean() - neg.mean()
            se = np.sqrt(pos.var(ddof=1) / len(pos) + neg.var(ddof=1) / len(neg))
            dfw = se ** 4 / ((pos.var(ddof=1) / len(pos)) ** 2 / (len(pos) - 1) +
                             (neg.var(ddof=1) / len(neg)) ** 2 / (len(neg) - 1))
            tq = stats.t.ppf(0.975, dfw)
            sp = np.sqrt(((len(pos) - 1) * pos.var(ddof=1) + (len(neg) - 1) * neg.var(ddof=1)) /
                         (len(pos) + len(neg) - 2))
            es.append(dict(cluster_key=key, roi=r, n_pos=len(pos), n_neg=len(neg),
                           mean_pos=pos.mean(), mean_neg=neg.mean(), diff_um2=diff,
                           ci_lo=diff - tq * se, ci_hi=diff + tq * se,
                           pct_of_neg=100 * diff / neg.mean(), cohens_d=diff / sp,
                           median_pos=np.median(pos), median_neg=np.median(neg)))
    es = pd.DataFrame(es)
    P("Welch mean difference with 95% CI; Cohen's d on pooled SD. n PHF1+ is 6-7 cells per cluster.")
    P(es.round(3).to_string(index=False))
    es.to_csv(f"{OUT}/stats_imc_phf1_nucleus_crops_effectsize.tsv", sep="\t", index=False)
    P("\n=== Displayed cells ===")
    P(sel[["cluster_key", "roi", "PHF1", "ObjectNumber", "area_um2", "area_pctile_in_group",
           "group_median_area", "selection"]].round(2).to_string(index=False))
    P("\n=== Session ===")
    P(f"python {platform.python_version()}  numpy {np.__version__}  pandas {pd.__version__}  "
      f"tifffile {tifffile.__version__}  matplotlib {matplotlib.__version__}")
    import scipy, skimage
    P(f"scipy {scipy.__version__}  skimage {skimage.__version__}  platform {platform.platform()}")

print("wrote", OUT)
