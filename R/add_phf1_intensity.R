#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# add_phf1_intensity.R
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
# add_phf1_intensity.R
#
# Step 4 of the PHF1 immunofluorescence intensity pipeline.
#
# Joins the per-cell PHF1/DAPI mask-intensity table produced by
# python/extract_phf1_intensity.py onto the canonical seu_PHF1.rds,
# adds per-sample-normalised columns, re-saves the object, and rebuilds sce.qs /
# sce_neuron.qs so the new colData propagates to the qs objects the DEG and
# distance scripts read. Structure mirrors R/add_PHF1.R.
#
# NOTE on units: phf1_intensity_mean/median/sd are raw 8-bit DN and are only
# comparable WITHIN a sample (confocal exposure/gain differ per acquisition).
# Use phf1_intensity_mean_z (per-sample z-score) or phf1_intensity_mean_pct
# (per-sample percentile) for any cross-sample / Braak comparison.
#
# Downstream per-celltype objects (celltype_sce_neighbours/*.qs) are regenerated
# by splitting_sce_cluster_celltype.r / label_phf1_neighbours.r and will pick up
# these columns when next rebuilt.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(SummarizedExperiment)
})

setwd("<PROJECT_ROOT>/phf1_v2")

FRAC_IN_BOUNDS_MIN <- 0.5   # cells whose mask is <50% inside the PHF1 image -> NA intensity

seu <- readRDS("seu_PHF1.rds")
intens <- read.csv("PHF1/phf1_intensity_per_cell.csv", stringsAsFactors = FALSE)

stopifnot("cell_id" %in% colnames(intens))
cat(sprintf("Loaded %d cells from extractor; seu has %d cells.\n",
            nrow(intens), ncol(seu)))

## --- columns to attach (keep only what we add to colData) ---
intensity_cols <- c("phf1_intensity_mean", "phf1_intensity_median", "phf1_intensity_sd",
                    "phf1_intensity_p90", "phf1_intensity_p95", "phf1_intensity_p99",
                    "phf1_intensity_max", "phf1_intensity_sum")
keep <- c("cell_id", intensity_cols, "phf1_n_mask_px", "phf1_frac_in_bounds",
          "phf1_frac_saturated", "phf1_frac_gt64", "phf1_frac_gt128", "dapi_intensity_mean")
missing <- setdiff(keep, colnames(intens))
if (length(missing))
  stop(sprintf("phf1_intensity_per_cell.csv is missing columns: %s.\nRe-run extract_phf1_intensity.py to regenerate it with the full stat set.",
               paste(missing, collapse = ", ")))
intens <- intens[, keep]

## --- NA-out poorly-registered cells (mask mostly outside the PHF1 image) ---
poorly <- !is.na(intens$phf1_frac_in_bounds) & intens$phf1_frac_in_bounds < FRAC_IN_BOUNDS_MIN
for (cc in intensity_cols) intens[[cc]][poorly] <- NA
cat(sprintf("Set %d cells with frac_in_bounds < %.2f to NA intensity.\n",
            sum(poorly), FRAC_IN_BOUNDS_MIN))

## --- idempotency: drop any PHF1-intensity columns from a previous run before re-joining ---
derived_cols <- c("phf1_intensity_mean_z", "phf1_intensity_mean_pct",
                  "phf1_intensity_p95_z", "phf1_intensity_p95_pct")
prev <- intersect(c(setdiff(keep, "cell_id"), derived_cols), colnames(seu@meta.data))
if (length(prev)) {
  seu@meta.data[prev] <- NULL
  cat("Dropped existing columns for idempotent re-run:", paste(prev, collapse = ", "), "\n")
}

## --- left join by rowname (== cell_id), preserving cell order (cf. add_metadata.R) ---
n_before <- ncol(seu)
order_before <- colnames(seu)

md <- seu@meta.data %>%
  rownames_to_column("._cell") %>%
  left_join(intens, by = c("._cell" = "cell_id")) %>%
  column_to_rownames("._cell")

stopifnot(identical(rownames(md), order_before))   # no reordering
seu@meta.data <- md
stopifnot(ncol(seu) == n_before)

## --- per-sample normalisation (cross-sample-comparable columns) ---
## robust z (median/MAD) so the heavy positive tail does not inflate the scale;
## provided for the two primary analysis statistics (mean and p95).
robust_z <- function(x) {
  m <- median(x, na.rm = TRUE); s <- mad(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) s <- sd(x, na.rm = TRUE)
  (x - m) / s
}
md <- seu@meta.data
md <- md %>%
  group_by(sample_id) %>%
  mutate(
    phf1_intensity_mean_z   = robust_z(phf1_intensity_mean),
    phf1_intensity_mean_pct = dplyr::percent_rank(phf1_intensity_mean),
    phf1_intensity_p95_z    = robust_z(phf1_intensity_p95),
    phf1_intensity_p95_pct  = dplyr::percent_rank(phf1_intensity_p95)
  ) %>%
  ungroup() %>%
  as.data.frame()
rownames(md) <- order_before
seu@meta.data <- md

## --- report ---
cat("\nNA summary for new columns:\n")
for (cc in c("phf1_intensity_mean", "phf1_intensity_p95", "phf1_intensity_max",
             "dapi_intensity_mean", "phf1_intensity_mean_z", "phf1_intensity_p95_z")) {
  cat(sprintf("  %-26s NA = %d / %d\n", cc, sum(is.na(seu@meta.data[[cc]])), ncol(seu)))
}
cat("\nMedian PHF1 stats by manual PHF1 status (mean / p95 / max):\n")
for (cc in c("phf1_intensity_mean", "phf1_intensity_p95", "phf1_intensity_max")) {
  cat(sprintf("  %s:\n", cc)); print(tapply(seu@meta.data[[cc]], seu@meta.data$PHF1, median, na.rm = TRUE))
}

saveRDS(seu, "seu_PHF1.rds")
cat("\nSaved seu_PHF1.rds\n")

## --- rebuild SCE objects (identical block to add_PHF1.R so columns propagate) ---
sce <- as.SingleCellExperiment(seu, assay = "RNA")
rna_features <- as.data.frame(seu[["RNA"]]@meta.features)
stopifnot(nrow(rna_features) == nrow(sce))
rowData(sce) <- cbind(
  data.frame(gene = rownames(sce), row.names = rownames(sce)),
  rna_features
)
qs::qsave(sce, "sce.qs")
cat("Saved sce.qs\n")

is_neuron <- grepl("Neuron|Exc|Inh", sce$celltype)
sce_neuron <- sce[, is_neuron]
qs::qsave(sce_neuron, "sce_neuron.qs")
cat(sprintf("Saved sce_neuron.qs (%d neuronal cells)\n", ncol(sce_neuron)))
