#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# add_phf1_threshold_calls.R
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
# add_phf1_threshold_calls.R
#
# Add two binary PHF1+ calls (per-sample p95 thresholds) to seu_PHF1.rds and the SCE
# objects the DEG scripts read, so they can be used as the DEG dependent variable:
#   phf1_p95_pos_5mad : phf1_intensity_p95 >= per-sample (median + 5*MAD)
#   phf1_p95_pos_3mad : phf1_intensity_p95 >= per-sample (median + 3*MAD)
#
# Encoded as logical TRUE/FALSE (matching the existing PHF1 column), so the category
# DEG works with --ref_class FALSE. Cells with NA p95 (mask outside the IF image) -> FALSE,
# mirroring how PHF1 treats non-positive cells.
#
# Run after add_phf1_intensity.R (needs phf1_intensity_p95). Then in run_pb_deg_category.sh
# just change  --dependent_var PHF1  ->  --dependent_var phf1_p95_pos_5mad  (keep --ref_class FALSE).

suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(dplyr); library(tibble); library(qs)
})
setwd("<PROJECT_ROOT>/phf1_v2")

COL5 <- "phf1_p95_pos_5mad"; COL3 <- "phf1_p95_pos_3mad"

seu <- readRDS("seu_PHF1.rds")
md <- seu@meta.data
stopifnot("phf1_intensity_p95" %in% colnames(md))

## per-sample threshold call (MAD, falling back to SD if MAD == 0)
pos_at_k <- function(x, k) {
  m <- median(x, na.rm = TRUE); s <- mad(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) s <- sd(x, na.rm = TRUE)
  !is.na(x) & x >= (m + k * s)
}
calls <- md %>%
  rownames_to_column("cell") %>%
  group_by(sample_id) %>%
  mutate(!!COL5 := pos_at_k(phf1_intensity_p95, 5),
         !!COL3 := pos_at_k(phf1_intensity_p95, 3)) %>%
  ungroup()

# named logical vectors keyed by cell barcode (== seu rowname == SCE colname)
b5 <- setNames(calls[[COL5]], calls$cell)
b3 <- setNames(calls[[COL3]], calls$cell)

cat("New positive calls (cohort):\n")
cat(sprintf("  %s : %d  (manual PHF1+ = %d)\n", COL5, sum(b5), sum(md$PHF1 %in% c(TRUE, "TRUE"))))
cat(sprintf("  %s : %d\n", COL3, sum(b3)))

## --- attach to seu_PHF1.rds (logical; $<- overwrites, so re-running is safe) ---
seu[[COL5]] <- unname(b5[colnames(seu)])
seu[[COL3]] <- unname(b3[colnames(seu)])
saveRDS(seu, "seu_PHF1.rds"); cat("Updated seu_PHF1.rds\n")

## --- attach to every SCE the DEG variants read, matched by colname ---
attach_to_sce <- function(path) {
  if (!file.exists(path)) return(invisible())
  sce <- qs::qread(path)
  v5 <- unname(b5[colnames(sce)]); v5[is.na(v5)] <- FALSE
  v3 <- unname(b3[colnames(sce)]); v3[is.na(v3)] <- FALSE
  sce[[COL5]] <- v5; sce[[COL3]] <- v3
  qs::qsave(sce, path)
  cat(sprintf("Updated %s : %s+=%d, %s+=%d\n", path, COL5, sum(v5), COL3, sum(v3)))
}
attach_to_sce("sce.qs")
attach_to_sce("sce_neuron.qs")
for (d in c("celltype_sce", "celltype_sce_neighbours")) {
  if (dir.exists(d))
    for (f in list.files(d, pattern = "_sce(_neighbours)?\\.qs$", full.names = TRUE)) attach_to_sce(f)
}
cat("Done. Set --dependent_var to", COL5, "or", COL3, "(--ref_class FALSE) in the DEG run script.\n")
