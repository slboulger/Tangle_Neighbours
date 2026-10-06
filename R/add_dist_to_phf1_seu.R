#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# add_dist_to_phf1_seu.R
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
# add_dist_to_phf1_seu.R
#
# Joins the per-cell distance-to-nearest-PHF1+-neuron (`dist_to_phf1_um`) onto the
# canonical seu_PHF1.rds, under the same column name it carries on the per-celltype
# SCE files.
#
# `dist_to_phf1_um` is computed by R/label_phf1_neighbours.r and stored on the
# per-celltype SingleCellExperiment files in celltype_sce_neighbours/*.qs, but it
# is NOT on the Seurat object that all plotting/loading defaults to. This script
# READS the already-computed values out of those SCEs (it does not recompute) so
# the Seurat column is byte-identical to what the SCE/DEG pipeline uses, and joins
# them onto seu_PHF1.rds. Structure mirrors R/add_phf1_intensity.R.
#
# NA for `dist_to_phf1_um` is by design: PHF1+ source neurons get NA (they are the
# source set, not query cells) in label_phf1_neighbours.r.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(qs)
  library(SummarizedExperiment)
})

setwd("<PROJECT_ROOT>/phf1_v2")

NEIGHBOUR_DIR <- "celltype_sce_neighbours"

seu <- readRDS("seu_PHF1.rds")
cat(sprintf("Loaded seu_PHF1.rds: %d cells.\n", ncol(seu)))

## --- gather dist_to_phf1_um from every per-celltype neighbour SCE ---
sce_files <- list.files(NEIGHBOUR_DIR, pattern = "_sce_neighbours\\.qs$", full.names = TRUE)
if (length(sce_files) == 0)
  stop(sprintf("No *_sce_neighbours.qs files found in %s", NEIGHBOUR_DIR))
cat(sprintf("Reading dist_to_phf1_um from %d neighbour SCE files.\n", length(sce_files)))

dist_list <- lapply(sce_files, function(f) {
  sce <- qread(f)
  if (!"dist_to_phf1_um" %in% colnames(colData(sce)))
    stop(sprintf("%s has no dist_to_phf1_um column.", basename(f)))
  d <- data.frame(
    cell_id         = colnames(sce),
    dist_to_phf1_um = sce$dist_to_phf1_um,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  rm(sce); gc()
  d
})
dist_tbl <- do.call(rbind, dist_list)
rm(dist_list); gc()

stopifnot(!any(duplicated(dist_tbl$cell_id)))   # each cell lives in exactly one celltype file
cat(sprintf("Gathered distances for %d cells.\n", nrow(dist_tbl)))

## --- idempotency: drop any dist column from a previous run before re-joining ---
if ("dist_to_phf1_um" %in% colnames(seu@meta.data)) {
  seu@meta.data[["dist_to_phf1_um"]] <- NULL
  cat("Dropped existing dist_to_phf1_um column for idempotent re-run.\n")
}

## --- left join by rowname (== cell_id), preserving cell order (cf. add_phf1_intensity.R) ---
n_before     <- ncol(seu)
order_before <- colnames(seu)

md <- seu@meta.data %>%
  rownames_to_column("._cell") %>%
  left_join(dist_tbl, by = c("._cell" = "cell_id")) %>%
  column_to_rownames("._cell")

stopifnot(identical(rownames(md), order_before))   # no reordering
seu@meta.data <- md
stopifnot(ncol(seu) == n_before)

## --- report ---
d  <- seu@meta.data$dist_to_phf1_um
n_na      <- sum(is.na(d))
n_matched <- sum(!is.na(d))
cat(sprintf("\nMatched %d / %d cells; %d NA.\n", n_matched, ncol(seu), n_na))

cat("\nis.na(dist_to_phf1_um) vs PHF1 (NA should be the PHF1+ source neurons):\n")
print(table(dist_NA = is.na(d), PHF1 = seu@meta.data$PHF1, useNA = "ifany"))

cat("\nsummary(dist_to_phf1_um), um:\n")
print(summary(d))

saveRDS(seu, "seu_PHF1.rds")
cat("\nSaved seu_PHF1.rds\n")
