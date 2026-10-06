#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# add_dist_to_full_sce.R
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
# add_dist_to_full_sce.R
#
# ONE-TIME: produce an all-cell SCE that carries the CANONICAL dist_to_phf1_um.
#
# dist_to_phf1_um is computed by label_phf1_neighbours.r and stored ONLY in the
# per-celltype celltype_sce_neighbours/*.qs files (which also EXCLUDE the PHF1+
# neurons). The all-cell sce.qs does not have the column. Downstream analyses need
# the PHF1+ neurons retained AND the canonical distances on the all-cell object --
# so here we map the precomputed distances onto sce.qs by cell_id (PHF1+ cells
# correctly receive NA) and save the result once.
#
# Using the precomputed values guarantees downstream distances are byte-identical
# to those behind the DEG pipeline -- no risk of drift from a different
# source-neuron list.
#
# Usage (base R args, so it runs without argparse):
#   Rscript R/add_dist_to_full_sce.R [in_sce] [neighbours_dir] [out_sce]
# Defaults: sce.qs  celltype_sce_neighbours  sce_phf1_dist.qs

suppressPackageStartupMessages({
  library(qs); library(SingleCellExperiment); library(SummarizedExperiment)
})

a <- commandArgs(trailingOnly = TRUE)
in_sce  <- if (length(a) >= 1) a[1] else "sce.qs"
nb_dir  <- if (length(a) >= 2) a[2] else "celltype_sce_neighbours"
out_sce <- if (length(a) >= 3) a[3] else "sce_phf1_dist.qs"

cat("Input SCE:        ", in_sce, "\n")
cat("Neighbours dir:   ", nb_dir, "\n")
cat("Output SCE:       ", out_sce, "\n\n")

## Build cell_id -> canonical dist map from all per-celltype neighbours files.
files <- list.files(nb_dir, pattern = "\\.qs$", full.names = TRUE)
if (length(files) == 0) stop("No .qs files in ", nb_dir)
map <- list()
for (f in files) {
  s  <- qread(f)
  id <- if ("cell_id" %in% colnames(colData(s))) s$cell_id else colnames(s)
  map[[f]] <- setNames(s$dist_to_phf1_um, id)
  rm(s); gc()
}
dist_map <- do.call(c, unname(map))                  # do NOT let unlist prefix names
dist_map <- dist_map[!duplicated(names(dist_map))]
cat("Cells in neighbours map:", length(dist_map),
    "| non-NA dist:", sum(!is.na(dist_map)), "\n")

## Map onto the all-cell SCE by cell_id.
sce <- qread(in_sce)
id  <- if ("cell_id" %in% colnames(colData(sce))) sce$cell_id else colnames(sce)
sce$dist_to_phf1_um <- dist_map[id]

phf1 <- as.character(sce$PHF1) == "TRUE"
cat("All-cell SCE cells:", ncol(sce), "\n")
cat("Matched to map:    ", sum(id %in% names(dist_map)), "\n")
cat("PHF1-negative cells with non-NA dist:", sum(!phf1 & !is.na(sce$dist_to_phf1_um)),
    "/", sum(!phf1), "\n")
cat("PHF1+ cells with NA dist (expected):  ", sum(phf1 & is.na(sce$dist_to_phf1_um)),
    "/", sum(phf1), "\n")
cat("Distance summary (um):\n"); print(summary(sce$dist_to_phf1_um))

qsave(sce, out_sce)
cat("\nWritten:", out_sce, "\n")
