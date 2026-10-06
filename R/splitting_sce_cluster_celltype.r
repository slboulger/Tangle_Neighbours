# ---------------------------------------------------------------------------
# splitting_sce_cluster_celltype.r
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
#
# WHAT THIS DOES
#   Splits the full SingleCellExperiment into one object per cell type, written to
#   celltype_sce/. Every differential-expression and distance script is stratified by
#   cell type and reads these per-type objects rather than the full matrix.
#
#   Run AFTER rename_celltypes.R: the split is keyed on the `celltype` column and the
#   downstream scripts expect the long excitatory labels (Exc-IT-L2-3-CBLN2-HOPX), not
#   the short ones carried by the upstream object.
#
# INPUTS   sce.qs
# OUTPUTS  celltype_sce/<celltype>.qs (17 files)
# NEXT     run_label_phf1_neighbours.sh, which adds dist_to_phf1_um to each
# ---------------------------------------------------------------------------
library(tidyverse)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(qs)


setwd("<PROJECT_ROOT>/phf1_v2")

sce <- qread("sce.qs")

outdir <- "celltype_sce"
dir.create(outdir)


sce <- sce[!is.na(SummarizedExperiment::rowData(sce)$gene), ]
sce <- sce[!duplicated(SummarizedExperiment::rowData(sce)$gene), ]

idx <- unique(sce$celltype)

# Optional: If celltype clusters below a certain number needs to be excluded 
# idx <- names(table(sce$cluster_celltype)[(table(sce$cluster_celltype) >= 1000)])

for (celltype in idx) {
  sce_subset <- sce[, sce$celltype == celltype]
  sce_name <- file.path(outdir, paste0(celltype, "_sce.qs"))
  qs::qsave(sce_subset, file = sce_name)
}
