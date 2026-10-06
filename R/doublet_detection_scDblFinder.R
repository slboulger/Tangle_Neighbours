# ---------------------------------------------------------------------------
# doublet_detection_scDblFinder.R
#
# Upstream pipeline - step 4 of 9
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
#   Doublet detection with scDblFinder, run per donor (samples = 'sample_id') at an
#   expected doublet rate of 0.03. Adds scDblFinder.score / .class and the derived
#   predicted_doublet flag.
#
#   NOTE: doublets are called with scDblFinder throughout (the per-sample table is
#   written as QC/scrublet_by_sample.tsv).
#
# INPUTS   seu_after_qc.RDS
# OUTPUTS  seu_after_doublet.RDS
# ---------------------------------------------------------------------------
library(scDblFinder)
library(Seurat)

setwd("<PROJECT_ROOT>/phf1_v2")

seu <- readRDS("seu_after_qc.RDS")


# Convert to SingleCellExperiment (in-memory, no files)
sce <- Seurat::as.SingleCellExperiment(seu)

# Run scDblFinder (tunes itself; works on counts)
sce <- scDblFinder(sce, samples = "sample_id", dbr = 0.03)

# Bring results back to Seurat
seu$scDblFinder.score <- SummarizedExperiment::colData(sce)$scDblFinder.score
seu$scDblFinder.class <- SummarizedExperiment::colData(sce)$scDblFinder.class
seu$predicted_doublet <- seu$scDblFinder.class != "singlet"

# Quick check and filtering
write.table(table(seu$predicted_doublet, seu$sample_id), file = "QC/scrublet_by_sample.tsv")

# Optional: remove predicted doublets
seu_doublet <- subset(seu, subset = !predicted_doublet)



saveRDS(seu_doublet,"seu_after_doublet.RDS")