# ---------------------------------------------------------------------------
# SCTransform.R
#
# Upstream pipeline - step 5 of 9
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
#   SCTransform normalisation with the glmGamPoi fitting backend.
#
# INPUTS   seu_after_doublet.RDS
# OUTPUTS  seu_transformed.RDS
# ---------------------------------------------------------------------------
library(Seurat)

setwd("<PROJECT_ROOT>/phf1_v2")

seu <- readRDS("seu_after_doublet.RDS")

min_cells <- ceiling(0.005 * ncol(seu)) # genes present in >0.5% cells
keep_genes <- rowSums(GetAssayData(seu, slot = "counts") > 0) >= min_cells
table(keep_genes)
seu <- subset(seu, features = rownames(seu)[keep_genes])

# Normalize and find variable features
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4)
seu <- FindVariableFeatures(seu, nfeatures = 2000, selection.method = "vst")

vars.to.regress <- intersect(c("nCount_RNA","percent.neg","percent.false"),
                             colnames(seu@meta.data))

seu <- SCTransform(
  seu,
  assay = "RNA",
  new.assay.name = "SCT",
  method = "glmGamPoi",                 # faster/more stable
  vars.to.regress = vars.to.regress,    # remove if you prefer minimal regression
  return.only.var.genes = FALSE,        # keep full panel
  residual.features = rownames(seu), # compute residuals for all genes
  verbose = TRUE
)

saveRDS(seu,"seu_transformed.RDS")