# ---------------------------------------------------------------------------
# rename_celltypes.R
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
# as.SingleCellExperiment() and rowData() below need these attached explicitly
# when the script is run in batch.
library(Seurat)
library(SummarizedExperiment)
library(dplyr)

setwd("<PROJECT_ROOT>/phf1_v2")
seu <- readRDS("seu_PHF1.rds")

exc_recode <- c(
  "Exc-IT-L3-5" = "Exc-IT-L3-5-CHGA-IL1RAPL2",
  "Exc-IT-L2-3" = "Exc-IT-L2-3-CBLN2-HOPX",
  "Exc-CT-L6"   = "Exc-CT-L6-SYNPO2-SEMA3E",
  "Exc-ET-L5"   = "Exc-ET-L5-SPON1-FGD4",
  "Exc-IT-L6"   = "Exc-IT-L6-CTXN1-ERC2"
)

seu$celltype <- ifelse(
  seu$celltype %in% names(exc_recode),
  exc_recode[as.character(seu$celltype)],
  as.character(seu$celltype)
)

# if celltype was a factor and you want to preserve that:
seu$celltype <- factor(seu$celltype)

table(seu$celltype)

saveRDS(seu, "seu_PHF1.rds")

sce <- as.SingleCellExperiment(seu, assay = "RNA")

# Pull Seurat feature meta.data
rna_features <- as.data.frame(seu[["RNA"]]@meta.features)

# Check that dimensions match
stopifnot(nrow(rna_features) == nrow(sce))

# Create a clean rowData: gene_id + Seurat feature stats
rowData(sce) <- cbind(
  data.frame(
    gene = rownames(sce),
    row.names = rownames(sce)
  ),
  rna_features
)

# Inspect
rowData(sce)[1:5, ]

qs::qsave(sce, "sce.qs")

# logical: TRUE for neuronal celltypes (Neuron / Exc / Inh in the string)
is_neuron <- grepl("Neuron|Exc|Inh", sce$celltype)

table(is_neuron, sce$celltype)  # optional: inspect how many per type

# subset SCE to neuronal cells only
sce_neuron <- sce[, is_neuron]
qs::qsave(sce_neuron, "sce_neuron.qs")
