# ---------------------------------------------------------------------------
# add_PHF1.R
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
#   Sets the `PHF1` boolean on the integrated Seurat object from the manual tangle
#   annotation (PHF1/PHF1_Score.xlsx, one `cell_id` column), then re-saves the Seurat
#   object and its SingleCellExperiment conversions.
#
#   This one assignment IS the analysis definition. `dist_to_phf1_um` is a
#   nearest-neighbour distance to the set of PHF1+ neurons, so editing the annotation
#   changes the distance for EVERY cell in the dataset, not only cells near an added
#   tangle, and every distance model downstream moves with it.
#
#   Samples lacking either PHF1+ or PHF1- cells are dropped. In the published run this
#   removes exactly one section, BBN_2932 (38,129 cells, 0 PHF1+), leaving the 9 analysed
#   donors. BBN_2932 was excluded from the study because the section was white matter,
#   not entorhinal cortex; it also has no PHF1+ cells, so this filter is where the
#   exclusion takes effect in the code.
#
# INPUTS   seu_reintegrated.rds, PHF1/PHF1_Score.xlsx
# OUTPUTS  seu_PHF1.rds, sce.qs, sce_neuron.qs
# NEXT     rename_celltypes.R (required - see README running order)
# ---------------------------------------------------------------------------
library(Seurat)
library(readxl)
library(SummarizedExperiment)
library(dplyr)

setwd("<PROJECT_ROOT>/phf1_v2")
seu <- readRDS("seu_reintegrated.rds")
PHF1 <- read_xlsx("PHF1/PHF1_Score.xlsx")

seu$PHF1 <- rownames(seu@meta.data) %in% PHF1$cell_id

table(seu$PHF1)

prop.table(table(seu$PHF1,seu$celltype))

table(seu$PHF1, useNA = "ifany")

# Build a sample × PHF1 table
phf1_by_sample <- seu@meta.data %>%
  group_by(sample_id, PHF1) %>%
  summarise(n = n(),.groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = PHF1,
    values_from = n,
    values_fill = 0,
    names_prefix = "PHF1_"
  )

phf1_by_sample

# Sample_ids with both PHF1 TRUE and FALSE. Drops BBN_2932, the excluded white-matter
# section (0 PHF1+ cells); the other 9 donors are retained.
valid_samples <- phf1_by_sample %>%
  filter(PHF1_FALSE > 0, PHF1_TRUE > 0) %>%
  pull(sample_id)

valid_samples

seu <- subset(seu, subset = sample_id %in% valid_samples)

table(seu$PHF1, seu$sample_id)

saveRDS(seu, "seu_PHF1.rds")
seu <- readRDS("seu_PHF1.rds")

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
