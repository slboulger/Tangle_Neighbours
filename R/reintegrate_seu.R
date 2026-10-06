# ---------------------------------------------------------------------------
# reintegrate_seu.R
#
# Upstream pipeline - step 9 of 9
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
#   Merges the round-1 and round-2 annotations into a single object carrying the final
#   celltype column, which is what add_PHF1.R then reads.
#
# INPUTS   AUCell_1/seu.rds, AUCell_Neuron/seu.rds
# OUTPUTS  seu_reintegrated.rds
# ---------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

setwd("<PROJECT_ROOT>/phf1_v2")

neuron_seu <- readRDS("AUCell_Neuron/seu.rds")
seu <- readRDS("AUCell_1/seu.rds")

head(neuron_seu)
head(seu)

# 1. Start from the coarse labels
seu$celltype <- seu$AUCell_label_SCT_propAdj

# 2. Get neuron labels, aligned to seu's cells
neuron_labels <- neuron_seu$AUCell_neuron_label_SCT_propAdj
names(neuron_labels) <- colnames(neuron_seu)

# 3. Restrict to cells that are actually in `seu`
common_cells <- intersect(colnames(seu), names(neuron_labels))
neuron_labels <- neuron_labels[common_cells]

# 4. Overwrite neuron cells in `seu$celltype` with neuron-specific labels
seu$celltype[common_cells] <- neuron_labels

# Extract metadata for Napari
meta <- seu@meta.data
meta <- meta %>% select(celltype)
meta$cell_ID <- row.names(meta) # adds cell_ID column
rownames(meta) <- NULL
meta <- meta %>% relocate(cell_ID) # moves cell_ID to first column position
write.table(meta, file="seu_metadata.csv", 
            sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

saveRDS(seu,"seu_reintegrated.rds")

# Separate slide metadata

meta <- seu@meta.data %>%
  select(SlideLabel, celltype)

# Add cell_ID column from rownames
meta$cell_ID <- rownames(meta)
rownames(meta) <- NULL

# Move cell_ID first
meta <- meta %>% relocate(cell_ID)

# Split by SlideLabel
meta_slide1 <- meta %>% filter(SlideLabel == 1) %>% select(-SlideLabel)
meta_slide2 <- meta %>% filter(SlideLabel == 2) %>% select(-SlideLabel)

# Write to separate files
write.table(
  meta_slide1,
  file = "stitched_slides/1b/_metadata.csv",
  sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE
)

write.table(
  meta_slide2,
  file = "stitched_slides/2a/_metadata.csv",
  sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE
)

meta <- seu@meta.data %>%
  mutate(cell_ID = rownames(seu@meta.data))

# Subset to Slide 2
meta_slide2 <- meta %>%
  filter(SlideLabel == 2)

ggplot(meta_slide2, aes(x = x_slide_mm, y = y_slide_mm, color = celltype)) +
  geom_point(size = 0.3, alpha = 0.8) +
  scale_color_discrete(name = "Cell type") +
  coord_fixed() +                      # preserve aspect ratio
  scale_y_reverse() +                  # often needed for image-like orientation
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank()
  ) +
  ggtitle("Cell types - Slide 2")
