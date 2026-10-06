# ---------------------------------------------------------------------------
# add_metadata.R
#
# Upstream pipeline - step 2 of 9
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
#   Joins the donor-level metadata (Braak stage, sex, age, post-mortem interval) onto the
#   merged object.
#
# INPUTS   seu_all.RDS
# OUTPUTS  seu_with_meta.RDS
# ---------------------------------------------------------------------------
library(readxl)
library(Seurat)
library(dplyr)
library(tibble)

setwd("<PROJECT_ROOT>/phf1_v2")

seu <- readRDS("seu_all.RDS")
meta <- read_xlsx("CosMx_Samples.xlsx")

head(meta)
head(seu@meta.data)

sum(seu$sample_id %in% meta$CaseID) # matching cells
setdiff(unique(seu$sample_id), meta$CaseID) # missing sample_ids

seu_meta <- seu@meta.data %>%
  rownames_to_column("cell") %>%
  left_join(meta, by = c("sample_id" = "CaseID")) %>%
  column_to_rownames("cell")

# Update Seurat metadata
seu@meta.data <- seu_meta

# Sanity checks
table(is.na(seu$BrainBank))  # FALSE means attached successfully

saveRDS(seu,"seu_with_meta.RDS")
