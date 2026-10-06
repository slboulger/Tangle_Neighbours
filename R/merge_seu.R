# ---------------------------------------------------------------------------
# merge_seu.R
#
# Upstream pipeline - step 1 of 9 (object construction)
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
#   Merges the two per-slide AtoMx Seurat exports into one object and assigns the donor
#   label to every cell from its FOV range (one contiguous FOV block per donor per slide,
#   plus one single-FOV exception on slide 1).
#
#   Donor labels here are BBN IDs. The FOV ranges are the mapping; if you re-run this on
#   a differently-ordered export the ranges will not transfer.
#
#   EXCLUDED DONOR: ten sections were run, so ten donors are labelled here. BBN_2932
#   (slide 2, FOVs 42-100) is NOT part of the analysed cohort: the section was found to
#   be white matter rather than entorhinal cortex, and the donor was excluded. Its FOVs
#   are labelled only because they are in the AtoMx export; the donor is removed in
#   add_PHF1.R. Every downstream object and result covers the remaining 9 donors.
#
# INPUTS   seuratObject_*_RNA_1b.RDS, seuratObject_*_RNA_2a.RDS (AtoMx exports)
# OUTPUTS  seu_all.RDS
# ---------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

setwd("<PROJECT_ROOT>/phf1_v2")

seu_s1 <- readRDS("seuratObject_IGFQ002102_matthews_13.1.2026_CosMx_RNA_1b.RDS")
seu_s2 <- readRDS("seuratObject_IGFQ002102_matthews_13.1.2026_CosMx_RNA_2a.RDS")

DefaultAssay(seu_s1)
DefaultAssay(seu_s2)

seu_s1$SlideLabel <- 1
seu_s2$SlideLabel <- 2

# Check for collisions
length(intersect(colnames(seu_s1), colnames(seu_s2)))  # 0 means no duplicates

# Add sample_id column based on FOV range

# Ranges for s1
ranges_s1 <- tribble(
  ~sample_id, ~start, ~end,
  "BBN_9928",        1,   38,
  "BBN_9889",       39,   89,
  "BBN_24895",           90,  133,
  "BBN_9931",      134,  200,
  "BBN00428931",        201,  253,
  "BBN_10231",        254,  321
)

# Exception for s1: FOV 322 should be BBN_9931
exceptions_s1 <- tribble(
  ~fov, ~sample_id,
  322, "BBN_9931"
)

# Ranges for s2. BBN_2932 is the excluded donor (white matter, not entorhinal cortex);
# it is removed in add_PHF1.R.
ranges_s2 <- tribble(
  ~sample_id, ~start, ~end,
  "BBN00635914",        1,   41,
  "BBN_2932",     42,  100,
  "BBN00628710",      101,  230,
  "BBN00638047",      231,  273
)

source("R/assign_sample_id_from_fov.R")

seu_s1  <- assign_sample_id_from_fov(
  seu_s1, fov_col = "fov", ranges = ranges_s1, exceptions = exceptions_s1
)

seu_s2  <- assign_sample_id_from_fov(
  seu_s2, fov_col = "fov", ranges = ranges_s2, exceptions = NULL
)

table(seu_s1$sample_id, useNA = "ifany")
table(seu_s2$sample_id, useNA = "ifany")

SaveSeuratRds(seu_s1,file = "seu_s1.RDS")
SaveSeuratRds(seu_s2,file = "seu_s2.RDS")

# 1) Make image-free copies
seu_s1_noimg <- seu_s1; seu_s1_noimg@images <- list()
seu_s2_noimg <- seu_s2; seu_s2_noimg@images <- list()

SaveSeuratRds(seu_s1_noimg,file = "seu_s1_noimg.RDS")
SaveSeuratRds(seu_s2_noimg,file = "seu_s2_noimg.RDS")

seu_all <- merge(seu_s1_noimg, seu_s2_noimg, project = "CosMxMerged")

SaveSeuratRds(seu_all,file = "seu_all.RDS")

table(seu_all$SlideLabel, seu_all$sample_id, useNA = "ifany")
