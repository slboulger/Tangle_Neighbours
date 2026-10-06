#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# export_seu_coords_for_phf1.R
#
# Upstream pipeline - writes PHF1/seu_coords.csv, the per-cell coordinate table read
# by the PHF1 image scripts (1E, 1F, 2C, 2D), the S1B/S1C sensitivity arms and the
# Table S3 builder
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# export_seu_coords_for_phf1.R
#
# Step 1 of the PHF1 immunofluorescence intensity extraction pipeline.
#
# Exports a slim per-cell coordinate table from the canonical seu_PHF1.rds so the
# Python extractor (python/extract_phf1_intensity.py) never has to load
# the ~8 GB Seurat object.
#
# seu_PHF1 carries BOTH FOV-local pixel coordinates (x_FOV_px / y_FOV_px) and slide
# millimetre coordinates (x_slide_mm / y_slide_mm) per cell, so the extractor can fit
# the mask-pixel -> slide-mm affine directly per FOV from real cell correspondences;
# no global-pixel calibration or flatFiles metadata are needed here.
#
# Output:
#   PHF1/seu_coords.csv  cell_id, sample_id, fov, slide_ID_numeric,
#                        x_FOV_px, y_FOV_px, x_slide_mm, y_slide_mm, PHF1
#   (cell_id encodes c_<slide_prefix>_<fov>_<cell_ID>; the extractor parses it.)
#
# Read-only on the object.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

setwd("<PROJECT_ROOT>/phf1_v2")

seu <- readRDS("seu_PHF1.rds")
md  <- seu@meta.data

pick <- function(md, candidates, what, required = TRUE) {
  hit <- candidates[candidates %in% colnames(md)]
  if (length(hit) == 0) {
    msg <- sprintf("Could not find a column for '%s'. Tried: %s.\nAvailable columns:\n%s",
                   what, paste(candidates, collapse = ", "),
                   paste(colnames(md), collapse = ", "))
    if (required) stop(msg) else { warning(msg); return(NA_character_) }
  }
  hit[1]
}

col_xmm  <- pick(md, c("x_slide_mm", "sdimx"), "x_slide_mm")
col_ymm  <- pick(md, c("y_slide_mm", "sdimy"), "y_slide_mm")
col_fov  <- pick(md, c("fov", "FOV"), "fov")
col_samp <- pick(md, c("sample_id", "sampleID"), "sample_id")
col_phf1 <- pick(md, c("PHF1", "phf1"), "PHF1 boolean")
col_slid <- pick(md, c("slide_ID_numeric", "SlideLabel", "slide_ID"), "slide id", required = FALSE)
col_xfov <- pick(md, c("x_FOV_px", "CenterX_local_px"), "x_FOV_px", required = FALSE)
col_yfov <- pick(md, c("y_FOV_px", "CenterY_local_px"), "y_FOV_px", required = FALSE)
col_ct   <- pick(md, c("celltype"), "celltype", required = FALSE)   # for the PHF1/celltype overlay figs

df <- data.frame(
  cell_id    = rownames(md),
  sample_id  = as.character(md[[col_samp]]),
  fov        = md[[col_fov]],
  x_slide_mm = as.numeric(md[[col_xmm]]),
  y_slide_mm = as.numeric(md[[col_ymm]]),
  PHF1       = as.logical(md[[col_phf1]] %in% c(TRUE, "TRUE", "True", 1, "1")),
  stringsAsFactors = FALSE
)
if (!is.na(col_slid)) df$slide_ID_numeric <- md[[col_slid]]
if (!is.na(col_xfov)) df$x_FOV_px <- as.numeric(md[[col_xfov]])
if (!is.na(col_yfov)) df$y_FOV_px <- as.numeric(md[[col_yfov]])
if (!is.na(col_ct))   df$celltype <- as.character(md[[col_ct]])

dir.create("PHF1", showWarnings = FALSE)
write.csv(df, "PHF1/seu_coords.csv", row.names = FALSE)

# celltype palette for the PHF1/celltype overlay figures (best-effort; needs palettes.R)
tryCatch({
  source("R/palettes.R")
  pal <- data.frame(celltype = names(celltype_palette), hex = unname(celltype_palette))
  write.csv(pal, "PHF1/celltype_palette.csv", row.names = FALSE)
  cat("Wrote PHF1/celltype_palette.csv\n")
}, error = function(e) warning("Could not write celltype_palette.csv: ", conditionMessage(e)))

cat(sprintf("Wrote PHF1/seu_coords.csv: %d cells, %d samples, PHF1+ = %d\n",
            nrow(df), length(unique(df$sample_id)), sum(df$PHF1)))
print(table(df$sample_id))
cat("Example cell_id:", paste(head(df$cell_id, 3), collapse = " | "), "\n")
