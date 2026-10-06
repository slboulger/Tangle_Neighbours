#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# label_phf1_neighbours.r
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
# label_phf1_neighbours.r
#
# ONE-TIME script: runs across all cells, computes per-cell distance to the
# nearest PHF1+ NEURON within the same sample_id, and writes back per-celltype
# SCE .qs files with two new colData columns:
#
#   dist_to_phf1_um   : distance in um to the nearest PHF1+ neuron in same sample
#                       (NA for PHF1+ neurons themselves)
#   phf1_neighbour    : character — "near" (<= 20 um), "distal" (> 20 um), or
#                       NA for PHF1+ neurons (excluded from DEG)
#
# Inputs
#   --all_sce         : path to a single SCE containing ALL cells (all celltypes)
#                       OR a directory of per-celltype _sce.qs files (see --mode)
#   --mode            : "single" (one big SCE) | "dir" (directory of per-celltype .qs)
#   --celltype_sce_dir: directory containing per-celltype _sce.qs files to annotate
#                       (can be same as --all_sce when --mode dir)
#   --output_dir      : where to write the annotated per-celltype SCE files
#   --neuron_celltypes: comma-separated list of neuron celltype labels that can
#                       be PHF1+ (e.g. "Exc-IT-L2-3,Exc-IT-L4-5,Inh-SST")
#   --distance_um     : neighbourhood radius in microns (default 20)
#   --coord_x         : colData column for x coordinate (default x_slide_mm)
#   --coord_y         : colData column for y coordinate (default y_slide_mm)
#   --sample_col      : colData column defining spatial sample / FOV grouping
#                       (default sample_id)

library(argparse)
library(SingleCellExperiment)
library(qs)
library(Matrix)
library(cli)

# Resolve this script's directory so the shared helper is found regardless of
# the working directory the job was submitted from.
.get_script_dir <- function() {
  ca <- commandArgs(FALSE)
  f  <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
source(file.path(.get_script_dir(), "phf1_distance_utils.r"))

## ---------------------------------------------------------------------------
## Arguments
## ---------------------------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("--all_sce",         required = TRUE,
  help = "Path to single all-cell SCE (.qs) or directory of per-celltype SCEs")
parser$add_argument("--mode",            default  = "single",
  help = "single | dir  [default: single]")
parser$add_argument("--celltype_sce_dir",required = TRUE,
  help = "Directory of per-celltype _sce.qs files to annotate and re-save")
parser$add_argument("--output_dir",      required = TRUE,
  help = "Output directory for annotated per-celltype SCE files")
parser$add_argument("--neuron_celltypes",required = TRUE,
  help = "Comma-separated neuron celltype labels eligible to be PHF1+ sources")
parser$add_argument("--distance_um",     default  = 20,  type = "double",
  help = "Neighbourhood radius in microns [default: 20]")
parser$add_argument("--coord_x",         default  = "x_slide_mm",
  help = "colData column for x [default: x_slide_mm]")
parser$add_argument("--coord_y",         default  = "y_slide_mm",
  help = "colData column for y [default: y_slide_mm]")
parser$add_argument("--sample_col",      default  = "sample_id",
  help = "colData column for spatial grouping [default: sample_id]")

args <- parser$parse_args()

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

neuron_celltypes <- trimws(strsplit(args$neuron_celltypes, ",")[[1]])
dist_threshold_mm <- args$distance_um / 1000   # convert um -> mm (same units as x_slide_mm)

## ---------------------------------------------------------------------------
## 1. Load all cells into a single colData data.frame
##    We only need cell_id, celltype, sample_id, x, y, PHF1
## ---------------------------------------------------------------------------
cli_h1("Loading cell metadata")

required_cols <- c("celltype", args$sample_col, args$coord_x, args$coord_y, "PHF1")

if (args$mode == "single") {

  cli_text("Reading single all-cell SCE from {.path {args$all_sce}}")
  sce_all <- qread(args$all_sce)
  cd_all  <- as.data.frame(colData(sce_all))
  cd_all$cell_id <- colnames(sce_all)
  rm(sce_all); gc()

} else {

  # directory mode: read colData only from each per-celltype file
  sce_files <- list.files(args$all_sce, pattern = "_sce\\.qs$", full.names = TRUE)
  cli_text("Reading colData from {length(sce_files)} per-celltype SCE files")
  cd_list <- lapply(sce_files, function(f) {
    s <- qread(f)
    d <- as.data.frame(colData(s))
    d$cell_id <- colnames(s)
    rm(s); gc()
    d
  })
  cd_all <- do.call(rbind, cd_list)
  rm(cd_list); gc()

}

missing_cols <- setdiff(required_cols, colnames(cd_all))
if (length(missing_cols) > 0) {
  stop("Missing columns in colData: ", paste(missing_cols, collapse = ", "))
}

# Coerce coordinates to numeric
cd_all[[args$coord_x]] <- as.numeric(cd_all[[args$coord_x]])
cd_all[[args$coord_y]] <- as.numeric(cd_all[[args$coord_y]])
cd_all$PHF1            <- as.character(cd_all$PHF1)

cat("Total cells loaded:", nrow(cd_all), "\n")
cat("PHF1+ cells:", sum(cd_all$PHF1 == "TRUE", na.rm = TRUE), "\n")
cat("PHF1+ neurons (source):",
    sum(cd_all$PHF1 == "TRUE" & cd_all$celltype %in% neuron_celltypes, na.rm = TRUE), "\n")

## ---------------------------------------------------------------------------
## 2. For each sample_id, compute distance from every PHF1-negative cell to
##    the nearest PHF1+ neuron using a KD-tree (RANN package).
##    PHF1+ neurons themselves get NA.
## ---------------------------------------------------------------------------
cli_h1("Computing spatial distances within each sample")

samples <- unique(cd_all[[args$sample_col]])
cat("Samples to process:", paste(samples, collapse = ", "), "\n")

# The source-neuron definition and the per-sample RANN::nn2 loop live in the
# shared helper phf1_distance_utils.r (single source of truth, also used by the
# label-permutation null engine).
dist_vec <- compute_dist_to_phf1_um(
  cd               = cd_all,
  sample_col       = args$sample_col,
  neuron_celltypes = neuron_celltypes,
  coord_x          = args$coord_x,
  coord_y          = args$coord_y,
  verbose          = TRUE
)

## ---------------------------------------------------------------------------
## 3. Assign neighbourhood label
## ---------------------------------------------------------------------------
cli_h1("Assigning neighbourhood labels")

label_vec <- rep(NA_character_, nrow(cd_all))
names(label_vec) <- cd_all$cell_id

phf1_neg_idx <- which(cd_all$PHF1 != "TRUE")
label_vec[phf1_neg_idx] <- ifelse(
  dist_vec[cd_all$cell_id[phf1_neg_idx]] <= args$distance_um,
  "near",
  "distal"
)

# Summary
cat("Label counts:\n")
print(table(label_vec, useNA = "always"))
cat("\nDistance summary (PHF1-negative cells, um):\n")
print(summary(dist_vec[phf1_neg_idx]))

## ---------------------------------------------------------------------------
## 4. Add labels back to per-celltype SCE files and re-save
## ---------------------------------------------------------------------------
cli_h1("Annotating and saving per-celltype SCE files")

ct_files <- list.files(args$celltype_sce_dir, pattern = "_sce\\.qs$", full.names = TRUE)
cat("Found", length(ct_files), "celltype SCE files to annotate.\n")

for (f in ct_files) {
  ct_name <- gsub("_sce\\.qs$", "", basename(f))
  cli_text("Annotating {ct_name}")

  sce_ct <- qread(f)
  rownames(sce_ct) <- rowData(sce_ct)$gene

  cids <- colnames(sce_ct)

  # Cells in this SCE that are PHF1+ neurons — will get NA label (excluded from DEG)
  # All others get near/distal
  sce_ct$dist_to_phf1_um  <- dist_vec[cids]
  sce_ct$phf1_neighbour   <- label_vec[cids]

  cat(sprintf("  %s: near=%d, distal=%d, PHF1+/NA=%d\n",
    ct_name,
    sum(sce_ct$phf1_neighbour == "near",   na.rm = TRUE),
    sum(sce_ct$phf1_neighbour == "distal", na.rm = TRUE),
    sum(is.na(sce_ct$phf1_neighbour))
  ))

  out_path <- file.path(args$output_dir, paste0(ct_name, "_sce_neighbours.qs"))
  qsave(sce_ct, out_path)
}

cli_alert_success("Done. Annotated SCE files written to {args$output_dir}")
