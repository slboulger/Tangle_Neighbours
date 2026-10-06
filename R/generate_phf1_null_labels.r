#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# generate_phf1_null_labels.r
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
# generate_phf1_null_labels.r
#
# LABEL-PERMUTATION ENGINE (runs ONCE, globally).
#
# Produces permuted ("null") PHF1 labels and the matching distance matrices. The
# same labellings are shared by the label-permutation nulls (docs/MODELS.md,
# Permutation nulls); each downstream analysis re-runs its own specification on
# the shuffled labels.
#
# Permutation scheme
#   The PHF1 TRUE/FALSE label is shuffled WITHIN each (celltype x sample) group,
#   holding the observed PHF1+ count of that group fixed. This preserves, exactly,
#   the per-celltype-per-sample PHF1+ composition (so the categorical DEG keeps an
#   identical group structure) while destroying which individual cells carry the
#   label.
#
# Two artifacts (from the SAME shuffled labels, so label-based and distance-based
# nulls are mutually consistent):
#   1. phf1_label_permutations.qs  — cell_id x (1 + n_perm) LOGICAL matrix.
#         col "observed" = real PHF1 (== "TRUE"); cols null0001..nullNNNN = shuffled.
#         Consumed by R/null_replica_nn3.R.
#   2. phf1_null_distance_matrix.qs — cell_id x (1 + n_perm) NUMERIC matrix.
#         Distance (um) to nearest PHF1+ neuron under each labelling; NA where the
#         cell is PHF1+ under that labelling (i.e. not in the modelling set).
#         Consumed by R/null_replica_nn3.R and R/null_glia_distance_within_celltype.R.
#
# Distances are computed through the SAME helper the canonical observed pipeline
# uses (compute_dist_to_phf1_um in phf1_distance_utils.r, source_idx_override = NULL,
# RANN::nn2 k=1, mm*1000 -> um): for each labelling we swap PHF1 to the shuffled
# column and call the observed path, so the null distances are byte-identical in
# METHOD to label_phf1_neighbours.r. The observed column reproduces the canonical
# dist_to_phf1_um exactly (asserted at runtime).
#
# Use run_generate_phf1_null_labels.sh

suppressPackageStartupMessages({
  library(argparse)
  library(SingleCellExperiment)
  library(qs)
  library(cli)
  library(BiocParallel)
})

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
  help = "Path to single all-cell SCE (.qs)")
parser$add_argument("--neuron_celltypes",required = TRUE,
  help = "Comma-separated neuron celltype labels eligible to be PHF1+ distance sources")
parser$add_argument("--coord_x",         default = "x_slide_mm",
  help = "colData column for x [default: x_slide_mm]")
parser$add_argument("--coord_y",         default = "y_slide_mm",
  help = "colData column for y [default: y_slide_mm]")
parser$add_argument("--sample_col",      default = "sample_id",
  help = "colData column for spatial grouping [default: sample_id]")
parser$add_argument("--n_perm",          default = 100, type = "integer",
  help = "Number of null labellings [default: 100]")
parser$add_argument("--seed",            default = 42,  type = "integer",
  help = "RNG seed [default: 42]")
parser$add_argument("--ncores",          default = 8,   type = "integer",
  help = "Cores for BiocParallel distance recomputation [default: 8]")
parser$add_argument("--output_dir",      required = TRUE,
  help = "Output dir (deg/null_label_deg/engine)")

args <- parser$parse_args()

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
neuron_celltypes <- trimws(strsplit(args$neuron_celltypes, ",")[[1]])
set.seed(args$seed)

## ---------------------------------------------------------------------------
## 1. Load all cells into a single colData data.frame (coords + labels only)
##    Mirrors label_phf1_neighbours.r exactly.
## ---------------------------------------------------------------------------
cli_h1("Loading cell metadata")
required_cols <- c("celltype", args$sample_col, args$coord_x, args$coord_y, "PHF1")

cli_text("Reading all-cell SCE from {.path {args$all_sce}}")
sce_all <- qread(args$all_sce)
cd_all  <- as.data.frame(colData(sce_all))
cd_all$cell_id <- colnames(sce_all)
rm(sce_all); gc()

missing_cols <- setdiff(required_cols, colnames(cd_all))
if (length(missing_cols) > 0) {
  stop("Missing columns in colData: ", paste(missing_cols, collapse = ", "))
}

cd_all[[args$coord_x]] <- as.numeric(cd_all[[args$coord_x]])
cd_all[[args$coord_y]] <- as.numeric(cd_all[[args$coord_y]])
cd_all$PHF1            <- as.character(cd_all$PHF1)

observed_pos <- cd_all$PHF1 == "TRUE"
cat("Total cells loaded:", nrow(cd_all), "\n")
cat("PHF1+ cells (all celltypes):", sum(observed_pos, na.rm = TRUE), "\n")
cat("PHF1+ neurons (distance source):",
    sum(observed_pos & cd_all$celltype %in% neuron_celltypes, na.rm = TRUE), "\n")

## ---------------------------------------------------------------------------
## 2. Shuffle PHF1 within each (celltype x sample) group, count preserved.
## ---------------------------------------------------------------------------
cli_h1(sprintf("Generating %d null labellings (within celltype x sample, seed %d)",
               args$n_perm, args$seed))

grp_key <- paste(cd_all$celltype, cd_all[[args$sample_col]], sep = "\r")
groups  <- split(seq_len(nrow(cd_all)), grp_key)
k_by_grp <- vapply(groups, function(ix) sum(observed_pos[ix], na.rm = TRUE), integer(1))

label_mat <- matrix(FALSE, nrow = nrow(cd_all), ncol = args$n_perm)
for (p in seq_len(args$n_perm)) {
  for (g in seq_along(groups)) {
    ix <- groups[[g]]
    k  <- k_by_grp[g]
    if (k > 0L) {
      # ix[sample.int(length(ix), k)] avoids the sample(length-1) gotcha
      label_mat[ix[sample.int(length(ix), k)], p] <- TRUE
    }
  }
}

# Assemble with the observed (real) column first.
label_full <- cbind(observed_pos, label_mat)
colnames(label_full) <- c("observed", sprintf("null%04d", seq_len(args$n_perm)))
rownames(label_full) <- cd_all$cell_id
storage.mode(label_full) <- "logical"

# Sanity: per-group count preserved in every column.
grp_counts <- t(vapply(groups, function(ix)
  colSums(label_full[ix, , drop = FALSE]), numeric(ncol(label_full))))
if (!all(grp_counts == k_by_grp)) {
  stop("Per-(celltype x sample) PHF1+ count not preserved across null labellings.")
}
if (!all(label_full[, "observed"] == observed_pos)) {
  stop("'observed' label column does not reproduce the real PHF1 labels.")
}
cat("Groups (celltype x sample):", length(groups),
    "| non-empty PHF1+ groups:", sum(k_by_grp > 0), "\n")
cat("Label matrix dims (cells x [1 + N]):",
    paste(dim(label_full), collapse = " x "), "\n")

lbl_path <- file.path(args$output_dir, "phf1_label_permutations.qs")
qsave(label_full, lbl_path)
cat("Written:", lbl_path, "\n")

## ---------------------------------------------------------------------------
## 3. Recompute distance-to-nearest-PHF1+-neuron for each labelling.
##    Swap PHF1 to the labelling column and call the OBSERVED path of the shared
##    helper (source_idx_override = NULL) => identical method to label_phf1_neighbours.r.
##    No RNG in this step (labels already fixed) -> fully deterministic.
## ---------------------------------------------------------------------------
cli_h1(sprintf("Recomputing distances for %d labellings (%d cores)",
               ncol(label_full), args$ncores))

do_dist <- function(j) {
  cd_j <- cd_all
  cd_j$PHF1 <- ifelse(label_full[, j], "TRUE", "FALSE")
  dv <- compute_dist_to_phf1_um(
    cd               = cd_j,
    sample_col       = args$sample_col,
    neuron_celltypes = neuron_celltypes,
    coord_x          = args$coord_x,
    coord_y          = args$coord_y,
    verbose          = FALSE
  )
  dv[cd_all$cell_id]   # align to master cell order
}

bpparam <- MulticoreParam(workers = args$ncores, progressbar = TRUE)
dist_list <- bplapply(seq_len(ncol(label_full)), do_dist, BPPARAM = bpparam)

dist_mat <- do.call(cbind, dist_list)
colnames(dist_mat) <- colnames(label_full)
rownames(dist_mat) <- cd_all$cell_id
cat("Distance matrix dims (cells x [1 + N]):",
    paste(dim(dist_mat), collapse = " x "), "\n")

## ---------------------------------------------------------------------------
## 4. Validate the observed distance column against canonical dist_to_phf1_um.
## ---------------------------------------------------------------------------
data_dir  <- dirname(dirname(dirname(normalizePath(args$output_dir))))  # engine -> null_label_deg -> deg -> data_dir
neigh_dir <- file.path(data_dir, "celltype_sce_neighbours")
if (dir.exists(neigh_dir)) {
  neigh_files <- list.files(neigh_dir, pattern = "_sce_neighbours\\.qs$", full.names = TRUE)
  if (length(neigh_files) > 0) {
    cli_text("Validating observed column against {.path {basename(neigh_files[1])}}")
    s_chk <- qread(neigh_files[1])
    chk_dist <- s_chk$dist_to_phf1_um
    names(chk_dist) <- colnames(s_chk)
    shared <- intersect(names(chk_dist), rownames(dist_mat))
    cmp <- all.equal(unname(dist_mat[shared, "observed"]), unname(chk_dist[shared]))
    if (!isTRUE(cmp)) {
      stop("Observed distances DO NOT match canonical dist_to_phf1_um: ", cmp)
    }
    cat(sprintf("Validation OK: observed distances match canonical for %d shared cells.\n",
                length(shared)))
    rm(s_chk); gc()
  } else {
    cat("WARNING: no *_sce_neighbours.qs found — skipping observed-distance validation.\n")
  }
} else {
  cat("WARNING: celltype_sce_neighbours/ not found — skipping observed-distance validation.\n")
}

dist_path <- file.path(args$output_dir, "phf1_null_distance_matrix.qs")
qsave(dist_mat, dist_path)
cat("Written:", dist_path, "\n")

## ---------------------------------------------------------------------------
## 5. Provenance
## ---------------------------------------------------------------------------
params <- list(
  scheme           = "negative control: PHF1 label shuffled WITHIN each (celltype x sample) group, per-group PHF1+ count held fixed. Primary DEG pipelines are run VERBATIM on these labels (coef/contrast vs zero, standard BH). Distances recomputed via compute_dist_to_phf1_um observed path (RANN::nn2 k=1, mm*1000 -> um).",
  n_perm           = args$n_perm,
  seed             = args$seed,
  neuron_celltypes = neuron_celltypes,
  coord_x          = args$coord_x,
  coord_y          = args$coord_y,
  sample_col       = args$sample_col,
  n_groups         = length(groups),
  n_cells          = nrow(cd_all),
  label_dims       = dim(label_full),
  distance_dims    = dim(dist_mat),
  all_sce          = normalizePath(args$all_sce)
)
saveRDS(params, file.path(args$output_dir, "params.rds"))

sink(file.path(args$output_dir, "params.txt"))
cat("PHF1 null-label negative control — scheme parameters\n")
cat("====================================================\n\n")
cat("Scheme:\n  ", params$scheme, "\n\n", sep = "")
cat("N null labellings: ", params$n_perm, "\n")
cat("Seed:              ", params$seed, "\n")
cat("Distance-source neuron celltypes:\n  ",
    paste(params$neuron_celltypes, collapse = "\n  "), "\n\n", sep = "")
cat("Coordinate cols:   ", params$coord_x, ",", params$coord_y, "\n")
cat("Sample col:        ", params$sample_col, "\n")
cat("All-cell object:   ", params$all_sce, "\n")
cat("Cells:             ", params$n_cells, "\n")
cat("Groups (ct x smp): ", params$n_groups, "\n")
cat("Label matrix dims: ", paste(params$label_dims, collapse = " x "), "\n")
cat("Distance matrix dims:", paste(params$distance_dims, collapse = " x "), "\n\n")
cat("Observed PHF1+ count per (celltype x sample) group (held fixed):\n")
kt <- data.frame(group = sub("\r", " | ", names(k_by_grp)),
                 phf1_pos = as.integer(k_by_grp), row.names = NULL)
print(kt[order(-kt$phf1_pos), ], row.names = FALSE)
sink()

sink(file.path(args$output_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

cli_alert_success("Done. Null labels + distance matrix + params written to {args$output_dir}")
