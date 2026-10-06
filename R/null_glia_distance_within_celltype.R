#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# null_glia_distance_within_celltype.R
#
# Figure panels: 6A
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# null_glia_distance_within_celltype.R
#
# PRECOMPUTE for the within-celltype observed-vs-null glia distance test.
#
# This script computes the statistic WITHIN each celltype: the mean distance to
# the nearest PHF1+ neuron for that glial type alone, observed vs the same
# quantity under each PHF1-label shuffle.
# No comparison to the other glia enters it.
#
# Two statistics per (celltype, draw) are written, because they weight donors
# differently:
#   mean_of_donor_means -- unweighted mean over the 9 donors' per-donor means.
#     PRIMARY: this is exactly what the crossbar in the canonical glia panel draws.
#   pooled_cell_mean    -- cell-level mean pooled across donors. Secondary; donors
#     contributing more cells of that type dominate it.
#
# Outputs (small TSVs consumed by R/plot_dist_to_phf1_by_glia_vs_null.R):
#   null_glia_within_celltype_per_draw.tsv    celltype x draw (incl. "observed")
#   null_glia_within_celltype_donor_draw.tsv  celltype x donor x draw
#
# Reads the 1.1 GB null distance matrix, so run this once; the plotting script
# then only touches the TSVs.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr)
  library(SingleCellExperiment); library(qs)
})

hpc <- "<PROJECT_ROOT>/phf1_v2"
loc <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(hpc)) hpc else loc)
source("R/palettes.R")

out_dir <- "plots/dist_to_phf1_glia_vs_null"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

f_mat <- "deg/null_label_deg/engine/phf1_null_distance_matrix.qs"
if (!file.exists(f_mat)) stop("Missing null distance matrix: ", f_mat)

## ---------------------------------------------------------------------------
## Glia cell-level metadata, taken from the per-celltype SCE neighbour objects
## (cheap: ~65 MB total) rather than the 8.4 GB Seurat object. Rownames of the
## SCE are the cell ids used by the null distance matrix.
## ---------------------------------------------------------------------------
glia_files <- file.path("celltype_sce_neighbours", paste0(glia_order, "_sce_neighbours.qs"))
names(glia_files) <- glia_order
missing <- glia_files[!file.exists(glia_files)]
if (length(missing)) stop("Missing glia SCE: ", paste(missing, collapse = ", "))

md <- bind_rows(lapply(glia_order, function(g) {
  s  <- qs::qread(glia_files[[g]], nthreads = 2)
  cd <- as.data.frame(SummarizedExperiment::colData(s))
  out <- tibble(cell_id   = rownames(cd),
                celltype  = as.character(cd$celltype),
                sample_id = as.character(cd$sample_id),
                PHF1      = cd$PHF1,
                dist_obs  = as.numeric(cd$dist_to_phf1_um))
  rm(s); gc()
  out
}))

# The canonical glia panel models PHF1-NEGATIVE glia only. Match that exactly.
md <- md |> filter(!PHF1, !is.na(dist_obs))
cat("Glia celltypes:", paste(sort(unique(md$celltype)), collapse = ", "), "\n")
cat("PHF1-negative glia cells:", nrow(md), "| donors:", length(unique(md$sample_id)), "\n")
if (nrow(md) != 63550L)
  warning("Expected 63550 PHF1-negative glia to match the canonical panel, got ", nrow(md))

## ---------------------------------------------------------------------------
## Null distance matrix, subset to those cells.
## ---------------------------------------------------------------------------
cat("Reading", f_mat, "...\n"); flush.console()
dist_full <- qs::qread(f_mat, nthreads = 4)
cat("Matrix:", nrow(dist_full), "x", ncol(dist_full), "\n")
if (!all(md$cell_id %in% rownames(dist_full)))
  stop("Some glia cells absent from the null distance matrix.")

dist_sub <- dist_full[md$cell_id, , drop = FALSE]
rm(dist_full); gc()

draw_names <- colnames(dist_sub)
if (!"observed" %in% draw_names) stop("No 'observed' column in the distance matrix.")
null_names <- grep("^null", draw_names, value = TRUE)
cat("Null draws:", length(null_names), "\n")

# Sanity: the matrix's observed column must reproduce the SCE distances.
obs_delta <- max(abs(as.numeric(dist_sub[, "observed"]) - md$dist_obs), na.rm = TRUE)
cat("Max |matrix observed - SCE dist_to_phf1_um|:", signif(obs_delta, 3), "um\n")
if (obs_delta > 1e-6)
  warning("Observed column disagrees with the SCE distances by up to ", obs_delta, " um.")

## ---------------------------------------------------------------------------
## Per (celltype x donor x draw) means, then the two celltype-level statistics.
## ---------------------------------------------------------------------------
glia_present <- sort(unique(md$celltype))
keys  <- paste(md$celltype, md$sample_id, sep = "|")
ukeys <- unique(keys)
idx_by_key <- split(seq_len(nrow(md)), factor(keys, levels = ukeys))

# NAs: PHF1+ labels are shuffled within (celltype x sample), so in a given null
# draw a glial cell that is PHF1-negative in the observed data can itself be
# labelled PHF1+, and the matrix then carries no distance for it. These cells are
# dropped per draw, and the rate is reported below.
na_by_draw <- colSums(is.na(dist_sub))
cat("Cell-draw entries with no distance: ",
    signif(100 * sum(na_by_draw) / length(dist_sub), 3), "% overall; ",
    "max per draw ", max(na_by_draw), " of ", nrow(dist_sub), " cells\n", sep = "")
na_by_ct <- tapply(seq_len(nrow(md)), md$celltype,
                   function(ii) sum(is.na(dist_sub[ii, null_names, drop = FALSE])) /
                                (length(ii) * length(null_names)))
cat("Per-celltype NA rate across null draws:\n")
print(round(100 * na_by_ct, 3))

donor_stats <- lapply(idx_by_key, function(ii) {
  sub <- dist_sub[ii, , drop = FALSE]
  list(m = colMeans(sub, na.rm = TRUE), n = colSums(!is.na(sub)))
})
donor_mat <- t(vapply(donor_stats, function(x) x$m, numeric(ncol(dist_sub))))
n_mat     <- t(vapply(donor_stats, function(x) x$n, numeric(ncol(dist_sub))))
colnames(donor_mat) <- colnames(n_mat) <- draw_names

long_mean <- as_tibble(donor_mat, rownames = "key") |>
  pivot_longer(cols = all_of(draw_names), names_to = "draw", values_to = "mean_dist")
long_n <- as_tibble(n_mat, rownames = "key") |>
  pivot_longer(cols = all_of(draw_names), names_to = "draw", values_to = "n_cells")

donor_draw <- long_mean |>
  left_join(long_n, by = c("key", "draw")) |>
  separate(key, into = c("celltype", "sample_id"), sep = "\\|") |>
  mutate(series = ifelse(draw == "observed", "observed", "null"))

if (any(is.na(donor_draw$mean_dist)))
  stop("A (celltype x donor x draw) cell had no usable cells at all.")

per_draw <- donor_draw |>
  group_by(celltype, draw, series) |>
  summarise(mean_of_donor_means = mean(mean_dist),
            pooled_cell_mean    = sum(mean_dist * n_cells) / sum(n_cells),
            n_donors            = n(),
            n_cells             = sum(n_cells),
            .groups = "drop")

write.table(per_draw, file.path(out_dir, "null_glia_within_celltype_per_draw.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(donor_draw |> select(celltype, sample_id, draw, series, mean_dist, n_cells),
            file.path(out_dir, "null_glia_within_celltype_donor_draw.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nObserved mean-of-donor-means vs null median, per glia:\n")
chk <- per_draw |>
  group_by(celltype) |>
  summarise(observed = mean_of_donor_means[series == "observed"],
            null_med = median(mean_of_donor_means[series == "null"]),
            null_sd  = sd(mean_of_donor_means[series == "null"]),
            .groups = "drop") |>
  mutate(delta = observed - null_med, z = delta / null_sd)
print(as.data.frame(chk), row.names = FALSE)

cat("\nWrote:\n  ", file.path(out_dir, "null_glia_within_celltype_per_draw.tsv"),
    "\n  ", file.path(out_dir, "null_glia_within_celltype_donor_draw.tsv"), "\n", sep = "")
