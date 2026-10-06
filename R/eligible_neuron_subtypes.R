# ---------------------------------------------------------------------------
# eligible_neuron_subtypes.R
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
# ---------------------------------------------------------------------
# eligible_neuron_subtypes(): neuronal subtypes well-represented for PHF1
# analysis -- those with >= min_phf1_cells PHF1+ cells in >= min_samples
# samples. Single source of truth for the "PHF1-filtered" figure variants
# (nDEG / nPathway by celltype), so every script agrees on the criterion.
#
# Reads the precomputed per-sample x per-celltype PHF1 count table written by
# R/plot_phf1_neuron_density_by_celltype.R (columns include celltype,
# sample_id, n_phf1_pos; a full zero-filled sample x neuron_order grid). This
# avoids an 8.4 GB readRDS(seu_PHF1.rds) in a lightweight plotting script. If
# seu_PHF1 changes, regenerate that plot first, then rerun these scripts.
#
# Returns a character vector of celltype labels (a subset of neuron_order when
# neuron_levels is supplied, preserving that order).
#
# Threshold sensitivity: >= 3 PHF1+ cells in >= 5 samples gives the same eligible
# set as the default >= 6 in >= 5 -- exactly Exc-IT-L2-3-CBLN2-HOPX and
# Exc-IT-L3-5-CHGA-IL1RAPL2.
# Per-subtype samples meeting each threshold (from the count table below):
#   Exc-IT-L2-3-CBLN2-HOPX     7 at >=3, 6 at >=6   (60,37,23,21,20,18,3,2,1)
#   Exc-IT-L3-5-CHGA-IL1RAPL2  7 at >=3, 5 at >=6   (25,11,11,6,6,4,4,2,1)
#   next best (Inh-VIP)        4 at >=3, 0 at >=6   (5,4,4,3,2,2,2,1,0)
# Exc-IT-L3-5 sits exactly on the boundary at 6 in 5 samples, so do not raise
# min_phf1_cells past 6 (or min_samples past 5) without expecting it to drop
# out and every PHF1-filtered figure to change.
# ---------------------------------------------------------------------
eligible_neuron_subtypes <- function(project_root,
                                     min_phf1_cells = 6,
                                     min_samples    = 5,
                                     counts_tsv = "plots/phf1_neuron_density/source_data_phf1_neuron_density_by_celltype.tsv",
                                     neuron_levels = NULL) {
  path <- file.path(project_root, counts_tsv)
  if (!file.exists(path)) stop("PHF1 count table not found: ", path)
  cnt <- read.delim(path, sep = "\t", header = TRUE, check.names = FALSE)
  if (!all(c("celltype", "n_phf1_pos") %in% names(cnt)))
    stop("Expected columns 'celltype' and 'n_phf1_pos' in ", path)

  # samples in which each subtype has >= min_phf1_cells PHF1+ cells
  ok   <- cnt[!is.na(cnt$n_phf1_pos) & cnt$n_phf1_pos >= min_phf1_cells, , drop = FALSE]
  tab  <- table(ok$celltype)
  elig <- names(tab)[tab >= min_samples]

  if (!is.null(neuron_levels)) elig <- neuron_levels[neuron_levels %in% elig]  # keep order
  elig
}
