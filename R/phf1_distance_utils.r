#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# phf1_distance_utils.r
#
# Shared utility - sourced by the scripts above, no panel of its own
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# phf1_distance_utils.r
#
# Shared helper: single source of truth for the PHF1+ source-neuron definition
# and the RANN::nn2 distance-to-nearest-PHF1+-neuron computation.
#
# Used by label_phf1_neighbours.r (observed labelling, canonical) and by the
# spatial label-permutation null engine, so that observed and permuted distances
# are computed through ONE implementation and can never drift apart.
#
# The observed path (source_idx_override = NULL): source = PHF1+ neurons in
# eligible celltypes, query = all PHF1-negative cells, k = 1, distances in
# mm * 1000 -> um, PHF1+ cells -> NA.
#
# The permutation path (source_idx_override supplied) keeps the query set fixed
# (observed PHF1-negative cells) but swaps in a caller-chosen source set; see
# compute_dist_to_phf1_um() for the self-match handling.

suppressPackageStartupMessages({
  library(cli)
})

if (!requireNamespace("RANN", quietly = TRUE)) {
  stop("Package 'RANN' is required. Install with: install.packages('RANN')")
}

## ---------------------------------------------------------------------------
## Source-neuron definition
##   PHF1+ neurons in eligible celltypes. `cd` must have character columns
##   `PHF1` (== "TRUE") and `celltype`. Returns integer row positions into `cd`.
## ---------------------------------------------------------------------------
phf1_source_idx <- function(cd, neuron_celltypes) {
  which(cd$PHF1 == "TRUE" & cd$celltype %in% neuron_celltypes)
}

## ---------------------------------------------------------------------------
## Per-sample distance-to-nearest-PHF1+-neuron
##
##   cd                 : data.frame with cell_id, PHF1 (chr), celltype,
##                        the coordinate columns and the sample column.
##   sample_col         : name of the sample grouping column.
##   neuron_celltypes   : character vector of eligible source celltypes.
##   coord_x, coord_y   : coordinate column names (mm; converted to um here).
##   source_idx_override: NULL for the observed source set (PHF1+ neurons), OR a
##                        named list keyed by sample value, each element an
##                        integer vector of row positions WITHIN that sample's
##                        subset (cd_smp, rows in their original order) to use as
##                        the source set instead. Used by the permutation engine.
##   verbose            : emit per-sample cli logging (off in the hot loop).
##
## Returns a numeric vector named by cell_id (NA for cells that are not in the
## query set, i.e. observed PHF1+ cells).
## ---------------------------------------------------------------------------
compute_dist_to_phf1_um <- function(cd, sample_col, neuron_celltypes,
                                     coord_x = "x_slide_mm",
                                     coord_y = "y_slide_mm",
                                     source_idx_override = NULL,
                                     verbose = TRUE) {

  stopifnot("cell_id" %in% colnames(cd))
  use_override <- !is.null(source_idx_override)

  samples  <- unique(cd[[sample_col]])
  dist_vec <- rep(NA_real_, nrow(cd))
  names(dist_vec) <- cd$cell_id

  for (smp in samples) {
    if (verbose) cli_text("Processing sample: {smp}")

    idx_smp <- which(cd[[sample_col]] == smp)
    cd_smp  <- cd[idx_smp, ]

    # Source set: observed PHF1+ neurons, OR caller-supplied (permuted) indices.
    if (use_override) {
      idx_source <- source_idx_override[[as.character(smp)]]
      if (is.null(idx_source)) idx_source <- integer(0)
    } else {
      idx_source <- phf1_source_idx(cd_smp, neuron_celltypes)
    }

    # Query set is ALWAYS the observed PHF1-negative cells in this sample.
    idx_query <- which(cd_smp$PHF1 != "TRUE")

    if (verbose) {
      cat(sprintf("  PHF1+ neuron sources: %d | PHF1-negative query cells: %d\n",
                  length(idx_source), length(idx_query)))
    }

    if (length(idx_source) == 0) {
      if (verbose) cat("  No PHF1+ neurons in this sample — all query cells will be NA.\n")
      next
    }
    if (length(idx_query) == 0) {
      if (verbose) cat("  No PHF1-negative cells in this sample — skipping.\n")
      next
    }

    coords_source <- as.matrix(cd_smp[idx_source, c(coord_x, coord_y)])
    coords_query  <- as.matrix(cd_smp[idx_query,  c(coord_x, coord_y)])
    query_cell_ids <- cd_smp$cell_id[idx_query]

    if (!use_override) {
      # ---- Observed path ------------------------------------------------
      nn_result <- RANN::nn2(data = coords_source, query = coords_query, k = 1)
      dists_um  <- nn_result$nn.dists[, 1] * 1000
    } else {
      # ---- Permutation path: query set fixed, source set permuted -------
      # A query cell (observed PHF1-negative) may coincide with a permuted
      # source cell (PHF1-negative eligible neuron chosen as source). Detect
      # such self-matches by cell_id and take the nearest NON-self source so
      # no spurious zero-distances and no NAs enter the matrix.
      k <- min(2L, nrow(coords_source))
      nn_result <- RANN::nn2(data = coords_source, query = coords_query, k = k)
      source_cell_ids <- cd_smp$cell_id[idx_source]
      nearest_ids <- source_cell_ids[nn_result$nn.idx[, 1]]
      is_self <- nearest_ids == query_cell_ids
      dnn <- nn_result$nn.dists[, 1]
      if (any(is_self) && k >= 2L) {
        dnn[is_self] <- nn_result$nn.dists[is_self, 2]
      }
      dists_um <- dnn * 1000
    }

    dist_vec[query_cell_ids] <- dists_um
  }

  dist_vec
}
