#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# imc_nn3_utils.R
#
# Figure panels: 4A - shared helper, no panel of its own
# Sourced by imc_cluster_loss_vs_tangle_distance.R, which uses nearest_dist_imc().
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# imc_nn3_utils.R
#
# Shared IMC spatial helpers:
#   * nearest_dist_imc()      -- per-ROI Euclidean distance from each target to the nearest anchor
#                                (the IMC distance-to-nearest-PHF1+ neuron).
#   * the IMC 3-NEAREST-NEURONAL-NEIGHBOUR SPACING metric (nn3_um) and the geometry it needs:
#     imc_neuron_pool(), compute_nn3_imc(), local_neuron_count_imc(), imc_edge_dist().
#   * small statistics helpers: roll_mean() / blank_sparse(), p_emp_two_tailed(), iv_meta().
# The nn3 metric is the IMC counterpart of R/nn3_utils.r on the CosMx side.
#
# ---------------------------------------------------------------------------
# WHY THIS METRIC
# ---------------------------------------------------------------------------
# Zwang et al. 2024 (Cell Rep., Hyman lab; PMC11441076) imaged rTg4510 cortex longitudinally and
# found tangle-FREE neurons died at >3x the rate of tangle-bearing neurons, and that before dying
# they became more distant from their neighbours -- quantified as the 3-nearest-neighbour distance
# among labelled neurons (33.9 +/- 1.9 um dying vs 24.0 +/- 0.5 um persisting, p = 0.02, mixed
# model). Neighbour spacing is therefore a LOCAL NEURON-LOSS readout, not merely a density
# descriptor.
#
# ---------------------------------------------------------------------------
# RELATIONSHIP TO R/nn3_utils.r -- A PORT, NOT A REWRITE
# ---------------------------------------------------------------------------
# compute_nn3_imc() below is a port of compute_nn3_um() (R/nn3_utils.r) and keeps every
# non-obvious part of it: the k + 2 neighbour slack, self-exclusion BY CELL INDEX rather than by
# dropping column 1, the colsel reordering, the NA rules, and both metric invariants. Three things
# differ, all forced by the platform:
#
#   1. GROUPING IS PER ROI (sample_id). IMC ROIs are separately ablated cortical strips with no
#      spatial relationship to each other, so a neighbour can never come from another ROI. Every
#      IMC distance script groups this way.
#   2. NO UNIT CONVERSION. CosMx coordinates are x_slide_mm / y_slide_mm and are multiplied by
#      1000; IMC ablation here is 1 um per pixel, so spatialCoords() Pos_X / Pos_Y are ALREADY
#      MICRONS.
#   3. POOL MEMBERSHIP IS PASSED IN AS A MASK, not matched against a fixed vector of labels. The
#      IMC neuron pool is defined by a regex on celltype_clusters (see imc_neuron_pool()), because
#      the cluster labels carry their marker names and are not a stable enumerated set.
#
# assert_nn3_reference_equivalence() runs both implementations on the same synthetic data and
# requires identical output. Run it with
#   Rscript R/imc_nn3_utils.R
# The CosMx helper R/nn3_utils.r is left unchanged rather than given extra arguments.
suppressPackageStartupMessages({
  library(cli)
})

if (!requireNamespace("RANN", quietly = TRUE)) {
  stop("Package 'RANN' is required. Install with: install.packages('RANN')")
}

## ---------------------------------------------------------------------------
## imc_neuron_pool()
##
## The neighbour pool: every neuron, REGARDLESS OF PHF1 STATUS.
##
## This is deliberately the SAME population that defines the PHF1+ anchors in the IMC distance
## analyses -- grepl("neuron", celltype_clusters, ignore.case = TRUE) -- so nn3_um and dist_um refer
## to the same underlying neuron population. As in nn3_utils.r, do not change one without changing
## the other.
##
## Note this is NOT the target set: targets are PHF1- EXCITATORY neurons. The pool is wider
## because spacing is a property of a neuron relative to all other neurons.
##
## Returns a logical vector over rows of `cd`.
## ---------------------------------------------------------------------------
imc_neuron_pool <- function(cd, celltype_col = "celltype_clusters", verbose = TRUE) {
  if (!celltype_col %in% colnames(cd))
    stop("imc_neuron_pool(): missing column '", celltype_col, "'")
  mask <- grepl("neuron", as.character(cd[[celltype_col]]), ignore.case = TRUE)
  if (!any(mask)) stop("imc_neuron_pool(): no cells matched the neuron regex.")
  if (verbose) {
    cat("Neuron pool (nn3 neighbour set; PHF1+ and PHF1- alike):\n")
    tb <- sort(table(as.character(cd[[celltype_col]][mask])), decreasing = TRUE)
    for (i in seq_along(tb))
      cat(sprintf("  %-45s %6d\n", names(tb)[i], tb[i]))
    cat(sprintf("  %-45s %6d\n", "TOTAL", sum(mask)))
  }
  mask
}

## ---------------------------------------------------------------------------
## compute_nn3_imc()
##
## Per-ROI k-nearest-neuron spacing among the neighbour pool.
##
##   cd         : data.frame with cell_id, the coordinate columns and the ROI column. May contain
##                non-pool cells; they are dropped (they are neither queries nor neighbours).
##   pool       : logical mask over rows of cd, from imc_neuron_pool().
##   sample_col : ROI grouping column.
##   coord_x/y  : coordinate columns, ALREADY IN MICRONS (IMC is 1 um/px).
##   k          : number of nearest neighbours (3 = the Zwang et al. metric).
##
## Both the QUERY set and the NEIGHBOUR set are the pool: spacing is a property of a neuron
## relative to other neurons, so it is defined for every pool cell and is INDEPENDENT OF PHF1
## STATUS, so nn3_um does not change when PHF1 labels are shuffled under a label-permutation null.
##
## Self-exclusion is done BY CELL INDEX, not by dropping column 1: cells sharing an exact centroid
## with another pool cell do exist, so the distance-0 match in column 1 is not always the cell
## itself and dropping column 1 blindly would silently measure the wrong neighbours.
##
## Returns a data.frame: cell_id, nn3_um (MEAN distance to the k nearest -- the headline metric),
## nn3rd_um (distance to the k-th nearest -- the sensitivity metric), n_pool_roi, n_coincident
## (pool neurons at distance 0, from segmentation; reported). NA where an ROI has <= k
## pool neurons, or where a cell has more than k + 1 coincident pool neurons.
##
## nn3_um is NEVER capped or winsorised.
## ---------------------------------------------------------------------------
compute_nn3_imc <- function(cd, pool, sample_col = "sample_id",
                            coord_x = "x", coord_y = "y", k = 3L, verbose = TRUE) {

  stopifnot("cell_id" %in% colnames(cd), k >= 1L)
  for (cl in c(sample_col, coord_x, coord_y))
    if (!cl %in% colnames(cd)) stop("compute_nn3_imc(): missing column '", cl, "'")
  stopifnot(is.logical(pool), length(pool) == nrow(cd), !anyNA(pool))

  p <- cd[pool, , drop = FALSE]
  if (nrow(p) == 0) stop("compute_nn3_imc(): the pool mask selected no cells.")

  out <- data.frame(
    cell_id      = p$cell_id,
    nn3_um       = NA_real_,
    nn3rd_um     = NA_real_,
    n_pool_roi   = NA_integer_,
    n_coincident = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (smp in unique(p[[sample_col]])) {
    idx <- which(p[[sample_col]] == smp)
    n_p <- length(idx)
    out$n_pool_roi[idx] <- n_p

    if (n_p <= k) {
      if (verbose)
        cat(sprintf("  %-12s pool neurons: %-6d  <= k = %d, nn3 left NA\n", smp, n_p, k))
      next
    }

    coords <- as.matrix(p[idx, c(coord_x, coord_y)])
    # k + 2 gives slack so the self-match can be dropped by index even when a coincident cell
    # also sits at distance 0.
    kk <- min(k + 2L, n_p)
    nn <- RANN::nn2(data = coords, query = coords, k = kk)

    self     <- nn$nn.idx == seq_len(n_p)             # TRUE where the match is the cell itself
    has_self <- rowSums(self) == 1L
    s_col    <- max.col(self, ties.method = "first")  # column holding the self-match

    # Columns are sorted ascending, so removing one element keeps the order: the j-th non-self
    # neighbour is column j if j < s_col, else column j + 1.
    jj     <- seq_len(k)
    colsel <- outer(s_col, jj, function(s_, j_) ifelse(j_ >= s_, j_ + 1L, j_))
    ok_col <- has_self & (colsel[, k] <= kk)

    d <- matrix(NA_real_, nrow = n_p, ncol = k)
    if (any(ok_col)) {
      r  <- which(ok_col)
      cs <- colsel[r, , drop = FALSE]
      d[r, ] <- nn$nn.dists[cbind(rep(r, k), as.vector(cs))]   # already um: IMC is 1 um/px
    }

    out$nn3_um[idx]       <- rowMeans(d)
    out$nn3rd_um[idx]     <- d[, k]
    out$n_coincident[idx] <- rowSums(nn$nn.dists == 0) - 1L

    n_drop <- sum(!ok_col)
    if (n_drop > 0L)
      cli_alert_warning(paste0("ROI {smp}: {n_drop} cell(s) have more than k + 1 coincident ",
                               "pool neurons; nn3 left NA for them."))
    if (verbose)
      cat(sprintf("  %-12s pool neurons: %-6d  median nn%d = %6.1f um  (median %d-th = %6.1f um)  coincident cells: %d\n",
                  smp, n_p, k, median(out$nn3_um[idx], na.rm = TRUE), k,
                  median(out$nn3rd_um[idx], na.rm = TRUE),
                  sum(out$n_coincident[idx] > 0L, na.rm = TRUE)))
  }

  # Metric invariants (cheap, and they catch coordinate-unit mistakes).
  ok <- !is.na(out$nn3_um)
  stopifnot(all(out$nn3_um[ok] > 0))
  stopifnot(all(out$nn3_um[ok] <= out$nn3rd_um[ok] + 1e-9))   # mean of k <= k-th

  out
}

## ---------------------------------------------------------------------------
## nearest_dist_imc()
##
## Euclidean distance from each target to the nearest anchor WITHIN THE SAME ROI, in microns.
## The same per-ROI nearest-anchor loop as the IMC distance analyses, in vector-argument form so it
## can be called repeatedly without rebuilding frames.
##
## A target in an ROI with no anchor stays NA and is dropped downstream -- never silently set to
## some large value.
##
## This is computed rather than read because the IMC object's own PHF1 proximity
## measure is BINARY and very local (colPair(spe, "neighborhood") is a steinbock 4-px
## mask-expansion graph, median edge 19.9 um). No continuous distance is stored anywhere.
## ---------------------------------------------------------------------------
nearest_dist_imc <- function(tx, ty, ts, ax, ay, as_) {
  out <- rep(NA_real_, length(tx))
  for (s in unique(ts)) {
    i <- which(ts == s); j <- which(as_ == s)
    if (!length(j)) next
    out[i] <- RANN::nn2(cbind(ax[j], ay[j]), cbind(tx[i], ty[i]), k = 1)$nn.dists[, 1]
  }
  out
}

## ---------------------------------------------------------------------------
## local_neuron_count_imc()
##
## Neurons within `radius_um`, self excluded -- CONSTRUCT-VALIDITY DIAGNOSTIC ONLY, never a
## headline model term. nn3_um should be strongly NEGATIVELY correlated with it; if it is not, the
## metric is wrong.
##
## The count SATURATES at k_max - 1 and is uncorrected for that cap, which is why it is a diagnostic.
##
## Queries are `q_idx` rows of cd; neighbours are the `pool` rows.
## ---------------------------------------------------------------------------
local_neuron_count_imc <- function(cd, pool, q_idx, sample_col = "sample_id",
                                   coord_x = "x", coord_y = "y",
                                   radius_um = 50, k_max = 60) {
  out <- rep(NA_real_, length(q_idx))
  pidx <- which(pool)
  for (s in unique(cd[[sample_col]][q_idx])) {
    qi <- which(cd[[sample_col]][q_idx] == s)
    pi <- pidx[cd[[sample_col]][pidx] == s]
    if (length(pi) < 2L) next
    q <- RANN::nn2(data  = cbind(cd[[coord_x]][pi], cd[[coord_y]][pi]),
                   query = cbind(cd[[coord_x]][q_idx[qi]], cd[[coord_y]][q_idx[qi]]),
                   k = min(k_max, length(pi)))
    # self excluded: every query here is itself a pool cell, so it matches at distance 0
    out[qi] <- rowSums(q$nn.dists <= radius_um) - 1
  }
  out
}

## ---------------------------------------------------------------------------
## imc_edge_dist()
##
## Distance from each cell to the edge of its own ROI, in microns, plus a capped copy.
##
## WHY THIS IS NEEDED, AND WHY IT MATTERS MORE HERE THAN IN COSMX. IMC ROIs are separately ablated
## cortical strips, median 478 x 2525 um. A neuron near the long
## edge of a 478 um-wide strip has some of its true nearest neighbours OUTSIDE the imaged region,
## so its nn3_um is inflated by edge censoring. With nn3 on the order of tens of microns that exposes
## a substantial share of every strip -- a far larger fraction than the CosMx FOV mosaic did.
##
## edge_dist quantifies the exposure and enters the model as a fixed effect. As in R/nn3_utils.r:
## ADJUST RATHER THAN EXCLUDE, so no cells are lost.
##
## METHOD. The plain per-ROI bounding-box distance, pmin over the four sides. No raster machinery is needed:
## an IMC ROI is a single ablated rectangle, unlike the CosMx FOV mosaic with interior holes that
## forced fov_edge_dist_um()'s occupancy raster.
##
## Values are capped at `cap_um` because edge censoring cannot affect an nn3_um smaller than the
## distance to the nearest unimaged region; beyond the cap the covariate would add leverage
## without meaning. Callers must check the cap sits above the upper tail of nn3rd_um.
##
## Returns the input frame with edge_dist_um and edge_dist_c added.
## ---------------------------------------------------------------------------
imc_edge_dist <- function(cd, sample_col = "sample_id", coord_x = "x", coord_y = "y",
                          cap_um = 150) {
  for (cl in c(sample_col, coord_x, coord_y))
    if (!cl %in% colnames(cd)) stop("imc_edge_dist(): missing column '", cl, "'")
  cd$edge_dist_um <- NA_real_
  for (s in unique(cd[[sample_col]])) {
    i  <- which(cd[[sample_col]] == s)
    xx <- cd[[coord_x]][i]; yy <- cd[[coord_y]][i]
    cd$edge_dist_um[i] <- pmin(xx - min(xx), max(xx) - xx,
                               yy - min(yy), max(yy) - yy)
  }
  cd$edge_dist_c <- pmin(cd$edge_dist_um, cap_um)
  cd
}

## ---------------------------------------------------------------------------
## roll_mean() / blank_sparse()
##
## The project's rolling-window mean, as in plot_modulescore_vs_phf1_distance_modelp.R, so
## distance curves are drawn the same way throughout. Model-free; SEM is over CELLS and is
## descriptive only.
## ---------------------------------------------------------------------------
roll_mean <- function(x, y, grid, half_window) {
  out <- lapply(grid, function(g) {
    idx <- which(x >= g - half_window & x <= g + half_window); n <- length(idx)
    if (n < 1) return(c(roll_mean = NA_real_, sem = NA_real_, n_window = 0))
    m <- mean(y[idx]); sv <- if (n > 1) stats::sd(y[idx]) / sqrt(n) else NA_real_
    c(roll_mean = m, sem = sv, n_window = n)
  })
  as.data.frame(do.call(rbind, out))
}

blank_sparse <- function(r, min_n = 10L) {
  r$roll_mean[r$n_window < min_n] <- NA
  r$sem[r$n_window < min_n]       <- NA
  r
}

## ---------------------------------------------------------------------------
## p_emp_two_tailed()
##
## Two-tailed empirical p-value of `obs` against the null's OWN distribution.
##
## This is NOT the (1 + #{|null| >= |obs|}) / (1 + n) form, which assumes a ZERO-CENTRED null.
## Distance-to-tangle is itself a nearest-neighbour distance to a subset of neurons, so an nn3 null
## need not be zero-centred, and comparing absolute magnitudes would ask the wrong question. What
## matters is where the observed coefficient sits inside the null's own distribution.
##
## Floor is 2 / (1 + n_null).
## ---------------------------------------------------------------------------
p_emp_two_tailed <- function(b, obs) {
  b <- b[is.finite(b)]
  if (!length(b) || !is.finite(obs)) return(NA_real_)
  n  <- length(b)
  lo <- (1 + sum(b <= obs)) / (1 + n)
  hi <- (1 + sum(b >= obs)) / (1 + n)
  min(1, 2 * min(lo, hi))
}

## ---------------------------------------------------------------------------
## iv_meta()
##
## Inverse-variance meta-analysis of per-donor estimates. The donor-level counterpart to a
## cell-level LMM coefficient, from the weighting block in R/nn3_phf1_vs_neg.R.
##
## WHY IT IS REPORTED. A cell-level fit on ~10^4 cells from 44 donors has p-values that are not
## magnitude cues. An UNWEIGHTED mean of donor slopes is dominated by donors carrying almost no
## PHF1+ anchors, while the LMM weights by information. The inverse-variance estimate is reported
## alongside the LMM so the donor-level and cell-level estimates can be compared.
##
## Returns estimate, se, ci_lo, ci_hi, p, n_donors and n_same_direction.
## ---------------------------------------------------------------------------
iv_meta <- function(est, se) {
  ok <- is.finite(est) & is.finite(se) & se > 0
  if (sum(ok) < 2L)
    return(list(est = NA_real_, se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                p = NA_real_, n_donors = sum(ok), n_same_dir = NA_integer_))
  e <- est[ok]; s <- se[ok]
  w  <- 1 / s^2
  mu <- sum(w * e) / sum(w)
  sm <- sqrt(1 / sum(w))
  list(est = mu, se = sm, ci_lo = mu - 1.96 * sm, ci_hi = mu + 1.96 * sm,
       p = 2 * stats::pnorm(-abs(mu / sm)), n_donors = length(e),
       n_same_dir = sum(sign(e) == sign(mu)))
}

## ---------------------------------------------------------------------------
## assert_nn3_reference_equivalence()
##
## Checks that compute_nn3_imc() is the same estimator as the CosMx compute_nn3_um(). Feeds the CosMx version mm coordinates so its internal * 1000
## cancels, and requires identical output on the same points.
##
## Also exercises the two cases the index-based self-exclusion exists for: exact coincident
## centroids, and an ROI with <= k pool neurons.
## ---------------------------------------------------------------------------
assert_nn3_reference_equivalence <- function(seed = 1L, k = 3L) {
  set.seed(seed)
  n <- 400L
  cd <- data.frame(
    cell_id   = sprintf("c%04d", seq_len(n)),
    sample_id = rep(c("R1", "R2", "R3"), length.out = n),
    x         = runif(n, 0, 500),
    y         = runif(n, 0, 2500),
    celltype_clusters = rep(c("Excitatory neuron cluster 1 (MAP2)",
                              "Astrocyte cluster 2 (S100b)"), length.out = n),
    stringsAsFactors = FALSE
  )
  # exact coincident centroids: the case that makes "drop column 1" wrong
  cd$x[2] <- cd$x[1]; cd$y[2] <- cd$y[1]
  cd$celltype_clusters[1:2] <- "Excitatory neuron cluster 1 (MAP2)"
  cd$sample_id[1:2] <- "R1"
  # a tiny ROI with <= k pool neurons, which must come back NA
  cd$sample_id[n - 2L] <- "R4"; cd$celltype_clusters[n - 2L] <- "Excitatory neuron cluster 1 (MAP2)"

  pool <- imc_neuron_pool(cd, verbose = FALSE)
  got  <- compute_nn3_imc(cd, pool, k = k, verbose = FALSE)

  ref_src <- "R/nn3_utils.r"
  if (!file.exists(ref_src)) stop("Cannot find ", ref_src, " to compare against.")
  ref_env <- new.env(parent = globalenv())
  sys.source(ref_src, envir = ref_env)

  # The CosMx helper matches celltype against a vector of labels and multiplies by 1000, so give
  # it mm coordinates and the pool labels actually present.
  cd_mm <- cd
  cd_mm$celltype    <- cd$celltype_clusters
  cd_mm$x_slide_mm  <- cd$x / 1000
  cd_mm$y_slide_mm  <- cd$y / 1000
  want <- ref_env$compute_nn3_um(cd_mm, sample_col = "sample_id",
                                 neuron_pool = unique(cd$celltype_clusters[pool]),
                                 k = k, verbose = FALSE)

  stopifnot(identical(got$cell_id, want$cell_id))
  for (cl in c("nn3_um", "nn3rd_um", "n_pool_roi", "n_coincident")) {
    ref_cl <- if (cl == "n_pool_roi") "n_pool_sample" else cl
    a <- got[[cl]]; b <- want[[ref_cl]]
    stopifnot(identical(is.na(a), is.na(b)))
    if (is.numeric(a))
      stopifnot(isTRUE(all.equal(a[!is.na(a)], b[!is.na(b)], tolerance = 1e-9)))
    else
      stopifnot(identical(a, b))
  }

  # the two exercised edge cases actually happened
  n_coin <- sum(got$n_coincident > 0L, na.rm = TRUE)
  n_na   <- sum(is.na(got$nn3_um))
  if (n_coin < 1L) stop("Self-test did not exercise the coincident-centroid path.")
  if (n_na  < 1L) stop("Self-test did not exercise the <= k pool neurons path.")

  cat(sprintf("compute_nn3_imc() == compute_nn3_um() on %d pool cells (%d coincident, %d NA).\n",
              nrow(got), n_coin, n_na))
  cat(sprintf("  median nn3 = %.2f um, median nn3rd = %.2f um\n",
              median(got$nn3_um, na.rm = TRUE), median(got$nn3rd_um, na.rm = TRUE)))
  invisible(TRUE)
}

# Run the self-test when executed directly (Rscript R/imc_nn3_utils.R), not when sourced.
# sys.nframe() is 0 only at the top level of a script; source() puts us inside a frame.
if (sys.nframe() == 0L && !interactive()) {
  setwd("<PROJECT_ROOT>/phf1_v2")
  cat("=== imc_nn3_utils.R self-test ===\n")
  assert_nn3_reference_equivalence()
  cat("PASS\n")
}
