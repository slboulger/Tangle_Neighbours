#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# phf1_hidden_anchors.R
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
# phf1_hidden_anchors.R
#
# ARM B (hidden-anchor augmentation) and ARM C (admixture bound) for the
# out-of-plane tangle correction. Unlike a sham-anchor null, which REPLACES the
# anchor set with an equal-sized matched set, augmentation ADDS to it.
#
# ============================================================================
# ARM B DILUTION
# ============================================================================
# Hidden anchors are a POSITION-WEIGHTED RANDOM DRAW, not an identification of
# which cells are actually mislabelled. Arm B therefore does NOT add back the
# missing anchors -- it adds ~133 anchors at plausible-but-wrong positions. The
# expected consequence is ATTENUATION OF THE DISTANCE SLOPE THROUGH ANCHOR-SET
# DILUTION, independent of whether the field effect is real.
#
# So the arm is read in ONE DIRECTION:
#
#     survival under B  -> robustness
#     weakening under B -> expected from dilution either way
#
# The size of that dilution is measured:
# R/test_tangle_3d_utils.R section 4 plants a known field, thins the anchors,
# adds them back by weighted draw and reports
#
#     ARM B DILUTION CONSTANT = slope_aug / slope_obs = 0.906 [0.896, 0.915]
#
# written to results/tangle_3d/armB_dilution_constant.tsv. Every Arm B slope is
# reported against that number; a ~9% weakening is what the PROCEDURE itself
# produces. arm_b_dilution_constant() reads it back so no caller re-derives it.
#
# ============================================================================
# WHY ARM C SELECTS ON INTENSITY AND NOT ON POSITION
# ============================================================================
# The same randomness breaks the obvious Arm C design. Dropping a
# position-selected RANDOM subset of near cells leaves the slope UNBIASED IN
# EXPECTATION, so it returns "no change" by construction and bounds nothing.
#
# To bound admixture the drop set must be ENRICHED FOR GENUINELY TANGLE-BEARING
# CELLS. select_intensity_drop() ranks PHF1-negative cells by per-cell PHF1 p95
# intensity RESIDUALISED WITHIN DISTANCE STRATA and drops the top H.
#
# The residualising is required. Per-cell PHF1 intensity varies with distance to
# the nearest tangle, so ranking on raw p95 would select near cells and make the
# arm a distance selection in disguise. Within-stratum residuals ask "bright FOR
# ITS DISTANCE".
#
# The paired PLACEBO drops H position-weighted cells at matched distances and
# should move nothing. THE ADMIXTURE BOUND IS THE DIFFERENCE BETWEEN THE TWO,
# not the intensity arm on its own.
#
# Exports:
#   draw_hidden_anchors()      one replicate -> source_idx_override (obs + hidden)
#   hidden_dist_matrix()       R replicates  -> n_cells x R distance matrix
#   select_intensity_drop()    Arm C drop set + matched placebo
#   arm_b_dilution_constant()  read the calibration constant from disk
#   hidden_balance_table()     hidden vs observed anchors, per covariate SMD

suppressPackageStartupMessages({
  library(stats)
})

if (!requireNamespace("RANN", quietly = TRUE))
  stop("Package 'RANN' is required. Install with: install.packages('RANN')")

HIDDEN_WEIGHTS   <- c("intensity", "uniform")
HIDDEN_SIGMA_UM  <- 200      # kernel sigma for the anchor-intensity weighting


## ---------------------------------------------------------------------------
## arm_b_dilution_constant()
##
## Read the calibration constant produced by R/test_tangle_3d_utils.R. Raises if
## it is missing rather than defaulting to 1, because a silent 1 would present
## the procedure's own attenuation as a finding about the field.
## ---------------------------------------------------------------------------
arm_b_dilution_constant <- function(
    path = file.path("results", "tangle_3d", "armB_dilution_constant.tsv")) {
  if (!file.exists(path))
    stop("Arm B dilution constant not found at ", path, ".\n",
         "Run: Rscript R/test_tangle_3d_utils.R --reps 40\n",
         "Arm B results are not interpretable without it -- the procedure ",
         "attenuates the slope on its own and the size of that attenuation is ",
         "exactly what this constant records.")
  tb <- utils::read.delim(path, stringsAsFactors = FALSE)
  setNames(tb$value, tb$quantity)
}


## ---------------------------------------------------------------------------
## .anchor_intensity()
##
## Kernel-smoothed observed PHF1+ anchor intensity at every candidate position.
## This is the weight that makes the augmented pattern inherit the observed
## clustering.
##
## MEASURED, AND IT IS THE OPPOSITE WAY ROUND FROM THE OBVIOUS GUESS. At
## p_detect = 0.75 on the real cohort (398 observed anchors, H = 131):
##
##     intensity-weighted   median d2D 304.6 -> 276.6 um   ( -9.2% )
##     uniform              median d2D 304.6 -> 259.0 um   ( -15.0% )
##
## Intensity weighting shortens the median LESS, not more. Hidden anchors placed
## near existing anchors are REDUNDANT -- they cover ground that is already
## covered -- whereas uniformly placed ones reach cells that had no anchor
## nearby. (Uniform lands within 2 um of the naive Poisson prediction of 264 um
## for a 33% denser anchor set, exactly as it should.)
##
## So uniform is the UPPER bound on the distance shortening and intensity is the
## realistic case, not the reverse. Since real missed tangles are clustered with
## the observed ones, INTENSITY IS THE HEADLINE and uniform is the optimistic
## bracket. This is the same mechanism as the clustering result in
## R/test_tangle_3d_utils.R section 3: clustering anchors leaves voids, and
## voids lengthen typical distances.
##
## Uses kernel_field_rann(score_mat = NULL), which returns the summed kernel
## weight and already guards RANN's silent radius truncation. Falls back to a
## direct truncated-Gaussian sum if kernel_field_utils.R was not sourced, so
## this file has no hard load-order dependency.
## ---------------------------------------------------------------------------
.anchor_intensity <- function(cand_xy_mm, anchor_xy_mm, sigma_um) {
  if (nrow(anchor_xy_mm) == 0L) return(rep(0, nrow(cand_xy_mm)))
  if (exists("kernel_field_rann", mode = "function")) {
    return(kernel_field_rann(cand_xy_mm, anchor_xy_mm, NULL, sigma_um,
                             verbose = FALSE)$dens)
  }
  k  <- min(nrow(anchor_xy_mm), 256L)
  nn <- RANN::nn2(anchor_xy_mm, cand_xy_mm, k = k)
  d  <- nn$nn.dists * 1000
  d[nn$nn.idx == 0] <- Inf
  rowSums(exp(-0.5 * (d / sigma_um)^2))
}


## ---------------------------------------------------------------------------
## draw_hidden_anchors()
##
## ONE replicate. Returns a named list keyed by sample value, each element the
## integer row positions WITHIN that sample's subset of `cd` -- precisely the
## `source_idx_override` contract of compute_dist_to_phf1_um() in
## R/phf1_distance_utils.r. Do not reimplement the nn2 call; pass this straight
## in. The existing caller of that contract is R/generate_phf1_null_labels.r.
##
## The returned index vector is c(OBSERVED anchors, HIDDEN anchors). THE
## OBSERVED PLANE STAYS GROUND TRUTH: no observed label is altered, only added
## to. That is what distinguishes this from a permutation null.
##
## Per sample:
##   n_obs     = observed PHF1+ neurons in the eligible celltypes
##   H         = n_obs (1 - p_detect)/p_detect, the axially-missed count
##   candidates= PHF1-NEGATIVE neurons of the SAME eligible subtypes, so the
##               augmented set keeps the observed celltype composition
##
## The self-match hazard is real and is handled downstream:
## compute_dist_to_phf1_um()'s override path takes k = 2 and drops a zero
## distance when the nearest "anchor" is the query cell itself. Hidden anchors
## come FROM the query pool, so without that branch every hidden anchor would
## report distance 0 to itself.
##
##   p_detect : from t3_geometry(); 1 means H = 0 and the observed set is
##              returned unchanged (the identity arm, used as a correctness
##              check against the canonical DEG table).
##
## attr "stats": per-sample n_obs, H, n_candidates, weight, mean/max weight.
## attr "hidden_ids": character vector of hidden cell_ids, for Arm C.
## ---------------------------------------------------------------------------
draw_hidden_anchors <- function(cd, neuron_celltypes, p_detect,
                                weight     = "intensity",
                                sigma_um   = HIDDEN_SIGMA_UM,
                                sample_col = "sample_id",
                                coord_x    = "x_slide_mm",
                                coord_y    = "y_slide_mm") {

  weight <- match.arg(weight, HIDDEN_WEIGHTS)
  stopifnot("cell_id" %in% colnames(cd))
  if (p_detect <= 0 || p_detect > 1)
    stop("draw_hidden_anchors(): p_detect must be in (0, 1].")

  samples <- unique(cd[[sample_col]])
  out  <- vector("list", length(samples)); names(out) <- as.character(samples)
  stat <- vector("list", length(samples))
  hid_ids <- character(0)

  for (si in seq_along(samples)) {
    smp    <- samples[si]
    in_smp <- which(cd[[sample_col]] == smp)
    cd_smp <- cd[in_smp, , drop = FALSE]

    idx_obs  <- phf1_source_idx(cd_smp, neuron_celltypes)
    idx_cand <- which(cd_smp$PHF1 != "TRUE" & cd_smp$celltype %in% neuron_celltypes)
    n_obs    <- length(idx_obs)
    H        <- as.integer(round(n_obs * (1 - p_detect) / p_detect))

    stat[[si]] <- data.frame(
      sample_id = as.character(smp), n_obs = n_obs, H = H,
      n_candidates = length(idx_cand), weight = weight,
      w_mean = NA_real_, w_max = NA_real_,
      stringsAsFactors = FALSE)

    if (n_obs == 0L || H == 0L || length(idx_cand) == 0L) {
      out[[si]] <- idx_obs
      next
    }
    if (length(idx_cand) <= H) {          # degenerate: take the whole pool
      out[[si]] <- c(idx_obs, idx_cand)
      hid_ids <- c(hid_ids, cd_smp$cell_id[idx_cand])
      next
    }

    if (weight == "uniform") {
      w <- rep(1, length(idx_cand))
    } else {
      w <- .anchor_intensity(
        as.matrix(cd_smp[idx_cand, c(coord_x, coord_y)]),
        as.matrix(cd_smp[idx_obs,  c(coord_x, coord_y)]),
        sigma_um)
      w[!is.finite(w) | w < 0] <- 0
      # A donor whose candidates all sit outside 5 sigma of every anchor has no
      # information to weight on. Falling back to uniform is the right default
      # but it must be VISIBLE in the log, not silent -- uniform gives a
      # materially larger shortening (see .anchor_intensity()'s note), so a
      # quiet fallback would inflate the correction for that donor.
      if (sum(w) <= 0) {
        w <- rep(1, length(idx_cand))
        stat[[si]]$weight <- "uniform_fallback"
      }
    }
    stat[[si]]$w_mean <- mean(w); stat[[si]]$w_max <- max(w)

    pick <- sample(idx_cand, H, prob = w)     # without replacement
    out[[si]] <- c(idx_obs, pick)
    hid_ids   <- c(hid_ids, cd_smp$cell_id[pick])
  }

  attr(out, "stats")      <- do.call(rbind, stat)
  attr(out, "hidden_ids") <- hid_ids
  attr(out, "p_detect")   <- p_detect
  attr(out, "weight")     <- weight
  out
}


## ---------------------------------------------------------------------------
## hidden_dist_matrix()
##
## R replicates of the augmented distance field. Rows aligned to `cell_ids` --
## the CALLER's order, not cd's -- because the analysis frame is one celltype
## while distances must be computed against anchors from the whole section.
## Getting that alignment wrong silently scrambles everything, so it is done
## here once and BY NAME.
##
## Returns an n_cells x R matrix with attrs "draw_stats" (rbind of per-replicate
## stats), "hidden_ids" (list of length R, for the paired Arm C), "p_detect",
## "weight".
## ---------------------------------------------------------------------------
hidden_dist_matrix <- function(cd, cell_ids, neuron_celltypes, p_detect,
                               R = 10L, weight = "intensity",
                               sigma_um = HIDDEN_SIGMA_UM,
                               sample_col = "sample_id",
                               coord_x = "x_slide_mm", coord_y = "y_slide_mm",
                               seed = 42L, verbose = TRUE) {

  M <- matrix(NA_real_, nrow = length(cell_ids), ncol = R,
              dimnames = list(cell_ids, NULL))
  st <- vector("list", R); hid <- vector("list", R)

  for (i in seq_len(R)) {
    set.seed(seed + i)
    dr <- draw_hidden_anchors(cd, neuron_celltypes, p_detect, weight,
                              sigma_um, sample_col, coord_x, coord_y)
    d  <- compute_dist_to_phf1_um(cd, sample_col, neuron_celltypes,
                                  coord_x, coord_y,
                                  source_idx_override = dr, verbose = FALSE)
    M[, i]   <- d[match(cell_ids, names(d))]
    s        <- attr(dr, "stats"); s$rep <- i
    st[[i]]  <- s
    hid[[i]] <- attr(dr, "hidden_ids")
    if (verbose && (i %% 5L == 0L || i == R))
      cat(sprintf("  hidden-anchor replicate %d/%d\n", i, R))
  }

  attr(M, "draw_stats") <- do.call(rbind, st)
  attr(M, "hidden_ids") <- hid
  attr(M, "p_detect")   <- p_detect
  attr(M, "weight")     <- weight
  M
}


## ---------------------------------------------------------------------------
## select_intensity_drop()
##
## ARM C. Returns the cells to drop from the comparator arm, by two selections
## that are only meaningful as a PAIR:
##
##   $intensity : top H by PHF1 p95 intensity RESIDUAL WITHIN DISTANCE STRATA.
##                Enriched for genuinely tangle-bearing cells.
##   $placebo   : H cells drawn at MATCHED distances (same strata, same counts)
##                weighted by anchor intensity but ignoring PHF1 intensity.
##                Should move nothing.
##
## THE ADMIXTURE BOUND IS THE DIFFERENCE. The intensity arm alone confounds
## "dropped cells that were really tangle-bearing" with "dropped H near cells",
## and the placebo is what removes the second.
##
## Stratification is on distance deciles by default. Residuals are taken within
## stratum AND within donor, because per-cell PHF1 intensity has a strong donor
## offset (staining, exposure) that would otherwise concentrate the whole drop
## set in one or two brightly-stained donors.
##
##   dat      : data.frame with cell_id, sample_id, dist_col, intensity_col
##   H        : total cells to drop (per cohort, apportioned across strata by
##              stratum size so the distance profile of the drop set matches)
## ---------------------------------------------------------------------------
select_intensity_drop <- function(dat, H, intensity_col = "phf1_intensity_p95",
                                  dist_col = "dist_to_phf1_um",
                                  donor_col = "sample_id",
                                  n_strata = 10L, seed = 42L) {

  for (cl in c("cell_id", intensity_col, dist_col, donor_col))
    if (!cl %in% colnames(dat)) stop("select_intensity_drop(): missing column '", cl, "'")
  if (H <= 0) return(list(intensity = character(0), placebo = character(0),
                          strata = NULL))

  d <- dat[is.finite(dat[[dist_col]]) & is.finite(dat[[intensity_col]]), , drop = FALSE]
  if (nrow(d) < 10L * n_strata)
    stop("select_intensity_drop(): only ", nrow(d), " usable cells for ",
         n_strata, " strata.")

  br <- unique(stats::quantile(d[[dist_col]], seq(0, 1, length.out = n_strata + 1),
                               na.rm = TRUE))
  d$stratum <- cut(d[[dist_col]], breaks = br, include.lowest = TRUE, labels = FALSE)

  # Residualise intensity within (stratum x donor). Median/MAD rather than
  # mean/SD: the tail we are selecting on is exactly what would inflate an SD
  # and shrink its own z-score.
  key <- paste(d$stratum, d[[donor_col]], sep = "|")
  d$resid <- unsplit(lapply(split(d[[intensity_col]], key), function(v) {
    m <- stats::median(v, na.rm = TRUE)
    s <- stats::mad(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) s <- stats::sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(0, length(v)))
    (v - m) / s
  }), key)

  # Apportion H across strata in proportion to stratum size, so the drop set's
  # DISTANCE PROFILE matches the population and the arm is not a near-cell
  # selection wearing an intensity badge.
  tab   <- table(d$stratum)
  quota <- as.integer(round(H * as.numeric(tab) / sum(tab)))
  names(quota) <- names(tab)
  # Rounding can lose or gain a cell or two; fix on the largest stratum.
  if (sum(quota) != H) quota[which.max(tab)] <- quota[which.max(tab)] + (H - sum(quota))

  set.seed(seed)
  int_ids <- pla_ids <- character(0)
  for (s in names(quota)) {
    k <- quota[[s]]
    if (k <= 0) next
    ds <- d[d$stratum == as.integer(s), , drop = FALSE]
    if (nrow(ds) <= k) { int_ids <- c(int_ids, ds$cell_id); next }
    ord <- order(ds$resid, decreasing = TRUE)
    int_ids <- c(int_ids, ds$cell_id[ord[seq_len(k)]])
    # Placebo: same stratum, same count, drawn at random from the cells NOT
    # taken by the intensity arm, so the two sets are disjoint and the
    # difference is attributable to the selection rather than to overlap.
    pool <- ds$cell_id[ord[-seq_len(k)]]
    pla_ids <- c(pla_ids, sample(pool, min(k, length(pool))))
  }

  list(intensity = int_ids, placebo = pla_ids,
       strata = data.frame(stratum = names(quota), n = as.integer(tab),
                           quota = as.integer(quota), row.names = NULL))
}


## ---------------------------------------------------------------------------
## hidden_balance_table()
##
## Standardised mean difference between OBSERVED and HIDDEN anchors on the
## context covariates. Unlike a sham null, |SMD| < 0.1 is NOT the target here:
## hidden anchors are supposed to sit where tangles are, so a systematic
## difference from a uniform draw is the point. What this table is for is
## showing that the INTENSITY weighting moved the draw toward the observed
## anchors relative to the UNIFORM draw -- i.e. that the weighting did
## something. Report both weightings side by side.
## ---------------------------------------------------------------------------
hidden_balance_table <- function(cd, hidden_ids, neuron_celltypes,
                                 match_cols = c("dens100", "edge_dist_um")) {
  match_cols <- match_cols[match_cols %in% colnames(cd)]
  if (!length(match_cols)) return(NULL)
  obs <- cd[cd$PHF1 == "TRUE" & cd$celltype %in% neuron_celltypes, , drop = FALSE]
  hid <- cd[cd$cell_id %in% hidden_ids, , drop = FALSE]
  do.call(rbind, lapply(match_cols, function(k) {
    mo <- mean(obs[[k]], na.rm = TRUE); so <- stats::sd(obs[[k]], na.rm = TRUE)
    mh <- mean(hid[[k]], na.rm = TRUE)
    data.frame(covariate = k, mean_observed = mo, mean_hidden = mh,
               smd = if (is.finite(so) && so > 0) (mh - mo) / so else NA_real_,
               row.names = NULL)
  }))
}
