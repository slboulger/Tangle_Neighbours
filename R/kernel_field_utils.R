#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# kernel_field_utils.R
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
# kernel_field_utils.R
#
# Shared helper: single source of truth for the GAUSSIAN KERNEL FIELD used to
# summarise the glial activation state surrounding a cell.
#
# Deliberately mirrors the design of R/phf1_distance_utils.r and R/nn3_utils.r
# -- the metric lives in ONE place so it cannot drift
# between the observed analysis and any permutation null built on top of it.
#
# ---------------------------------------------------------------------------
# THE TWO QUANTITIES, AND WHY THEY MUST NOT BE COLLAPSED
# ---------------------------------------------------------------------------
# For a target cell i and each glial source cell j of a given celltype in the
# SAME sample, with w_ij = exp(-d_ij^2 / (2 * sigma^2)):
#
#   1. kernel-weighted MEAN activation   sum_j(w_ij * score_j) / sum_j(w_ij)
#      "how activated is the glia around this cell", independent of how much
#      glia there is. NA when no source falls inside the truncation.
#
#   2. kernel DENSITY of that population sum_j(w_ij)
#      "how much glia is around this cell", independent of its state. 0 (not
#      NA) when no source falls inside the truncation -- an empty neighbourhood
#      is a real zero, not a missing observation.
#
# These are biologically different questions. A mean that rises because the
# glia are more activated is not the same finding as a mean that rises because
# there are more glia, and the only way to tell them apart is to carry both.
#
# ---------------------------------------------------------------------------
# TRUNCATION: 5 SIGMA, NOT 3
# ---------------------------------------------------------------------------
# TRUNC_MULT defaults to 5. At 3 sigma the kernel leaves ~21% of neurons with no
# microglial source inside the truncation, and those uncovered neurons are
# systematically FURTHER from a tangle, so dropping them would be a selection on
# the exposure. At 5 sigma the furthest contributing weight is
# exp(-250^2 / (2*50^2)) = 4e-6, i.e. numerically irrelevant but not selective.
#
# Callers may pass trunc_mult = 3, but kernel_coverage_audit() should be written
# out alongside so the cost of that choice is visible rather than assumed.
#
# ---------------------------------------------------------------------------
# CONVENTIONS
# ---------------------------------------------------------------------------
# Coordinates are x_slide_mm / y_slide_mm (slide-global centroids, MM). They are
# slide-global and NOT sample-global -- two samples can share a slide -- so every
# computation here is done WITHIN one sample by the caller, exactly as
# compute_dist_to_phf1_um() and compute_nn3_um() do. Distances are converted to
# MICRONS (mm * 1000) internally, the same mm*1000 convention as
# R/phf1_distance_utils.r. sigma is always in microns.
#
# TWO IMPLEMENTATIONS, PINNED TOGETHER
#   kernel_lag()        dense chunked outer-product, kept as the
#                       reference implementation.
#   kernel_field_rann() RANN radius search. Same maths, but it never
#                       materialises the full n_query x n_source matrix and it
#                       does all modules in one search, so it is the one to use
#                       for a multi-module x multi-sigma sweep.
# The two are tested to agree to 1e-8. If you change one, change the other.

suppressPackageStartupMessages({
  library(cli)
})

if (!requireNamespace("RANN", quietly = TRUE)) {
  stop("Package 'RANN' is required. Install with: install.packages('RANN')")
}

## ---------------------------------------------------------------------------
## Default truncation multiplier. See the header. Exported so callers report the
## value they used rather than hard-coding a second copy of it.
## ---------------------------------------------------------------------------
KERNEL_TRUNC_MULT_DEFAULT <- 5

## ---------------------------------------------------------------------------
## kernel_lag()  --  REFERENCE IMPLEMENTATION
##
## Gaussian kernel weighted mean of `gscore` from source points (gx, gy) onto
## query points (nx, ny). Coordinates are in mm and converted to um here.
## Returns the weighted mean (NA where no source falls inside the truncation) and
## the kernel mass sum, which is the local density readout.
## Chunked over query points so the dense matrix stays bounded regardless of how
## many source cells the sample has.
## ---------------------------------------------------------------------------
kernel_lag <- function(nx, ny, gx, gy, gscore, h_um, trunc_um) {
  n <- length(nx); m <- length(gx)
  out_mean <- rep(NA_real_, n); out_dens <- rep(0, n)
  if (m == 0L || n == 0L) return(list(mean = out_mean, dens = out_dens))
  h2      <- 2 * h_um^2
  trunc2  <- trunc_um^2
  chunk   <- max(200L, as.integer(4e6 / m))          # ~4M cells per dense block
  for (s in seq(1L, n, by = chunk)) {
    ii <- s:min(s + chunk - 1L, n)
    dx <- outer(nx[ii], gx, "-") * 1000
    dy <- outer(ny[ii], gy, "-") * 1000
    d2 <- dx^2 + dy^2
    rm(dx, dy)
    w  <- exp(-d2 / h2)
    w[d2 > trunc2] <- 0
    rm(d2)
    sw <- rowSums(w)
    out_dens[ii] <- sw
    ok <- sw > 0
    if (any(ok)) out_mean[ii][ok] <- as.numeric(w[ok, , drop = FALSE] %*% gscore) / sw[ok]
    rm(w)
  }
  list(mean = out_mean, dens = out_dens)
}

## ---------------------------------------------------------------------------
## kernel_field_rann()  --  the one to actually call
##
## Same kernel as kernel_lag(), evaluated through a RANN radius search so that
## only the neighbours inside trunc_mult * sigma are ever touched.
##
##   query_xy_mm  : n_query x 2 numeric matrix (mm). The target cells.
##   source_xy_mm : n_source x 2 numeric matrix (mm). The glial source cells.
##   score_mat    : n_source x M numeric matrix of source scores, one column per
##                  module. Column names become the names of the returned means.
##                  Pass NULL (or a 0-column matrix) for a density-only call.
##   sigma_um     : kernel sigma, MICRONS.
##   trunc_mult   : kernel truncated at trunc_mult * sigma. Default 5, see header.
##   k_probe      : initial neighbour-list depth. NOT a cap -- see below.
##
## ALL MODULES IN ONE SEARCH. The weights w_ij depend only on geometry, so they
## are computed once per (sample, source celltype, sigma) and multiplied through
## every module column at once. An M-module x S-sigma sweep therefore costs one
## neighbour search per (sample, sigma), not M x S of them.
##
## NO SILENT TRUNCATION. RANN's searchtype = "radius" returns at most k
## neighbours even when more lie inside the radius, padding the rest with
## nn.idx == 0 and nn.dists == 1.340781e+154. A quietly truncated neighbour list
## would bias every weight downward, so if the k-th slot is still a real
## neighbour the depth is DOUBLED and the search repeated, up to n_source, at
## which point truncation is impossible by construction.
##
## NA SOURCE SCORES are not tolerated: a source cell whose score is NA would have
## to be dropped from the denominator for that module but kept for the others,
## which makes the density and the mean rest on different neighbour sets. Filter
## upstream; this function raises instead of guessing.
##
## Returns list(mean = n_query x M matrix (NA where the neighbourhood is empty),
##              dens = length-n_query numeric (0 where the neighbourhood is
##                     empty -- an empty neighbourhood is a real zero),
##              k_used = final probe depth,
##              n_uncovered = number of query cells with an empty neighbourhood)
## ---------------------------------------------------------------------------
kernel_field_rann <- function(query_xy_mm, source_xy_mm, score_mat, sigma_um,
                              trunc_mult = KERNEL_TRUNC_MULT_DEFAULT,
                              k_probe = 256L, verbose = FALSE) {

  query_xy_mm  <- as.matrix(query_xy_mm)
  n_query      <- nrow(query_xy_mm)
  source_xy_mm <- if (is.null(source_xy_mm)) matrix(numeric(0), 0, 2) else as.matrix(source_xy_mm)
  n_source     <- nrow(source_xy_mm)

  if (is.null(score_mat)) score_mat <- matrix(numeric(0), nrow = n_source, ncol = 0)
  score_mat <- as.matrix(score_mat)
  M         <- ncol(score_mat)
  mod_names <- if (M > 0) colnames(score_mat) else character(0)

  stopifnot(
    "sigma_um must be a single positive number"   = length(sigma_um) == 1L && sigma_um > 0,
    "trunc_mult must be a single positive number" = length(trunc_mult) == 1L && trunc_mult > 0,
    "score_mat must have one row per source cell" = nrow(score_mat) == n_source
  )
  if (M > 0 && any(!is.finite(score_mat))) {
    stop("kernel_field_rann(): non-finite source scores. Drop or impute them before ",
         "calling -- the weighted mean and the density must rest on the same neighbour set.")
  }

  out_mean <- matrix(NA_real_, nrow = n_query, ncol = M, dimnames = list(NULL, mod_names))
  out_dens <- rep(0, n_query)
  if (n_query == 0L || n_source == 0L) {
    return(list(mean = out_mean, dens = out_dens, k_used = 0L, n_uncovered = n_query))
  }

  radius_mm <- trunc_mult * sigma_um / 1000        # mm: same units as the coordinates
  two_s2    <- 2 * sigma_um^2

  ## Probe depth. Doubled until the k-th slot is empty for every query cell, i.e.
  ## until no neighbour list can have been cut short. Bounded by n_source.
  k <- min(as.integer(k_probe), n_source)
  repeat {
    nn <- RANN::nn2(data = source_xy_mm, query = query_xy_mm, k = k,
                    searchtype = "radius", radius = radius_mm)
    if (k >= n_source) break
    if (!any(nn$nn.idx[, k] > 0L)) break
    k_new <- min(k * 2L, n_source)
    if (verbose)
      cli_alert_info(paste0("kernel_field_rann(): probe depth {k} was reached inside the ",
                            "radius; re-probing at {k_new}."))
    k <- k_new
  }

  idx   <- nn$nn.idx                                # n_query x k, 0 = no neighbour
  found <- idx > 0L
  ## RANN pads unfilled slots with 1.340781e+154, which would overflow when
  ## squared. Never let a padded distance reach the exponential.
  d_um  <- nn$nn.dists * 1000
  d_um[!found] <- 0

  w <- exp(-(d_um^2) / two_s2)
  w[!found] <- 0
  rm(d_um)

  out_dens <- rowSums(w)
  covered  <- out_dens > 0

  if (M > 0 && any(covered)) {
    idx_safe <- idx
    idx_safe[!found] <- 1L                          # any valid index; weight is 0
    for (m in seq_len(M)) {
      vals <- matrix(score_mat[idx_safe, m], nrow = n_query)
      num  <- rowSums(w * vals)
      out_mean[covered, m] <- num[covered] / out_dens[covered]
    }
  }

  list(mean = out_mean, dens = out_dens, k_used = k, n_uncovered = sum(!covered))
}

## ---------------------------------------------------------------------------
## kernel_coverage_audit()
##
## The coverage check for the truncation radius. A kernel that leaves cells with
## an empty neighbourhood has silently dropped them from the analysis, and the
## question that matters is not HOW MANY were dropped but WHETHER THE DROPPED
## ONES DIFFER -- specifically, whether they sit further from a tangle, because
## then the exclusion is a selection on the exposure and the resulting estimate
## is biased no matter how small the percentage looks.
##
##   covered  : logical, one per query cell (dens > 0)
##   dist_um  : dist_to_phf1_um for the same cells
##
## Returns a one-row data.frame. Warns when >1% are uncovered AND the covered /
## uncovered distance distributions differ, which is the combination that
## actually matters -- either alone is tolerable.
## ---------------------------------------------------------------------------
kernel_coverage_audit <- function(covered, dist_um, sigma_um, trunc_mult,
                                  source_celltype = NA_character_,
                                  sample_id = NA_character_) {

  covered <- as.logical(covered)
  n_tot   <- length(covered)
  n_unc   <- sum(!covered)
  pct_unc <- if (n_tot > 0) 100 * n_unc / n_tot else NA_real_

  d_cov <- dist_um[covered  & is.finite(dist_um)]
  d_unc <- dist_um[!covered & is.finite(dist_um)]

  wp <- if (length(d_cov) > 1L && length(d_unc) > 1L) {
    suppressWarnings(stats::wilcox.test(d_cov, d_unc)$p.value)
  } else NA_real_

  out <- data.frame(
    sample_id           = as.character(sample_id),
    source_celltype     = as.character(source_celltype),
    sigma_um            = sigma_um,
    trunc_mult          = trunc_mult,
    trunc_um            = sigma_um * trunc_mult,
    n_total             = n_tot,
    n_uncovered         = n_unc,
    pct_uncovered       = pct_unc,
    median_dist_covered   = if (length(d_cov)) stats::median(d_cov) else NA_real_,
    median_dist_uncovered = if (length(d_unc)) stats::median(d_unc) else NA_real_,
    wilcox_p            = wp,
    stringsAsFactors = FALSE
  )

  if (isTRUE(pct_unc > 1) && isTRUE(wp < 0.05)) {
    cli_alert_warning(paste0(
      "Kernel coverage: {round(pct_unc, 1)}% of cells have no {source_celltype} source ",
      "inside {round(sigma_um * trunc_mult)} um, and the uncovered cells differ in ",
      "distance-to-tangle (median {round(out$median_dist_uncovered)} vs ",
      "{round(out$median_dist_covered)} um, Wilcoxon p = {signif(wp, 3)}). ",
      "That is a selection on the exposure -- raise --trunc_mult."))
  }
  out
}

## ---------------------------------------------------------------------------
## within_donor_z()
##
## Z-score a per-cell value WITHIN each donor.
##
## This is not a stylistic choice. Between-donor differences in level ride on top
## of between-donor differences in how much tissue sits at each distance, so a
## field built from un-normalised scores would carry donor identity (see the
## within-donor z-scoring rationale in R/imc_phf1_glia_distance.R).
##
## Donors whose values have zero or non-finite SD return 0 (a constant carries no
## within-donor information) rather than NaN, so one degenerate donor cannot
## delete an entire sample from the fit.
## ---------------------------------------------------------------------------
within_donor_z <- function(x, donor) {
  x <- as.numeric(x)
  donor <- as.character(donor)
  out <- rep(NA_real_, length(x))
  for (d in unique(donor)) {
    ii <- which(donor == d)
    v  <- x[ii]
    ok <- is.finite(v)
    if (!any(ok)) next
    mu <- mean(v[ok]); sdv <- stats::sd(v[ok])
    out[ii] <- if (!is.finite(sdv) || sdv == 0) 0 else (v - mu) / sdv
  }
  out
}

## ---------------------------------------------------------------------------
## resolve_moderator_genes()
##
## One definition of "which genes are in the moderator module", shared by the
## script that BUILDS the score and the DEG script that needs the membership list
## to decide which genes get a leave-one-out moderator.
##
## If these two ever disagreed, member genes would be tested against a moderator
## that still contains them while the output column claimed otherwise -- a silent
## error with no symptom. Hence one function, called by both.
##
## Three accepted forms:
##   directory  per-celltype gene sets, <dir>/<ct_safe>/<ct_safe>_phf1_geneset.txt.
##              THE DEFAULT. This is the PHF1 module from
##              find_phf1_markers_by_celltype.R, the same one scored by
##              plot_phf1_module_vs_phf1_distance_modelp.R (its read_geneset()).
##              Those genes are ALREADY the up direction --
##              "upregulated in PHF1+ vs PHF1- within that celltype" -- so there
##              is no logFC filter to apply here and applying one would be wrong.
##              The signature is DIFFERENT FOR EVERY CELLTYPE, which is the point
##              of it, so `celltype` is required.
##   .rds       a named list of sets; `set` picks one (e.g. otero_L23_up). These
##              are GLOBAL: the same genes in every celltype.
##   .txt       one gene symbol per line.
##
## The celltype is sanitised with the SAME gsub the plotting script uses, so
## "Unassigned Neuron" resolves to the Unassigned_Neuron directory that exists on
## disk rather than to a path with a space in it.
##
## Returns a character vector; raises if the celltype has no file. The plotting
## script returns NULL there because a missing signature just drops a panel; here
## it would silently become an empty moderator and the whole model would be
## meaningless, so it is an error instead.
## ---------------------------------------------------------------------------
resolve_moderator_genes <- function(path, set = NA_character_, celltype = NA_character_) {

  if (dir.exists(path)) {
    if (is.na(celltype))
      stop("resolve_moderator_genes(): '", path, "' is a directory of per-celltype ",
           "gene sets, so a celltype must be supplied.")
    ## Same sanitisation as plot_phf1_module_vs_phf1_distance_modelp.R.
    ct_safe <- gsub("[^A-Za-z0-9_-]", "_", celltype)
    f <- file.path(path, ct_safe, paste0(ct_safe, "_phf1_geneset.txt"))
    if (!file.exists(f)) {
      alt <- file.path(path, paste0(ct_safe, "_phf1_geneset.txt"))
      if (file.exists(alt)) f <- alt else
        stop("No PHF1 gene set for celltype '", celltype, "' under ", path,
             "\nExpected: ", f,
             "\nThe PHF1 module only exists for celltypes eligible for PHF1 ",
             "analysis (enough PHF1+ cells across enough donors); there is ",
             "nothing to moderate on in the others.")
    }
    g <- readLines(f, warn = FALSE)

  } else if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    gl <- readRDS(path)
    if (is.na(set) || !set %in% names(gl))
      stop("Set '", set, "' not found in ", path,
           ". Available: ", paste(names(gl), collapse = ", "))
    g <- as.character(gl[[set]])

  } else {
    if (!file.exists(path)) stop("Moderator gene set not found: ", path)
    g <- readLines(path, warn = FALSE)
  }

  unique(trimws(g[nzchar(trimws(g))]))
}

## ---------------------------------------------------------------------------
## moderator_pipeline()
##
## The single definition of how a raw module score becomes the moderator A_resid.
## It lives here, not in the DEG script, because the leave-one-out moderator and
## the standard moderator MUST go through byte-identical processing -- if they
## drift, the member genes silently become a different analysis from the rest of
## the panel, and nothing downstream would reveal it.
##
##   S        numeric vector, or an n x m matrix of columns to process
##            IDENTICALLY (one per leave-one-out gene)
##   qr_tech  QR of the technical design (~ nUMI_log + percent_neg)
##   qr_dist  QR of technical + log(distance); NULL to skip distance adjustment
##   donor    donor label per cell, for the within-donor z
##
## Residualisation is a JOINT projection onto technical + distance, not the
## technical residual subsequently regressed on distance. Sequential projections
## onto correlated covariates do not equal one joint projection, and "adjusted
## for both" means the joint one.
##
## Returns list(z = n x m matrix, r2_tech = per column, r2_dist = per column).
## r2_dist is the INCREMENTAL fraction of technical-adjusted variance that
## distance removes -- the number that says how much of a putative tau signature
## was really proximity to a tangle.
## ---------------------------------------------------------------------------
moderator_pipeline <- function(S, qr_tech, qr_dist = NULL, donor) {
  S <- as.matrix(S)
  tss <- colSums(sweep(S, 2, colMeans(S))^2)

  resid_tech <- qr.resid(qr_tech, S)
  r2_tech    <- 1 - colSums(resid_tech^2) / tss

  if (!is.null(qr_dist)) {
    S2      <- qr.resid(qr_dist, S)
    r2_dist <- 1 - colSums(S2^2) / colSums(resid_tech^2)
  } else {
    S2      <- resid_tech
    r2_dist <- rep(NA_real_, ncol(S))
  }

  Z <- apply(S2, 2, function(v) within_donor_z(v, as.character(donor)))
  Z <- matrix(as.numeric(Z), nrow = nrow(S), ncol = ncol(S),
              dimnames = list(NULL, colnames(S)))

  list(z = Z, r2_tech = r2_tech, r2_dist = r2_dist)
}

## ---------------------------------------------------------------------------
## loo_score_matrix()
##
## The leave-one-out module score for every member gene at once:
##
##     S_i_LOO(g) = (k * M_i - x_ig) / (k - 1) - C_i
##
## which is exactly the module score rebuilt from the OTHER k - 1 genes, against
## the same control mean. Vectorised over genes -- AddModuleScore is never called
## again, and the control set is deliberately NOT recomputed: the control mean is
## an expression-matched background, not a member of the module, so leaving out a
## module gene does not change it.
##
##   M  length-n raw module mean
##   C  length-n control mean
##   X  n x m matrix of member-gene expression on the SAME scale the score was
##      built on (SCT data slot)
##   k  full module size (>= 2), NOT ncol(X) -- the module may contain genes that
##      are not on the tested panel, and the algebra needs the true k.
## ---------------------------------------------------------------------------
loo_score_matrix <- function(M, C, X, k) {
  X <- as.matrix(X)
  stopifnot(
    "k must be at least 2 for leave-one-out" = k >= 2,
    "M must have one value per cell" = length(M) == nrow(X),
    "C must have one value per cell" = length(C) == nrow(X)
  )
  (k * M - X) / (k - 1) - C
}

## ---------------------------------------------------------------------------
## vif_from_design()
##
## Variance inflation factors from a fixed-effects design matrix, so collinearity
## between the glial field, the moderator, distance and the interaction can be
## inspected without re-fitting the mixed model. VIF is a property of the design,
## so dropping the random intercept does not change it materially.
## Returns Inf rather than a divide-by-zero for an exactly aliased term.
## ---------------------------------------------------------------------------
vif_from_design <- function(X) {
  X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
  vapply(colnames(X), function(v) {
    y  <- X[, v]
    Xo <- X[, setdiff(colnames(X), v), drop = FALSE]
    if (ncol(Xo) == 0) return(1)
    r2 <- summary(stats::lm(y ~ Xo))$r.squared
    if (!is.finite(r2) || r2 >= 1) return(Inf)
    1 / (1 - r2)
  }, numeric(1))
}

## ---------------------------------------------------------------------------
## kernel_col_names()
##
## Single definition of the output column naming, so script 1 (which writes the
## columns) and scripts 2/3 (which read them) can never disagree about the
## spelling. sigma is formatted without a decimal point where it is a whole
## number, so 50 -> "s50" rather than "s50.0".
## ---------------------------------------------------------------------------
.fmt_sigma <- function(sigma_um) {
  gsub("\\.", "p", format(sigma_um, trim = TRUE, scientific = FALSE))
}

kernel_mean_col <- function(module, sigma_um) {
  sprintf("kern_%s_mean_s%s", gsub("[^A-Za-z0-9]", "_", module), .fmt_sigma(sigma_um))
}

kernel_dens_col <- function(celltype, sigma_um) {
  sprintf("kern_%s_dens_s%s", gsub("[^A-Za-z0-9]", "_", celltype), .fmt_sigma(sigma_um))
}
