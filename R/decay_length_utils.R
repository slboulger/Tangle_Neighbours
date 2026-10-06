# ---------------------------------------------------------------------------
# decay_length_utils.R
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
# decay_length_utils.R
#
# Length constant (lambda, um) of a module score's decay away from the nearest
# PHF1+ neuron. Sourced by plot_phf1_module_vs_phf1_distance_modelp.R and
# plot_modulescore_vs_phf1_distance_modelp.R, which call it as the "decay" model
# type alongside linear / log / gam.
#
# THE MODEL. For donor j, cell i at distance d_i (um):
#
#   S_ij = P_j + A_j * exp(-d_i / lambda_j) + b1*nUMI_log_i + b2*percent_neg_i + e_ij
#
#   P_j      donor plateau (the background level far from any tangle)
#   A_j      SIGNED amplitude at d = 0 (the excess at the tangle)
#   lambda_j length constant: the distance over which the excess falls by 1/e
#
# A is signed ON PURPOSE. Some modules rise away from the tangle (e.g. Cameron
# Neuroprotective3 has a positive log-distance slope); the identical functional
# form describes that with A < 0 and the same lambda. Never constrain A > 0.
# A single exponential is monotone by construction, so no monotone-spline
# machinery (scam) is needed to keep the fit interpretable.
#
# Sex / Age / PMI are DROPPED from the per-donor fits: they are donor-level and
# exactly aliased with P_j. nUMI_log and percent_neg are retained (they vary
# within donor and are part of the canonical adjustment in the sibling models).
#
# WHY PROFILE LEAST SQUARES AND NOT nls(). The model is LINEAR in P, A and the
# covariate coefficients once lambda is fixed. So lambda is found by 1-D profile
# least squares -- for each candidate lambda fit a plain linear model and record
# RSS, then minimise over lambda. This cannot fail to converge, which is what
# matters when the estimator is run separately inside each of 9 donors, and it
# yields a profile-likelihood CI directly:
#
#   { lambda : RSS(lambda) <= RSS(lambda_hat) * (1 + F_{1,n-p}(0.95) / (n-p)) }
#
# decay_nls_check() runs nls(y ~ SSasymp(...)) as a cross-check only.
#
# DISTANCE IS RAW MICRONS HERE -- a deliberate, documented departure from the
# canonical log(dist)/sd(log(dist)) transform (see docs/MODELS.md). lambda has units
# of microns and is undefined on a log/SD-rescaled axis; that is the whole point
# of the analysis. The existing "gam" model type is the precedent (it also
# smooths raw distance). The canonical transform is unchanged everywhere else.
#
# HEADLINE ESTIMATOR: per-donor fits, then a random-effects (DerSimonian-Laird)
# meta-analysis of log lambda with a Knapp-Hartung t-adjusted CI -- with k <= 9
# donors the unadjusted DL interval is far too narrow. fit_decay_shared() gives
# a one-shared-lambda / donor-specific-P-and-A anchor as the consistency check.
# The meta-analysis is hand-rolled (no metafor dependency).
#
# GATES. A single exponential returns a number even when the truth is log-linear,
# in which case lambda is determined by where the distance cap sits.
# decay_analyse() therefore always computes, and the callers always print:
#   1. model comparison  -- decay must beat log-linear on AIC (the log model is
#      the primary model; if it also wins on fit, the reading is
#      "there is a gradient but no identifiable length constant")
#   2. cap sensitivity   -- lambda refitted at half the cap; a lambda that tracks
#      the cap is not a biological length constant
#   3. identifiability   -- lambda or its upper CI beyond cap/2 means the curve
#      has not plateaued inside the observation window
# Their conjunction sets `quotable`. Nothing here decides significance: that
# stays with the canonical log model, which also GATES which modules get fitted.
#
# No side effects: constants and functions only.

DECAY_LAMBDA_MIN_UM   <- 5      # below ~a cell diameter lambda is not resolvable
DECAY_LAMBDA_MAX_MULT <- 3      # search up to MULT * distance cap
DECAY_N_GRID          <- 200    # log-spaced profile grid points
DECAY_MIN_CELLS_DONOR <- 100    # min cells for a donor to contribute a lambda
DECAY_MIN_SPAN_UM     <- 200    # donor needs a plateau region to identify lambda
DECAY_CONF            <- 0.95
DECAY_CAP_RATIO_LO    <- 0.67   # lambda(half cap) / lambda(full cap) tolerance
DECAY_CAP_RATIO_HI    <- 1.50
DECAY_MIN_DAIC        <- 10     # decay must beat log-linear by this many AIC units
DECAY_MAX_CI_RATIO    <- 5      # lambda_hi / lambda_lo above this is not a usable number
# The last two exist because, in simulation on ~10^4 cells with a realistically
# right-skewed distance distribution, a purely LOG-LINEAR (plateau-free) truth
# is fitted by an exponential that wins on AIC by ~3 units and survives the
# cap-sensitivity check -- but its profile CI spans a 27-fold range. So `quotable`
# requires a MARGIN over log-linear and a CI narrow enough to place lambda.
#
# The 5-fold CI limit is set by what the number has to be able to DO, not by what
# the data happen to give: lambda is worth quoting if it can be placed in one of
# the biologically distinct spatial bins -- juxtacellular (~10-30 um), local
# microenvironment (~50-150 um), columnar/laminar (~300-500 um), region-wide
# (>1000 um). Discriminating adjacent bins needs a CI inside roughly 5-fold.
# A CI of, say, [37, 127] um (3.5-fold) still says "tens, not hundreds";
# [15, 400] um (27-fold) does not discriminate the bins.

# Gaussian AIC from an RSS, up to an additive constant that is IDENTICAL across
# models compared on the same rows. k = estimated parameters + 1 for sigma.
# Vectorised over rss / k so a whole comparison table can be filled in one call.
.gauss_aic <- function(rss, n, k) {
  ifelse(is.finite(rss) & rss > 0, n * log(rss / n) + 2 * k, NA_real_)
}

# Profile least squares for lambda given a design-matrix builder.
#   build_X(lambda) -> model matrix whose columns are linear in the parameters
#   n_par_extra     -> parameters estimated outside build_X (just lambda: 1)
# Returns lambda_hat, the profile CI, and RSS / AIC at the optimum.
.profile_lambda <- function(build_X, y, lambda_min, lambda_max,
                           n_grid = DECAY_N_GRID, conf = DECAY_CONF) {
  n <- length(y)
  rss_at <- function(lam) {
    X <- build_X(lam)
    if (anyNA(X) || !all(is.finite(X))) return(Inf)
    f <- .lm.fit(X, y)
    # A rank drop means this lambda makes the exponential term collinear with the
    # intercept (or a donor dummy) -- i.e. not identified. Inf RSS pushes the
    # profile away from it and, if the whole upper range is like this, leaves the
    # CI open at the bound, which is exactly the flag we want.
    if (f$rank < ncol(X)) return(Inf)
    sum(f$residuals^2)
  }

  grid  <- exp(seq(log(lambda_min), log(lambda_max), length.out = n_grid))
  rss_g <- vapply(grid, rss_at, numeric(1))
  if (!any(is.finite(rss_g)))
    return(list(lambda = NA_real_, lambda_lo = NA_real_, lambda_hi = NA_real_,
                rss = NA_real_, n = n, n_par = NA_integer_, aic = NA_real_,
                converged = FALSE, at_bound_lo = NA, at_bound_hi = NA,
                ci_open_lo = NA, ci_open_hi = NA))

  i0 <- which.min(rss_g)
  # Refine inside the bracketing grid cells (RSS is smooth in log lambda here).
  lo_b <- grid[max(1, i0 - 1)]; hi_b <- grid[min(n_grid, i0 + 1)]
  opt  <- stats::optimize(function(l) rss_at(exp(l)), lower = log(lo_b), upper = log(hi_b))
  lam_hat <- exp(opt$minimum); rss_hat <- opt$objective
  if (!is.finite(rss_hat) || rss_hat > rss_g[i0]) { lam_hat <- grid[i0]; rss_hat <- rss_g[i0] }

  n_par <- ncol(build_X(lam_hat)) + 1L                  # linear params + lambda
  df_res <- n - n_par
  if (df_res <= 1)
    return(list(lambda = lam_hat, lambda_lo = NA_real_, lambda_hi = NA_real_,
                rss = rss_hat, n = n, n_par = n_par, aic = NA_real_,
                converged = FALSE, at_bound_lo = NA, at_bound_hi = NA,
                ci_open_lo = NA, ci_open_hi = NA))

  thresh <- rss_hat * (1 + stats::qf(conf, 1, df_res) / df_res)
  gfun   <- function(l) rss_at(exp(l)) - thresh

  # Walk the grid outwards from the optimum to the first crossing, then uniroot.
  cross <- function(idx_seq) {
    prev <- NA_real_
    for (i in idx_seq) {
      if (is.finite(rss_g[i]) && rss_g[i] > thresh) return(c(grid[i], prev))
      if (!is.finite(rss_g[i])) return(c(grid[i], prev))   # unidentified beyond here
      prev <- grid[i]
    }
    c(NA_real_, prev)
  }
  up <- cross(if (i0 < n_grid) seq(i0 + 1L, n_grid) else integer(0))
  dn <- cross(if (i0 > 1L)     seq(i0 - 1L, 1L)     else integer(0))

  # Bracket each root between lambda_hat (where RSS - thresh is negative by
  # construction) and the first grid point beyond the threshold. Bracketing on the
  # last sub-threshold grid point instead would be NA when the crossing lies in the
  # grid cell immediately next to the optimum, and the resulting uniroot failure
  # would wrongly report the CI as open at the bound -- which silently drops the
  # donor from the meta-analysis.
  lam_hi <- lambda_max; ci_open_hi <- TRUE
  if (is.finite(up[1]) && up[1] > lam_hat) {
    r <- try(stats::uniroot(gfun, lower = log(lam_hat), upper = log(up[1])), silent = TRUE)
    if (!inherits(r, "try-error")) { lam_hi <- exp(r$root); ci_open_hi <- FALSE }
  }
  lam_lo <- lambda_min; ci_open_lo <- TRUE
  if (is.finite(dn[1]) && dn[1] < lam_hat) {
    r <- try(stats::uniroot(gfun, lower = log(dn[1]), upper = log(lam_hat)), silent = TRUE)
    if (!inherits(r, "try-error")) { lam_lo <- exp(r$root); ci_open_lo <- FALSE }
  }

  at_lo <- lam_hat <= lambda_min * 1.01
  at_hi <- lam_hat >= lambda_max * 0.99
  list(lambda = lam_hat, lambda_lo = lam_lo, lambda_hi = lam_hi,
       rss = rss_hat, n = n, n_par = n_par,
       aic = .gauss_aic(rss_hat, n, n_par + 1L),
       converged = !at_lo && !at_hi && !ci_open_lo && !ci_open_hi,
       at_bound_lo = at_lo, at_bound_hi = at_hi,
       ci_open_lo = ci_open_lo, ci_open_hi = ci_open_hi)
}

# Single-population decay fit: one plateau, one amplitude, one lambda.
#   d  distance in um (>= 0)
#   y  module score
#   Z  optional numeric matrix / data.frame of within-donor covariates
# Plateau / amplitude estimates and their SEs are CONDITIONAL on lambda_hat
# (they ignore the uncertainty in lambda) -- stated wherever they are reported.
fit_decay_profile <- function(d, y, Z = NULL,
                              lambda_min = DECAY_LAMBDA_MIN_UM, lambda_max = 3000,
                              n_grid = DECAY_N_GRID, conf = DECAY_CONF,
                              center = TRUE) {
  stopifnot(length(d) == length(y))
  if (!is.null(Z)) Z <- as.matrix(Z)
  ok <- is.finite(d) & is.finite(y) & (if (is.null(Z)) TRUE else stats::complete.cases(Z))
  d <- d[ok]; y <- y[ok]; if (!is.null(Z)) Z <- Z[ok, , drop = FALSE]
  if (any(d < 0)) stop("Negative distance passed to fit_decay_profile().")
  # Covariates are CENTRED so the intercept is the background score at TYPICAL
  # nUMI_log / percent_neg -- i.e. an actual plateau level, comparable across
  # donors. Uncentred, the intercept is the score extrapolated to nUMI_log = 0,
  # which is meaningless, varies wildly between donors, and makes any
  # plateau-relative quantity nonsense. Centring shifts only the intercept: it
  # leaves lambda and the amplitude untouched (the same argument that justifies
  # dropping scale()'s centring from the distance transform).
  # decay_analyse() centres on COHORT means instead and passes center = FALSE, so
  # every donor's plateau is referenced to the same covariate values.
  if (center && !is.null(Z)) Z <- sweep(Z, 2, colMeans(Z), "-")

  build_X <- function(lam) {
    X <- cbind(`(Intercept)` = 1, expterm = exp(-d / lam))
    if (!is.null(Z)) X <- cbind(X, Z)
    X
  }
  pr <- .profile_lambda(build_X, y, lambda_min, lambda_max, n_grid, conf)
  out <- c(pr, list(plateau = NA_real_, amplitude = NA_real_, amp_se = NA_real_,
                    amp_ci_lo = NA_real_, amp_ci_hi = NA_real_, amp_p = NA_real_,
                    dist_span_um = if (length(d)) diff(range(d)) else NA_real_))
  if (!is.finite(pr$lambda)) return(out)

  dd <- data.frame(y = y, expterm = exp(-d / pr$lambda))
  if (!is.null(Z)) dd <- cbind(dd, as.data.frame(Z))
  m  <- stats::lm(y ~ ., data = dd)
  cf <- stats::coef(summary(m))
  if ("expterm" %in% rownames(cf)) {
    ci <- try(stats::confint(m, "expterm", level = conf), silent = TRUE)
    out$plateau   <- unname(stats::coef(m)[["(Intercept)"]])
    out$amplitude <- unname(cf["expterm", "Estimate"])
    out$amp_se    <- unname(cf["expterm", "Std. Error"])
    out$amp_p     <- unname(cf["expterm", "Pr(>|t|)"])
    if (!inherits(ci, "try-error")) { out$amp_ci_lo <- ci[1, 1]; out$amp_ci_hi <- ci[1, 2] }
  }
  out
}

# One shared lambda across donors, donor-specific plateau and amplitude:
#   score ~ 0 + donor + donor:exp(-d/lambda) + covariates
# Still linear given lambda, so the same profile machinery applies. This is the
# consistency anchor for the per-donor meta-analysis: it uses every cell and
# cannot drop a donor for having a flat gradient.
fit_decay_shared <- function(d, y, donor, Z = NULL,
                             lambda_min = DECAY_LAMBDA_MIN_UM, lambda_max = 3000,
                             n_grid = DECAY_N_GRID, conf = DECAY_CONF,
                             center = TRUE) {
  donor <- droplevels(factor(donor))
  if (!is.null(Z)) Z <- as.matrix(Z)
  ok <- is.finite(d) & is.finite(y) & !is.na(donor) &
        (if (is.null(Z)) TRUE else stats::complete.cases(Z))
  d <- d[ok]; y <- y[ok]; donor <- droplevels(donor[ok])
  if (!is.null(Z)) Z <- Z[ok, , drop = FALSE]
  if (center && !is.null(Z)) Z <- sweep(Z, 2, colMeans(Z), "-")   # see fit_decay_profile()
  D <- stats::model.matrix(~ 0 + donor)

  build_X <- function(lam) {
    X <- cbind(D, D * exp(-d / lam))
    colnames(X) <- c(paste0("P_", levels(donor)), paste0("A_", levels(donor)))
    if (!is.null(Z)) X <- cbind(X, Z)
    X
  }
  pr <- .profile_lambda(build_X, y, lambda_min, lambda_max, n_grid, conf)
  pr$n_donors <- nlevels(donor)
  pr$plateau_mean <- NA_real_; pr$amplitude_mean <- NA_real_
  if (is.finite(pr$lambda)) {
    f  <- stats::lm.fit(build_X(pr$lambda), y)
    cf <- f$coefficients
    pr$plateau_mean   <- mean(cf[paste0("P_", levels(donor))], na.rm = TRUE)
    pr$amplitude_mean <- mean(cf[paste0("A_", levels(donor))], na.rm = TRUE)
  }
  pr
}

# AIC comparison of the decay against the models it has to beat, on IDENTICAL
# rows. Requires d > 0 (the log-linear competitor needs it) -- the callers apply
# the canonical stop() on non-positive distance before getting here.
#
# Pass `donor` to make the comparison DONOR-STRATIFIED, which is what you want
# whenever the estimator itself is: every shape then gets a donor-specific
# intercept AND a donor-specific slope/amplitude, so the AIC difference reflects
# the SHAPE of the distance response and not between-donor level differences that
# happen to correlate with the distance distribution. Donor-blind (donor = NULL)
# is kept for single-donor calls, where there is nothing to stratify by.
decay_model_comparison <- function(d, y, Z = NULL,
                                   lambda_min = DECAY_LAMBDA_MIN_UM, lambda_max = 3000,
                                   n_grid = DECAY_N_GRID, donor = NULL) {
  if (!is.null(Z)) Z <- as.matrix(Z)
  ok <- is.finite(d) & is.finite(y) & (if (is.null(Z)) TRUE else stats::complete.cases(Z)) &
        (if (is.null(donor)) TRUE else !is.na(donor))
  d <- d[ok]; y <- y[ok]; if (!is.null(Z)) Z <- Z[ok, , drop = FALSE]
  if (!is.null(donor)) donor <- droplevels(factor(donor[ok]))
  if (any(d <= 0)) stop("Non-positive dist_to_phf1_um in decay_model_comparison(); ",
                        "log-linear comparison requires dist > 0.")
  n  <- length(y)
  nz <- if (is.null(Z)) 0L else ncol(Z)

  lin_rss <- function(X) {
    f <- .lm.fit(X, y)
    if (f$rank < ncol(X)) return(NA_real_)
    sum(f$residuals^2)
  }

  if (is.null(donor)) {
    base <- if (is.null(Z)) matrix(1, n, 1, dimnames = list(NULL, "(Intercept)")) else cbind(1, Z)
    dec   <- fit_decay_profile(d, y, Z, lambda_min, lambda_max, n_grid)
    r_log <- lin_rss(cbind(base, log(d)))
    r_lin <- lin_rss(cbind(base, d))
    r_nul <- lin_rss(base)
    npar  <- c(2L + nz + 1L, 2L + nz, 2L + nz, 1L + nz)
    lam   <- dec$lambda; r_dec <- dec$rss
  } else {
    D <- stats::model.matrix(~ 0 + donor)
    J <- ncol(D)
    base <- if (is.null(Z)) D else cbind(D, Z)
    dec   <- fit_decay_shared(d, y, donor, Z, lambda_min, lambda_max, n_grid)
    r_log <- lin_rss(cbind(base, D * log(d)))   # donor-specific log-linear slopes
    r_lin <- lin_rss(cbind(base, D * d))        # donor-specific linear slopes
    r_nul <- lin_rss(base)                      # donor-specific intercepts only
    npar  <- c(2L * J + nz + 1L, 2L * J + nz, 2L * J + nz, J + nz)
    lam   <- dec$lambda; r_dec <- dec$rss
  }

  out <- data.frame(
    model = c("decay", "log_linear", "linear", "plateau_only"),
    n_par = npar, rss = c(r_dec, r_log, r_lin, r_nul), stringsAsFactors = FALSE)
  out$aic   <- .gauss_aic(out$rss, n, out$n_par + 1L)
  out$delta <- out$aic - min(out$aic, na.rm = TRUE)
  out$n     <- n
  attr(out, "best")      <- out$model[which.min(out$aic)]
  attr(out, "lambda")    <- lam
  attr(out, "stratified") <- !is.null(donor)
  out
}

# nls cross-check on the covariate-adjusted score: does a general-purpose
# Gauss-Newton optimiser land on the same optimum the profile scan found?
# Started FROM the profile estimate (SSasymp's self-start assumes a response
# rising to its asymptote and fails outright on decaying data). Cross-check
# only: never the reported estimate, and allowed to fail.
decay_nls_check <- function(d, y, Z = NULL, start_lambda = NULL,
                            start_plateau = NULL, start_amp = NULL) {
  if (!is.null(Z)) {
    dd <- cbind(data.frame(y = y), as.data.frame(as.matrix(Z)))
    y  <- stats::residuals(stats::lm(y ~ ., data = dd)) + mean(y, na.rm = TRUE)
  }
  fit_with <- function(st) suppressWarnings(try(stats::nls(
    y ~ plateau + amp * exp(-d / lam), data = data.frame(d = d, y = y), start = st,
    control = stats::nls.control(warnOnly = FALSE, maxiter = 200)), silent = TRUE))

  m <- NULL
  if (all(vapply(list(start_lambda, start_plateau, start_amp), is.numeric, logical(1))) &&
      all(is.finite(c(start_lambda, start_plateau, start_amp))))
    m <- fit_with(list(plateau = start_plateau, amp = start_amp, lam = start_lambda))
  if (is.null(m) || inherits(m, "try-error"))                       # last resort: self-start
    m <- suppressWarnings(try(stats::nls(y ~ stats::SSasymp(d, Asym, R0, lrc),
                                         data = data.frame(d = d, y = y)), silent = TRUE))
  if (inherits(m, "try-error") || is.null(m))
    return(list(lambda = NA_real_, ok = FALSE, message = paste(as.character(m), collapse = " ")))

  cf <- stats::coef(m)
  if ("lrc" %in% names(cf))                                          # SSasymp parameterisation
    return(list(lambda = unname(exp(-cf[["lrc"]])), ok = TRUE,
                plateau = unname(cf[["Asym"]]),
                amplitude = unname(cf[["R0"]] - cf[["Asym"]]), message = ""))
  list(lambda = unname(cf[["lam"]]), ok = TRUE, plateau = unname(cf[["plateau"]]),
       amplitude = unname(cf[["amp"]]), message = "")
}

# Random-effects meta-analysis of log lambda. DerSimonian-Laird tau^2 with a
# Knapp-Hartung t-adjusted CI on k-1 df (essential at k <= 9); the HK scale
# factor is truncated below at 1 so the interval is never narrower than DL.
meta_log_lambda <- function(y, se, conf = DECAY_CONF) {
  keep <- is.finite(y) & is.finite(se) & se > 0
  y <- y[keep]; se <- se[keep]; k <- length(y)
  na <- list(k = k, mu = NA_real_, se_dl = NA_real_, se_hk = NA_real_,
             ci_lo = NA_real_, ci_hi = NA_real_, tau2 = NA_real_, Q = NA_real_,
             Q_df = NA_integer_, Q_p = NA_real_, I2 = NA_real_,
             lambda = NA_real_, lambda_lo = NA_real_, lambda_hi = NA_real_,
             lambda_geomean = NA_real_, lambda_median = NA_real_, hk_factor = NA_real_)
  if (k < 1) return(na)
  if (k == 1) {
    crit <- stats::qnorm(1 - (1 - conf) / 2)
    return(utils::modifyList(na, list(
      mu = y, se_dl = se, se_hk = se, tau2 = 0, Q = 0, Q_df = 0L, I2 = 0,
      ci_lo = y - crit * se, ci_hi = y + crit * se, hk_factor = 1,
      lambda = exp(y), lambda_lo = exp(y - crit * se), lambda_hi = exp(y + crit * se),
      lambda_geomean = exp(y), lambda_median = exp(y))))
  }
  v <- se^2; w <- 1 / v
  mu_fe <- sum(w * y) / sum(w)
  Q     <- sum(w * (y - mu_fe)^2)
  C     <- sum(w) - sum(w^2) / sum(w)
  tau2  <- max(0, (Q - (k - 1)) / C)
  ws    <- 1 / (v + tau2)
  mu    <- sum(ws * y) / sum(ws)
  se_dl <- sqrt(1 / sum(ws))
  hk    <- max(sum(ws * (y - mu)^2) / (k - 1), 1)   # truncated Knapp-Hartung
  se_hk <- se_dl * sqrt(hk)
  crit  <- stats::qt(1 - (1 - conf) / 2, df = k - 1)
  list(k = k, mu = mu, se_dl = se_dl, se_hk = se_hk, hk_factor = hk,
       ci_lo = mu - crit * se_hk, ci_hi = mu + crit * se_hk,
       tau2 = tau2, Q = Q, Q_df = k - 1L,
       Q_p = stats::pchisq(Q, df = k - 1, lower.tail = FALSE),
       I2 = if (Q > k - 1) 100 * (Q - (k - 1)) / Q else 0,
       lambda = exp(mu), lambda_lo = exp(mu - crit * se_hk),
       lambda_hi = exp(mu + crit * se_hk),
       lambda_geomean = exp(mean(y)), lambda_median = stats::median(exp(y)))
}

# Systematic distance profile of a null (SEG) score, per donor.
#
# WHY THIS EXISTS. Using the RAW per-cell SEG score as the baseline covariate does
# not work: a module score is a noisy per-cell measurement, and the systematic
# distance-dependent part of it is small next to that noise (in this project the
# SEG profile spans ~0.005 while per-cell SEG scatter is an order of magnitude
# larger). Regressing on a proxy that noisy attenuates its coefficient towards
# zero -- classical regression dilution -- so the raw covariate silently removes
# almost NONE of the distance drift it was added to remove. Averaging SEG within
# distance bins inside each donor estimates the same systematic profile with the
# per-cell noise divided out, which is what "use the null as the baseline" needs.
#
# Quantile bins (not fixed-width) so every bin carries cells despite the heavily
# right-skewed distance distribution; bin count is capped so each holds >= min_per_bin.
# Returns one value per input cell: its donor's mean null score at that distance.
seg_distance_profile <- function(dist, seg, donor, n_bins = 10, min_per_bin = 50) {
  out <- rep(NA_real_, length(seg))
  for (s in unique(as.character(donor))) {
    i <- which(as.character(donor) == s & is.finite(dist) & is.finite(seg))
    if (!length(i)) next
    nb <- max(1L, min(n_bins, floor(length(i) / min_per_bin)))
    if (nb <= 1L) { out[i] <- mean(seg[i]); next }
    br <- unique(stats::quantile(dist[i], probs = seq(0, 1, length.out = nb + 1), na.rm = TRUE))
    if (length(br) < 3L) { out[i] <- mean(seg[i]); next }
    b  <- cut(dist[i], breaks = br, include.lowest = TRUE)
    mu <- tapply(seg[i], b, mean, na.rm = TRUE)
    out[i] <- as.numeric(mu[as.character(b)])
  }
  # any cell the binning could not place keeps the donor mean, else the global mean
  if (anyNA(out)) {
    gm <- mean(seg, na.rm = TRUE)
    for (s in unique(as.character(donor))) {
      i <- which(as.character(donor) == s & is.na(out))
      if (length(i)) out[i] <- mean(seg[as.character(donor) == s], na.rm = TRUE)
    }
    out[is.na(out)] <- gm
  }
  out
}

# Gradient magnitude of a signature against the null, per donor, in SD units.
#
# WHY THIS EXISTS, and why it is the RIGHT use of a stably-expressed-gene null.
# Putting the null in as a covariate cannot answer "is the signature's gradient
# just the null's gradient?" when the two share a shape: at the fitted lambda the
# exponential term and the null profile are near-collinear, so least squares splits
# the shared variance arbitrarily (the exact functional form outcompetes the binned
# null, and the amplitude barely moves even when the null IS a scaled copy of the
# signal). Collinear predictors cannot arbitrate ownership of a gradient.
#
# What the null CAN do is calibrate MAGNITUDE. Both scores are put on their own SD
# scale and the near-to-far drop is measured WITHIN each donor (so donor
# composition cannot drive it), then averaged across donors with a t CI. The ratio
# says how many times steeper the signature's gradient is than the null
# expectation: ~1 means the signature is doing nothing a stably expressed gene set
# does not also do; >>1 is a genuine signature-specific gradient. It is an
# enrichment-over-null statement, not a length constant.
null_benchmark_gradient <- function(dat, score_col, seg_col,
                                    dist_col = "dist_to_phf1_um",
                                    donor_col = "sample_id",
                                    near_um = 100, far_um = 500,
                                    min_per_side = 25, conf = DECAY_CONF) {
  donors <- unique(as.character(dat[[donor_col]]))
  rows <- lapply(donors, function(s) {
    i <- as.character(dat[[donor_col]]) == s
    d <- as.numeric(dat[[dist_col]])[i]
    one <- function(col) {
      y <- as.numeric(dat[[col]])[i]
      nr <- y[d <= near_um & is.finite(y)]; fr <- y[d >= far_um & is.finite(y)]
      if (length(nr) < min_per_side || length(fr) < min_per_side) return(NA_real_)
      sdy <- stats::sd(y, na.rm = TRUE)
      if (!is.finite(sdy) || sdy <= 0) return(NA_real_)
      (mean(nr) - mean(fr)) / sdy          # positive = higher NEAR tangles
    }
    data.frame(sample_id = s, sig_drop_sd = one(score_col),
               seg_drop_sd = if (is.null(seg_col)) NA_real_ else one(seg_col),
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  est <- function(v) {
    v <- v[is.finite(v)]; k <- length(v)
    if (k < 2) return(c(mean = if (k) v else NA_real_, lo = NA_real_, hi = NA_real_, k = k))
    tc <- stats::qt(1 - (1 - conf) / 2, k - 1); se <- stats::sd(v) / sqrt(k)
    c(mean = mean(v), lo = mean(v) - tc * se, hi = mean(v) + tc * se, k = k)
  }
  s <- est(tab$sig_drop_sd); g <- est(tab$seg_drop_sd)
  paired <- tab[is.finite(tab$sig_drop_sd) & is.finite(tab$seg_drop_sd), , drop = FALSE]
  diff_est <- est(paired$sig_drop_sd - paired$seg_drop_sd)
  list(per_donor = tab, near_um = near_um, far_um = far_um,
       sig = s, seg = g, excess = diff_est,
       ratio = if (is.finite(g[["mean"]]) && abs(g[["mean"]]) > 1e-12)
                 s[["mean"]] / g[["mean"]] else NA_real_,
       n_donors_same_dir = sum(paired$sig_drop_sd > paired$seg_drop_sd, na.rm = TRUE),
       n_donors_paired = nrow(paired))
}

# Human-readable block for the null benchmark.
decay_write_null_benchmark <- function(nb, score_label = "signature") {
  cat(sprintf("Near = within %g um, far = beyond %g um; drop = (mean near - mean far) / SD of\n",
              nb$near_um, nb$far_um))
  cat("that score, computed WITHIN each donor then averaged across donors (t CI on k-1 df).\n")
  cat("Positive = higher near tangles. This is a magnitude benchmark, NOT a length constant.\n")
  f <- function(e) if (is.finite(e[["mean"]]))
    sprintf("%+.3f SD [%+.3f, %+.3f] (k = %d donors)", e[["mean"]], e[["lo"]], e[["hi"]], e[["k"]])
    else "not computable"
  cat(sprintf("  %-24s %s\n", paste0(score_label, ":"), f(nb$sig)))
  cat(sprintf("  %-24s %s\n", "SEG null:", f(nb$seg)))
  cat(sprintf("  %-24s %s\n", "excess (signature - null):", f(nb$excess)))
  cat(sprintf("  ratio signature/null = %s\n",
              if (is.finite(nb$ratio)) sprintf("%.2fx", nb$ratio) else "NA"))
  cat(sprintf("  donors where the signature's gradient exceeds the null's: %d / %d\n",
              nb$n_donors_same_dir, nb$n_donors_paired))
  cat("  A ratio near 1, or an excess CI covering 0, means the signature's distance gradient\n")
  cat("  is no steeper than that of the SEG null.\n")
  invisible(NULL)
}

# Distance at which a fraction `frac` of the excess has resolved.
# d_frac(lambda, 0.95) = lambda * log(20) = 2.996 * lambda.
d_frac <- function(lambda, frac = 0.95) lambda * log(1 / (1 - frac))

decay_curve <- function(grid, plateau, amplitude, lambda) {
  plateau + amplitude * exp(-grid / lambda)
}

# Per-donor fits with the pre-registered inclusion rules. Returns one row per
# donor INCLUDING the excluded ones, with the reason -- dropping donors whose
# gradient is flat biases the pooled lambda toward the steepest donors, so the
# exclusions must be visible in the source data, not silently filtered.
decay_fit_donors <- function(dat, score_col, dist_col = "dist_to_phf1_um",
                             donor_col = "sample_id",
                             covar_cols = c("nUMI_log", "percent_neg"),
                             lambda_min = DECAY_LAMBDA_MIN_UM, lambda_max = 3000,
                             min_cells = DECAY_MIN_CELLS_DONOR,
                             min_span = DECAY_MIN_SPAN_UM,
                             n_grid = DECAY_N_GRID, conf = DECAY_CONF,
                             center = TRUE) {
  covar_cols <- intersect(covar_cols, colnames(dat))
  donors <- sort(unique(as.character(dat[[donor_col]])))
  rows <- lapply(donors, function(s) {
    sub <- dat[as.character(dat[[donor_col]]) == s, , drop = FALSE]
    d <- as.numeric(sub[[dist_col]]); y <- as.numeric(sub[[score_col]])
    Z <- if (length(covar_cols)) as.matrix(sub[, covar_cols, drop = FALSE]) else NULL
    ok <- is.finite(d) & is.finite(y) & (if (is.null(Z)) TRUE else stats::complete.cases(Z))
    d <- d[ok]; y <- y[ok]; if (!is.null(Z)) Z <- Z[ok, , drop = FALSE]
    span <- if (length(d)) diff(range(d)) else NA_real_
    blank <- data.frame(
      sample_id = s, n_cells = length(d), dist_span_um = span,
      lambda_um = NA_real_, lambda_lo = NA_real_, lambda_hi = NA_real_,
      log_lambda = NA_real_, se_log_lambda = NA_real_,
      plateau = NA_real_, amplitude = NA_real_, amp_se = NA_real_,
      converged = NA, at_bound = NA, ci_open = NA,
      included = FALSE, exclusion_reason = NA_character_,
      stringsAsFactors = FALSE)

    if (length(d) < min_cells) {
      blank$exclusion_reason <- sprintf("n_cells %d < %d", length(d), min_cells); return(blank)
    }
    if (!is.finite(span) || span < min_span) {
      blank$exclusion_reason <- sprintf("distance span %.0f um < %.0f um", span, min_span)
      return(blank)
    }
    fit <- fit_decay_profile(d, y, Z, lambda_min, lambda_max, n_grid, conf, center)
    blank$lambda_um  <- fit$lambda
    blank$lambda_lo  <- fit$lambda_lo
    blank$lambda_hi  <- fit$lambda_hi
    blank$plateau    <- fit$plateau
    blank$amplitude  <- fit$amplitude
    blank$amp_se     <- fit$amp_se
    blank$converged  <- isTRUE(fit$converged)
    blank$at_bound   <- isTRUE(fit$at_bound_lo) || isTRUE(fit$at_bound_hi)
    blank$ci_open    <- isTRUE(fit$ci_open_lo)  || isTRUE(fit$ci_open_hi)
    if (!is.finite(fit$lambda)) { blank$exclusion_reason <- "profile fit failed"; return(blank) }
    if (blank$at_bound) {
      blank$exclusion_reason <- "lambda at a search bound (not identified)"; return(blank)
    }
    if (!is.finite(fit$lambda_lo) || !is.finite(fit$lambda_hi) || blank$ci_open) {
      blank$exclusion_reason <- "profile CI open at a bound (not identified)"; return(blank)
    }
    blank$log_lambda    <- log(fit$lambda)
    blank$se_log_lambda <- (log(fit$lambda_hi) - log(fit$lambda_lo)) /
                           (2 * stats::qnorm(1 - (1 - conf) / 2))
    if (!is.finite(blank$se_log_lambda) || blank$se_log_lambda <= 0) {
      blank$exclusion_reason <- "non-finite SE(log lambda)"; return(blank)
    }
    blank$included <- TRUE
    blank
  })
  do.call(rbind, rows)
}

# Full pipeline for one module: per-donor fits -> meta-analysis -> shared-lambda
# anchor -> model comparison -> cap sensitivity -> gates -> plot curve.
#
# `cap_um` is the distance cap the caller has ALREADY applied to `dat`; it is
# used for the half-cap refit and the identifiability rule, not to filter again.
decay_analyse <- function(dat, score_col, cap_um,
                          dist_col = "dist_to_phf1_um", donor_col = "sample_id",
                          covar_cols = c("nUMI_log", "percent_neg"),
                          lambda_min = DECAY_LAMBDA_MIN_UM,
                          lambda_max_mult = DECAY_LAMBDA_MAX_MULT,
                          min_cells = DECAY_MIN_CELLS_DONOR,
                          min_span = DECAY_MIN_SPAN_UM,
                          n_grid = DECAY_N_GRID, conf = DECAY_CONF,
                          grid_n = 200) {
  lambda_max <- lambda_max_mult * cap_um
  covar_cols <- intersect(covar_cols, colnames(dat))

  # Centre the covariates ONCE, on cohort means, so every donor's plateau is the
  # background score at the same typical nUMI_log / percent_neg and the plateaus
  # are comparable with each other and with the cohort fit. Everything downstream
  # is then told not to re-centre per donor. Lambda and the amplitudes are
  # unaffected by this -- only the intercept moves.
  for (cc in covar_cols) dat[[cc]] <- as.numeric(dat[[cc]]) - mean(as.numeric(dat[[cc]]), na.rm = TRUE)

  donors <- decay_fit_donors(dat, score_col, dist_col, donor_col, covar_cols,
                             lambda_min, lambda_max, min_cells, min_span, n_grid, conf,
                             center = FALSE)
  inc  <- donors[donors$included, , drop = FALSE]
  meta <- meta_log_lambda(inc$log_lambda, inc$se_log_lambda, conf)

  d <- as.numeric(dat[[dist_col]]); y <- as.numeric(dat[[score_col]])
  Z <- if (length(covar_cols)) as.matrix(dat[, covar_cols, drop = FALSE]) else NULL
  score_sd <- stats::sd(y, na.rm = TRUE)
  shared  <- fit_decay_shared(d, y, dat[[donor_col]], Z, lambda_min, lambda_max, n_grid,
                              conf, center = FALSE)
  # Donor-stratified, to match the donor-stratified estimator (see the function's
  # comment): otherwise the comparison can crown `decay` on between-donor level
  # differences that track the distance distribution rather than on curve shape.
  compare <- decay_model_comparison(d, y, Z, lambda_min, lambda_max, n_grid,
                                    donor = dat[[donor_col]])
  nls_chk <- decay_nls_check(d, y, Z, start_lambda = shared$lambda,
                             start_plateau = shared$plateau_mean,
                             start_amp = shared$amplitude_mean)

  # Per-donor model comparison: in how many donors does the decay actually win?
  dwin <- vapply(inc$sample_id, function(s) {
    sub <- dat[as.character(dat[[donor_col]]) == s, , drop = FALSE]
    cmp <- try(decay_model_comparison(as.numeric(sub[[dist_col]]),
                                      as.numeric(sub[[score_col]]),
                                      if (length(covar_cols)) as.matrix(sub[, covar_cols, drop = FALSE]) else NULL,
                                      lambda_min, lambda_max, n_grid), silent = TRUE)
    if (inherits(cmp, "try-error")) return(NA_character_)
    attr(cmp, "best")
  }, character(1))

  # Cap sensitivity: same estimator on the inner half of the distance window.
  half <- dat[as.numeric(dat[[dist_col]]) <= cap_um / 2, , drop = FALSE]
  meta_half <- list(lambda = NA_real_, k = 0)
  if (nrow(half) >= min_cells) {
    dh <- decay_fit_donors(half, score_col, dist_col, donor_col, covar_cols,
                           lambda_min, lambda_max_mult * cap_um / 2,
                           min_cells, min_span / 2, n_grid, conf)
    ih <- dh[dh$included, , drop = FALSE]
    if (nrow(ih) >= 1) meta_half <- meta_log_lambda(ih$log_lambda, ih$se_log_lambda, conf)
  }
  cap_ratio <- meta_half$lambda / meta$lambda

  aic_dec  <- compare$aic[compare$model == "decay"]
  aic_log  <- compare$aic[compare$model == "log_linear"]
  aic_lin  <- compare$aic[compare$model == "linear"]
  ci_ratio <- meta$lambda_hi / meta$lambda_lo
  gates <- list(
    delta_aic_vs_log      = aic_dec - aic_log,          # negative = decay wins
    beats_log             = isTRUE(aic_dec < aic_log - DECAY_MIN_DAIC),
    # The decay must also beat a PLAIN LINEAR decline, because a straight line over
    # the window is exactly the lambda -> infinity limit of the exponential: a slope
    # with no plateau, and therefore no length constant.
    delta_aic_vs_linear   = aic_dec - aic_lin,
    beats_linear          = isTRUE(aic_dec < aic_lin - DECAY_MIN_DAIC),
    best_model            = attr(compare, "best"),
    n_donors_decay_wins   = sum(dwin == "decay", na.rm = TRUE),
    n_donors_included     = nrow(inc),
    n_donors_total        = nrow(donors),
    lambda_ratio_half_cap = cap_ratio,
    cap_sensitive         = !is.finite(cap_ratio) ||
                            cap_ratio < DECAY_CAP_RATIO_LO || cap_ratio > DECAY_CAP_RATIO_HI,
    ci_ratio              = ci_ratio,
    precise               = isTRUE(is.finite(ci_ratio) && ci_ratio <= DECAY_MAX_CI_RATIO),
    identifiable          = isTRUE(is.finite(meta$lambda) && meta$lambda <= cap_um / 2 &&
                                   is.finite(meta$lambda_hi) && meta$lambda_hi <= cap_um / 2))
  gates$quotable <- isTRUE(gates$beats_log && gates$beats_linear && gates$identifiable &&
                           gates$precise && !gates$cap_sensitive &&
                           gates$n_donors_included >= 3)
  gates$verdict <- if (!is.finite(meta$lambda)) "lambda not estimable"
    else if (!gates$beats_log) sprintf("NOT quotable: decay does not beat log-linear by %g AIC (dAIC %+.1f) -- gradient without an identifiable length constant", DECAY_MIN_DAIC, gates$delta_aic_vs_log)
    else if (!gates$beats_linear) sprintf("NOT quotable: decay does not beat a plain linear decline by %g AIC (dAIC %+.1f) -- a slope across the whole window, i.e. no plateau and no length constant", DECAY_MIN_DAIC, gates$delta_aic_vs_linear)
    else if (!gates$identifiable) sprintf("NOT quotable: lambda (%.0f um, CI to %.0f um) not resolved inside the %.0f um window", meta$lambda, meta$lambda_hi, cap_um)
    else if (!gates$precise) sprintf("NOT quotable: CI spans %.1f-fold [%.0f, %.0f] um -- too wide to quote", ci_ratio, meta$lambda_lo, meta$lambda_hi)
    else if (gates$cap_sensitive) sprintf("NOT quotable: lambda tracks the cap (half-cap ratio %.2f)", cap_ratio)
    else if (gates$n_donors_included < 3) sprintf("NOT quotable: only %d donors identified lambda", gates$n_donors_included)
    else sprintf("quotable length constant: lambda = %.0f um [%.0f, %.0f]", meta$lambda, meta$lambda_lo, meta$lambda_hi)

  # Plot curve at the pooled lambda, donor-mean plateau and amplitude. Ribbon =
  # the lambda-CI envelope widened by the between-donor SEM of plateau and
  # amplitude (so it carries both the length-constant and the level uncertainty).
  curve <- NULL; plateau <- NA_real_; amplitude <- NA_real_
  d_bg <- NA_real_
  amp_ci <- c(NA_real_, NA_real_); plat_ci <- c(NA_real_, NA_real_)
  if (is.finite(meta$lambda) && nrow(inc) >= 1) {
    plateau   <- mean(inc$plateau,   na.rm = TRUE)
    amplitude <- mean(inc$amplitude, na.rm = TRUE)
    sem_p <- if (nrow(inc) > 1) stats::sd(inc$plateau,   na.rm = TRUE) / sqrt(nrow(inc)) else 0
    sem_a <- if (nrow(inc) > 1) stats::sd(inc$amplitude, na.rm = TRUE) / sqrt(nrow(inc)) else 0
    # Amplitude / plateau CIs are BETWEEN-DONOR t-intervals on k-1 df, matching how
    # lambda is treated. The per-donor CIs in `donors` are within-donor and
    # conditional on that donor's lambda_hat, so they are not pooled here.
    if (nrow(inc) > 1) {
      tc <- stats::qt(1 - (1 - conf) / 2, df = nrow(inc) - 1)
      amp_ci  <- amplitude + c(-1, 1) * tc * sem_a
      plat_ci <- plateau   + c(-1, 1) * tc * sem_p
    }
    g   <- seq(0, cap_um, length.out = grid_n)
    fit <- decay_curve(g, plateau, amplitude, meta$lambda)
    f1  <- decay_curve(g, plateau, amplitude, meta$lambda_lo)
    f2  <- decay_curve(g, plateau, amplitude, meta$lambda_hi)
    half_w <- 1.96 * sqrt(sem_p^2 + (sem_a * exp(-g / meta$lambda))^2)
    curve <- data.frame(dist_to_phf1_um = g, fitted_score = fit,
                        ci_lo = pmin(fit, f1, f2) - half_w,
                        ci_hi = pmax(fit, f1, f2) + half_w,
                        plateau = plateau, amplitude = amplitude,
                        lambda_um = meta$lambda)
    # Descriptive "reaches background": first grid distance at which the curve's
    # band overlaps the plateau's band. Noise-dependent, hence secondary to d95.
    p_lo <- plateau - 1.96 * sem_p; p_hi <- plateau + 1.96 * sem_p
    hit  <- which(curve$ci_lo <= p_hi & curve$ci_hi >= p_lo)
    if (length(hit)) d_bg <- g[hit[1]]
  }

  list(donors = donors, meta = meta, shared = shared, compare = compare,
       nls_check = nls_chk, donor_best = data.frame(sample_id = inc$sample_id,
                                                    best_model = unname(dwin),
                                                    stringsAsFactors = FALSE),
       meta_half = meta_half, gates = gates, curve = curve,
       plateau = plateau, plateau_ci_lo = plat_ci[1], plateau_ci_hi = plat_ci[2],
       amplitude = amplitude, amp_ci_lo = amp_ci[1], amp_ci_hi = amp_ci[2],
       # Standardised amplitude: the excess at the tangle in SD units of the module
       # score. Robust where the plateau sits near zero -- which it does, because
       # AddModuleScore is centred on control genes, so "% of plateau" is unstable
       # and is deliberately NOT reported.
       score_sd = score_sd, amp_in_sd = amplitude / score_sd,
       amp_in_sd_lo = amp_ci[1] / score_sd, amp_in_sd_hi = amp_ci[2] / score_sd,
       d_background_um = d_bg,
       cap_um = cap_um, lambda_min = lambda_min, lambda_max = lambda_max,
       covar_cols = covar_cols,
       d95_um = d_frac(meta$lambda), d95_lo = d_frac(meta$lambda_lo),
       d95_hi = d_frac(meta$lambda_hi))
}

# Human-readable block for a stats_*.txt log. `res` is a decay_analyse() result.
decay_write_stats_block <- function(res, label = "module score") {
  m <- res$meta; g <- res$gates
  cat("MODEL: S = P_j + A_j * exp(-d / lambda_j) + b1*nUMI_log + b2*percent_neg\n")
  cat("  fitted PER DONOR by profile least squares (the model is linear in P, A and the\n")
  cat("  covariate coefficients once lambda is fixed, so lambda comes from a 1-D profile\n")
  cat("  and its CI is a profile-likelihood CI -- no nls convergence failures).\n")
  cat("  Sex/Age/PMI are omitted: donor-level, exactly aliased with P_j.\n")
  cat("  A is SIGNED -- a module that RISES away from the tangle has A < 0 and the same lambda.\n\n")
  cat("DISTANCE TRANSFORM: raw microns. This is a deliberate, documented departure from the\n")
  cat("  canonical log(dist)/sd(log(dist)) transform: lambda has units of microns\n")
  cat("  and is undefined on a log/SD-rescaled axis. The 'gam' model type is the precedent.\n")
  cat("  The canonical transform is unchanged in every other analysis.\n\n")
  cat(sprintf("Cap: %.0f um | lambda search range [%.0f, %.0f] um | covariates: %s\n\n",
              res$cap_um, res$lambda_min, res$lambda_max,
              if (length(res$covar_cols)) paste(res$covar_cols, collapse = ", ") else "none"))

  cat("=== Per-donor fits (inclusion rules applied; excluded donors shown) ===\n")
  print(res$donors, row.names = FALSE, digits = 4)
  cat(sprintf("\nIncluded %d of %d donors.\n", g$n_donors_included, g$n_donors_total))
  if (any(!res$donors$included))
    cat("  Excluded donors can bias the pooled lambda toward the steepest gradients; see the\n",
        "  exclusion_reason column.\n", sep = "")

  cat("\n=== Random-effects meta-analysis of log lambda (DerSimonian-Laird + Knapp-Hartung) ===\n")
  cat(sprintf("k = %d donors\n", m$k))
  cat(sprintf("tau^2 = %.4f | Q = %.2f on %s df (p = %s) | I^2 = %.1f%%\n",
              m$tau2, m$Q, format(m$Q_df), format.pval(m$Q_p, digits = 3), m$I2))
  cat(sprintf("Knapp-Hartung scale factor = %.3f (truncated below at 1)\n", m$hk_factor))
  cat(sprintf("SE(log lambda): DL %.4f -> HK %.4f\n", m$se_dl, m$se_hk))
  cat(sprintf("Unweighted geometric mean lambda = %.1f um | median = %.1f um\n",
              m$lambda_geomean, m$lambda_median))

  cat("\n=== Effect sizes ===\n")
  cat(sprintf("LENGTH CONSTANT  lambda = %.1f um   95%% CI [%.1f, %.1f]\n",
              m$lambda, m$lambda_lo, m$lambda_hi))
  cat(sprintf("DISTANCE TO BACKGROUND  d95 = 2.996*lambda = %.0f um   95%% CI [%.0f, %.0f]\n",
              res$d95_um, res$d95_lo, res$d95_hi))
  cat("  (d95 = distance at which 95% of the excess has resolved; the CI is the exact\n")
  cat("   monotone transform of the pooled lambda CI.)\n")
  cat(sprintf("  descriptive band-overlap distance = %s um (secondary, noise-dependent)\n",
              if (is.finite(res$d_background_um)) sprintf("%.0f", res$d_background_um) else "NA"))
  cat(sprintf("AMPLITUDE at d=0  A = %+.4f   95%% CI [%+.4f, %+.4f]  (between-donor, %d donors)\n",
              res$amplitude, res$amp_ci_lo, res$amp_ci_hi, m$k))
  cat(sprintf("  standardised: %+.3f SD of the module score (SD = %.4f)   95%% CI [%+.3f, %+.3f]\n",
              res$amp_in_sd, res$score_sd, res$amp_in_sd_lo, res$amp_in_sd_hi))
  cat(sprintf("BACKGROUND PLATEAU  P = %.4f   95%% CI [%.4f, %.4f]\n",
              res$plateau, res$plateau_ci_lo, res$plateau_ci_hi))
  cat("  Covariates are centred on cohort means, so P is the background score at typical\n")
  cat("  nUMI_log / percent_neg. Amplitude is NOT expressed as a % of the plateau:\n")
  cat("  AddModuleScore is centred on control genes so the plateau sits near zero and the\n")
  cat("  ratio is unstable -- the SD-standardised amplitude above is the scale-free measure.\n")
  cat("  A > 0 means the score is ELEVATED near tangles; A < 0 means it rises with distance.\n")
  cat("  Per-donor amplitude CIs in the table above are within-donor and CONDITIONAL on\n")
  cat("  that donor's lambda_hat; the CI quoted here is the between-donor one.\n")

  cat("\n=== Gates (these decide whether lambda is quotable) ===\n")
  cat(sprintf("1. Model comparison on identical rows (%s):\n",
              if (isTRUE(attr(res$compare, "stratified")))
                "donor-stratified: every shape gets a donor-specific intercept AND slope/amplitude,\n   so the AIC gap reflects curve SHAPE, not between-donor levels"
              else "donor-pooled"))
  print(res$compare, row.names = FALSE, digits = 6)
  cat(sprintf("   best = %s | dAIC(decay - log_linear) = %+.1f -> decay %s (needs a %g-unit margin)\n",
              g$best_model, g$delta_aic_vs_log,
              if (g$beats_log) "wins" else "LOSES", DECAY_MIN_DAIC))
  cat(sprintf("   dAIC(decay - linear) = %+.1f -> decay %s. A straight line over the window IS\n",
              g$delta_aic_vs_linear, if (g$beats_linear) "wins" else "LOSES"))
  cat("   the lambda -> infinity limit, so losing here means a slope with no plateau.\n")
  cat(sprintf("   decay wins in %d of %d included donors individually\n",
              g$n_donors_decay_wins, g$n_donors_included))
  cat(sprintf("1b. CI precision: lambda_hi/lambda_lo = %s -> %s (max %g-fold)\n",
              if (is.finite(g$ci_ratio)) sprintf("%.1f", g$ci_ratio) else "NA",
              if (g$precise) "usable" else "TOO WIDE", DECAY_MAX_CI_RATIO))
  cat(sprintf("2. Cap sensitivity: lambda(half cap) / lambda(full cap) = %s -> %s\n",
              if (is.finite(g$lambda_ratio_half_cap)) sprintf("%.2f", g$lambda_ratio_half_cap) else "NA",
              if (g$cap_sensitive) "CAP-SENSITIVE" else sprintf("stable (tolerance %.2f-%.2f)",
                                                                DECAY_CAP_RATIO_LO, DECAY_CAP_RATIO_HI)))
  cat(sprintf("   lambda at half cap = %s um (k = %d donors)\n",
              if (is.finite(res$meta_half$lambda)) sprintf("%.1f", res$meta_half$lambda) else "NA",
              res$meta_half$k))
  cat(sprintf("3. Identifiability: lambda %.0f um and CI upper %.0f um vs cap/2 = %.0f um -> %s\n",
              m$lambda, m$lambda_hi, res$cap_um / 2,
              if (g$identifiable) "resolved inside the window" else "NOT RESOLVED"))
  cat(sprintf("\nVERDICT: %s\n", g$verdict))

  cat("\n=== Consistency checks (not the reported estimate) ===\n")
  cat(sprintf("Shared-lambda fit (one lambda, donor-specific P and A, all %d donors, no donor dropped):\n",
              res$shared$n_donors))
  cat(sprintf("  lambda = %s um, profile CI [%s, %s]\n",
              sprintf("%.1f", res$shared$lambda),
              if (is.finite(res$shared$lambda_lo)) sprintf("%.1f", res$shared$lambda_lo) else "open",
              if (is.finite(res$shared$lambda_hi)) sprintf("%.1f", res$shared$lambda_hi) else "open"))
  cat(sprintf("nls(SSasymp) cross-check on covariate-adjusted score: lambda = %s\n",
              if (isTRUE(res$nls_check$ok)) sprintf("%.1f um", res$nls_check$lambda)
              else "failed to converge (expected; the profile estimator is the reported one)"))

  cat("\nNOTE (CI): the reported CI is the BETWEEN-donor meta-analytic one (k = ", m$k,
      " donors);\n", sep = "")
  cat("  lambda itself is the effect size, in microns.\n")
  cat("NOTE (gating): only modules with a significant LOG-model distance effect are fitted.\n")
  cat("  The log model stays primary for 'is there a gradient?'; this answers 'how far?'.\n")
  invisible(NULL)
}

# One-row summary for the written-once headline table.
decay_summary_row <- function(res, ...) {
  m <- res$meta; g <- res$gates
  data.frame(..., n_donors_included = g$n_donors_included,
             n_donors_total = g$n_donors_total,
             lambda_um = m$lambda, lambda_lo = m$lambda_lo, lambda_hi = m$lambda_hi,
             se_log_lambda_hk = m$se_hk, tau2 = m$tau2, Q = m$Q, I2 = m$I2,
             lambda_geomean_unweighted = m$lambda_geomean, lambda_median = m$lambda_median,
             d95_um = res$d95_um, d95_lo = res$d95_lo, d95_hi = res$d95_hi,
             d_background_um = res$d_background_um,
             plateau = res$plateau, plateau_ci_lo = res$plateau_ci_lo,
             plateau_ci_hi = res$plateau_ci_hi,
             amplitude = res$amplitude, amp_ci_lo = res$amp_ci_lo,
             amp_ci_hi = res$amp_ci_hi, score_sd = res$score_sd,
             amp_in_sd = res$amp_in_sd, amp_in_sd_lo = res$amp_in_sd_lo,
             amp_in_sd_hi = res$amp_in_sd_hi,
             lambda_shared = res$shared$lambda,
             lambda_shared_lo = res$shared$lambda_lo,
             lambda_shared_hi = res$shared$lambda_hi,
             lambda_nls_check = res$nls_check$lambda,
             aic_decay = res$compare$aic[res$compare$model == "decay"],
             aic_log = res$compare$aic[res$compare$model == "log_linear"],
             aic_linear = res$compare$aic[res$compare$model == "linear"],
             aic_plateau_only = res$compare$aic[res$compare$model == "plateau_only"],
             delta_aic_vs_log = g$delta_aic_vs_log, best_model = g$best_model,
             beats_log = g$beats_log,
             delta_aic_vs_linear = g$delta_aic_vs_linear, beats_linear = g$beats_linear,
             n_donors_decay_wins = g$n_donors_decay_wins,
             lambda_ratio_half_cap = g$lambda_ratio_half_cap,
             lambda_ci_ratio = g$ci_ratio, precise = g$precise,
             cap_sensitive = g$cap_sensitive, identifiable = g$identifiable,
             quotable = g$quotable, verdict = g$verdict,
             stringsAsFactors = FALSE)
}

# Forest plot of lambda (log x-axis): pooled diamond per row, per-donor points
# behind it. `tab` needs row_label, lambda_um, lambda_lo, lambda_hi, quotable;
# `donor_tab` (optional) needs row_label, lambda_um.
decay_lambda_forest <- function(tab, donor_tab = NULL, cap_um = NULL, colours = NULL) {
  tab <- tab[is.finite(tab$lambda_um), , drop = FALSE]
  if (!nrow(tab)) return(NULL)
  lev <- rev(as.character(tab$row_label))
  tab$row_label <- factor(as.character(tab$row_label), levels = lev)
  tab$shape_q   <- factor(ifelse(tab$quotable, "quotable", "flagged"),
                          levels = c("quotable", "flagged"))
  if (is.null(tab$series)) tab$series <- "lambda"

  p <- ggplot2::ggplot()
  if (!is.null(cap_um))
    p <- p + ggplot2::geom_vline(xintercept = cap_um / 2, linetype = "dotted",
                                 linewidth = 0.3, colour = "grey50")
  if (!is.null(donor_tab) && nrow(donor_tab)) {
    donor_tab <- donor_tab[is.finite(donor_tab$lambda_um) &
                             as.character(donor_tab$row_label) %in% lev, , drop = FALSE]
    if (nrow(donor_tab)) {
      donor_tab$row_label <- factor(as.character(donor_tab$row_label), levels = lev)
      p <- p + ggplot2::geom_point(
        data = donor_tab, ggplot2::aes(x = lambda_um, y = row_label),
        colour = "grey65", size = 0.7, alpha = 0.8,
        position = ggplot2::position_jitter(height = 0.14, width = 0, seed = 42))
    }
  }
  p <- p +
    # geom_errorbar(orientation = "y"), not geom_errorbarh(): the latter is
    # deprecated from ggplot2 4.0.
    ggplot2::geom_errorbar(data = tab,
      ggplot2::aes(y = row_label, xmin = lambda_lo, xmax = lambda_hi, colour = series),
      orientation = "y", width = 0, linewidth = 0.4) +
    ggplot2::geom_point(data = tab,
      ggplot2::aes(x = lambda_um, y = row_label, colour = series, shape = shape_q),
      size = 1.5) +
    ggplot2::scale_shape_manual(values = c(quotable = 18, flagged = 4),
                                limits = c("quotable", "flagged"), drop = FALSE, name = NULL) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = expression("Length constant " * lambda * " (" * mu * "m)"), y = NULL)
  p <- p + if (!is.null(colours))
    ggplot2::scale_colour_manual(values = colours, guide = "none")
  else ggplot2::scale_colour_discrete(guide = "none")
  p
}
