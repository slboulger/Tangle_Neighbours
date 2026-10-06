# ---------------------------------------------------------------------------
# imc_utils.R
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
# imc_utils.R
#
# Shared helpers for the IMC (EC project) analysis scripts. Sourced by every imc_*.R and
# null_replica_imc_*.R script; holds nothing analysis-specific.
#
# Every IMC stats log prints the donor random-effect variance and a singularity flag. Where the
# outcome is z-scored within donor (docs/MODELS.md Set 3, within-donor z-scoring), the outcome has
# no between-donor variance, the donor variance is estimated at zero and the fit is equivalent to
# pooled least squares on the standardised outcome; printing the variance makes that explicit.

# -------------------------------------------------------------------
# Random-effect variance reporting
# -------------------------------------------------------------------

#' Tabulate the variance components of a fitted merMod.
#'
#' @param m       A model from lmer/lmerTest::lmer/glmer, or NULL if the fit failed.
#' @param model   Optional label for the model, carried into the returned frame.
#' @param outcome Optional label for the outcome, carried into the returned frame.
#' @return data.frame with one row per variance component: group, term, vcov, sdcor,
#'         pct_of_total, plus is_singular / n_groups / note. One NA row if `m` is NULL.
re_variance_summary <- function(m, model = NA_character_, outcome = NA_character_) {
  na_row <- function(note) data.frame(
    outcome = outcome, model = model, group = NA_character_, term = NA_character_,
    vcov = NA_real_, sdcor = NA_real_, pct_of_total = NA_real_,
    is_singular = NA, n_groups = NA_integer_, note = note, stringsAsFactors = FALSE)

  if (is.null(m)) return(na_row("model not fitted"))
  if (!inherits(m, "merMod")) return(na_row(paste("not a merMod:", class(m)[1])))

  vc <- tryCatch(as.data.frame(lme4::VarCorr(m)), error = function(e) NULL)
  if (is.null(vc) || !nrow(vc)) return(na_row("VarCorr unavailable"))

  # var2 is non-NA only on covariance rows; those are not variances and must not be summed.
  vc <- vc[is.na(vc$var2), , drop = FALSE]
  total <- sum(vc$vcov, na.rm = TRUE)

  # glmer with a binomial/poisson family has no Residual row, so "% of total" is a share of the
  # random-effect variances alone. Flagged in the note rather than silently rescaled.
  has_resid <- "Residual" %in% vc$grp
  ng <- tryCatch(vapply(lme4::getME(m, "flist"), function(f) length(levels(f)), integer(1)),
                 error = function(e) integer(0))

  data.frame(
    outcome = outcome, model = model,
    group = vc$grp,
    term  = ifelse(is.na(vc$var1), "", vc$var1),
    vcov = vc$vcov, sdcor = vc$sdcor,
    pct_of_total = if (isTRUE(total > 0)) 100 * vc$vcov / total else NA_real_,
    is_singular = isTRUE(lme4::isSingular(m, tol = 1e-4)),
    n_groups = ifelse(vc$grp %in% names(ng), ng[match(vc$grp, names(ng))], NA_integer_),
    note = if (has_resid) "" else "no residual variance (non-Gaussian family): pct is of RE variance only",
    stringsAsFactors = FALSE)
}

#' Print the random-effect variance block into a sink()'d stats log.
#'
#' Emits the variance and SD of every component, each group's share of total variance, and -- when
#' the fit is singular -- an explicit warning that the model does not in fact adjust for that
#' grouping factor. Safe to call on a NULL model.
#'
#' @param m      A merMod, or NULL.
#' @param label  Heading text.
#' @param indent Leading whitespace for the body lines.
cat_re_variance <- function(m, label = "Random effects", indent = "  ") {
  cat("\n=== ", label, " ===\n", sep = "")
  s <- re_variance_summary(m)
  if (all(is.na(s$vcov))) {
    cat(indent, "not available: ", s$note[1], "\n", sep = "")
    return(invisible(s))
  }
  for (i in seq_len(nrow(s))) {
    # var1 is already parenthesised for an intercept term, so do not add another pair
    trm <- s$term[i]
    grp <- if (!nzchar(trm)) s$group[i]
           else if (startsWith(trm, "(")) paste(s$group[i], trm)
           else paste0(s$group[i], " (", trm, ")")
    cat(sprintf("%s%-34s variance %11.4g   SD %11.4g", indent, grp, s$vcov[i], s$sdcor[i]))
    if (is.finite(s$pct_of_total[i])) cat(sprintf("   %6.2f%% of total", s$pct_of_total[i]))
    if (!is.na(s$n_groups[i])) cat(sprintf("   n groups %d", s$n_groups[i]))
    cat("\n")
  }
  if (nzchar(s$note[1])) cat(indent, "note: ", s$note[1], "\n", sep = "")
  if (isTRUE(s$is_singular[1])) {
    cat(indent, "*** SINGULAR FIT (lme4::isSingular, tol = 1e-4) ***\n", sep = "")
    cat(indent, "  A variance component is estimated at (effectively) zero, so the fit is\n", sep = "")
    cat(indent, "  equivalent to pooling over that grouping factor.\n", sep = "")
  }
  invisible(s)
}

#' Total variance of a fitted merMod (all variance components, residual included).
#'
#' The denominator for a model-based Cohen's d, kept in one place so every caller pairs it with
#' the estimate from the same model.
re_total_sd <- function(m) {
  if (is.null(m) || !inherits(m, "merMod")) return(NA_real_)
  vc <- as.data.frame(lme4::VarCorr(m))
  sqrt(sum(vc$vcov[is.na(vc$var2)], na.rm = TRUE))
}
