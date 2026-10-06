#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# phf1_morphology_filters.R
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
# phf1_morphology_filters.R
#
# Shared constants and helpers for the CosMx morphology / channel-intensity family, so the
# analysis scripts, their figures and their label nulls cannot drift apart on the two things
# that silently change a result: the pixel scale and the nucleus-validity filter.
#
# Sourced by:
#   R/phf1_morphology_exc.R
#   R/phf1_distance_morphology_exc.R
#   R/phf1_morphology_dist_phf1_panel.R
#
# PIXEL SCALE. Area.um2 / Area = 0.014467506 for every cell, exactly, so 1 px = 0.120281 um.
# NucArea and Perimeter ship in pixels with no .um2 companion, unlike Area.
#
# NUCLEUS VALIDITY. Two distinct problems:
#   1. NucArea == 0  -- no nucleus was segmented at all. 26.1% of CBLN2 cells.
#   2. NucArea > 0 but negligible -- a segmentation fragment, not a nucleus. The smallest
#      non-zero value in this dataset is 4 px = 0.058 um^2, and NucArea is quantised in
#      4-px steps. 132 CBLN2 cells sit below 1 um^2 and 1020 below 10 um^2.
# A neuronal nucleus 4 um across has a cross-section of ~12.6 um^2; 10 um^2 corresponds to a
# 3.57 um equivalent diameter, below anything biologically plausible. NUC_MIN_UM2 is therefore
# a HARD floor on implausibility, not a distributional trim, and it is applied to the two
# NUCLEAR measures only (nucarea, nucaspect) -- the cell-level shape metrics are defined for
# cells with no segmented nucleus at all and must not inherit a nucleus filter.
#
# The filter is mildly PHF1-dependent, in the CONSERVATIVE direction: 5.20% of PHF1- vs 3.52% of
# PHF1+ nucleated cells fall below 10 um^2, so removing them slightly SHRINKS the PHF1+
# nucleus-area excess rather than manufacturing it.

PX_UM       <- 0.120281      # microns per pixel
PX2_UM2     <- 0.014467506   # square microns per square pixel (Area.um2 / Area)
NUC_MIN_UM2 <- 10            # hard floor for a plausible nucleus (~3.57 um diameter)
NUC_MIN_PX  <- NUC_MIN_UM2 / PX2_UM2   # ~691 px, for reporting in pixel terms

## Assert the pixel scale on the object rather than trusting the constant.
assert_px_scale <- function(cd) {
  err <- max(abs(cd$Area.um2 / cd$Area - PX2_UM2), na.rm = TRUE)
  if (!is.finite(err) || err > 1e-6)
    stop(sprintf(paste0("Area.um2/Area is not the constant %.9f (max deviation %.3g). The ",
                        "px->um conversion used for NucArea and Perimeter is invalid."),
                 PX2_UM2, err), call. = FALSE)
  invisible(err)
}

## TRUE where the cell has a nucleus large enough to be believable.
nucleus_ok <- function(nuc_area_px) {
  !is.na(nuc_area_px) & nuc_area_px > 0 & (nuc_area_px * PX2_UM2) >= NUC_MIN_UM2
}

## One-row-per-group audit of what the filter removes, for the stats logs.
nucleus_filter_audit <- function(nuc_area_px, group = NULL) {
  a <- nuc_area_px * PX2_UM2
  d <- data.frame(
    group      = if (is.null(group)) "all" else as.character(group),
    n          = 1L,
    n_zero     = as.integer(!is.na(nuc_area_px) & nuc_area_px == 0),
    n_tiny     = as.integer(!is.na(nuc_area_px) & nuc_area_px > 0 & a < NUC_MIN_UM2),
    n_kept     = as.integer(nucleus_ok(nuc_area_px))
  )
  out <- stats::aggregate(cbind(n, n_zero, n_tiny, n_kept) ~ group, data = d, FUN = sum)
  out$pct_zero <- 100 * out$n_zero / out$n
  out$pct_tiny <- 100 * out$n_tiny / out$n
  out$pct_kept <- 100 * out$n_kept / out$n
  out
}

## Standard block for a stats log, so the wording does not drift between scripts.
cat_nucleus_filter_note <- function(audit = NULL) {
  cat("\n*** NUCLEUS VALIDITY FILTER ***\n")
  cat(sprintf("  Nuclear measures require NucArea > 0 AND >= %g um^2 (%.0f px).\n",
              NUC_MIN_UM2, NUC_MIN_PX))
  cat("  Two separate problems: cells with NO segmented nucleus (NucArea == 0), and cells\n")
  cat("  with a non-zero but negligible one -- the smallest value in this dataset is 4 px =\n")
  cat("  0.058 um^2. A 4 um nucleus has a cross-section of ~12.6 um^2, so the floor is a\n")
  cat("  hard implausibility bound, not a distributional trim.\n")
  cat("  Applied to nucarea / nucaspect ONLY: the cell-level shape metrics are defined for\n")
  cat("  cells with no nucleus and must not inherit a nucleus filter.\n")
  cat("  The filter is mildly PHF1-dependent in the CONSERVATIVE direction (5.20% of\n")
  cat("  PHF1- vs 3.52% of PHF1+ nucleated cells fall below the floor), so it shrinks\n")
  cat("  rather than creates any PHF1+ nucleus-area excess.\n")
  if (!is.null(audit)) { cat("\n"); print(audit, row.names = FALSE, digits = 4) }
}
