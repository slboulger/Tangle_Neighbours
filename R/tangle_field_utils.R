#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# tangle_field_utils.R
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
# tangle_field_utils.R
#
# Shared helper: single source of truth for the EDGE-CORRECTED KERNEL-WEIGHTED
# TANGLE FIELD -- a smoothed intensity surface of PHF1+ ("tangle-bearing")
# neurons, evaluated at every tangle-free cell coordinate.
#
# Deliberately mirrors the design of R/phf1_distance_utils.r, R/nn3_utils.r and
# R/kernel_field_utils.R -- the metric lives in ONE place so it cannot drift
# between the observed analysis and the permutation null built on top of it.
#
# ---------------------------------------------------------------------------
# WHY NOT kernel_field_rann() -- THE EXISTING KERNEL IS A DIFFERENT QUANTITY
# ---------------------------------------------------------------------------
# R/kernel_field_utils.R already computes sum_j exp(-d_ij^2 / (2 sigma^2)) via a
# RANN radius search, and that is the right tool for the GLIAL STATE field it was
# written for. It cannot serve here, for two reasons that are not stylistic:
#
#   1. IT IS NOT AN INTENSITY. It is an unnormalised kernel mass: a dimensionless
#      weighted count, not points per unit area. It is missing the 1/(2 pi sigma^2)
#      normalisation, so its magnitude is not comparable across sigma and it is not
#      the quantity an intensity surface ratio D/N is defined on.
#   2. IT HAS NO EDGE CORRECTION. Tangles outside the imaged footprint contribute
#      nothing, so an uncorrected sum is biased LOW exactly at tissue and FOV
#      borders. That is not a small effect here: measured on this cohort's geometry
#      the median edge-corrected / uncorrected ratio runs 1.00 (sigma = 25 um) ->
#      1.18 (100) -> 3.99 (500), max 6.9, and integral.im recovers 13.97 of 14
#      points with edge = TRUE against 10.92 without -- a 22% deficit.
#
# Border cells assigned a spuriously low tangle burden would produce a monotone
# gradient from edge geometry alone; edge correction removes that mechanism.
#
# What IS reused from kernel_field_utils.R: vif_from_design() and within_donor_z().
# Sourced, not copied.
#
# ---------------------------------------------------------------------------
# densityfun(), NEVER density.ppp() + interp.im()
# ---------------------------------------------------------------------------
# The field must be evaluated at arbitrary cell coordinates, which are not the
# anchor locations, so density.ppp(at = "points") does not apply. The obvious
# route -- rasterise with density.ppp() then interpolate with interp.im() -- is
# WRONG here and quietly so: interp.im() returns NA for a point whose bilinear
# stencil touches an outside-mask pixel, i.e. for cells near the mask boundary.
# Measured on a test fixture: 4 of 83 query points lost, all of them at the edge
# -- exactly the cells edge correction exists to fix. Silently dropping them
# reintroduces the bias by a different door.
#
# densityfun() evaluates the same edge-corrected estimator at arbitrary
# coordinates with no grid and no interpolation. On the same fixture it returned
# zero NAs and agreed with interp.im() at r = 0.9996. It is also not slower:
# ~0.1 s per donor per sigma for 4,300 query points.
#
# ---------------------------------------------------------------------------
# THE OBSERVATION WINDOW IS THE FOV UNION, NOT A BOUNDING BOX OR A HULL
# ---------------------------------------------------------------------------
# CosMx FOVs here are ~0.504 mm tiles on a ~0.51 mm pitch placed along the
# cortical ribbon; they cover roughly two thirds of the sample bounding box in
# offset strips with genuine interior holes (R/nn3_utils.r:192-197). Edge
# correction divides by the kernel mass falling inside the window, so a window
# that claims the gaps are imaged under-corrects across them and produces
# spurious low-intensity bands. A convex hull cannot represent holes at all
# (R/nn3_utils.r:321-330).
#
# fov_rectangles() is factored out of fov_edge_dist_um() (R/nn3_utils.r:252-271)
# and keeps its key property: the pixel -> mm affine is FITTED PER FOV off the
# data, so no nominal FOV size or um/px constant is assumed.
#
# ---------------------------------------------------------------------------
# THE RASTER RESOLUTION *IS* THE ACCURACY OF THE EDGE CORRECTION
# ---------------------------------------------------------------------------
# Because the window is a MASK, spatstat computes the edge-correction denominator
# e(u) = integral over W of K_sigma(u - v) dv on that mask's own grid, so a
# coarser raster means a coarser correction. Measured against an independent
# midpoint-rule reference, relative error at sigma = 150 um:
#
#   raster (um)    20      10       5       2        1
#   max rel      0.077   0.022   0.015   0.005    0.002
#   median rel   0.010   0.004   0.003   0.001    0.0004
#
# i.e. clean O(h) convergence -- the two implementations agree in the limit.
#
# TF_RASTER_UM = 10 is the default because 0.4% median error is negligible against
# the 1.0-1.5x edge factors these sections actually show, and because the raster
# drives the cost of the permutation null: one sigma x 9 donors x 1000 replicates
# runs ~54 min at 10 um against ~3 h at 5 um (measured on donor BBN00628710, the largest
# at 34 mm^2). The null MUST use the SAME raster as the observed run -- a different
# raster is a different transform and the coefficients would not be comparable --
# so this is one dial for both. --raster_um 5 is available for a sensitivity run.
#
# ---------------------------------------------------------------------------
# IDENTIFIABILITY OF SMALL-SIGMA FIELDS
# ---------------------------------------------------------------------------
# There are 398 PHF1+ pool-neuron anchors across the 9 donors, ranging 4
# (BBN00428931) to 120 (BBN_9931), over 10-34 mm^2 imaged each: densities of
# 0.29-6.72 anchors/mm^2, median 3.24. The expected number of anchors inside the
# kernel's effective area 2 pi sigma^2 is therefore:
#
#   sigma (um)        50    100    200    350    500    750   1000
#   median donor    0.05   0.20   0.81   2.49   5.09   11.4   20.4
#   donors with >=3    0      0      0      3      5      7      8
#
# No donor reaches 3 expected anchors below sigma = 350 um. The sigma needed for
# >=3 is ~384 um at the median donor, 890 um in BBN_9928 and 1288 um in BBN00428931;
# the two sparsest donors (4 and 6 anchors) set the floor.
#
# With ~26,400 cells a noisy exposure is ATTENUATED, not biased, but at
# sigma <= 100 um the field approaches a smoothed "is there a tangle within
# ~2 sigma" indicator, i.e. close to a thresholded nearest distance.
# field_diagnostics() measures all of this per (donor, sigma) so the gate is
# data-driven rather than asserted here.
#
# ---------------------------------------------------------------------------
# CONVENTIONS
# ---------------------------------------------------------------------------
# Coordinates are x_slide_mm / y_slide_mm (slide-global centroids, MM). They are
# slide-global and NOT sample-global -- two samples can share a slide -- so every
# computation here is done WITHIN one sample by the caller, exactly as
# compute_dist_to_phf1_um(), compute_nn3_um() and kernel_field_rann() do.
# sigma is always in MICRONS; the returned intensity is always per mm^2;
# distances are microns (mm * 1000).

suppressPackageStartupMessages({
  library(cli)
})

## SPATSTAT AVAILABILITY -- NOTHING NEEDS INSTALLING.
## Only the DOTTED sub-packages are used, never library(spatstat): the umbrella
## package is not installed in this project's environments, so library(spatstat)
## would fail.
##   local R 4.4      spatstat.geom 3.8.2, spatstat.explore 3.8.2
##   HPC dgenv        spatstat.geom 3.6-0, spatstat.explore 3.5-3
## On the HPC they arrive as SEURAT TRANSITIVE DEPENDENCIES rather than a
## deliberate pin, so the version follows Seurat. If a version must change,
## install with upgrade = "never" so other packages are not upgraded alongside.
##
## The one API detail this file depends on is densityfun()'s `drop` argument;
## tangle_field() asserts the returned length, so a version change that removed it
## would stop rather than silently corrupt a join.
for (p in c("spatstat.geom", "spatstat.explore", "RANN")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("Package '", p, "' is required by tangle_field_utils.R. Use the dotted ",
         "sub-packages (spatstat.geom / spatstat.explore); the umbrella 'spatstat' ",
         "package is not installed in this project and library(spatstat) will fail.")
}

## Nominal FOV side in pixels. Same value as R/nn3_utils.r:226. Used only to
## reconstruct the rectangle from the FITTED affine and as an area cross-check --
## the mm/px scale itself is never assumed.
TF_FOV_PX      <- 4256L
## Nominal FOV side in mm, for the area assertion only. 4256 px * 0.120281 um/px.
TF_FOV_MM      <- 0.51192
## Default raster for the window mask. 10 um matches fov_edge_dist_um().
TF_RASTER_UM   <- 10
## A (donor, sigma) cell is "not resolved" above this floor fraction.
TF_MAX_FLOOR_FRAC <- 0.20
## "geometry dominated" thresholds -- see field_diagnostics().
TF_MAX_EDGE_FACTOR      <- 5
TF_MIN_WITHIN_VAR_SHARE <- 0.25


## ---------------------------------------------------------------------------
## robust_slope()
##
## VERBATIM from R/nn3_utils.r:229-241. Copied rather than sourced because
## nn3_utils.r is a library but fov_rectangles() below is a factoring-out of
## fov_edge_dist_um()'s internals, and the two must not drift: the rectangle this
## produces has to be the same rectangle the edge-distance covariate is measured
## against.
##
## Robust mm-per-pixel slope, estimated per sample from every FOV with enough
## pixel spread, then median-pooled. Sign carries the axis orientation (y is
## flipped between FOV-local pixels and slide mm).
## ---------------------------------------------------------------------------
robust_slope <- function(px, mm, grp) {
  sl <- vapply(split(seq_along(px), grp), function(i) {
    if (length(i) < 20L) return(NA_real_)
    vx <- var(px[i])
    if (!is.finite(vx) || vx < 1e6) return(NA_real_)   # need a real pixel spread
    unname(stats::coef(stats::lm(mm[i] ~ px[i]))[2])
  }, numeric(1))
  sl <- sl[is.finite(sl)]
  if (!length(sl)) stop("Could not estimate a pixel->mm slope for a sample.")
  stats::median(sl)
}


## ---------------------------------------------------------------------------
## fov_rectangles()
##
## Per-FOV imaged rectangle in mm, for ONE sample's rows.
##
##   cd  : data.frame for a single sample, with fov_col, coord_x/y (mm) and
##         px_x/px_y (FOV-local pixels).
##
## Method, from fov_edge_dist_um() (R/nn3_utils.r:252-271): within a FOV,
## x_slide_mm is affine in x_FOV_px, so a robust per-sample slope plus a per-FOV
## MEDIAN intercept maps pixel 0 and pixel fov_px - 1 onto the rectangle edges.
## The intercept is a median so a handful of stray centroids cannot shift a
## rectangle.
##
## Returns data.frame(fov, x0, x1, y0, y1) in mm, one row per FOV.
## ---------------------------------------------------------------------------
fov_rectangles <- function(cd, fov_col = "fov",
                           coord_x = "x_slide_mm", coord_y = "y_slide_mm",
                           px_x = "x_FOV_px", px_y = "y_FOV_px",
                           fov_px = TF_FOV_PX) {

  for (cl in c(fov_col, coord_x, coord_y, px_x, px_y))
    if (!cl %in% colnames(cd)) stop("fov_rectangles(): missing column '", cl, "'")

  fv <- as.character(cd[[fov_col]])
  sx <- robust_slope(cd[[px_x]], cd[[coord_x]], fv)
  sy <- robust_slope(cd[[px_y]], cd[[coord_y]], fv)

  ix <- vapply(split(seq_len(nrow(cd)), fv),
               function(i) stats::median(cd[[coord_x]][i] - sx * cd[[px_x]][i]), numeric(1))
  iy <- vapply(split(seq_len(nrow(cd)), fv),
               function(i) stats::median(cd[[coord_y]][i] - sy * cd[[px_y]][i]), numeric(1))

  ex <- cbind(ix, ix + sx * (fov_px - 1))
  ey <- cbind(iy, iy + sy * (fov_px - 1))

  data.frame(
    fov = names(ix),
    x0  = pmin(ex[, 1], ex[, 2]), x1 = pmax(ex[, 1], ex[, 2]),
    y0  = pmin(ey[, 1], ey[, 2]), y1 = pmax(ey[, 1], ey[, 2]),
    stringsAsFactors = FALSE, row.names = NULL
  )
}


## ---------------------------------------------------------------------------
## fov_owin()
##
## Rasterise the union of per-FOV rectangles into a spatstat owin MASK.
##
## Three implementation details that are load-bearing:
##
##  1. PIXEL-CENTRE CONTAINMENT, not floor/ceiling expansion. fov_edge_dist_um()
##     expands each rectangle outward to whole raster pixels, which is right for
##     its purpose (a boundary set must not fall inside the tissue) but overstates
##     AREA by ~2.5%, and area is the edge-correction denominator. Centre
##     containment measured 0.9994 of the analytic area on a 68-FOV fixture.
##
##  2. owin(mask = ) IS INDEXED [row = y, col = x]. The occupancy matrix here is
##     built [x, y] to match fov_edge_dist_um()'s `occ`, so it is TRANSPOSED on
##     the way in. Getting this wrong yields a window that is the transpose of
##     the tissue -- which for these elongated cortical strips is a plausible-
##     looking but completely wrong mask.
##
##  3. The area is ASSERTED against n_fov * TF_FOV_MM^2. FOVs abut on a pitch
##     equal to their side, so overlap is possible but small; the tolerance is
##     one-sided-generous (the union can be SMALLER than the sum if FOVs overlap,
##     but never larger). A gross failure means the affine fit collapsed.
##
## Returns a spatstat.geom owin.
## ---------------------------------------------------------------------------
fov_owin <- function(rect, raster_um = TF_RASTER_UM, area_tol = 0.05,
                     verbose = FALSE) {

  stopifnot(nrow(rect) > 0)
  raster_mm <- raster_um / 1000
  pad <- 2L

  ox <- min(rect$x0) - pad * raster_mm
  oy <- min(rect$y0) - pad * raster_mm
  nx <- as.integer(ceiling((max(rect$x1) - ox) / raster_mm)) + pad
  ny <- as.integer(ceiling((max(rect$y1) - oy) / raster_mm)) + pad

  xcol <- ox + (seq_len(nx) - 0.5) * raster_mm
  yrow <- oy + (seq_len(ny) - 0.5) * raster_mm

  occ <- matrix(FALSE, nrow = nx, ncol = ny)
  for (r in seq_len(nrow(rect))) {
    ii <- which(xcol >= rect$x0[r] & xcol <= rect$x1[r])
    jj <- which(yrow >= rect$y0[r] & yrow <= rect$y1[r])
    if (length(ii) && length(jj)) occ[ii, jj] <- TRUE
  }
  if (!any(occ)) stop("fov_owin(): the rasterised FOV union is empty.")

  w <- spatstat.geom::owin(
    xrange = range(xcol) + c(-1, 1) * raster_mm / 2,
    yrange = range(yrow) + c(-1, 1) * raster_mm / 2,
    mask   = t(occ),                     # see note 2
    xy     = list(x = xcol, y = yrow)
  )

  a_obs <- spatstat.geom::area(w)
  a_nom <- nrow(rect) * TF_FOV_MM^2
  if (!is.finite(a_obs) || a_obs <= 0)
    stop("fov_owin(): non-positive window area.")
  if (a_obs > a_nom * (1 + area_tol))
    stop(sprintf(paste0("fov_owin(): window area %.3f mm2 exceeds %d FOVs x %.5f mm2 ",
                        "= %.3f mm2 by more than %.0f%%. The pixel->mm affine has ",
                        "probably collapsed; do not proceed."),
                 a_obs, nrow(rect), TF_FOV_MM^2, a_nom, 100 * area_tol))
  if (verbose)
    cli_alert_info(paste0("owin: {nrow(rect)} FOVs, area {round(a_obs, 3)} mm2 ",
                          "({round(100 * a_obs / a_nom)}% of nominal)"))
  w
}


## ---------------------------------------------------------------------------
## snap_to_window()
##
## Move cells that fall marginally OUTSIDE the rasterised window onto the nearest
## imaged pixel centre, and raise if any cell is further out than `tol_um`.
##
## WHY THIS EXISTS -- SILENT ROW LOSS FROM densityfun()'s DEFAULT, MEASURED
## The function densityfun() returns DROPS points outside the window by default
## (drop = TRUE) rather than returning NA for them. On sample BBN_9931 that gives
## 4,805 values for 4,807 CBLN2 query cells: the two cells sitting within one
## 10 um raster pixel of the mask edge simply vanished. A short vector joined back
## to a cell table by position does not error -- it misaligns every row after the
## first dropped one, or is silently recycled. This is the same class of failure
## as interp.im() returning NA at the boundary (see the header), and it strikes
## the same cells: the ones at the tissue edge, which is where a spurious gradient
## would come from.
##
## The two defences are both mandatory and neither is sufficient alone:
##   * tangle_field() always evaluates with drop = FALSE and hard-asserts the
##     returned length, so loss can never be silent.
##   * this function removes the cause, ONCE PER SAMPLE, so the same cell gets the
##     same coordinate at every sigma and in every null replicate. Snapping
##     per-sigma would make the exposure depend on the bandwidth through the
##     coordinates as well as the kernel.
##
## A centroid half a raster pixel outside the mask is a rasterisation artefact,
## not a cell outside the tissue, so snapping it (a <= tol_um shift, and only for
## cells already at the boundary) is the right repair. A cell further out than
## tol_um means the coordinates and the FOV pixel columns disagree -- a
## coordinate-frame error -- and that is an error, not something to snap away.
##
## Returns list(x, y, n_snapped, max_snap_um, inside).
## ---------------------------------------------------------------------------
snap_to_window <- function(x_mm, y_mm, w, tol_um = 15, what = "cells",
                           sample_id = NA_character_, verbose = TRUE) {

  x_mm <- as.numeric(x_mm); y_mm <- as.numeric(y_mm)
  stopifnot(length(x_mm) == length(y_mm))
  if (!length(x_mm))
    return(list(x = x_mm, y = y_mm, n_snapped = 0L, max_snap_um = 0,
                inside = logical(0)))

  ok <- spatstat.geom::inside.owin(x_mm, y_mm, w)
  if (all(ok))
    return(list(x = x_mm, y = y_mm, n_snapped = 0L, max_snap_um = 0, inside = ok))

  ## Pixel centres of the imaged part of the mask. w$m is [row = y, col = x].
  m  <- w$m
  jj <- which(m, arr.ind = TRUE)                     # col 1 = row (y), col 2 = col (x)
  ctr <- cbind(w$xcol[jj[, 2]], w$yrow[jj[, 1]])

  bad <- which(!ok)
  nn  <- RANN::nn2(data = ctr, query = cbind(x_mm[bad], y_mm[bad]), k = 1)
  snap_um <- nn$nn.dists[, 1] * 1000

  if (max(snap_um) > tol_um)
    stop(sprintf(paste0("snap_to_window(): %d of %d %s in sample %s fall outside the ",
                        "imaged FOV window, the worst by %.1f um (tolerance %.1f um). ",
                        "That is a coordinate-frame error, not a rasterisation one -- ",
                        "the coordinates and the FOV pixel columns disagree. Do not ",
                        "proceed."),
                 length(bad), length(x_mm), what, as.character(sample_id),
                 max(snap_um), tol_um))

  x_mm[bad] <- ctr[nn$nn.idx[, 1], 1]
  y_mm[bad] <- ctr[nn$nn.idx[, 1], 2]

  if (verbose)
    cli_alert_info(paste0("{length(bad)} of {length(ok)} {what} in sample ",
                          "{sample_id} snapped onto the window edge ",
                          "(max {round(max(snap_um), 1)} um)."))

  list(x = x_mm, y = y_mm, n_snapped = length(bad),
       max_snap_um = max(snap_um), inside = ok)
}


## ---------------------------------------------------------------------------
## tangle_field()
##
## Edge-corrected kernel intensity of `anchor` evaluated at `query`, per mm^2.
##
##   anchor_xy_mm : n_anchor x 2 (mm). The tangle-bearing neurons (or, under the
##                  null, the permuted anchor set).
##   query_xy_mm  : n_query x 2 (mm). The cells to evaluate at.
##   w            : owin from fov_owin().
##   sigma_um     : kernel sigma, MICRONS.
##
## Returns list(D = edge-corrected intensity per mm^2,
##              D_uncorr = the same without edge correction,
##              n_anchor, sigma_um).
## D_uncorr is returned so the per-cell edge factor D / D_uncorr is available as a
## diagnostic -- at large sigma the correction does most of the work and that has
## to be visible rather than buried inside the estimator.
##
## An empty anchor set returns all-zero (a section with no tangles has zero tangle
## intensity -- that is a real zero, not a missing observation, the same
## convention kernel_field_utils.R uses for its density).
##
## drop = FALSE IS NOT OPTIONAL. densityfun()'s returned function DROPS points
## outside the window by default, so f(x, y) can be SHORTER than the query and
## every downstream row silently misaligns -- measured at 4,805 returned for 4,807
## query cells on BBN_9931. With drop = FALSE those points come back NA instead,
## the length is guaranteed, and the assert below turns a would-be silent
## corruption into a stop(). Run snap_to_window() first to remove the cause.
## ---------------------------------------------------------------------------
tangle_field <- function(anchor_xy_mm, query_xy_mm, w, sigma_um) {

  anchor_xy_mm <- if (is.null(anchor_xy_mm)) matrix(numeric(0), 0, 2) else as.matrix(anchor_xy_mm)
  query_xy_mm  <- as.matrix(query_xy_mm)
  n_a <- nrow(anchor_xy_mm); n_q <- nrow(query_xy_mm)

  stopifnot("sigma_um must be a single positive number" =
              length(sigma_um) == 1L && is.finite(sigma_um) && sigma_um > 0)

  if (n_q == 0L)
    return(list(D = numeric(0), D_uncorr = numeric(0), n_anchor = n_a, sigma_um = sigma_um))
  if (n_a == 0L)
    return(list(D = rep(0, n_q), D_uncorr = rep(0, n_q), n_anchor = 0L, sigma_um = sigma_um))

  sigma_mm <- sigma_um / 1000
  T_ppp <- spatstat.geom::ppp(anchor_xy_mm[, 1], anchor_xy_mm[, 2], window = w,
                              check = FALSE)

  ## "data contain duplicated points" is EXPECTED and is deliberately muffled.
  ## A handful of cells share a centroid exactly with another cell -- the
  ## coincident-centroid segmentation artefact that compute_nn3_um() reports as
  ## n_coincident (R/nn3_utils.r:99-102). For a kernel density this is not a
  ## problem needing a fix: two separately called PHF1+ cells at the same
  ## coordinate legitimately contribute two units of tangle intensity there, and
  ## unlike a nearest-neighbour metric there is no self-match to get wrong. It is
  ## COUNTED in field_diagnostics() as n_anchor_dup rather than silenced outright,
  ## so a sudden rise in duplicates cannot hide behind a muffled warning.
  quiet_dup <- function(expr) {
    withCallingHandlers(expr, warning = function(v) {
      if (grepl("duplicated points", conditionMessage(v))) invokeRestart("muffleWarning")
    })
  }
  f_on  <- quiet_dup(spatstat.explore::densityfun(T_ppp, sigma = sigma_mm, edge = TRUE))
  f_off <- quiet_dup(spatstat.explore::densityfun(T_ppp, sigma = sigma_mm, edge = FALSE))

  D  <- as.numeric(f_on(query_xy_mm[, 1], query_xy_mm[, 2], drop = FALSE))
  Du <- as.numeric(f_off(query_xy_mm[, 1], query_xy_mm[, 2], drop = FALSE))

  ## Length guard first: this is the assert that makes silent row loss impossible.
  if (length(D) != n_q || length(Du) != n_q)
    stop(sprintf(paste0("tangle_field(): densityfun returned %d/%d values for %d query ",
                        "points at sigma = %g um. drop = FALSE should make this ",
                        "impossible; the spatstat API has changed."),
                 length(D), length(Du), n_q, sigma_um))
  if (anyNA(D) || anyNA(Du))
    stop(sprintf(paste0("tangle_field(): %d query point(s) lie outside the window at ",
                        "sigma = %g um. Run snap_to_window() on the query and anchor ",
                        "coordinates first -- do NOT drop them, they are edge cells and ",
                        "dropping them reintroduces exactly the edge bias the correction ",
                        "exists to remove."),
                 sum(is.na(D) | is.na(Du)), sigma_um))
  if (any(D < 0) || any(Du < 0))
    stop("tangle_field(): negative intensity returned.")

  list(D = D, D_uncorr = Du, n_anchor = n_a, sigma_um = sigma_um)
}


## ---------------------------------------------------------------------------
## dist_to_phf1_k_um()
##
## Distance to the k-th nearest anchor, in microns. The k-generalisation of
## compute_dist_to_phf1_um() (R/phf1_distance_utils.r:56), which is k = 1 only.
##
## Written here rather than reusing compute_nn3_um(): that function measures
## NEURON SPACING within a pool (query set == neighbour set == the pool), which is
## a different quantity from the distance to a k-th ANCHOR. Substituting it would
## silently change the exposure.
##
## k = 1 reproduces compute_dist_to_phf1_um() exactly, which is what pins this to
## the canonical metric. Same conventions: coordinates in mm, distances
## mm * 1000 -> um.
##
##   anchor_xy_mm : n_anchor x 2 (mm)
##   query_xy_mm  : n_query x 2 (mm)
##   k            : neighbour rank(s), a vector
##
## Returns an n_query x length(k) matrix, columns named "d1", "d3", ...
## NA where the sample has fewer than k anchors.
## ---------------------------------------------------------------------------
dist_to_phf1_k_um <- function(anchor_xy_mm, query_xy_mm, k = c(1L, 3L, 5L)) {

  anchor_xy_mm <- if (is.null(anchor_xy_mm)) matrix(numeric(0), 0, 2) else as.matrix(anchor_xy_mm)
  query_xy_mm  <- as.matrix(query_xy_mm)
  k <- as.integer(k)
  stopifnot(all(k >= 1L))

  n_a <- nrow(anchor_xy_mm); n_q <- nrow(query_xy_mm)
  out <- matrix(NA_real_, nrow = n_q, ncol = length(k),
                dimnames = list(NULL, paste0("d", k)))
  if (n_q == 0L || n_a == 0L) return(out)

  kk <- min(max(k), n_a)
  nn <- RANN::nn2(data = anchor_xy_mm, query = query_xy_mm, k = kk)
  for (j in seq_along(k)) {
    if (k[j] <= kk) out[, j] <- nn$nn.dists[, k[j]] * 1000
  }
  out
}


## ---------------------------------------------------------------------------
## field_transform()
##
## log(D + eps), then divide by SD. NOT centred, matching the canonical distance
## transform (docs/MODELS.md, Distance transform).
##
##   eps   : NULL -> half the minimum NON-ZERO value of D. Passed explicitly when
##           the null must use the observed run's eps (a null replicate with a
##           different eps is a different transform and its coefficients are not
##           comparable).
##   scale : "celltype" divides by ONE SD over all cells -- the default, and the
##           only choice that makes beta comparable to the existing
##           per-SD-log-distance coefficients, which use a single within-celltype
##           SD. "donor" divides by each donor's own SD.
##
## WHY "celltype" IS THE DEFAULT. Per-donor scaling is actively harmful on this
## cohort: BBN00428931 has 4 anchors and BBN_9928 has 6, so their within-donor field is
## nearly noise, and rescaling it to unit variance would present that noise as
## real exposure spread with the same weight as BBN_9931's 118 anchors. "donor" is
## available as a sensitivity, not as the headline.
##
## Returns list(z = transformed vector, eps = eps used, sd = SD(s) used).
## ---------------------------------------------------------------------------
field_transform <- function(D, donor = NULL, eps = NULL,
                            scale = c("celltype", "donor")) {
  scale <- match.arg(scale)
  D <- as.numeric(D)
  if (any(!is.finite(D))) stop("field_transform(): non-finite field value.")
  if (any(D < 0))         stop("field_transform(): negative field value.")

  if (is.null(eps)) {
    pos <- D[D > 0]
    if (!length(pos))
      stop("field_transform(): the field is zero everywhere; eps is undefined.")
    eps <- min(pos) / 2
  }
  stopifnot("eps must be a single positive number" =
              length(eps) == 1L && is.finite(eps) && eps > 0)

  lg <- log(D + eps)

  if (scale == "celltype") {
    s <- stats::sd(lg)
    if (!is.finite(s) || s == 0)
      stop("field_transform(): the transformed field has zero SD; nothing to scale.")
    return(list(z = lg / s, eps = eps, sd = s))
  }

  if (is.null(donor)) stop("field_transform(scale = 'donor'): `donor` is required.")
  donor <- as.character(donor)
  z  <- rep(NA_real_, length(lg))
  sd_by <- setNames(rep(NA_real_, length(unique(donor))), unique(donor))
  for (d in unique(donor)) {
    ii <- which(donor == d)
    s  <- stats::sd(lg[ii])
    sd_by[d] <- s
    ## A donor with a constant field carries no within-donor information; 0 rather
    ## than NaN so one degenerate donor cannot delete a whole sample from the fit
    ## (same convention as within_donor_z() in kernel_field_utils.R).
    z[ii] <- if (!is.finite(s) || s == 0) 0 else lg[ii] / s
  }
  list(z = z, eps = eps, sd = sd_by)
}


## ---------------------------------------------------------------------------
## field_diagnostics()
##
## The diagnostic check for a bandwidth, per (donor, sigma). This is the gate that
## decides which sigma may be a HEADLINE winner; every sigma is still plotted,
## because the profile is the output.
##
## Five things are measured, and each answers a specific way the field can fail:
##
##   frac_at_floor          the field underflowed. At sigma = 25 um, 84% of cells
##                          sit on the numerical floor; the exposure is a binary
##                          indicator, not a continuous field. -> not_resolved
##   expected_anchors       2 pi sigma^2 * lambda, the anchor count the kernel
##                          actually averages over. Below ~3 the per-cell
##                          intensity is dominated by Poisson noise in the anchor
##                          count. Reported, not gated on, because attenuation is
##                          not invalidity -- but it is the number that says
##                          whether "local" is measurable at all.
##   edge_factor_median/max how much of the value is the edge correction rather
##                          than the data. -> geometry_dominated
##   within_donor_var_share within-donor variance of the TRANSFORMED field as a
##                          share of its total. The model carries (1|sample_id),
##                          so only within-donor variation identifies beta; once
##                          the kernel spans a large part of the section the field
##                          approaches a per-donor constant and there is nothing
##                          left to fit. -> geometry_dominated
##   cv, sd_log             plain spread, for the profile figure.
##
##   D, D_uncorr : per-cell vectors for ONE sigma, all donors
##   donor       : donor label per cell
##   area_mm2    : named numeric, imaged area per donor (from area(owin))
##   n_anchor    : named integer, anchors per donor
##
## Returns a data.frame, one row per donor plus one "ALL" row.
## ---------------------------------------------------------------------------
field_diagnostics <- function(D, D_uncorr, donor, sigma_um, area_mm2, n_anchor,
                              z = NULL, n_anchor_dup = NULL, label = "field",
                              max_floor_frac = TF_MAX_FLOOR_FRAC,
                              max_edge_factor = TF_MAX_EDGE_FACTOR,
                              min_within_var_share = TF_MIN_WITHIN_VAR_SHARE) {

  donor <- as.character(donor)
  stopifnot(length(D) == length(donor), length(D_uncorr) == length(D))

  ## Within-donor variance share of the transformed field. Computed on z when the
  ## caller supplies it (the quantity the model actually sees); otherwise on
  ## log(D + eps) with the default eps.
  if (is.null(z)) z <- field_transform(D)$z
  gm <- tapply(z, donor, mean)
  within_ss <- sum((z - gm[donor])^2)
  total_ss  <- sum((z - mean(z))^2)
  share_all <- if (total_ss > 0) within_ss / total_ss else NA_real_

  one <- function(ii, lab) {
    d  <- D[ii]; du <- D_uncorr[ii]
    pos <- d[d > 0]
    mn  <- if (length(pos)) min(pos) else NA_real_
    ef  <- if (length(du)) d[du > 0] / du[du > 0] else numeric(0)
    ## Anchor density. On the ALL row this is the COHORT density (total anchors
    ## over total imaged area), not NA -- otherwise the summary row silently loses
    ## expected_anchors, which is the number that says whether "local" is
    ## measurable at this bandwidth at all.
    lam <- if (lab %in% names(area_mm2) && lab %in% names(n_anchor)) {
      n_anchor[[lab]] / area_mm2[[lab]]
    } else {
      sum(n_anchor, na.rm = TRUE) / sum(area_mm2, na.rm = TRUE)
    }
    zz <- z[ii]
    data.frame(
      sample_id             = lab,
      sigma_um              = sigma_um,
      n                     = length(ii),
      n_anchor              = if (lab %in% names(n_anchor)) n_anchor[[lab]] else sum(n_anchor),
      ## Anchors sharing a centroid exactly -- the coincident-centroid
      ## segmentation artefact (R/nn3_utils.r:99-102). Harmless for a kernel
      ## density but counted so a rise cannot hide behind a muffled warning.
      n_anchor_dup          = if (is.null(n_anchor_dup)) NA_integer_
                              else if (lab %in% names(n_anchor_dup)) n_anchor_dup[[lab]]
                              else sum(n_anchor_dup),
      area_mm2              = if (lab %in% names(area_mm2)) area_mm2[[lab]] else sum(area_mm2),
      anchor_density_mm2    = lam,
      expected_anchors      = 2 * pi * (sigma_um / 1000)^2 * lam,
      min_pos               = mn,
      frac_at_floor         = if (is.na(mn)) 1 else mean(abs(d - mn) < 1e-6 * max(1, mn)),
      cv                    = if (mean(d) > 0) stats::sd(d) / mean(d) else NA_real_,
      sd_log                = stats::sd(zz),
      edge_factor_median    = if (length(ef)) stats::median(ef) else NA_real_,
      edge_factor_max       = if (length(ef)) max(ef) else NA_real_,
      within_donor_var_share = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  rows <- lapply(unique(donor), function(d) one(which(donor == d), d))
  all_row <- one(seq_along(D), "ALL")
  out <- do.call(rbind, c(rows, list(all_row)))

  ## The within-donor share is a cohort-level quantity; per-donor it is 1 by
  ## construction, so it is only meaningful on the ALL row.
  out$within_donor_var_share[out$sample_id == "ALL"] <- share_all

  out$not_resolved       <- out$frac_at_floor > max_floor_frac
  out$geometry_dominated <- (!is.na(out$edge_factor_median) &
                               out$edge_factor_median > max_edge_factor) |
    (out$sample_id == "ALL" & !is.na(share_all) & share_all < min_within_var_share)

  n_bad <- sum(out$not_resolved[out$sample_id != "ALL"])
  if (n_bad > 0)
    cli_alert_warning(paste0("[{label}] sigma = {sigma_um} um: {n_bad} donor(s) have >",
                             "{round(100 * max_floor_frac)}% of cells on the field ",
                             "floor -- the exposure is a binary indicator there, ",
                             "not a continuous field."))
  if (isTRUE(out$geometry_dominated[out$sample_id == "ALL"])) {
    ar <- out[out$sample_id == "ALL", ]
    cli_alert_warning(paste0("[{label}] sigma = {sigma_um} um: geometry-dominated ",
                             "(edge factor {round(ar$edge_factor_median, 3)}, ",
                             "within-donor variance share {round(share_all, 3)}). ",
                             "It cannot be a headline bandwidth."))
  }
  out$label <- label
  rownames(out) <- NULL
  out
}


## ---------------------------------------------------------------------------
## r2_marginal_lmm()
##
## Nakagawa & Schielzeth MARGINAL R^2: fixed-effects variance / total variance.
##
## IDENTICAL to r2m() in R/phf1_module_distance_intensity_adjust.R and
## R/phf1_intensity_gradient_neurons.R; copied here so this file is self-contained.
## ---------------------------------------------------------------------------
r2_marginal_lmm <- function(m) {
  if (is.null(m)) return(NA_real_)
  vf <- stats::var(as.numeric(stats::model.matrix(m) %*% lme4::fixef(m)))
  vc <- as.data.frame(lme4::VarCorr(m))
  vr <- sum(vc$vcov[vc$grp != "Residual"])
  ve <- vc$vcov[vc$grp == "Residual"]
  as.numeric(vf / (vf + vr + ve))
}


## ---------------------------------------------------------------------------
## loo_donor_r2()
##
## Leave-one-DONOR-out predictive R^2 -- the arbiter for comparing non-nested
## exposures. Fit on 8 donors, predict the 9th.
##
## TWO NUMBERS, AND THE SECOND IS THE ONE THAT MATTERS.
## A held-out donor has no estimable random intercept, so prediction must use
## re.form = NA (population level). That means the donor's unknown offset lands in
## the residual and dominates it -- with 9 donors the marginal R^2 is mostly a
## measure of how far that donor's mean sits from the cohort's, which is the same
## for every exposure and therefore cannot discriminate between them.
##
##   r2_pred_marginal  1 - SS_res / SS_tot about the TRAINING mean. Reported for
##                     completeness; do not arbitrate on it.
##   r2_pred_within    the same after removing the held-out donor's own mean from
##                     both y and yhat. This isolates the WITHIN-donor gradient,
##                     which is what the exposure is supposed to explain and what
##                     the (1|sample_id) model identifies beta from. THE ARBITER.
##
## POOLED BY SUMMING SS ACROSS FOLDS, not by averaging per-fold R^2: donors range
## from ~1,500 to ~4,900 CBLN2 cells, and a mean of ratios would let the smallest
## donor count as much as the largest. Per-fold values are returned as well, so
## the influence of any single donor is visible.
##
##   dat       : model frame, complete cases only
##   formula   : an lmer formula; must contain (1 | <donor_col>)
##   donor_col : donor column name
##
## A fold that fails to converge or cannot be fitted is SKIPPED with a message and
## counted, never silently imputed.
##
## Returns list(r2_pred_marginal, r2_pred_within, per_fold, n_folds_used,
##              n_folds_skipped).
## ---------------------------------------------------------------------------
loo_donor_r2 <- function(dat, formula, donor_col = "sample_id", verbose = FALSE) {

  stopifnot(donor_col %in% colnames(dat))
  y_col <- all.vars(formula)[1]
  donors <- unique(as.character(dat[[donor_col]]))

  ss_res_m <- 0; ss_tot_m <- 0
  ss_res_w <- 0; ss_tot_w <- 0
  per <- list(); skipped <- character(0)

  for (d in donors) {
    te <- which(as.character(dat[[donor_col]]) == d)
    tr <- setdiff(seq_len(nrow(dat)), te)
    if (!length(te) || !length(tr)) { skipped <- c(skipped, d); next }

    train <- droplevels(dat[tr, , drop = FALSE])
    test  <- dat[te, , drop = FALSE]

    ## A covariate that is constant in the training set (e.g. Sex, when the
    ## held-out donor is the only one of that sex) makes the fit rank-deficient.
    ## Skip the fold rather than silently changing the model for one fold.
    m <- tryCatch(
      suppressMessages(suppressWarnings(
        lme4::lmer(formula, data = train, REML = TRUE,
                   control = lme4::lmerControl(optimizer = "bobyqa",
                                               optCtrl = list(maxfun = 2e5)))
      )), error = function(e) NULL)
    if (is.null(m)) { skipped <- c(skipped, d); next }

    yhat <- tryCatch(
      as.numeric(stats::predict(m, newdata = test, re.form = NA,
                                allow.new.levels = TRUE)),
      error = function(e) NULL)
    if (is.null(yhat) || anyNA(yhat)) { skipped <- c(skipped, d); next }

    y  <- as.numeric(test[[y_col]])
    mu <- mean(as.numeric(train[[y_col]]))

    r_m <- sum((y - yhat)^2);        t_m <- sum((y - mu)^2)
    yc  <- y - mean(y); hc <- yhat - mean(yhat)
    r_w <- sum((yc - hc)^2);         t_w <- sum(yc^2)

    ss_res_m <- ss_res_m + r_m; ss_tot_m <- ss_tot_m + t_m
    ss_res_w <- ss_res_w + r_w; ss_tot_w <- ss_tot_w + t_w

    per[[d]] <- data.frame(
      held_out = d, n_cells = length(te),
      r2_marginal = 1 - r_m / t_m,
      r2_within   = if (t_w > 0) 1 - r_w / t_w else NA_real_,
      stringsAsFactors = FALSE)
    if (verbose)
      cat(sprintf("  fold %-10s n = %6d  r2_within = %+.4f\n",
                  d, length(te), 1 - r_w / t_w))
  }

  if (length(skipped))
    cli_alert_info("loo_donor_r2(): {length(skipped)} fold(s) skipped: {paste(skipped, collapse = ', ')}")

  list(
    r2_pred_marginal = if (ss_tot_m > 0) 1 - ss_res_m / ss_tot_m else NA_real_,
    r2_pred_within   = if (ss_tot_w > 0) 1 - ss_res_w / ss_tot_w else NA_real_,
    per_fold         = if (length(per)) do.call(rbind, per) else NULL,
    n_folds_used     = length(per),
    n_folds_skipped  = length(skipped)
  )
}


## ---------------------------------------------------------------------------
## Column naming. Single definition so the builder (stage 1) and the model
## scripts (stage 2) cannot disagree about the spelling. sigma is formatted
## without a decimal point where it is a whole number, so 50 -> "s50", matching
## the .fmt_sigma()/kernel_mean_col() idiom in R/kernel_field_utils.R:507-517.
## ---------------------------------------------------------------------------
.tf_fmt_sigma <- function(sigma_um) {
  gsub("\\.", "p", format(sigma_um, trim = TRUE, scientific = FALSE))
}

tf_abs_col      <- function(sigma_um) sprintf("field_abs_s%s",       .tf_fmt_sigma(sigma_um))
tf_abs_unc_col  <- function(sigma_um) sprintf("field_absunc_s%s",    .tf_fmt_sigma(sigma_um))
tf_den_col      <- function(sigma_um) sprintf("field_den_s%s",       .tf_fmt_sigma(sigma_um))
tf_den_unc_col  <- function(sigma_um) sprintf("field_denunc_s%s",    .tf_fmt_sigma(sigma_um))
tf_rel_col      <- function(sigma_um) sprintf("field_rel_s%s",       .tf_fmt_sigma(sigma_um))


## ---------------------------------------------------------------------------
## load_cell_table()
##
## ONE definition of how the per-cell table is read, shared by stage 1 (which
## needs coordinates + FOV pixels + PHF1 + celltype) and stage 2 (which needs the
## Set 1 covariates). Both stages must see the SAME cells with the SAME labels; a
## stage 1 built from one source and a stage 2 keyed to another would join on
## cell_id and silently model a different cohort.
##
## THE DEFAULT SOURCE IS sce_phf1_dist.qs, the canonical all-cell object. It carries
## every column both stages need in one place: cell_id, sample_id, celltype, PHF1,
## x_slide_mm, y_slide_mm, fov, x_FOV_px, y_FOV_px, dist_to_phf1_um, SlideLabel,
## nCount_RNA, percent.neg, Sex, Age, PMI, Braak.
##
## Expected counts on the published object: 601 PHF1+ cells of any celltype, 398 of
## them in NEURON_POOL_10 (CBLN2 188, Exc-IT-L3-5 71). Only the 398 neuron-pool cells
## are distance anchors.
##
## Paths here are relative to the working directory; no project root is hard-coded.
##
##   sce  : path to a SingleCellExperiment .qs (preferred)
##   csv  : path to a per-cell CSV, used only if `sce` is NULL/missing
##   need : required column names; missing ones raise
##
## PHF1 is normalised to the STRING "TRUE"/"FALSE" that phf1_source_idx() compares
## against, from either R logicals or the Python-style "True"/"False" the CSVs
## carry, and the resulting count is asserted -- a silent normalisation failure
## would make every anchor set empty and every field identically zero.
## ---------------------------------------------------------------------------
load_cell_table <- function(sce = NULL, csv = NULL, need = character(0),
                            min_phf1 = 100L, verbose = TRUE) {

  cd <- NULL; src <- NA_character_
  if (!is.null(sce) && nzchar(sce) && file.exists(sce)) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE))
      stop("SingleCellExperiment is required to read '", sce, "'.")
    if (!requireNamespace("qs", quietly = TRUE))
      stop("qs is required to read '", sce, "'.")
    s  <- qs::qread(sce)
    cd <- as.data.frame(SummarizedExperiment::colData(s))
    ## cell_id resolution idiom from R/add_dist_to_full_sce.R:41-44.
    if (!"cell_id" %in% colnames(cd)) cd$cell_id <- colnames(s)
    rm(s); invisible(gc())
    src <- sce
  } else if (!is.null(csv) && nzchar(csv) && file.exists(csv)) {
    cd  <- utils::read.csv(csv, stringsAsFactors = FALSE)
    src <- csv
  } else {
    stop("No per-cell table found. Tried sce = '", sce %||% "", "' and csv = '",
         csv %||% "", "'. Use sce_phf1_dist.qs: PHF1/tangle_cells.csv ",
         "is a derived intermediate and may not exist.")
  }

  cd$cell_id   <- as.character(cd$cell_id)
  cd$sample_id <- as.character(cd$sample_id)
  if ("celltype" %in% colnames(cd)) cd$celltype <- as.character(cd$celltype)
  if ("PHF1" %in% colnames(cd))
    cd$PHF1 <- ifelse(cd$PHF1 %in% c(TRUE, "TRUE", "True", "true"), "TRUE", "FALSE")
  for (nm in c("x_slide_mm", "y_slide_mm", "x_FOV_px", "y_FOV_px"))
    if (nm %in% colnames(cd)) cd[[nm]] <- as.numeric(cd[[nm]])

  miss <- setdiff(need, colnames(cd))
  if (length(miss))
    stop("Missing column(s) in ", src, ": ", paste(miss, collapse = ", "),
         "\nx_FOV_px/y_FOV_px are what make the FOV observation window ",
         "reconstructible; sce_phf1_dist.qs carries them.")

  if ("PHF1" %in% colnames(cd)) {
    n_pos <- sum(cd$PHF1 == "TRUE")
    if (n_pos < min_phf1)
      stop("Only ", n_pos, " PHF1-positive cells after normalising the PHF1 column ",
           "(expected ~601). The label ",
           "encoding in ", src, " has changed.")
    if (verbose)
      cli_alert_info("{basename(src)}: {nrow(cd)} cells, {n_pos} PHF1+, {length(unique(cd$sample_id))} samples")
  }
  attr(cd, "source") <- src
  cd
}
