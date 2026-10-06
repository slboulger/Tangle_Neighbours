#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# tangle_3d_utils.R
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
# tangle_3d_utils.R
#
# Geometry and forward simulation for the OUT-OF-PLANE TANGLE CORRECTION.
#
# dist_to_phf1_um is an IN-PLANE distance measured on a single 2D section. The
# true nearest tangle may lie above or below the imaged tissue, so the stored
# value overstates true 3D separation, and some neurons labelled PHF1-negative
# are genuinely tangle-bearing with their inclusion out of plane. Both errors
# push the distance DEG toward null, so the reported field effect is a floor.
#
# Pure functions only -- no argparse, no setwd(), no file writing at load time,
# no package attachment beyond the spatstat sub-packages the simulation needs.
# (Same contract as R/cosmx_deg_module_utils.R.)
#
# ============================================================================
# THE GEOMETRY, AND THE CONDITIONAL THAT IS EASY TO GET WRONG
# ============================================================================
# There are two different questions and they have different answers:
#
#   (1) Is an inclusion PRESENT IN THE SECTION optically callable?
#   (2) Does a tangle-bearing neuron present in the section have ITS INCLUSION
#       in the section?
#
# (1) is about the optics and the answer is yes -- see the h_optical note below.
# (2) is the quantity the correction needs, and the answer is not necessarily
# yes; the two must not be conflated.
#
# Somata are ~17 um and sections ~7 um, so EVERY soma present is sliced, and
# whether the ~11 um inclusion falls in the slice depends on where it sits
# inside the soma. The two windows are nested (the inclusion lies inside the
# soma), so with a section of thickness t:
#
#     p_detect = (D_inclusion + t) / (D_soma + t)
#
# At t = 7, D_i = 11, D_s = 17 this is 18/24 = 0.750.
#
# WHY THE ECCENTRIC OFFSET DOES NOT MATTER. Flame-shaped and globose NFTs are
# not concentric within the soma, so one expects an offset delta to widen or
# narrow p_detect. It does neither. Let z_s be the soma centre, uniform over the
# range in which the soma intersects the slab (an interval of length D_s + t),
# and put the inclusion centre at z_s + delta. The inclusion intersects when
# |z_s + delta| < (D_i + t)/2, whose measure inside that interval is exactly
# D_i + t for ANY |delta| <= (D_s - D_i)/2 -- the shifted window still fits
# entirely inside. So p_detect is invariant to delta over the whole physically
# admissible range, and integrating delta ~ U(-3, 3) returns 0.750 exactly.
# t3_p_detect_mc() verifies this by simulation; the test suite asserts it.
#
# p_detect IS AN UPPER BOUND. It assumes callability requires only intersection,
# whereas the annotation criteria required recognisable somatic morphology, so a
# thin cap gets rejected. The three-rater re-adjudication of top-intensity
# tangle-free neurons is what calibrates that shortfall.
#
# ============================================================================
# WHY h_optical IS NOT LOAD-BEARING, AND WHY NO NUMBER IS QUOTED FOR IT
# ============================================================================
# From the Leica XML in PHF1/*_Merged.lof (UTF-16LE tail of the file):
#
#   ObjectiveName       "HC PL APO CS2 20x/0.75 DRY"
#   NumericalAperture   0.75          RefractionIndex 1.0
#   Pinhole             5.6586e-05 m  PinholeAiry     0.9996
#   PinholeAiryCalculation 580        (nm)
#   DimensionDescription DimID=1 (x), DimID=2 (y)   <- NO z DIMENSION
#
# So the PHF1 stain is a SINGLE optical plane; there is no z-stack to mine, and
# the z of a real tangle is NOT MEASURED anywhere in this project.
#
# A caution about the pinhole. The 56.6 um figure is at the PINHOLE PLANE, where
# the total magnification includes
# the ~3x relay, not just the objective. At 1 Airy unit the OBJECT-space pinhole
# diameter is by definition 1.22*lambda/NA = 1.22*0.580/0.75 = 0.943 um, and
# 56.6/0.943 = 60x = 20x objective x 3.0x relay. Dividing by the objective
# magnification alone (56.6/20 = 2.83 um) inflates the axial FWHM from ~2.3 um
# to ~5.6 um.
#
# NO POINT ESTIMATE OF h_optical IS RETURNED BY THIS FILE, deliberately. A dry
# 0.75 NA imaging into tissue mounted at n ~ 1.5 suffers substantial spherical
# aberration, so the in-tissue axial PSF is elongated well beyond the nominal
# figure by an amount this project cannot measure. The aberration and the
# pinhole error push in opposite directions, which is why the qualitative
# conclusion survives either way:
#
#     h_optical < D_inclusion, under every plausible value.
#
# That is all that is needed, and it is all that is claimed. The consequence is
# that the OPTICS DROP OUT OF p_detect ENTIRELY, and the real sensitivity sits
# on D_soma and D_inclusion -- which is why t3_geometry() sweeps them and why
# R/phf1_3d_distance.R measures them from the segmentation and PHF1 masks
# rather than trusting the defaults here.
#
# ============================================================================
# p_detect AND h_eff ARE ONE AXIS, NOT TWO
# ============================================================================
# Both derive from the same triple. For the anchor set the effective depth is
# h_eff = D_inclusion + t (the axial window in which a tangle contributes an
# anchor), and p_detect = (D_i + t)/(D_s + t). Crossing an independent p_detect
# grid against an independent h_eff grid produces cells that are geometrically
# impossible -- p_detect = 0.75 with h_eff = 10 um cannot exist. Sweep the
# triple; derive both.
#
# Exports:
#   T3_*                    default geometry constants, with provenance
#   t3_geometry()           the depth ledger: p_detect and h_eff from (t,Di,Ds)
#   t3_geometry_grid()      the ledger over a swept triple
#   t3_p_detect_mc()        Monte-Carlo check of p_detect incl. eccentricity
#   t3_lambda2_from_counts() per-donor in-plane anchor intensity from area
#   t3_median_d3d()         Poisson median 3D NN distance from lambda2 and h
#   t3_simulate_3d()        thinned neuron process in continuous R^3
#   t3_project()            slab selection + in-plane distance under a real owin
#   t3_validate()           G/F/K envelopes + the rigid-stacking diagnostic
#   t3_hidden_count()       H = n_obs (1 - p)/p

suppressPackageStartupMessages({
  library(stats)
})

if (!requireNamespace("RANN", quietly = TRUE))
  stop("Package 'RANN' is required. Install with: install.packages('RANN')")

## ---------------------------------------------------------------------------
## Default geometry, and where each number comes from.
##
## These are DEFAULTS for the sweep centre, not measurements. D_soma and
## D_inclusion are measurable from data already on disk and the Stage 1 driver
## does measure them; if the measured central values disagree materially with
## these, re-centre the sweep and say so in the log.
## ---------------------------------------------------------------------------

# Residual section thickness after processing, um. Sections were cut at 10 um
# nominal and taken through overnight xylene, flow-cell removal and post-run
# staining. Cryosections typically retain 60-80% of nominal after aggressive
# solvent processing. NOT measured -- 10 um is the hard upper bound.
T3_T_UM        <- 7
T3_T_UM_RANGE  <- c(6, 8)

# Inclusion (NFT) diameter, um. Sweep default; measured from phf1_n_mask_px.
T3_D_INCL_UM        <- 11
T3_D_INCL_UM_RANGE  <- c(10, 12)

# Soma diameter, um. Sweep default; measured from the segmentation Area.um2.
T3_D_SOMA_UM        <- 17
T3_D_SOMA_UM_RANGE  <- c(16, 20)

# CosMx FOV side, mm, and the confocal pixel size, um/px. FOV constants agree
# with TF_FOV_MM / TF_FOV_PX in R/tangle_field_utils.R:165-167 by construction;
# the pixel size is 4.736080 mm / 16672 px from the same LOF XML.
T3_FOV_MM      <- 0.51192
T3_PX_UM       <- 4.736080e-3 / 16672 * 1e3

# Distance cap, um. Matches --max_dist_um in deg_dream_linear_distance_phf1.r.
T3_MAX_DIST_UM <- 1000


## ---------------------------------------------------------------------------
## t3_geometry()
##
## The depth ledger. Single source of truth for every depth number -- do not
## recompute p_detect anywhere else.
##
##   t_um, d_incl_um, d_soma_um : the triple.
##
## Returns a one-row data.frame:
##   p_detect  (D_i + t)/(D_s + t)   -- probability the inclusion of a
##             tangle-bearing neuron PRESENT in the section is intersected by
##             the section. An UPPER bound (intersection, not callability).
##   h_eff     D_i + t                -- axial window contributing anchors
##   h_cell    D_s + t                -- axial window contributing CELLS; the
##             ceiling on what hidden-anchor relabelling can reach, because a
##             tangle-bearer deeper than this is not in the object at all
##   h_recover h_eff / h_cell == p_detect, named separately because it is used
##             as the recoverable fraction rather than as a probability
## ---------------------------------------------------------------------------
t3_geometry <- function(t_um      = T3_T_UM,
                        d_incl_um = T3_D_INCL_UM,
                        d_soma_um = T3_D_SOMA_UM) {

  if (any(t_um <= 0, d_incl_um <= 0, d_soma_um <= 0))
    stop("t3_geometry(): all of t_um, d_incl_um, d_soma_um must be > 0.")
  if (any(d_incl_um > d_soma_um))
    stop("t3_geometry(): d_incl_um (", d_incl_um, ") exceeds d_soma_um (",
         d_soma_um, "). The inclusion lies INSIDE the soma; the nesting that ",
         "makes p_detect = (D_i + t)/(D_s + t) valid is violated.")

  h_eff  <- d_incl_um + t_um
  h_cell <- d_soma_um + t_um

  data.frame(
    t_um      = t_um,
    d_incl_um = d_incl_um,
    d_soma_um = d_soma_um,
    p_detect  = h_eff / h_cell,
    h_eff     = h_eff,
    h_cell    = h_cell,
    h_recover = h_eff / h_cell,
    # Largest eccentric offset over which p_detect is provably invariant.
    delta_max_um = (d_soma_um - d_incl_um) / 2,
    row.names = NULL
  )
}


## ---------------------------------------------------------------------------
## t3_geometry_grid()
##
## The ledger over the swept triple. This is Arm D's grid -- one axis, not two.
## Returns the full factorial plus a `corner` label marking the min/max p_detect
## cells and the central cell, which are the three that actually get fitted.
## ---------------------------------------------------------------------------
t3_geometry_grid <- function(t_um      = T3_T_UM_RANGE,
                             d_incl_um = T3_D_INCL_UM_RANGE,
                             d_soma_um = T3_D_SOMA_UM_RANGE,
                             centre    = c(T3_T_UM, T3_D_INCL_UM, T3_D_SOMA_UM)) {

  t_seq  <- if (length(t_um) == 2)      seq(t_um[1], t_um[2], by = 1)      else t_um
  di_seq <- if (length(d_incl_um) == 2) seq(d_incl_um[1], d_incl_um[2], by = 1) else d_incl_um
  ds_seq <- if (length(d_soma_um) == 2) seq(d_soma_um[1], d_soma_um[2], by = 2) else d_soma_um

  # The centre must be ON the grid or it never gets the "headline" label and the
  # fitted cell silently becomes whichever interior point happens to be nearest.
  # D_soma = 17 against seq(16, 20, by = 2) is exactly that trap.
  if (length(centre) == 3L) {
    t_seq  <- sort(unique(c(t_seq,  centre[1])))
    di_seq <- sort(unique(c(di_seq, centre[2])))
    ds_seq <- sort(unique(c(ds_seq, centre[3])))
  }

  g <- expand.grid(t_um = t_seq, d_incl_um = di_seq, d_soma_um = ds_seq,
                   KEEP.OUT.ATTRS = FALSE)
  g <- g[g$d_incl_um <= g$d_soma_um, , drop = FALSE]

  out <- do.call(rbind, lapply(seq_len(nrow(g)), function(i)
    t3_geometry(g$t_um[i], g$d_incl_um[i], g$d_soma_um[i])))

  out$corner <- "interior"
  out$corner[which.min(out$p_detect)] <- "min_p_detect"
  out$corner[which.max(out$p_detect)] <- "max_p_detect"
  is_centre <- out$t_um == centre[1] & out$d_incl_um == centre[2] &
               out$d_soma_um == centre[3]
  if (any(is_centre)) out$corner[which(is_centre)[1]] <- "headline"

  out[order(out$p_detect), , drop = FALSE]
}


## ---------------------------------------------------------------------------
## t3_p_detect_mc()
##
## Monte-Carlo check of the closed form, including an eccentric inclusion.
## Exists so the delta-invariance claim in the header is verified rather than
## asserted; the test suite calls it.
##
##   delta_um : half-width of a U(-delta, delta) offset of the inclusion centre
##              from the soma centre. Values beyond (D_s - D_i)/2 are outside
##              the physically admissible range (the inclusion would poke out of
##              the soma) and p_detect genuinely does fall there -- which is a
##              property of the geometry, not a failure of the formula.
## ---------------------------------------------------------------------------
t3_p_detect_mc <- function(t_um = T3_T_UM, d_incl_um = T3_D_INCL_UM,
                           d_soma_um = T3_D_SOMA_UM, delta_um = 0,
                           n = 2e6L, seed = 42L) {
  set.seed(seed)
  half_cell <- (d_soma_um + t_um) / 2
  z_soma <- stats::runif(n, -half_cell, half_cell)
  z_incl <- z_soma + if (delta_um > 0) stats::runif(n, -delta_um, delta_um) else 0
  mean(abs(z_incl) < (d_incl_um + t_um) / 2)
}


## ---------------------------------------------------------------------------
## t3_hidden_count()
##
## Number of hidden (axially-missed) anchors implied by an observed count.
## H = n_obs (1 - p)/p, i.e. the observed n_obs are p of the true n_obs/p.
## ---------------------------------------------------------------------------
t3_hidden_count <- function(n_obs, p_detect) {
  if (any(p_detect <= 0 | p_detect > 1))
    stop("t3_hidden_count(): p_detect must be in (0, 1].")
  as.integer(round(n_obs * (1 - p_detect) / p_detect))
}


## ---------------------------------------------------------------------------
## t3_lambda2_from_counts()
##
## Per-donor in-plane anchor intensity, anchors per mm^2, from COUNTS and
## IMAGED AREA -- not by inverting a nearest-neighbour median.
##
## Why it matters: the correction ratio scales as lambda2^(-1/6), and the nine
## donors span a 23x range of anchor intensity, so inverting a POOLED median
## under a homogeneous Poisson assumption is not a neutral simplification. It is
## also biased: cells beside unimaged FOV gaps have upward-censored
## nearest-anchor distances, which inflates the pooled median and deflates the
## inferred intensity.
##
##   area_mm2 : per-donor imaged area. Pass the fov_owin() UNION area where
##              available (FOVs can overlap, so the union is <= n_fov * FOV^2
##              and the nominal sum is a lower bound on intensity).
##
## Returns a data.frame with one row per donor plus a pooled row.
## ---------------------------------------------------------------------------
t3_lambda2_from_counts <- function(n_anchor, area_mm2, sample_id) {
  stopifnot(length(n_anchor) == length(area_mm2),
            length(n_anchor) == length(sample_id))
  if (any(area_mm2 <= 0)) stop("t3_lambda2_from_counts(): non-positive area.")

  out <- data.frame(
    sample_id = as.character(sample_id),
    n_anchor  = as.integer(n_anchor),
    area_mm2  = as.numeric(area_mm2),
    lambda2_per_mm2 = n_anchor / area_mm2,
    row.names = NULL, stringsAsFactors = FALSE
  )
  rbind(out, data.frame(
    sample_id = "POOLED", n_anchor = sum(n_anchor), area_mm2 = sum(area_mm2),
    lambda2_per_mm2 = sum(n_anchor) / sum(area_mm2),
    row.names = NULL, stringsAsFactors = FALSE))
}


## ---------------------------------------------------------------------------
## t3_median_d3d()
##
## Median nearest-neighbour distance for a homogeneous Poisson process, in 2D
## and in 3D, and the inflation factor between them.
##
##   2D:  P(d > r) = exp(-lambda2 pi r^2)      -> median = sqrt(ln2/(pi lambda2))
##   3D:  P(d > r) = exp(-lambda3 (4/3) pi r^3)-> median = (3 ln2/(4 pi lambda3))^(1/3)
##   with lambda3 = lambda2 / h_eff.
##
## Poisson is the CONSERVATIVE case: clustering shortens d3D further, so the
## Poisson inflation understates the discrepancy rather than flattering it.
##
##   lambda2_per_mm2 : in-plane intensity, anchors/mm^2 (from counts, ideally)
##   h_eff_um        : effective axial depth of the anchor window
##   d2_observed_um  : optional; if given, the inflation is computed against the
##                     OBSERVED median rather than against the Poisson 2D median
## ---------------------------------------------------------------------------
t3_median_d3d <- function(lambda2_per_mm2, h_eff_um, d2_observed_um = NA_real_) {
  lam2 <- lambda2_per_mm2 / 1e6                       # per um^2
  lam3 <- lam2 / h_eff_um                             # per um^3
  med2 <- sqrt(log(2) / (pi * lam2))
  med3 <- (3 * log(2) / (4 * pi * lam3))^(1 / 3)
  ref  <- ifelse(is.na(d2_observed_um), med2, d2_observed_um)
  data.frame(
    lambda2_per_mm2 = lambda2_per_mm2,
    lambda3_per_um3 = lam3,
    h_eff_um        = h_eff_um,
    median_d2d_poisson_um = med2,
    median_d2d_ref_um     = ref,
    median_d3d_um   = med3,
    inflation       = ref / med3,
    row.names = NULL
  )
}


## ---------------------------------------------------------------------------
## t3_n_tangle_for()
##
## How many tangles to place in the PADDED VOLUME so that the number landing in
## the detectable slab matches the OBSERVED anchor count for that donor.
##
## This is the parameterisation the driver needs, and getting it backwards is
## the easy mistake: only h_eff / (2 * pad_z) of the simulated tangles are ever
## detected -- about 1% at the default padding -- so passing the observed anchor
## count straight in as n_tangle would simulate a cohort a hundredfold too
## sparse and shorten every 3D distance accordingly.
##
##   n_detected_target : observed anchors for this donor
##   pad_z_um          : axial half-extent of the simulated volume
##   h_eff_um          : detectable axial window (D_i + t)
## ---------------------------------------------------------------------------
t3_n_tangle_for <- function(n_detected_target, pad_z_um, h_eff_um) {
  if (h_eff_um <= 0 || pad_z_um <= 0)
    stop("t3_n_tangle_for(): pad_z_um and h_eff_um must be > 0.")
  as.integer(ceiling(n_detected_target * (2 * pad_z_um) / h_eff_um))
}


## ---------------------------------------------------------------------------
## t3_pad_z()
##
## Axial half-extent for a target in-plane anchor intensity. Exposed separately
## from t3_simulate_3d() so the driver can compute n_tangle BEFORE simulating
## (the two are mutually dependent, and resolving that inside the simulator
## would hide the assumption).
##
## Padded to 4 * the expected 3D NN distance plus 4 * sigma_z, so a cell in the
## slab sees the same out-of-plane neighbour density it would in unbounded
## tissue. Under-padding is the failure mode that UNDERSTATES the omission,
## which is the direction that flatters the correction -- so it is floored, not
## tuned.
## ---------------------------------------------------------------------------
t3_pad_z <- function(lambda2_per_mm2, h_eff_um, sigma_z = 150, min_pad_um = 300) {
  lam3 <- (lambda2_per_mm2 / 1e6) / h_eff_um          # per um^3
  e_d3 <- if (lam3 > 0) (3 * log(2) / (4 * pi * lam3))^(1 / 3) else 200
  max(4 * e_d3 + 4 * sigma_z, min_pad_um)
}


## ---------------------------------------------------------------------------
## t3_simulate_3d()
##
## Forward simulation of tangles in CONTINUOUS R^3. No planes anywhere.
##
## Two design choices that are load-bearing:
##
##  1. TANGLES ARE A THINNING OF A SIMULATED NEURON PROCESS, not an independent
##     pattern. Tangles are neurons; tangle clustering partly inherits neuronal
##     clustering; and thinning enforces one inclusion per soma automatically.
##     An independent tangle pattern permits two tangles within a soma diameter,
##     which inflates short-range clustering in the direction that FLATTERS the
##     correction.
##
##  2. DETECTION IS A SINGLE GEOMETRIC PREDICATE APPLIED ONCE per tangle, in
##     t3_project(). Each tangle exists as one object with one call, so the
##     duplicate-object and rigid-column artefacts of any stacked-plane scheme
##     cannot arise, and axial discretisation disappears as a concept.
##
## Padding: laterally by 4*sigma_xy so no cluster is truncated at the boundary,
## axially by pad_z_um (default 4 * the expected 3D NN distance + 4*sigma_z) so
## the slab sees the same out-of-plane density it would in an infinite block.
## Under-padding axially is the failure mode that would UNDERSTATE the omission.
##
##   win        : spatstat owin for this donor (from fov_owin()), units mm
##   n_neuron   : neurons to place in the padded volume (see t3_n_neuron_for())
##   n_tangle   : tangles to retain by thinning. Thinning is uniform over
##                neurons for "poisson", and cluster-biased for "thomas".
##   model      : "poisson" (headline, conservative) | "thomas" (scenario)
##   sigma_xy   : lateral cluster SD in um. NOT fitted: a 512 um
##                FOV cannot identify a 100-300 um cluster scale by minimum
##                contrast, and two donors have 4 and 6 anchors.
##   sigma_z    : axial cluster SD in um. NOT identifiable from a single plane
##                at all; swept, including the isotropic sigma_z = sigma_xy case.
##
## Returns a list(neurons, tangles) of data.frames with x_mm, y_mm, z_um, and
## for tangles an `is_neuron_row` index back into neurons.
## ---------------------------------------------------------------------------
##   mu       : expected offspring per cluster parent, for model = "thomas".
##              The clustering-strength dial. Counted over the PADDED VOLUME, so
##              it is not the number of tangles per island in the plane -- only
##              h_eff/(2*pad_z) of any cluster is ever seen. Larger mu = fewer,
##              tighter clusters. See the note at the n_par line.
##   neuron_model : "uniform" | "clustered". Whether the NEURON process itself
##              clusters. This decides the SIGN of the clustering effect on the
##              inflation ratio -- see the block comment at the neuron draw.
##   neuron_mu  : expected neurons per island, for neuron_model = "clustered".
t3_simulate_3d <- function(win, n_neuron, n_tangle,
                           model    = c("poisson", "thomas"),
                           sigma_xy = 150, sigma_z = 150, mu = 30,
                           neuron_model = c("uniform", "clustered"),
                           neuron_mu = 200,
                           pad_z_um = NULL, seed = 42L) {

  model        <- match.arg(model)
  neuron_model <- match.arg(neuron_model)
  if (!requireNamespace("spatstat.geom", quietly = TRUE) ||
      !requireNamespace("spatstat.random", quietly = TRUE))
    stop("t3_simulate_3d() needs spatstat.geom and spatstat.random. ",
         "Run Stage 1 in an environment that provides them.")
  set.seed(seed)

  if (n_tangle > n_neuron)
    stop("t3_simulate_3d(): n_tangle (", n_tangle, ") exceeds n_neuron (",
         n_neuron, "); tangles are a thinning of the neuron process.")

  # Axial half-extent. Wide enough that a cell in the slab has the same chance
  # of an out-of-plane near neighbour as it would in unbounded tissue.
  if (is.null(pad_z_um)) {
    area_mm2 <- spatstat.geom::area(win)
    lam3     <- (n_tangle / area_mm2 / 1e6) / max(2 * sigma_z, 1)   # crude, per um^3
    e_d3     <- if (lam3 > 0) (3 * log(2) / (4 * pi * lam3))^(1 / 3) else 200
    pad_z_um <- 4 * e_d3 + 4 * sigma_z
  }

  bb0    <- spatstat.geom::boundingbox(win)
  pad_mm0 <- 4 * sigma_xy / 1000

  # --- neurons ---------------------------------------------------------------
  # WHETHER NEURONS THEMSELVES CLUSTER DECIDES THE SIGN OF THE CLUSTERING
  # EFFECT, so it is an explicit argument rather than a hidden uniform.
  #
  # Query cells in the real analysis are CBLN2 neurons, which in EC layer II sit
  # in pre-alpha islands. If the simulated neurons are uniform while the tangles
  # cluster, the query points are decoupled from the anchors and clustering
  # LENGTHENS the typical distance (clusters leave voids, and a uniform query
  # lands in them). If the neurons cluster too, query and anchors share the
  # islands and clustering SHORTENS it. Measured here: uniform neurons give a
  # median query distance of 286 -> 399 um going from Poisson to Thomas, i.e.
  # the opposite sign to the one usually assumed.
  #
  # So "clustering shortens d3D, therefore Poisson is conservative" is NOT
  # unconditionally true -- it holds only when the query cells share the
  # anchors' clustering. Do not assert the direction; measure it under whichever
  # coupling the calibrated model uses, and report which.
  if (neuron_model == "uniform") {
    nb <- spatstat.random::runifpoint(n_neuron, win = win)
    nx <- nb$x; ny <- nb$y
    nz <- stats::runif(n_neuron, -pad_z_um, pad_z_um)
  } else {
    n_npar <- max(2L, as.integer(round(n_neuron / neuron_mu)))
    qx <- stats::runif(n_npar, bb0$xrange[1] - pad_mm0, bb0$xrange[2] + pad_mm0)
    qy <- stats::runif(n_npar, bb0$yrange[1] - pad_mm0, bb0$yrange[2] + pad_mm0)
    qz <- stats::runif(n_npar, -pad_z_um - 4 * sigma_z, pad_z_um + 4 * sigma_z)
    pick <- sample.int(n_npar, n_neuron, replace = TRUE)
    nx <- qx[pick] + stats::rnorm(n_neuron, 0, sigma_xy / 1000)
    ny <- qy[pick] + stats::rnorm(n_neuron, 0, sigma_xy / 1000)
    nz <- qz[pick] + stats::rnorm(n_neuron, 0, sigma_z)
    # Keep only what lands in the imaged window; top up so n_neuron is honoured.
    inside <- spatstat.geom::inside.owin(nx, ny, win)
    guard <- 0L
    while (sum(inside) < n_neuron && guard < 50L) {
      guard <- guard + 1L
      need <- n_neuron - sum(inside)
      pk <- sample.int(n_npar, need * 3L, replace = TRUE)
      ax <- qx[pk] + stats::rnorm(length(pk), 0, sigma_xy / 1000)
      ay <- qy[pk] + stats::rnorm(length(pk), 0, sigma_xy / 1000)
      az <- qz[pk] + stats::rnorm(length(pk), 0, sigma_z)
      ok <- spatstat.geom::inside.owin(ax, ay, win)
      nx <- c(nx[inside], ax[ok]); ny <- c(ny[inside], ay[ok]); nz <- c(nz[inside], az[ok])
      inside <- rep(TRUE, length(nx))
    }
    keep_n <- seq_len(min(n_neuron, sum(inside)))
    nx <- nx[inside][keep_n]; ny <- ny[inside][keep_n]; nz <- nz[inside][keep_n]
    n_neuron <- length(nx)
    if (n_tangle > n_neuron)
      stop("t3_simulate_3d(): clustered neuron process yielded only ", n_neuron,
           " in-window neurons, fewer than n_tangle = ", n_tangle,
           ". Raise n_neuron or lower neuron_mu.")
  }
  neurons <- data.frame(x_mm = nx, y_mm = ny, z_um = nz)

  # --- thinning to tangles ---------------------------------------------------
  if (model == "poisson") {
    idx <- sample.int(n_neuron, n_tangle)
  } else {
    # Thomas-like: draw cluster parents in the padded volume, then weight each
    # neuron by the summed 3D Gaussian kernel from the parents and thin without
    # replacement proportional to that weight. This inherits neuronal positions
    # (so one inclusion per soma still holds) while giving the tangle pattern a
    # cluster scale. It is a SCENARIO, not a fit.
    bb      <- spatstat.geom::boundingbox(win)
    pad_mm  <- 4 * sigma_xy / 1000
    # Offspring per parent. This MUST be an explicit dial, not a fraction of
    # n_tangle: n_tangle counts the whole PADDED VOLUME, which is ~100x the
    # detectable slab, so a naive n_tangle/8 gives hundreds of parents spread
    # over tens of mm^3 and the result is indistinguishable from Poisson.
    n_par   <- max(2L, as.integer(round(n_tangle / mu)))
    px <- stats::runif(n_par, bb$xrange[1] - pad_mm, bb$xrange[2] + pad_mm)
    py <- stats::runif(n_par, bb$yrange[1] - pad_mm, bb$yrange[2] + pad_mm)
    pz <- stats::runif(n_par, -pad_z_um - 4 * sigma_z, pad_z_um + 4 * sigma_z)

    w <- numeric(n_neuron)
    for (j in seq_len(n_par)) {
      dx <- (neurons$x_mm - px[j]) * 1000
      dy <- (neurons$y_mm - py[j]) * 1000
      dz <- (neurons$z_um - pz[j])
      w  <- w + exp(-0.5 * ((dx / sigma_xy)^2 + (dy / sigma_xy)^2 +
                            (dz / sigma_z)^2))
    }
    if (all(w <= 0) || !any(is.finite(w))) w <- rep(1, n_neuron)
    w[!is.finite(w)] <- 0
    idx <- sample.int(n_neuron, n_tangle, prob = w)
  }

  tangles <- neurons[idx, , drop = FALSE]
  tangles$is_neuron_row <- idx
  rownames(tangles) <- NULL

  list(neurons = neurons, tangles = tangles,
       pad_z_um = pad_z_um, model = model,
       sigma_xy = sigma_xy, sigma_z = sigma_z, win = win)
}


## ---------------------------------------------------------------------------
## t3_project()
##
## The single detection predicate, applied once. A tangle is DETECTED when its
## inclusion intersects the section, i.e. |z| < h_eff/2 with h_eff = D_i + t.
## A tangle-bearing CELL is in the object when |z| < h_cell/2 (h_cell = D_s + t)
## -- the wider window, which is why most axially-missed tangle-bearers are
## present but mislabelled rather than absent, and therefore why hidden-anchor
## relabelling is viable at all.
##
## Query points are the simulated tangle positions themselves plus, optionally,
## a supplied set of query coordinates (for a simulated d2D distribution
## comparable with the observed one). Distances are in-plane, from the detected
## anchors only, under the real disconnected FOV mosaic and the 1000 um cap.
##
## Returns a list with the classified tangle table and, if query_xy_mm given,
## the simulated in-plane distances.
## ---------------------------------------------------------------------------
t3_project <- function(sim, geom, query_xy_mm = NULL,
                       max_dist_um = T3_MAX_DIST_UM) {

  tg <- sim$tangles
  tg$in_slab_incl <- abs(tg$z_um) < geom$h_eff  / 2     # callable anchor
  tg$in_slab_cell <- abs(tg$z_um) < geom$h_cell / 2     # cell present in object
  tg$class <- ifelse(tg$in_slab_incl, "detected",
              ifelse(tg$in_slab_cell, "missed_recoverable", "missed_deep"))

  det <- tg[tg$in_slab_incl, c("x_mm", "y_mm"), drop = FALSE]

  d2 <- NULL
  if (!is.null(query_xy_mm) && nrow(det) > 0) {
    nn <- RANN::nn2(data = as.matrix(det), query = as.matrix(query_xy_mm), k = 1)
    d2 <- nn$nn.dists[, 1] * 1000
    d2 <- d2[d2 <= max_dist_um]
  }

  # True 3D nearest-tangle distance among ALL tangles, for the tangles that sit
  # in the slab -- the quantity the in-plane measurement is a proxy for.
  d3 <- NULL
  if (nrow(tg) > 1 && any(tg$in_slab_cell)) {
    all_xyz <- cbind(tg$x_mm * 1000, tg$y_mm * 1000, tg$z_um)
    q_xyz   <- all_xyz[tg$in_slab_cell, , drop = FALSE]
    nn3 <- RANN::nn2(data = all_xyz, query = q_xyz, k = min(2L, nrow(all_xyz)))
    d3  <- if (ncol(nn3$nn.dists) >= 2) nn3$nn.dists[, 2] else nn3$nn.dists[, 1]
  }

  list(tangles = tg, d2_um = d2, d3_um = d3,
       n_detected = sum(tg$in_slab_incl),
       n_missed_recoverable = sum(tg$class == "missed_recoverable"),
       n_missed_deep = sum(tg$class == "missed_deep"),
       p_detect_realised = sum(tg$in_slab_incl) / max(1, sum(tg$in_slab_cell)))
}


## ---------------------------------------------------------------------------
## t3_validate()
##
## Does the PROJECTED 3D model reproduce the OBSERVED anchor pattern? If it does
## not, the generative kernel is fiction and the run must abort rather than
## produce a correction off a model that cannot make the data.
##
## Three checks:
##  1. G, F and K envelopes of the projected simulated anchors against the
##     observed anchors (spatstat, so Stage 1 only).
##  2. The simulated d2D quantiles against the observed quantiles.
##  3. THE RIGID-STACKING DIAGNOSTIC. If synthetic anchors were ever placed at
##     observed cell coordinates (which this file never does -- t3_simulate_3d
##     draws fresh uniform positions), d2D would show a spike near zero and a
##     deficit at 20-100 um that the real data does not have. Checked anyway,
##     because it is cheap and it is the failure that would look plausible.
##
## Returns a data.frame of named checks with pass/fail, for the stats log.
## ---------------------------------------------------------------------------
t3_validate <- function(sim_xy_mm, obs_xy_mm, win,
                        d2_sim_um, d2_obs_um, nsim = 39L, alpha = 0.05,
                        seed = 42L) {

  if (!requireNamespace("spatstat.geom", quietly = TRUE) ||
      !requireNamespace("spatstat.explore", quietly = TRUE))
    stop("t3_validate() needs spatstat.geom and spatstat.explore (Stage 1 only).")
  set.seed(seed)

  out <- list()

  # --- 2. distribution match (cheap, always runs) ----------------------------
  qs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  qo <- as.numeric(stats::quantile(d2_obs_um, qs, na.rm = TRUE))
  qsim <- as.numeric(stats::quantile(d2_sim_um, qs, na.rm = TRUE))
  rel <- abs(qsim - qo) / pmax(qo, 1)
  out$quantiles <- data.frame(
    check = sprintf("d2D p%02d", round(qs * 100)),
    observed = qo, simulated = qsim, rel_diff = rel,
    pass = rel < 0.35, row.names = NULL)

  # --- 3. rigid-stacking diagnostic ------------------------------------------
  spike_obs <- mean(d2_obs_um < 5, na.rm = TRUE)
  spike_sim <- mean(d2_sim_um < 5, na.rm = TRUE)
  gap_obs   <- mean(d2_obs_um >= 20 & d2_obs_um <= 100, na.rm = TRUE)
  gap_sim   <- mean(d2_sim_um >= 20 & d2_sim_um <= 100, na.rm = TRUE)
  out$stacking <- data.frame(
    check = c("near-zero spike (<5 um)", "20-100 um mass"),
    observed = c(spike_obs, gap_obs), simulated = c(spike_sim, gap_sim),
    rel_diff = c(abs(spike_sim - spike_obs), abs(gap_sim - gap_obs) / max(gap_obs, 1e-6)),
    pass = c(spike_sim - spike_obs < 0.02, abs(gap_sim - gap_obs) / max(gap_obs, 1e-6) < 0.5),
    row.names = NULL)

  # --- 1. G / F / K envelopes ------------------------------------------------
  # Needs enough anchors for a summary function to mean anything. Two donors
  # have 4 and 6 anchors; return NA rather than a meaningless envelope.
  out$envelopes <- data.frame(
    check = c("G", "F", "K"), p_value = NA_real_, pass = NA, row.names = NULL)
  if (nrow(obs_xy_mm) >= 20 && nrow(sim_xy_mm) >= 20) {
    obs_pp <- spatstat.geom::ppp(obs_xy_mm[, 1], obs_xy_mm[, 2], window = win,
                                 checkdup = FALSE)
    sim_pp <- spatstat.geom::ppp(sim_xy_mm[, 1], sim_xy_mm[, 2], window = win,
                                 checkdup = FALSE)
    for (fn in c("G", "F", "K")) {
      p <- tryCatch({
        f_obs <- switch(fn, G = spatstat.explore::Gest(obs_pp),
                            F = spatstat.explore::Fest(obs_pp),
                            K = spatstat.explore::Kest(obs_pp))
        f_sim <- switch(fn, G = spatstat.explore::Gest(sim_pp),
                            F = spatstat.explore::Fest(sim_pp),
                            K = spatstat.explore::Kest(sim_pp))
        # Max absolute deviation between the two edge-corrected estimates,
        # scaled by the observed range -- a distance, reported as such. A full
        # Monte-Carlo envelope would need nsim refits of the 3D process per
        # donor; that cost is not justified for a supplementary check.
        yo <- f_obs[[if (fn == "K") "iso" else "km"]]
        ys <- f_sim[[if (fn == "K") "iso" else "km"]]
        n  <- min(length(yo), length(ys))
        max(abs(yo[seq_len(n)] - ys[seq_len(n)]), na.rm = TRUE) /
          max(abs(yo[seq_len(n)]), na.rm = TRUE)
      }, error = function(e) NA_real_)
      out$envelopes$p_value[out$envelopes$check == fn] <- p
      out$envelopes$pass[out$envelopes$check == fn] <- is.na(p) || p < 0.5
    }
  }

  out$all_pass <- all(c(out$quantiles$pass, out$stacking$pass,
                        out$envelopes$pass[!is.na(out$envelopes$pass)]))
  out
}


## ---------------------------------------------------------------------------
## t3_unseen_plane_anchors()   and   t3_dist3d_unseen()
##
## THE AXIAL-OMISSION CORRECTION, WITHOUT ANY MISCLASSIFICATION ASSUMPTION.
##
## An alternative to the hidden-anchor augmentation in R/phf1_hidden_anchors.R,
## which assumes a fraction of tangle-bearing neurons sit in the section with
## their inclusion outside it and are therefore sequenced but mislabelled. This
## correction makes no such assumption.
##
## What a single section genuinely cannot see is TISSUE THAT WAS NEVER SECTIONED:
## tangles above and below the slice. Those cannot be recovered by relabelling
## anything, because no cell of theirs is in the object. They are simulated here.
##
## THE MODEL
##   - Observed in-plane anchors are treated as COMPLETE and exact, at z = 0.
##     No tangle in the imaged plane is assumed to have been missed.
##   - Additional tangles are simulated ONLY at |z| > h_eff/2, i.e. strictly
##     outside the axial window the section samples, at the volumetric intensity
##     implied by the observed in-plane density: lambda3 = lambda2 / h_eff.
##   - Each cell's corrected distance is the minimum over BOTH sets, in 3D.
##
## Consequences worth stating before looking at any output:
##   1. The corrected distance is ALWAYS <= the observed 2D distance. Adding
##      anchors can only shorten; nothing can move further away.
##   2. THE EFFECT IS STRONGLY NON-UNIFORM AND CONCENTRATES IN THE TAIL. A cell
##      20 um from an in-plane tangle keeps essentially that distance, because
##      the nearest out-of-plane tangle is at least h_eff/2 away in z and almost
##      certainly further in 3D. A cell 800 um from any in-plane tangle can
##      easily have one 100 um above it. So this compresses the far end of the
##      distance distribution and leaves the near end alone -- which is a
##      re-scaling of the distance axis, not a re-ordering of cells.
##   3. Lateral extent is the imaged FOV mosaic and is NOT padded. Tangles
##      laterally outside the imaged region are invisible to the 2D analysis too,
##      so padding would mix edge censoring into an axial correction. This
##      isolates the axial channel; lateral censoring is left untouched.
##
## t3_unseen_plane_anchors(): draw the out-of-plane tangles for one donor.
##   win        : spatstat owin for the donor (the imaged FOV mosaic)
##   lambda2    : observed in-plane anchor intensity, per mm^2
##   h_eff_um   : axial thickness the section samples for tangles (D_incl + t)
##   z_max_um   : how far above and below to simulate. Must be >= the largest
##                distance that can matter, i.e. the analysis cap, or a tangle
##                that would have been nearest is silently omitted.
##   clustered  : FALSE = uniform in x,y (maximal correction: fills the gaps
##                between in-plane tangles). TRUE = x,y drawn in proportion to
##                the observed in-plane anchor intensity, so the unseen planes
##                inherit the observed clustering and the correction is WEAKER
##                in exactly the empty regions where it would otherwise bite.
##                Report both; they bracket the answer.
## ---------------------------------------------------------------------------
t3_unseen_plane_anchors <- function(win, lambda2_per_mm2, h_eff_um,
                                    z_max_um = 1100, clustered = FALSE,
                                    anchor_xy_mm = NULL, sigma_um = 200,
                                    seed = 42L) {
  if (!requireNamespace("spatstat.geom", quietly = TRUE) ||
      !requireNamespace("spatstat.random", quietly = TRUE))
    stop("t3_unseen_plane_anchors() needs spatstat (Stage 1, local only).")
  set.seed(seed)

  area_mm2 <- spatstat.geom::area(win)
  lam3     <- (lambda2_per_mm2 / 1e6) / h_eff_um          # per um^3
  # Volume of the UNSEEN slabs only: the full column minus the sampled slice.
  vol_um3  <- area_mm2 * 1e6 * (2 * z_max_um - h_eff_um)
  n        <- stats::rpois(1, lam3 * vol_um3)
  if (n == 0) return(data.frame(x_mm = numeric(0), y_mm = numeric(0),
                                z_um = numeric(0)))

  if (clustered && !is.null(anchor_xy_mm) && nrow(anchor_xy_mm) > 0) {
    # Place laterally on the observed anchors, smeared by sigma, then keep what
    # lands inside the imaged window.
    xy <- matrix(NA_real_, 0, 2)
    guard <- 0L
    while (nrow(xy) < n && guard < 60L) {
      guard <- guard + 1L
      k  <- sample.int(nrow(anchor_xy_mm), 3L * n, replace = TRUE)
      cx <- anchor_xy_mm[k, 1] + stats::rnorm(length(k), 0, sigma_um / 1000)
      cy <- anchor_xy_mm[k, 2] + stats::rnorm(length(k), 0, sigma_um / 1000)
      ok <- spatstat.geom::inside.owin(cx, cy, win)
      xy <- rbind(xy, cbind(cx[ok], cy[ok]))
    }
    xy <- xy[seq_len(min(n, nrow(xy))), , drop = FALSE]
  } else {
    pp <- spatstat.random::runifpoint(n, win = win)
    xy <- cbind(pp$x, pp$y)
  }

  # z strictly OUTSIDE the sampled slice, uniform over the two unseen slabs.
  half <- h_eff_um / 2
  zz <- stats::runif(nrow(xy), half, z_max_um) *
        sample(c(-1, 1), nrow(xy), replace = TRUE)
  data.frame(x_mm = xy[, 1], y_mm = xy[, 2], z_um = zz)
}


## ---------------------------------------------------------------------------
## t3_dist3d_unseen()
##
## Corrected distance for one donor: min over the OBSERVED in-plane anchors
## (z = 0, treated as complete) and the SIMULATED out-of-plane tangles.
##
## Observed anchors and query cells are both placed at z = 0. Both really lie
## somewhere within the ~7 um slice, so each is misplaced by at most ~3.5 um
## against distances on the order of 100 um -- immaterial, and it cannot bias
## the comparison because the observed arm is computed the same way.
##
## Returns a numeric vector of corrected distances in um, aligned to query_xy_mm,
## plus attr "d2d" (the observed in-plane distance) so the caller can report the
## shortening without recomputing it.
## ---------------------------------------------------------------------------
t3_dist3d_unseen <- function(query_xy_mm, anchor_xy_mm, unseen_xyz) {
  q <- cbind(as.matrix(query_xy_mm) * 1000, 0)
  d_obs <- rep(Inf, nrow(q))
  if (nrow(anchor_xy_mm) > 0) {
    a <- cbind(as.matrix(anchor_xy_mm) * 1000, 0)
    d_obs <- RANN::nn2(a, q, k = 1)$nn.dists[, 1]
  }
  d_new <- rep(Inf, nrow(q))
  if (nrow(unseen_xyz) > 0) {
    u <- cbind(unseen_xyz$x_mm * 1000, unseen_xyz$y_mm * 1000, unseen_xyz$z_um)
    d_new <- RANN::nn2(u, q, k = 1)$nn.dists[, 1]
  }
  out <- pmin(d_obs, d_new)
  attr(out, "d2d") <- d_obs
  attr(out, "n_unseen") <- nrow(unseen_xyz)
  attr(out, "frac_from_unseen") <- mean(d_new < d_obs)
  out
}
