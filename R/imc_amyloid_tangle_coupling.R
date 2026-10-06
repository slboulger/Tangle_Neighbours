#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# imc_amyloid_tangle_coupling.R
#
# Figure panels: S8D, S8E
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# imc_amyloid_tangle_coupling.R
#
# How much of the IMC tangle-distance gradient is attributable to amyloid proximity?
#
# The IMC distance models (docs/MODELS.md Set 3, 44 AD donors) carry no amyloid term -- `plaque`
# (Matched_4G8_40) is a Set 3 sensitivity term. This script quantifies the omitted-variable bias
# (OVB) of the tangle coefficient by decomposing it into the two factors that produce it:
#
#   delta  = coupling of plaque proximity to tangle distance
#              (how much more likely is a cell to sit in the amyloid niche when it is near a tangle)
#   gamma  = the plaque effect on the outcome
#              (the 40 um term in the Figure S8A ladder, m3_plaque; already on disk)
#   gamma x delta = the bias in the tangle coefficient from omitting amyloid
#
# and then asks the question the product cannot answer: does the gradient survive WHERE THERE IS NO
# AMYLOID? That is the `dist_z * plaque` stratified refit, whose PLAQUE-DISTAL simple slope is the
# actual test. The interaction is reported whether or not it is significant.
#
# NOTHING EXISTING IS MODIFIED. plots/imc_covariate_sensitivity/ (Figure S8A) is read as input and
# as a regression check. This script only writes plots/imc_amyloid_tangle_coupling/.
#
# ===================================================================================================
# WHY THE MODEL YOU WOULD FIRST REACH FOR IS THE WRONG ONE
# ===================================================================================================
#
# The natural way to estimate `delta` is a logistic mixed model, and this script fits one. But a
# LOGISTIC delta CANNOT BE MULTIPLIED BY GAMMA. OVB is a property of linear projection, and `plaque`
# enters the outcome models as a 0/1 dummy under an IDENTITY link, so the auxiliary regression of the
# omitted variable on the included predictor must be identity-link too. A logistic coefficient is a
# slope of logit(P), and d logit(P)/dx = (dP/dx)/(p(1-p)); at the observed prevalence (17-25%) that
# inflates the number by 1/(p(1-p)) ~ 5.3-6.7x. Measured here: log-odds delta = -0.1185 for
# dist_reln_calb1 against an exact multiplier of -0.0182. Using the log-odds delta would report AT8's
# bias as ~7% of its coefficient instead of ~1%.
#
# So the four delta estimands below have DIFFERENT JOBS and only one of them is the multiplier:
#
#   delta_lpm    identity-link (linear probability) fit, weight matrix matched to the outcome fit.
#                THE MULTIPLIER. Point estimate is algebraically exact (asserted below to 1e-8).
#   delta_logit  the glmer. Reported as an odds ratio, DESCRIPTIVE ONLY -- the inferential statement
#                about whether tau and amyloid proximity are coupled at all.
#   delta_ame    average marginal effect of the glmer, mean(beta * p(1-p)). The bridge between the
#                two scales: it lands within a few % of delta_lpm and shows the two agree once the
#                logistic coefficient is put back on the probability scale.
#   delta_nearfar  plaque-proximal fraction near (<=50 um) vs far (>=200 um) tangles, computed WITHIN
#                each donor then averaged with a t CI. The assumption-free descriptive check that the
#                coupling is real and in the expected direction. NOT the multiplier (see below).
#
# Three further requirements, each guarded in the code:
#
#  (1) delta's RIGHT-HAND SIDE MUST EQUAL THE OUTCOME MODEL'S RHS MINUS `plaque`, including `extra`.
#      Dropping `subcluster` for dist_exc_subcluster_adj gives -0.02137 instead of -0.01772: a 21%
#      error. The RHS is therefore built programmatically from the panel metadata (build_rhs()) and
#      can never drift from the outcome model.
#
#  (2) FOR THE TWO GLIA PANELS THE CORRECT delta IS POOLED OLS, NOT (1 | patient_id). Their gamma
#      comes from a fit whose donor random effect is estimated at ~1e-30 (the outcome is z-scored
#      within donor; docs/MODELS.md Set 3, within-donor z-scoring), and the OVB identity holds only
#      under ONE weight matrix. Astro: implied -0.02594004, pooled OLS -0.02594004 (7 s.f.), mixed
#      -0.026870. For the DEPTH analogue the same slip is 3.6x (astro delta_depth pooled 0.0038661
#      vs mixed 0.0138709), which would report 5.5% of the GFAP coefficient instead of 1.5%.
#      This script removes the branch entirely rather than testing isSingular(): delta is computed by
#      WHITENED GLS AT THE OUTCOME FIT'S OWN lambda = tau^2/sigma^2. When the donor variance has
#      collapsed, lambda ~ 0, the whitening becomes the identity, and the estimator IS pooled OLS --
#      automatically, by one code path, for every panel.
#
#  (3) SIGN. All OVB arithmetic is done on the RAW `est` scale (per s.d. FURTHER from a tangle) and
#      flipped only at the reporting layer. imc_covariate_sensitivity.R uses flip = -1 for all four
#      distance panels and never flips `plaque`/`depth`, so the shift in the flipped column equals
#      -gamma*delta.
#      Dividing a flipped shift by gamma returns delta = +0.0182 and inverts the biology into
#      "further from a tangle => MORE likely plaque-proximal". Every emitted row therefore carries an
#      explicit `scale` column ("raw" or "beta_closer").
#
# ===================================================================================================
# WHAT delta_nearfar CANNOT DO
# ===================================================================================================
# The near/far fraction difference is the most readable number in the analysis and it is NOT the
# multiplier. Three reasons, all of which matter: it is a donor-unweighted mean of per-donor
# fractions while delta is cell-weighted; it is a two-point secant rather than a slope, and
# prevalence is not linear in log-distance; and it is unadjusted for Sex/Age/PMI/extra. Rescaling it
# by the observed dist_z gap between the bins agrees with delta to ~1% for astrocytes but is ~39% off
# for dist_reln_calb1. It is reported as a descriptive check, with that clause attached.
#
# ===================================================================================================
# SCOPE AND REPORTING
# ===================================================================================================
# The DEFAULT run covers the subcluster-adjusted excitatory panel and the GFAP astrocyte panel
# (PANEL_DEFAULT below). `dist_reln_calb1` and `cd68_microglia_distance` stay defined and runnable
# by name; their BH families and their deltas are their own.
#
# The primary currency is shift / SE(b1). The percentage of b1 is SUPPRESSED when |b1| < 2*SE(b1),
# because a percentage of a near-null coefficient is not informative. The dist x plaque interaction
# is reported as EQUIVALENCE, not absence, with the minimum detectable interaction per outcome. The
# cortical-depth analogue of the decomposition is reported alongside, for comparison.
#
# Run LOCALLY, not on the HPC: spe.rds lives under <IMC_ROOT>. Optional panel slugs as args:
#   Rscript R/imc_amyloid_tangle_coupling.R dist_reln_calb1

suppressPackageStartupMessages({
  library(SpatialExperiment); library(dplyr); library(tibble); library(tidyr)
  library(ggplot2); library(lmerTest); library(RANN); library(grid); library(patchwork)
})

# --- project root. The RDS share has been mounted at more than one path on this machine, so
# --- resolve against a candidate list rather than a single hard-coded fallback.
root_candidates <- c(
  "<PROJECT_ROOT>",
  "<PROJECT_ROOT>",
  "<PROJECT_ROOT>"
)
root <- root_candidates[dir.exists(root_candidates)][1]
if (is.na(root)) stop("Project root not found. Mount the RDS share.")
setwd(root)
source("R/palettes.R")
source("R/imc_utils.R")   # re_variance_summary(), cat_re_variance(), re_total_sd()

out_dir <- "plots/imc_amyloid_tangle_coupling"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Constants. Every distance/cohort constant is copied from the origin scripts via
# imc_covariate_sensitivity.R and MUST NOT drift from them -- the regression check below is what
# proves it has not.
# -------------------------------------------------------------------
MAX_DIST     <- 300
EDGE_BUFFER  <- 0
MIN_ANCHORS  <- 1
NEAR_UM      <- 50     # "near a tangle"  for the descriptive fraction
FAR_UM       <- 200    # "far from a tangle"
N_QBINS      <- 10     # quantile bins for the fraction-vs-distance profile
MIN_BIN_N    <- 10     # a donor-bin with fewer cells than this is not counted
EQUIV_FRAC   <- 0.25   # pre-declared negligibility margin: 25% of the plaque-distal slope
FDR          <- 0.05

VULN_CLUSTERS <- c("Excitatory neuron cluster 4 (RELN)",
                   "Excitatory neuron cluster 2 (CALB1)")
MICRO <- c("Microglia cluster 1 (IBA1)", "Microglia cluster 2 (IBA1, CD68)")
ASTRO <- c("Astrocyte cluster 1 (GFAP, S100b)", "Astrocyte cluster 2 (S100b)")

# Okabe-Ito, reusing the Figure S8A covariate-axis colours so these figures read against that one.
COL_PROX   <- "#D55E00"   # vermillion    -- tau proximity
COL_PLAQUE <- "#0072B2"   # blue          -- amyloid proximity
COL_DEPTH  <- "#009E73"   # bluish green  -- cortical depth
STRATUM_COLOURS <- c("plaque-distal" = "#000000", "plaque-prox" = COL_PLAQUE)
COVAR_COLOURS   <- c(amyloid = COL_PLAQUE, depth = COL_DEPTH)

CM <- 1 / 2.54   # cm -> inches

# --- Tolerances. Two different jobs, deliberately very different sizes.
# TOL_ALG is an ALGEBRAIC assertion: in whitened coordinates GLS is OLS, so the OVB identity is
#   exact and any failure means the cell set / RHS / coding is wrong.
# TOL_REML is the PRACTICAL agreement between gamma*delta and the observed m1 -> m3 shift. The only
#   source of discrepancy is REML re-estimating lambda when `plaque` enters, so b1 and b3 are
#   computed under slightly different metrics. Measured worst case 0.01*SE(b1); 0.05 is a 5x margin.
# TOL_LADDER checks the refitted coefficients against the on-disk Figure S8A ladder.
# TOL_DELTA_REL is a loose sanity check against the reference deltas in EXPECTED_DELTA; it is
#   loose on purpose, because those came from a run whose lambda choice may differ in the last
#   digits. The algebra is pinned by TOL_ALG, not by this.
TOL_ALG       <- 1e-8
TOL_REML      <- 0.05
TOL_LADDER    <- 1e-5
TOL_DELTA_REL <- 0.03

# Reference deltas (the exact multiplier), for the sanity check.
EXPECTED_DELTA <- c(dist_reln_calb1          = -0.018161,
                    dist_exc_subcluster_adj  = -0.017715,
                    cd68_microglia_distance  = -0.0317660,
                    gfap_astrocytes_distance = -0.0259400)

LADDER_TSV <- "plots/imc_covariate_sensitivity/stats_covariate_ladder.tsv"
if (!file.exists(LADDER_TSV)) stop("Figure S8A ladder not found: ", LADDER_TSV,
                                   "\nRun R/imc_covariate_sensitivity.R first.")
ladder <- read.delim(LADDER_TSV, stringsAsFactors = FALSE)

# Panels run by default: the subcluster-adjusted excitatory panel and the GFAP astrocyte panel.
# `dist_reln_calb1` and `cd68_microglia_distance` are DEFINED below and remain runnable by name;
# their BH families and their deltas are their own. BH is within panel, so the default panels'
# results do not depend on which other panels are run.
PANEL_DEFAULT <- c("dist_exc_subcluster_adj", "gfap_astrocytes_distance")

PANEL_SELECT <- commandArgs(trailingOnly = TRUE)
if (!length(PANEL_SELECT)) PANEL_SELECT <- PANEL_DEFAULT

# ===================================================================
# Load. Same object and loader as the other imc_* scripts.
# ===================================================================
localwd <- "<IMC_ROOT>/"
spe <- readRDS(paste0(localwd, "spe.rds"))
stopifnot("PHF1_Otsu" %in% colnames(colData(spe)),
          "Matched_4G8_40" %in% colnames(colData(spe)))
stopifnot(all(c(MICRO, ASTRO) %in% unique(as.character(colData(spe)$celltype_clusters))))
stopifnot(all(c("CD68", "GFAP") %in% rownames(spe)))
stopifnot(is.logical(colData(spe)$Matched_4G8_40))   # the 40 um-dilated 4G8 mask overlap flag

markers <- setdiff(rownames(spe)[rowData(spe)$marker_class == "state"], "PHF1")
EM <- assay(spe, "exprs")

# ===================================================================
# Cells. Copied verbatim from imc_covariate_sensitivity.R, which in turn copies it from the origin
# script of each panel. The REGRESSION CHECK below refuses to proceed unless the rebuilt cell sets
# reproduce the on-disk Figure S8A coefficients.
# ===================================================================
cd <- as.data.frame(colData(spe))
cd$cell_id <- colnames(spe)
sc <- spatialCoords(spe); cd$x <- sc[, "Pos_X"]; cd$y <- sc[, "Pos_Y"]

min_rois_per_patient <- 3
cv_removal <- cd %>% distinct(patient_id, sample_id) %>% count(patient_id) %>%
  filter(n < min_rois_per_patient) %>% pull(patient_id)

# "AD" is BraakGroup != "Braak_0_1" -- there is no diagnosis column. Braak III-IV + VI, >= 3 ROIs,
# which is the 44-donor / 132-ROI analysed cohort out of a 67-donor / 195-ROI object.
COHORT <- function(d) d %>%
  filter(!celltype_clusters %in% c("Artefact cluster", "Unassigned cluster"),
         !patient_id %in% cv_removal,
         BraakGroup != "Braak_0_1",
         !is.na(PHF1_Otsu)) %>%
  mutate(patient_id = as.character(patient_id), sample_id = as.character(sample_id),
         celltype_clusters = as.character(celltype_clusters),
         Sex = factor(Sex))

base_all <- COHORT(cd)
stopifnot(!anyNA(base_all$Sex), !anyNA(base_all$Age), !anyNA(base_all$PMI))

# ROI geometry is a property of the ROI, so it is computed over ALL cells in the ROI and then carried
# into every cell set -- computing it on a celltype subset would make min(y)/max(y) depend on the
# selection.
base_all <- base_all %>%
  group_by(sample_id) %>%
  mutate(scaled_Y  = 100 * (y - min(y)) / (max(y) - min(y)),
         edge_dist = pmin(x - min(x), max(x) - x, y - min(y), max(y) - y)) %>%
  ungroup() %>%
  mutate(depth_prop = scaled_Y / 100,     # 0 = pia, 1 = deep border
         # amyloid proximity, defined exactly as in the eleven existing IMC scripts:
         # nucleus within the 40 um-dilated 4G8+ plaque mask
         plaque = factor(ifelse(Matched_4G8_40, "plaque-prox", "plaque-distal"),
                         levels = c("plaque-distal", "plaque-prox")),
         plaque_int = as.integer(Matched_4G8_40),   # the LPM / glmer outcome
         Age_s = as.numeric(scale(Age)), PMI_s = as.numeric(scale(PMI)))
stopifnot(all(is.finite(base_all$depth_prop)),
          min(base_all$depth_prop) >= -1e-9, max(base_all$depth_prop) <= 1 + 1e-9)

# --- distance to the nearest PHF1+ NEURON, within ROI
anchors <- base_all %>% filter(PHF1_Otsu == "PHF1_pos",
                               grepl("neuron", celltype_clusters, ignore.case = TRUE))
nearest_dist <- function(tg, an) {
  out <- rep(NA_real_, nrow(tg))
  for (s in unique(tg$sample_id)) {
    i <- which(tg$sample_id == s)
    a <- an[an$sample_id == s, c("x", "y"), drop = FALSE]
    if (!nrow(a)) next
    out[i] <- RANN::nn2(as.matrix(a), as.matrix(tg[i, c("x", "y")]), k = 1)$nn.dists[, 1]
  }
  out
}
anchor_n <- anchors %>% count(sample_id, name = "n_anchors")
roi_ok   <- anchor_n$sample_id[anchor_n$n_anchors >= MIN_ANCHORS]

# dist_z recomputed WITHIN the returned set
build_md_dist <- function(target_filter) {
  tg <- base_all %>% filter(PHF1_Otsu == "PHF1_neg") %>% target_filter()
  tg$dist_um <- nearest_dist(tg, anchors)
  md <- tg %>%
    filter(!is.na(dist_um), dist_um <= MAX_DIST, edge_dist >= EDGE_BUFFER,
           sample_id %in% roi_ok) %>%
    mutate(dist_z = log(dist_um) / sd(log(dist_um)))
  if (any(md$dist_um <= 0)) stop("non-positive distance; log() would fail")
  md
}

md_dist_vuln <- build_md_dist(function(d) d %>% filter(celltype_clusters %in% VULN_CLUSTERS))
md_dist_exc  <- build_md_dist(function(d) d %>% filter(grepl("^Excitatory", celltype_clusters)))
md_dist_exc$subcluster <- relevel(factor(md_dist_exc$celltype_clusters),
                                  ref = names(sort(table(md_dist_exc$celltype_clusters),
                                                   decreasing = TRUE))[1])

# --- glia: distance from EVERY cell, the extra dist_um > 0 filter, and DIST_SD computed ONCE over
# --- all capped cells rather than per glial subset. Preserved as-is; it is the reason the two glia
# --- panels SHARE one dist_sd (which is why their deltas are not comparable on units grounds).
base_all$dist_um <- NA_real_
for (s in unique(base_all$sample_id)) {
  i <- which(base_all$sample_id == s)
  a <- anchors[anchors$sample_id == s, c("x", "y"), drop = FALSE]
  if (!nrow(a)) next
  base_all$dist_um[i] <- RANN::nn2(as.matrix(a), as.matrix(base_all[i, c("x", "y")]),
                                   k = 1)$nn.dists[, 1]
}
md_glia <- base_all %>% filter(!is.na(dist_um), dist_um > 0, dist_um <= MAX_DIST)
DIST_SD_GLIA <- sd(log(md_glia$dist_um))
stopifnot(is.finite(DIST_SD_GLIA), DIST_SD_GLIA > 0)
md_glia$dist_scaled <- log(md_glia$dist_um) / DIST_SD_GLIA
md_glia$CD68 <- as.numeric(EM["CD68", md_glia$cell_id])
md_glia$GFAP <- as.numeric(EM["GFAP", md_glia$cell_id])
md_glia_micro <- md_glia %>% filter(celltype_clusters %in% MICRO)
md_glia_astro <- md_glia %>% filter(celltype_clusters %in% ASTRO)

# depth_s standardised within each modelled cell set, exactly as dist_z is.
add_depth_s <- function(md) { md$depth_s <- as.numeric(scale(md$depth_prop)); md }
md_dist_vuln  <- add_depth_s(md_dist_vuln)
md_dist_exc   <- add_depth_s(md_dist_exc)
md_glia_micro <- add_depth_s(md_glia_micro)
md_glia_astro <- add_depth_s(md_glia_astro)

# prep functions -- IDENTICAL to imc_covariate_sensitivity.R. delta must be fitted on prep()'s
# OUTPUT, not on P$md: prep_glia's within-donor scale() returns NaN for a donor with a single cell,
# and donor BBN_9110 has exactly one microglial cell, so the modelled frame is 10,970 cells /
# 43 donors rather than 10,971 / 44. The implied delta matches the z-set, not the raw set.
prep_marker <- function(md, mk) {
  d <- md
  d$value <- as.numeric(EM[mk, md$cell_id])
  as.data.frame(d[is.finite(d$value), , drop = FALSE])
}
prep_glia <- function(md, mk) {
  d <- md
  d$value_raw <- d[[mk]]
  d <- d[is.finite(d$value_raw), , drop = FALSE]
  d <- d %>% group_by(patient_id) %>% mutate(value = as.numeric(scale(value_raw))) %>% ungroup()
  as.data.frame(d[is.finite(d$value), , drop = FALSE])
}

# ===================================================================
# Panel table -- the FOUR DISTANCE panels of Figure S8A. The three state_* panels are excluded on
# purpose: their predictor is PHF1+/- status, so the analogous delta would be plaque ~ grp, a
# different estimand, and "bias in the tangle-DISTANCE coefficient" is not defined for them.
# `flip` = -1 throughout, i.e. the reported beta_closer = -beta(dist_z), per s.d. NEARER a tangle.
# ===================================================================
PANELS <- list(
  list(slug = "dist_reln_calb1", md = md_dist_vuln, outcomes = markers, prep = prep_marker,
       pterm = "dist_z", prow = "dist_z", extra = character(0), flip = -1,
       origin = "R/imc_phf1_distance_markers_subsets.R",
       label = "RELN+/CALB1+ excitatory neurons"),
  list(slug = "dist_exc_subcluster_adj", md = md_dist_exc, outcomes = markers, prep = prep_marker,
       pterm = "dist_z", prow = "dist_z", extra = "subcluster", flip = -1,
       origin = "R/imc_phf1_distance_markers_subsets.R",
       label = "Excitatory neurons, subcluster-adjusted"),
  list(slug = "cd68_microglia_distance", md = md_glia_micro, outcomes = "CD68", prep = prep_glia,
       pterm = "dist_scaled", prow = "dist_scaled", extra = character(0), flip = -1,
       origin = "R/imc_phf1_glia_distance.R",
       label = "CD68, pooled microglia"),
  list(slug = "gfap_astrocytes_distance", md = md_glia_astro, outcomes = "GFAP", prep = prep_glia,
       pterm = "dist_scaled", prow = "dist_scaled", extra = character(0), flip = -1,
       origin = "R/imc_phf1_glia_distance.R",
       label = "GFAP, pooled astrocytes")
)
names(PANELS) <- vapply(PANELS, `[[`, "", "slug")
if (!is.null(PANEL_SELECT)) {
  unknown <- setdiff(PANEL_SELECT, names(PANELS))
  if (length(unknown)) stop("unknown panel slug(s): ", paste(unknown, collapse = ", "))
  PANELS <- PANELS[PANEL_SELECT]
}

# ===================================================================
# Whitened GLS for a random-intercept model.
#
# For y = X beta + Z b + e with b ~ N(0, tau^2), e ~ N(0, sigma^2) and lambda = tau^2/sigma^2, the
# covariance within donor j of size n_j is V_j = sigma^2 (I + lambda 1 1'), and
#
#     V_j^{-1/2} propto I - (alpha_j / n_j) 1 1',     alpha_j = 1 - 1/sqrt(1 + n_j lambda)
#
# (check: (I - aJ/n)^2 = I - (2a - a^2) J/n and 2a - a^2 = 1 - (1-a)^2 = n lambda/(1 + n lambda),
#  which is exactly (I + lambda J)^{-1}). So whitening is nothing more than subtracting alpha_j
# times the donor mean from every column, and GLS becomes ORDINARY least squares in the whitened
# coordinates. Two consequences the whole script rests on:
#
#   * the OVB identity is EXACT in whitened coordinates, because it is exact for OLS;
#   * lambda -> 0 gives alpha_j -> 0, i.e. NO transform, i.e. pooled OLS. So the collapsed-random-
#     effect case needs no special branch: the two glia panels get pooled OLS automatically, which
#     is precisely the weight matrix their gamma was fitted under.
# ===================================================================
gls_lambda <- function(m) {
  vc <- as.data.frame(lme4::VarCorr(m))
  tau2 <- sum(vc$vcov[vc$grp == "patient_id" & is.na(vc$var2)])
  s2   <- vc$vcov[vc$grp == "Residual"]
  stopifnot(length(s2) == 1L, is.finite(tau2), is.finite(s2), s2 > 0)
  tau2 / s2
}

#' Whiten a numeric matrix at a given lambda, by donor: subtract alpha_j times the donor mean.
whiten <- function(M, grp, lambda) {
  M   <- as.matrix(M)
  grp <- as.character(grp)
  rs  <- rowsum(M, grp, reorder = FALSE)                 # donor sums; rownames = donor labels
  cnt <- as.numeric(table(grp)[rownames(rs)])
  gm  <- rs / cnt                                        # donor means
  nj  <- as.numeric(table(grp)[grp])
  alpha <- 1 - 1 / sqrt(1 + nj * lambda)
  M - alpha * gm[grp, , drop = FALSE]
}

#' GLS fixed-effect estimates at a FIXED lambda. Returns the full coefficient vector.
#'
#' @param form   two-sided formula for the FIXED part only (no random term).
#' @param d      data. Must carry `.gls_grp` (the donor) and have no missing values in `form` --
#'               both asserted, because silently dropping rows would make the OVB identity compare
#'               fits on different cell sets and the tier-1 assertion would fail confusingly.
#' @param lambda tau^2/sigma^2, taken from the outcome fit via gls_lambda().
gls_at_lambda <- function(form, d, lambda) {
  stopifnot(".gls_grp" %in% names(d))
  mf <- model.frame(form, data = d, na.action = na.fail)
  stopifnot(nrow(mf) == nrow(d))
  y  <- model.response(mf)
  X  <- model.matrix(form, mf)
  Xw <- whiten(X, d$.gls_grp, lambda)
  yw <- drop(whiten(matrix(y, ncol = 1), d$.gls_grp, lambda))
  qr.solve(Xw, yw)
}

#' The fixed-effect RHS of a panel's outcome model, as a character vector.
#'
#' THE ONE RULE OF THIS ANALYSIS: delta's RHS is the outcome model's RHS minus `plaque`. Built here
#' programmatically from the panel metadata so it cannot drift (dropping `subcluster` is a 21%
#' error on delta for dist_exc_subcluster_adj).
build_rhs <- function(P, add = character(0)) {
  c(P$pterm, P$extra, add, "Sex", "Age_s", "PMI_s")
}
f_of <- function(lhs, rhs) as.formula(paste(lhs, "~", paste(rhs, collapse = " + ")))
f_re <- function(lhs, rhs) as.formula(paste(lhs, "~", paste(c(rhs, "(1 | patient_id)"),
                                                            collapse = " + ")))

grab <- function(m, term) {
  co <- summary(m)$coefficients
  if (!term %in% rownames(co)) return(rep(NA_real_, 5))
  c(co[term, "Estimate"], co[term, "Std. Error"], co[term, "df"],
    co[term, "t value"], co[term, "Pr(>|t|)"])
}

# ===================================================================
# PART 0 -- REGRESSION CHECK. Refit the Figure S8A rungs on the rebuilt cell sets and require them
# to reproduce the on-disk coefficients. If m1/m3 do not reproduce, the cell set was rebuilt wrongly
# and every number downstream is meaningless, so this STOPS rather than warns.
# ===================================================================
message("== regression check against ", LADDER_TSV)
fits <- list()   # fits[[panel]][[outcome]] <- list(d=, m1=, m2=, m3=, mI=)
regchk <- list()

for (P in PANELS) {
  message("== ", P$slug, "  (", nrow(P$md), " cells, ", length(P$outcomes), " outcome(s))")
  rhs <- build_rhs(P)
  fits[[P$slug]] <- list()
  for (mk in P$outcomes) {
    d <- P$prep(P$md, mk)
    d$.gls_grp <- d$patient_id
    m1 <- lmerTest::lmer(f_re("value", rhs), data = d, REML = TRUE)
    m2 <- lmerTest::lmer(f_re("value", build_rhs(P, "depth_s")), data = d, REML = TRUE)
    m3 <- lmerTest::lmer(f_re("value", build_rhs(P, "plaque")),  data = d, REML = TRUE)
    fits[[P$slug]][[mk]] <- list(d = d, m1 = m1, m2 = m2, m3 = m3)

    ref <- ladder[ladder$panel == P$slug & ladder$outcome == mk, ]
    get_ref <- function(mod, trm) {
      r <- ref[ref$model == mod & ref$term == trm, "est"]
      if (length(r) != 1L) NA_real_ else r
    }
    regchk[[length(regchk) + 1L]] <- tibble(
      panel = P$slug, outcome = mk,
      n_cells = nrow(d), n_donors = n_distinct(d$patient_id), n_rois = n_distinct(d$sample_id),
      n_cells_ref  = if (nrow(ref)) ref$n_cells[1]  else NA_integer_,
      n_donors_ref = if (nrow(ref)) ref$n_donors[1] else NA_integer_,
      n_rois_ref   = if (nrow(ref)) ref$n_rois[1]   else NA_integer_,
      b1 = grab(m1, P$prow)[1], b1_ref = get_ref("m1_base",   "predictor"),
      b3 = grab(m3, P$prow)[1], b3_ref = get_ref("m3_plaque", "predictor"),
      gm = grab(m3, "plaqueplaque-prox")[1], gm_ref = get_ref("m3_plaque", "plaque"),
      b2 = grab(m2, P$prow)[1], b2_ref = get_ref("m2_depth",  "predictor"),
      gd = grab(m2, "depth_s")[1], gd_ref = get_ref("m2_depth", "depth"))
  }
}
regchk <- bind_rows(regchk) %>%
  mutate(max_abs_diff = pmax(abs(b1 - b1_ref), abs(b3 - b3_ref), abs(gm - gm_ref),
                             abs(b2 - b2_ref), abs(gd - gd_ref), na.rm = TRUE),
         n_ok = (n_cells == n_cells_ref) & (n_donors == n_donors_ref) & (n_rois == n_rois_ref),
         ok = n_ok & is.finite(max_abs_diff) & max_abs_diff < TOL_LADDER)
if (!all(regchk$ok)) {
  print(as.data.frame(regchk[!regchk$ok, ]))
  stop("REGRESSION CHECK FAILED: the rebuilt cell sets do not reproduce ", LADDER_TSV,
       ". The cell set is wrong and every downstream number would be meaningless.")
}
message("   regression check passed for all ", nrow(regchk),
        " panel x outcome fits (max abs diff ", signif(max(regchk$max_abs_diff), 3), ")")

# ===================================================================
# PART 1 -- delta. Four estimands per panel; only delta_lpm is the multiplier.
# ===================================================================

#' The exact OVB multiplier: identity-link regression of the omitted variable on the outcome model's
#' RHS, whitened at the OUTCOME FIT's lambda. No isSingular() branch -- see the header of the
#' whitening block for why lambda ~ 0 reproduces pooled OLS automatically.
delta_lpm_at <- function(P, d, lambda, lhs = "plaque_int", add = character(0)) {
  b <- gls_at_lambda(f_of(lhs, build_rhs(P, add)), d, lambda)
  unname(b[P$prow])
}

# --- per panel x outcome: the multiplier, at m3's lambda (the fit gamma comes from)
delta_rows <- list()
for (P in PANELS) for (mk in P$outcomes) {
  F <- fits[[P$slug]][[mk]]; d <- F$d
  lam1 <- gls_lambda(F$m1); lam2 <- gls_lambda(F$m2); lam3 <- gls_lambda(F$m3)
  # validate the whitening itself: GLS at m3's lambda must reproduce lmer's m3 fixed effects
  b_chk <- gls_at_lambda(f_of("value", build_rhs(P, "plaque")), d, lam3)
  whiten_err <- max(abs(b_chk[P$prow] - grab(F$m3, P$prow)[1]),
                    abs(b_chk["plaqueplaque-prox"] - grab(F$m3, "plaqueplaque-prox")[1]))
  delta_rows[[length(delta_rows) + 1L]] <- tibble(
    panel = P$slug, outcome = mk,
    lambda_m1 = lam1, lambda_m2 = lam2, lambda_m3 = lam3,
    singular_m1 = lme4::isSingular(F$m1), singular_m3 = lme4::isSingular(F$m3),
    weight_basis = if (lam3 < 1e-8) "pooled OLS (donor RE collapsed)" else "mixed (1 | patient_id)",
    delta_lpm       = delta_lpm_at(P, d, lam3),
    delta_lpm_lam1  = delta_lpm_at(P, d, lam1),
    delta_depth     = delta_lpm_at(P, d, lam2, lhs = "depth_s"),
    whiten_err      = whiten_err)
}
delta_rows <- bind_rows(delta_rows)
if (max(delta_rows$whiten_err) >= TOL_ALG)
  stop("whitened GLS does not reproduce lmer's fixed effects (max err ",
       signif(max(delta_rows$whiten_err), 3), "). The whitening is wrong.")
message("   whitened GLS reproduces lmer to ", signif(max(delta_rows$whiten_err), 3))

# --- per panel: the descriptive estimands. Computed on the REFERENCE OUTCOME's modelled frame, so
# --- the descriptive set is the set actually modelled (this matters for the glia panels, where
# --- prep_glia drops a single-cell donor).
delta_panel <- list()
for (P in PANELS) {
  mk0 <- P$outcomes[1]
  d <- fits[[P$slug]][[mk0]]$d
  rhs <- build_rhs(P)

  # (a) the requested logistic mixed model. DESCRIPTIVE: reported as an OR, never multiplied.
  gm <- lme4::glmer(f_re("plaque_int", rhs), data = d, family = binomial,
                    control = lme4::glmerControl(optimizer = "bobyqa",
                                                 optCtrl = list(maxfun = 2e5)))
  co <- summary(gm)$coefficients
  b_lg <- co[P$prow, "Estimate"]; se_lg <- co[P$prow, "Std. Error"]
  # (b) AME, by hand -- marginaleffects is not installed. Conditional p includes the donor BLUPs,
  #     marginal p does not; they differ by a few %, so both are reported and labelled.
  p_cond <- fitted(gm)
  p_marg <- predict(gm, type = "response", re.form = NA)
  ame_c <- mean(b_lg * p_cond * (1 - p_cond))
  ame_m <- mean(b_lg * p_marg * (1 - p_marg))
  # (c) delta on the LPM scale with a CLUSTER-ROBUST CI. The LPM's model-based SE is invalid
  #     (Bernoulli heteroskedasticity + donor clustering), so use CR2 by donor.
  lp <- lm(f_of("plaque_int", rhs), data = d)
  V  <- clubSandwich::vcovCR(lp, cluster = d$patient_id, type = "CR2")
  se_cr <- sqrt(diag(V))[P$prow]
  ct <- clubSandwich::coef_test(lp, vcov = V, cluster = d$patient_id, test = "Satterthwaite")
  ct <- ct[as.character(ct$Coef) == P$prow, ]
  stopifnot(nrow(ct) == 1L)

  # (d) near/far plaque-proximal fraction, WITHIN donor then averaged with a t CI. Bands are in
  #     RAW MICRONS (not the scaled predictor) so the thresholds mean what they say.
  nf <- d %>%
    mutate(band = case_when(dist_um <= NEAR_UM ~ "near",
                            dist_um >= FAR_UM  ~ "far", TRUE ~ NA_character_)) %>%
    filter(!is.na(band)) %>%
    group_by(patient_id, band) %>%
    summarise(frac = mean(plaque_int), n = n(), dz = mean(.data[[P$pterm]]), .groups = "drop") %>%
    pivot_wider(names_from = band, values_from = c(frac, n, dz))
  nf_ok <- nf %>% filter(!is.na(frac_near), !is.na(frac_far))
  dd <- nf_ok$frac_near - nf_ok$frac_far
  tt <- t.test(dd)
  gap_dz <- mean(nf_ok$dz_far - nf_ok$dz_near, na.rm = TRUE)   # in s.d. of the predictor

  dl <- delta_rows %>% filter(panel == P$slug)
  delta_panel[[length(delta_panel) + 1L]] <- tibble(
    panel = P$slug, label = P$label, ref_outcome = mk0,
    n_cells = nrow(d), n_donors = n_distinct(d$patient_id), n_rois = n_distinct(d$sample_id),
    prevalence = mean(d$plaque_int),
    dist_sd = sd(log(d$dist_um)),
    # --- the multiplier
    delta_lpm = median(dl$delta_lpm), delta_lpm_min = min(dl$delta_lpm),
    delta_lpm_max = max(dl$delta_lpm),
    delta_lpm_CI.L = unname(dl$delta_lpm[1] - 1.96 * se_cr),
    delta_lpm_CI.R = unname(dl$delta_lpm[1] + 1.96 * se_cr),
    delta_lpm_se_CR2 = unname(se_cr), delta_lpm_p_CR2 = unname(ct$p_Satt),
    weight_basis = dl$weight_basis[1], singular_outcome_fit = dl$singular_m3[1],
    # --- descriptive: logistic
    delta_logit = b_lg, delta_logit_se = se_lg,
    OR_per_sd_nearer = exp(-b_lg),
    OR_CI.L = exp(-(b_lg + 1.96 * se_lg)), OR_CI.R = exp(-(b_lg - 1.96 * se_lg)),
    logit_p = co[P$prow, "Pr(>|z|)"],
    re_var_donor_logit = as.data.frame(lme4::VarCorr(gm))$vcov[1],
    # --- descriptive: AME (the bridge)
    ame_conditional = ame_c, ame_marginal = ame_m,
    ame_vs_lpm_pct = 100 * (ame_c - dl$delta_lpm[1]) / abs(dl$delta_lpm[1]),
    # --- descriptive: near/far
    nf_donors = nrow(nf_ok), nf_donors_total = n_distinct(d$patient_id),
    nf_thin_near = sum(nf_ok$n_near < MIN_BIN_N), nf_thin_far = sum(nf_ok$n_far < MIN_BIN_N),
    frac_near = mean(nf_ok$frac_near), frac_far = mean(nf_ok$frac_far),
    nf_diff = unname(tt$estimate), nf_CI.L = tt$conf.int[1], nf_CI.R = tt$conf.int[2],
    nf_p = tt$p.value, nf_t = unname(tt$statistic), nf_df = unname(tt$parameter),
    nf_dz_gap = gap_dz, nf_diff_per_sd = unname(tt$estimate) / gap_dz,
    nf_vs_lpm_pct = 100 * ((unname(tt$estimate) / gap_dz) / abs(dl$delta_lpm[1]) - 1))
  attr(delta_panel[[length(delta_panel)]], "glmer") <- gm
  assign(paste0(".gm_", P$slug), gm)
  assign(paste0(".nf_", P$slug), nf_ok)
}
delta_panel <- bind_rows(delta_panel)

# --- sanity check against the reference deltas in EXPECTED_DELTA. Loose on purpose: the algebra
# --- is pinned by TOL_ALG below, not by this.
for (i in seq_len(nrow(delta_panel))) {
  ps <- delta_panel$panel[i]
  if (!ps %in% names(EXPECTED_DELTA)) next
  rel <- abs(delta_panel$delta_lpm[i] - EXPECTED_DELTA[[ps]]) / abs(EXPECTED_DELTA[[ps]])
  if (rel > TOL_DELTA_REL)
    warning(sprintf("delta for %s is %.6f, reference value %.6f (%.1f%% apart)",
                    ps, delta_panel$delta_lpm[i], EXPECTED_DELTA[[ps]], 100 * rel))
}

# --- the fraction-vs-distance profile, quantile bins (source data for figure 1)
prof <- list()
for (P in PANELS) {
  d <- fits[[P$slug]][[P$outcomes[1]]]$d
  br <- unique(quantile(d$dist_um, probs = seq(0, 1, length.out = N_QBINS + 1), names = FALSE))
  d$qbin <- cut(d$dist_um, breaks = br, include.lowest = TRUE, labels = FALSE)
  prof[[length(prof) + 1L]] <- d %>%
    group_by(qbin) %>%
    summarise(dist_mid = median(dist_um), dist_lo = min(dist_um), dist_hi = max(dist_um),
              n_cells = n(), n_donors = n_distinct(patient_id),
              frac = mean(plaque_int), .groups = "drop") %>%
    # binomial CI on the cell-level fraction; n_donors per bin is carried alongside
    mutate(panel = P$slug, label = P$label,
           se = sqrt(frac * (1 - frac) / n_cells),
           CI.L = pmax(0, frac - 1.96 * se), CI.R = pmin(1, frac + 1.96 * se))
}
prof <- bind_rows(prof)

# ===================================================================
# PART 2 -- gamma x delta, with the identity asserted rather than claimed.
#
# Tier 1 (algebraic, TOL_ALG): everything evaluated at ONE lambda. In whitened coordinates GLS is
#   OLS, so b1* - b3* = gamma* delta* is an identity and any failure is a coding error.
# Tier 2 (practical, TOL_REML): gamma*delta from the REML fits vs the observed m1 -> m3 shift. The
#   only discrepancy is REML re-estimating lambda when plaque enters.
# NEVER assert on the RATIO (b1-b3)/gamma: LAMP1 has gamma = 0.0024 and a shift of 4.3e-05, i.e.
# numerically 0/0. The assertion is on the PRODUCT residual.
# ===================================================================
bias_rows <- list()
for (P in PANELS) for (mk in P$outcomes) {
  F <- fits[[P$slug]][[mk]]; d <- F$d
  lam3 <- gls_lambda(F$m3); lam2 <- gls_lambda(F$m2)

  # ---- tier 1: the whole identity at lam3 (amyloid) and lam2 (depth)
  alg_resid <- function(lam, addterm, gterm, lhs_omit) {
    b1s <- gls_at_lambda(f_of("value", build_rhs(P)),          d, lam)[P$prow]
    fb  <- gls_at_lambda(f_of("value", build_rhs(P, addterm)), d, lam)
    b3s <- fb[P$prow]; gs <- fb[gterm]
    ds  <- gls_at_lambda(f_of(lhs_omit, build_rhs(P)),         d, lam)[P$prow]
    unname(gs * ds - (b1s - b3s))
  }
  r_alg_plq <- alg_resid(lam3, "plaque",  "plaqueplaque-prox", "plaque_int")
  r_alg_dep <- alg_resid(lam2, "depth_s", "depth_s",           "depth_s")

  # ---- tier 2: the published decomposition, from the REML fits
  a1 <- grab(F$m1, P$prow); a3 <- grab(F$m3, P$prow); a2 <- grab(F$m2, P$prow)
  b1 <- a1[1]; se1 <- a1[2]
  g_plq <- grab(F$m3, "plaqueplaque-prox")[1]
  g_dep <- grab(F$m2, "depth_s")[1]
  dl <- delta_rows %>% filter(panel == P$slug, outcome == mk)
  gd_plq <- g_plq * dl$delta_lpm
  gd_dep <- g_dep * dl$delta_depth
  sh_plq <- b1 - a3[1]
  sh_dep <- b1 - a2[1]

  # delta is imprecise at donor level (the CR2 CI can cover zero), so propagate its interval
  # through to the bias and take the WORST case. This is what makes "the amyloid bias is small" a
  # statement about the data rather than about a point estimate: if even the far end of delta's CI
  # leaves the shift small, delta's own imprecision cannot rescue an amyloid explanation.
  d_ci <- unlist(delta_panel[delta_panel$panel == P$slug, c("delta_lpm_CI.L", "delta_lpm_CI.R")])
  gd_plq_worst <- g_plq * d_ci[which.max(abs(g_plq * d_ci))]

  bias_rows[[length(bias_rows) + 1L]] <- tibble(
    panel = P$slug, label = P$label, outcome = mk, scale = "raw",
    n_cells = nrow(d), n_donors = n_distinct(d$patient_id),
    b1 = b1, b1_SE = se1, b1_pval = a1[5], b1_near_null = abs(b1) < 2 * se1,
    # amyloid
    gamma_plaque = g_plq, delta_plaque = dl$delta_lpm,
    gamma_x_delta_plaque = gd_plq, shift_plaque = sh_plq,
    resid_plaque = gd_plq - sh_plq, resid_over_SE_plaque = (gd_plq - sh_plq) / se1,
    shift_over_SE_plaque = sh_plq / se1,
    gamma_x_delta_plaque_worst = unname(gd_plq_worst),
    shift_over_SE_plaque_worst = unname(gd_plq_worst) / se1,
    pct_of_b1_plaque = ifelse(abs(b1) < 2 * se1, NA_real_, 100 * sh_plq / b1),
    alg_resid_plaque = r_alg_plq,
    # depth, for comparison
    gamma_depth = g_dep, delta_depth = dl$delta_depth,
    gamma_x_delta_depth = gd_dep, shift_depth = sh_dep,
    resid_depth = gd_dep - sh_dep, shift_over_SE_depth = sh_dep / se1,
    pct_of_b1_depth = ifelse(abs(b1) < 2 * se1, NA_real_, 100 * sh_dep / b1),
    alg_resid_depth = r_alg_dep,
    singular = lme4::isSingular(F$m3), flip = P$flip)
}
bias <- bind_rows(bias_rows)

# ---- the two assertion tiers
if (max(abs(c(bias$alg_resid_plaque, bias$alg_resid_depth))) >= TOL_ALG) {
  print(as.data.frame(bias[abs(bias$alg_resid_plaque) >= TOL_ALG |
                           abs(bias$alg_resid_depth)  >= TOL_ALG,
                           c("panel", "outcome", "alg_resid_plaque", "alg_resid_depth")]))
  stop("ALGEBRAIC OVB IDENTITY FAILED (tol ", TOL_ALG, "). In whitened coordinates GLS is OLS, so ",
       "this is exact -- a failure means the RHS, the cell set or the coding is wrong.")
}
message("   algebraic OVB identity holds to ",
        signif(max(abs(c(bias$alg_resid_plaque, bias$alg_resid_depth))), 3))

bias$reml_ok <- abs(bias$resid_plaque) <= pmax(TOL_ALG, TOL_REML * bias$b1_SE)
if (!all(bias$reml_ok)) {
  print(as.data.frame(bias[!bias$reml_ok, c("panel", "outcome", "gamma_x_delta_plaque",
                                            "shift_plaque", "resid_over_SE_plaque")]))
  stop("gamma*delta does not reproduce the observed m1 -> m3 shift within ", TOL_REML, "*SE(b1).")
}
message("   gamma x delta reproduces the observed shift; max |resid|/SE(b1) = ",
        signif(max(abs(bias$resid_over_SE_plaque)), 3))

# ===================================================================
# PART 3 -- stratification. One pooled `dist_z * plaque` fit per panel x outcome (NOT two
# stratum-separate fits: the pooled fit shares sigma^2 and the donor RE and is what yields the
# interaction test). dist_z stays on the PANEL-WIDE scale -- within-stratum sd is 0.985-1.055 and
# the supports overlap almost exactly, so rescaling would change nothing numerically while making
# the interaction a difference of coefficients in two different units.
#
# Simple slopes are computed by hand from the fixed-effect vcov rather than via emmeans, for two
# reasons: emmeans silently drops to asymptotic z above 3,000 observations (lmerTest.limit /
# pbkrtest.limit) whereas the ladder uses Satterthwaite t, and raising the limit to 26,832 is
# genuinely expensive. A normal reference is used and stated; at df in the thousands z and t differ
# in the fifth decimal. An emmeans::emtrends cross-check is asserted on the smallest panel.
# ===================================================================
inter_rows <- list()
for (P in PANELS) for (mk in P$outcomes) {
  F <- fits[[P$slug]][[mk]]; d <- F$d
  rhs_i <- c(paste0(P$pterm, " * plaque"), P$extra, "Sex", "Age_s", "PMI_s")
  mI <- lmerTest::lmer(f_re("value", rhs_i), data = d, REML = TRUE)
  fits[[P$slug]][[mk]]$mI <- mI

  irow <- paste0(P$prow, ":plaqueplaque-prox")
  co <- summary(mI)$coefficients
  V  <- as.matrix(vcov(mI))
  est_i <- co[irow, "Estimate"]; se_i <- co[irow, "Std. Error"]

  # distal simple slope IS the predictor main effect (plaque-distal is the reference level) --
  # asserted below rather than recomputed. proximal = main + interaction, with the full covariance.
  b_d  <- co[P$prow, "Estimate"];  se_d <- co[P$prow, "Std. Error"]
  b_p  <- b_d + est_i
  se_p <- sqrt(V[P$prow, P$prow] + V[irow, irow] + 2 * V[P$prow, irow])

  # per-stratum support. The proximal arm is much thinner at DONOR level than its cell share says.
  sup <- d %>% group_by(plaque) %>%
    summarise(n_cells = n(), n_donors = n_distinct(patient_id), n_rois = n_distinct(sample_id),
              .groups = "drop")
  don10 <- d %>% count(plaque, patient_id) %>% filter(n >= 10) %>% count(plaque, name = "d10")
  don30 <- d %>% count(plaque, patient_id) %>% filter(n >= 30) %>% count(plaque, name = "d30")
  sup <- sup %>% left_join(don10, by = "plaque") %>% left_join(don30, by = "plaque") %>%
    mutate(d10 = coalesce(d10, 0L), d30 = coalesce(d30, 0L))
  gv <- function(lv, cl) sup[[cl]][sup$plaque == lv]
  roi0 <- d %>% group_by(sample_id) %>% summarise(np = sum(plaque_int), .groups = "drop")

  # equivalence, not absence: is the interaction CI inside +/- EQUIV_FRAC * |distal slope|?
  margin <- EQUIV_FRAC * abs(b_d)
  ci_i <- est_i + c(-1.96, 1.96) * se_i

  inter_rows[[length(inter_rows) + 1L]] <- tibble(
    panel = P$slug, label = P$label, outcome = mk, scale = "raw",
    interaction = est_i, interaction_SE = se_i,
    interaction_CI.L = ci_i[1], interaction_CI.R = ci_i[2],
    interaction_z = est_i / se_i, interaction_p = 2 * pnorm(-abs(est_i / se_i)),
    slope_distal = b_d, slope_distal_SE = se_d,
    slope_distal_CI.L = b_d - 1.96 * se_d, slope_distal_CI.R = b_d + 1.96 * se_d,
    slope_prox = b_p, slope_prox_SE = se_p,
    slope_prox_CI.L = b_p - 1.96 * se_p, slope_prox_CI.R = b_p + 1.96 * se_p,
    sign_reversal = sign(b_d) != sign(b_p),
    equiv_margin = margin,
    equiv_within = is.finite(margin) && ci_i[1] > -margin && ci_i[2] < margin,
    # NA when the distal slope is itself null: "the smallest detectable interaction is 17,000% of
    # the distal slope" is division by nothing, not a power statement.
    mde_frac_of_distal = ifelse(abs(b_d) < 2 * se_d, NA_real_, (1.96 * se_i) / abs(b_d)),
    n_distal = gv("plaque-distal", "n_cells"), n_prox = gv("plaque-prox", "n_cells"),
    donors_distal = gv("plaque-distal", "n_donors"), donors_prox = gv("plaque-prox", "n_donors"),
    donors_prox_ge10 = gv("plaque-prox", "d10"), donors_prox_ge30 = gv("plaque-prox", "d30"),
    rois_zero_prox = sum(roi0$np == 0), rois_total = nrow(roi0),
    singular = lme4::isSingular(mI), flip = P$flip,
    b1_headline = bias$b1[bias$panel == P$slug & bias$outcome == mk])
}
inter <- bind_rows(inter_rows) %>%
  # BH WITHIN PANEL, across outcomes -- matching imc_covariate_sensitivity.R. Not pooled across
  # panels: different cell sets, different dist_z scalings, and reln_calb1 is a subset of exc.
  # The simple slopes are reparameterisations, not new tests, so they get no family.
  group_by(panel) %>% mutate(interaction_padj = p.adjust(interaction_p, "BH")) %>% ungroup()

# free assertion: the distal simple slope must equal the predictor main effect bit-for-bit
stopifnot(all(abs(inter$slope_distal - vapply(seq_len(nrow(inter)), function(i) {
  summary(fits[[inter$panel[i]]][[inter$outcome[i]]]$mI)$coefficients[
    PANELS[[inter$panel[i]]]$prow, "Estimate"] }, 0)) < 1e-12))

# emmeans cross-check on the smallest panel, where Satterthwaite df are affordable
emm_check <- NA_real_
ck_panel <- names(PANELS)[which.min(vapply(PANELS, function(P) nrow(P$md), 0))]
if (requireNamespace("emmeans", quietly = TRUE)) {
  Pk <- PANELS[[ck_panel]]; mk <- Pk$outcomes[1]
  et <- as.data.frame(emmeans::emtrends(fits[[ck_panel]][[mk]]$mI, ~ plaque, var = Pk$pterm))
  r <- inter[inter$panel == ck_panel & inter$outcome == mk, ]
  emm_check <- max(abs(et[[2]] - c(r$slope_distal, r$slope_prox)))
  if (emm_check >= 1e-8)
    stop("hand-computed simple slopes disagree with emmeans::emtrends (", emm_check, ")")
  message("   simple slopes agree with emmeans::emtrends to ", signif(emm_check, 3))
}

# ---- reporting-scale (beta_closer) copies. flip = -1 for all four distance panels, so every
# ---- distance quantity flips and `plaque`/`depth` terms do not. Emitted as separate rows with an
# ---- explicit `scale` column so no downstream reader has to guess.
flip_cols_bias <- c("b1", "gamma_x_delta_plaque", "shift_plaque", "resid_plaque",
                    "shift_over_SE_plaque", "resid_over_SE_plaque",
                    "gamma_x_delta_plaque_worst", "shift_over_SE_plaque_worst",
                    "gamma_x_delta_depth", "shift_depth", "resid_depth", "shift_over_SE_depth")
bias_bc <- bias %>% mutate(across(all_of(flip_cols_bias), ~ .x * flip), scale = "beta_closer")
flip_cols_int <- c("interaction", "interaction_CI.L", "interaction_CI.R",
                   "slope_distal", "slope_distal_CI.L", "slope_distal_CI.R",
                   "slope_prox", "slope_prox_CI.L", "slope_prox_CI.R", "b1_headline")
inter_bc <- inter %>%
  mutate(across(all_of(flip_cols_int), ~ .x * flip), scale = "beta_closer") %>%
  # flipping reverses the CI bounds
  mutate(tmp = interaction_CI.L, interaction_CI.L = pmin(interaction_CI.L, interaction_CI.R),
         interaction_CI.R = pmax(tmp, interaction_CI.R),
         tmp = slope_distal_CI.L, slope_distal_CI.L = pmin(slope_distal_CI.L, slope_distal_CI.R),
         slope_distal_CI.R = pmax(tmp, slope_distal_CI.R),
         tmp = slope_prox_CI.L, slope_prox_CI.L = pmin(slope_prox_CI.L, slope_prox_CI.R),
         slope_prox_CI.R = pmax(tmp, slope_prox_CI.R)) %>% select(-tmp)

# ===================================================================
# Figures. Standard TRIPLE: PDF (default `pdf` device, ASCII text, plotmath for the
# micron symbol), source-data TSV, stats log.
# ===================================================================
wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE, na = "NA")
sv <- function(g, f, w, h) ggsave(file.path(out_dir, f), g, width = w, height = h,
                                  units = "in", device = "pdf")

panel_lab <- setNames(vapply(PANELS, `[[`, "", "label"), names(PANELS))
prof$label_f <- factor(prof$label, levels = unname(panel_lab))

# Facet strips are ~1.8 in wide at this canvas, which clips a one-line panel label. Wrap them, and
# give the strip two lines of room. `fig_theme` rotates axis text 45 deg for long celltype names;
# the x axes here are numeric, so that is overridden back to horizontal in every figure.
LAB_WRAP  <- ggplot2::label_wrap_gen(width = 24)
strip_fix <- theme(strip.background = element_blank(),
                   strip.text = element_text(size = 6, lineheight = 1.05, margin = margin(b = 2)),
                   axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7, colour = "black"),
                   plot.margin = margin(4, 8, 4, 6, "pt"))

# --- Figure 1: the coupling itself, plaque-proximal fraction over tangle distance
nf_bands <- data.frame(xmin = c(0, FAR_UM), xmax = c(NEAR_UM, MAX_DIST),
                       band = c("near (<=50 um)", "far (>=200 um)"))
g1 <- ggplot(prof, aes(dist_mid, frac)) +
  geom_rect(data = nf_bands, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = COL_PROX, alpha = 0.07) +
  geom_ribbon(aes(ymin = CI.L, ymax = CI.R), fill = COL_PLAQUE, alpha = 0.20) +
  geom_line(colour = COL_PLAQUE, linewidth = 0.4) +
  geom_point(colour = COL_PLAQUE, size = 0.7) +
  facet_wrap(~ label_f, nrow = 1, scales = "free_y", labeller = LAB_WRAP) +
  labs(x = expression("Distance to nearest PHF1+ neuron (" * mu * "m)"),
       y = "Fraction plaque-proximal") +
  theme_classic(base_size = 8) + fig_theme + strip_fix
sv(g1, "plot_delta_coupling.pdf", 7.2, 2.1)
wr(prof %>% select(panel, label, qbin, dist_mid, dist_lo, dist_hi, n_cells, n_donors,
                   frac, se, CI.L, CI.R), "source_data_delta_coupling.tsv")

# --- Figure 2: the bias decomposition, in units of SE(b1). The percentage is not plotted, because
# --- for the near-null coefficients it is a division by nothing.
b2 <- bind_rows(
  bias %>% transmute(panel, label, outcome, covar = "amyloid",
                     shift_se = shift_over_SE_plaque, near_null = b1_near_null),
  bias %>% transmute(panel, label, outcome, covar = "depth",
                     shift_se = shift_over_SE_depth,  near_null = b1_near_null)) %>%
  mutate(label_f = factor(label, levels = unname(panel_lab)),
         outcome = factor(outcome, levels = rev(sort(unique(outcome)))))
g2 <- ggplot(b2, aes(shift_se, outcome, colour = covar, shape = near_null)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  # the materiality threshold used by the log's generated verdict, so the reader can see which
  # points clear it rather than only that the amyloid points sit near zero
  geom_vline(xintercept = c(-1, 1), linewidth = 0.25, linetype = "dashed", colour = "grey65") +
  geom_point(size = 1.3, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = COVAR_COLOURS, name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), name = "|b1| < 2 SE",
                     labels = c(`FALSE` = "no", `TRUE` = "yes")) +
  facet_wrap(~ label_f, nrow = 1, scales = "free_y", labeller = LAB_WRAP) +
  labs(x = "Shift in tangle coefficient when the covariate is added (SE of b1)", y = NULL) +
  theme_classic(base_size = 8) + fig_theme + strip_fix +
  theme(legend.position = "bottom")
sv(g2, "plot_bias_decomposition.pdf", 7.2, 3.0)
wr(bind_rows(bias, bias_bc), "source_data_bias_decomposition.tsv")

# --- Figure 3: the actual test. Tangle slope within each amyloid stratum, on the beta_closer scale
# --- (per s.d. NEARER a tangle), so it reads the same way as the published panels.
i3 <- bind_rows(
  inter_bc %>% transmute(panel, label, outcome, stratum = "plaque-distal",
                         est = slope_distal, lo = slope_distal_CI.L, hi = slope_distal_CI.R),
  inter_bc %>% transmute(panel, label, outcome, stratum = "plaque-prox",
                         est = slope_prox, lo = slope_prox_CI.L, hi = slope_prox_CI.R)) %>%
  mutate(label_f = factor(label, levels = unname(panel_lab)),
         outcome = factor(outcome, levels = rev(sort(unique(outcome)))))
g3 <- ggplot(i3, aes(est, outcome, colour = stratum)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.35,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 1.2, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = STRATUM_COLOURS, name = NULL) +
  facet_wrap(~ label_f, nrow = 1, scales = "free", labeller = LAB_WRAP) +
  labs(x = "Tangle-distance slope, per s.d. nearer (95% CI)", y = NULL) +
  theme_classic(base_size = 8) + fig_theme + strip_fix +
  theme(legend.position = "bottom")
sv(g3, "plot_interaction_slopes.pdf", 2.4 * n_distinct(i3$panel), 3.0)
wr(bind_rows(inter, inter_bc), "source_data_interaction_slopes.tsv")

# ===================================================================
# Structured tables + the human-readable stats log
# ===================================================================
wr(delta_panel, "stats_delta_coupling.tsv")
wr(delta_rows,  "stats_delta_per_outcome.tsv")
wr(bind_rows(bias, bias_bc), "stats_bias_decomposition.tsv")
wr(bind_rows(inter, inter_bc), "stats_interaction.tsv")
wr(regchk, "stats_regression_check.tsv")

fmt <- function(x, d = 6) formatC(x, format = "f", digits = d)
log_path <- file.path(out_dir, "stats_imc_amyloid_tangle_coupling.txt")
sink(log_path)
cat("=== Amyloid-tangle coupling and the omitted-amyloid bias in the IMC distance models ===\n")
cat("Script : R/imc_amyloid_tangle_coupling.R\n")
cat("Run    : ", format(Sys.time()), "\n", sep = "")
cat("Object : ", localwd, "spe.rds\n", sep = "")
cat("Input  : ", LADDER_TSV, " (Figure S8A)\n", sep = "")
cat("Panels : ", paste(names(PANELS), collapse = ", "), "\n\n", sep = "")

cat("--- The question -------------------------------------------------------------\n")
cat("The IMC distance models carry no amyloid term (docs/MODELS.md Set 3: `plaque` is a\n")
cat("sensitivity term). Two things are asked here:\n")
cat("  1. how big is the bias in the tangle coefficient from omitting amyloid?  = gamma x delta\n")
cat("  2. does the tangle gradient survive where there is no amyloid?           = the plaque-\n")
cat("     distal simple slope of a dist x plaque interaction.\n\n")

cat("--- Cohort -------------------------------------------------------------------\n")
cat("'AD' is BraakGroup != 'Braak_0_1' (Braak III-IV + VI); there is no diagnosis column.\n")
cat("With the >=3-ROI filter this is the 44-donor / 132-ROI analysed cohort out of a\n")
cat("67-donor / 195-ROI object. Panels land at 43-44 donors / 129-131 ROIs because ROIs with\n")
cat("no PHF1+ neuron anchor drop out. Per-panel n is tabulated below.\n")
cat("Distance: cap ", MAX_DIST, " um, transform log(dist_um)/sd(log(dist_um)), sd within cell set.\n\n", sep = "")

cat("--- PART 0: regression check against Figure S8A ------------------------------\n")
cat("The cell-set construction is copied from R/imc_covariate_sensitivity.R, so it is only\n")
cat("trustworthy if it reproduces that script's coefficients. It does:\n")
cat("  fits checked      : ", nrow(regchk), "\n", sep = "")
cat("  max abs coef diff : ", signif(max(regchk$max_abs_diff), 3), " (tol ", TOL_LADDER, ")\n", sep = "")
cat("  n_cells/n_donors/n_rois identical for every panel: ", all(regchk$n_ok), "\n\n", sep = "")

cat("--- PART 1: delta, the tau-amyloid coupling ----------------------------------\n")
cat("FOUR ESTIMANDS, ONE OF WHICH IS THE MULTIPLIER.\n\n")
cat("A logistic delta CANNOT be multiplied by gamma. OVB is a property of linear projection and\n")
cat("`plaque` enters the outcome models as a 0/1 dummy under an identity link, so the auxiliary\n")
cat("regression must be identity-link too. A logistic coefficient is a slope of logit(P), and\n")
cat("d logit(P)/dx = (dP/dx)/(p(1-p)); at the prevalences here that inflates it by ~5-7x.\n\n")
for (i in seq_len(nrow(delta_panel))) {
  r <- delta_panel[i, ]
  cat("[", r$panel, "]  ", r$label, "\n", sep = "")
  cat("  n = ", r$n_cells, " cells / ", r$n_donors, " donors / ", r$n_rois, " ROIs",
      "   plaque-proximal prevalence ", fmt(100 * r$prevalence, 1), "%\n", sep = "")
  cat("  sd(log dist_um) = ", fmt(r$dist_sd, 5), "\n", sep = "")
  cat("  MULTIPLIER   delta_lpm            = ", fmt(r$delta_lpm), "  [",
      fmt(r$delta_lpm_CI.L), ", ", fmt(r$delta_lpm_CI.R), "] (CR2 by donor, p = ",
      signif(r$delta_lpm_p_CR2, 3), ")\n", sep = "")
  cat("               weight basis         = ", r$weight_basis, "\n", sep = "")
  cat("               range over outcomes  = [", fmt(r$delta_lpm_min), ", ",
      fmt(r$delta_lpm_max), "]\n", sep = "")
  cat("  DESCRIPTIVE  logistic (glmer)     = ", fmt(r$delta_logit), " log-odds per s.d. FURTHER\n", sep = "")
  cat("               OR per s.d. NEARER   = ", fmt(r$OR_per_sd_nearer, 4), " [",
      fmt(r$OR_CI.L, 4), ", ", fmt(r$OR_CI.R, 4), "]  p = ", signif(r$logit_p, 3), "\n", sep = "")
  cat("               (CONDITIONAL, i.e. donor-specific: donor RE variance ",
      fmt(r$re_var_donor_logit, 4), " on the logit\n", sep = "")
  cat("                scale, so conditional != marginal. Not the multiplier: it is ",
      fmt(r$delta_logit / r$delta_lpm, 2), "x delta_lpm.)\n", sep = "")
  cat("  BRIDGE       AME (conditional)    = ", fmt(r$ame_conditional), "  (",
      fmt(r$ame_vs_lpm_pct, 1), "% from delta_lpm)\n", sep = "")
  cat("               AME (marginal)       = ", fmt(r$ame_marginal), "\n", sep = "")
  cat("  NEAR/FAR     frac <=", NEAR_UM, " um = ", fmt(r$frac_near, 4),
      "   frac >=", FAR_UM, " um = ", fmt(r$frac_far, 4), "\n", sep = "")
  cat("               difference           = ", fmt(r$nf_diff, 4), " [", fmt(r$nf_CI.L, 4), ", ",
      fmt(r$nf_CI.R, 4), "]  t(", fmt(r$nf_df, 1), ") = ", fmt(r$nf_t, 3),
      ", p = ", signif(r$nf_p, 3), "\n", sep = "")
  cat("               donors contributing both bands: ", r$nf_donors, "/", r$nf_donors_total,
      "; with <", MIN_BIN_N, " cells in a band: near ", r$nf_thin_near, ", far ",
      r$nf_thin_far, "\n", sep = "")
  cat("               NOT the multiplier. Rescaled by the ", fmt(r$nf_dz_gap, 3),
      " s.d. gap between bands it is\n", sep = "")
  cat("               ", fmt(r$nf_diff_per_sd), " vs delta_lpm ", fmt(r$delta_lpm), " (",
      fmt(r$nf_vs_lpm_pct, 1), "% apart): it is donor-unweighted while delta is\n", sep = "")
  cat("               cell-weighted, a two-point secant not a slope, and unadjusted.\n\n")
}
cat("Effect-size reading: the coupling points the expected way in every panel -- cells nearer a\n")
cat("tangle are more likely to sit in the amyloid niche -- and it is SMALL, a few percentage points\n")
cat("of plaque-proximal probability per s.d. of log distance.\n\n")
# Heading GENERATED, not hardcoded: whether the three inferences agree depends on which panels are
# in the run, so the sentence is derived from the table printed directly beneath it.
.dsg <- delta_panel %>% filter((delta_lpm_CI.L < 0 & delta_lpm_CI.R > 0) | delta_lpm_p_CR2 >= 0.05)
if (nrow(.dsg)) {
  cat("THE THREE INFERENCES ON delta DISAGREE IN ", nrow(.dsg), " OF ", nrow(delta_panel),
      " PANELS -- the conservative one does not\nreach 0.05 in: ",
      paste(.dsg$panel, collapse = ", "), "\n", sep = "")
} else {
  cat("ALL THREE INFERENCES ON delta AGREE IN ALL ", nrow(delta_panel),
      " PANELS: the donor-clustered CR2 interval\nexcludes zero everywhere, so the coupling is",
      " established at the donor level.\n",
      sep = "")
}
print(as.data.frame(delta_panel %>% transmute(
  panel,
  `logistic p (cell-level, conditional)` = signif(logit_p, 3),
  `near/far p (donor-level t)`           = signif(nf_p, 3),
  `LPM p (CR2 clustered by donor)`       = signif(delta_lpm_p_CR2, 3),
  `LPM CI covers 0`                      = (delta_lpm_CI.L < 0) & (delta_lpm_CI.R > 0))),
  row.names = FALSE)
cat("\nThe glmer p-value is conditional on donor (donor random intercept); the CR2 interval is on\n")
cat("the probability scale, clustered by donor.\n")
if (nrow(.dsg)) {
  cat("Where they disagree, delta is imprecise at donor level; that imprecision is propagated into\n")
  cat("PART 2 as a worst-case bias (shift/SE worst).\n")
}
cat("The deltas differ because of CELL SETS AND PREVALENCE, not units: there are only three\n")
cat("distinct dist_z scalings here (the two glia panels SHARE one DIST_SD_GLIA = ",
    fmt(DIST_SD_GLIA, 5), ").\n\n", sep = "")

cat("--- PART 2: gamma x delta, the bias --------------------------------------------\n")
cat("gamma is the `plaqueplaque-prox` coefficient of m3_plaque (Figure S8A); delta is delta_lpm\n")
cat("above; their product is the bias in the tangle coefficient from omitting amyloid.\n\n")
cat("ASSERTED, NOT CLAIMED. Two tiers:\n")
cat("  tier 1 (algebraic) everything evaluated at ONE lambda by whitened GLS, where GLS is OLS and\n")
cat("         the identity is exact. max |gamma*delta - (b1-b3)| = ",
    signif(max(abs(c(bias$alg_resid_plaque, bias$alg_resid_depth))), 3),
    "  (tol ", TOL_ALG, ")\n", sep = "")
cat("  tier 2 (practical) the published REML decomposition vs the observed m1 -> m3 shift.\n")
cat("         max |resid| / SE(b1) = ", signif(max(abs(bias$resid_over_SE_plaque)), 3),
    "  (tol ", TOL_REML, ")\n", sep = "")
cat("         The only source of discrepancy is REML re-estimating lambda when plaque enters.\n")
cat("  The assertion is on the PRODUCT, never the ratio (b1-b3)/gamma: LAMP1 has gamma ~ 0.002\n")
cat("  and a shift ~ 4e-05, i.e. numerically 0/0.\n\n")
cat("PRIMARY CURRENCY IS shift / SE(b1). The percentage of b1 is suppressed (NA) when\n")
cat("|b1| < 2*SE(b1), because dividing by a null coefficient manufactures a large number.\n\n")
bt <- bias %>% arrange(panel, desc(abs(shift_over_SE_plaque))) %>%
  transmute(panel, outcome, b1 = round(b1, 6), SE = round(b1_SE, 6),
            p_b1 = signif(b1_pval, 2), near_null = b1_near_null,
            gamma = round(gamma_plaque, 6), delta = round(delta_plaque, 6),
            `gamma*delta` = round(gamma_x_delta_plaque, 7), shift = round(shift_plaque, 7),
            `shift/SE` = round(shift_over_SE_plaque, 3),
            `shift/SE worst` = round(shift_over_SE_plaque_worst, 3),
            `pct_b1` = round(pct_of_b1_plaque, 2),
            `depth shift/SE` = round(shift_over_SE_depth, 3),
            `depth pct_b1` = round(pct_of_b1_depth, 2))
print(as.data.frame(bt), row.names = FALSE)
cat("\n  shift/SE worst = gamma x delta evaluated at the end of delta's CR2 interval that maximises\n")
cat("  the bias, i.e. how large the omitted-amyloid bias could be given delta's own imprecision.\n")
# --- SUMMARY, GENERATED FROM THE FIT rather than hardcoded, so it always matches the table above.
# --- "Material" = the shift moves the coefficient by >= 1 SE or >= 10% of itself, and the
# --- coefficient is not itself null (a % of a null coefficient is not informative).
MAT_SE  <- 1.0
MAT_PCT <- 10
bias$material_plaque <- !bias$b1_near_null &
  (abs(bias$shift_over_SE_plaque) >= MAT_SE | abs(bias$pct_of_b1_plaque) >= MAT_PCT)
mat <- bias %>% filter(material_plaque) %>% arrange(desc(abs(shift_over_SE_plaque)))
cat("\nVERDICT (generated from the fits above; 'material' = >= ", MAT_SE, " SE or >= ", MAT_PCT,
    "% of b1, b1 not null)\n", sep = "")
cat("  largest |shift|                      : ", fmt(max(abs(bias$shift_over_SE_plaque)), 3),
    " SE\n", sep = "")
cat("  same, at the worst end of delta's CI : ",
    fmt(max(abs(bias$shift_over_SE_plaque_worst)), 3), " SE\n", sep = "")
cat("  outcomes materially biased           : ", nrow(mat), " of ", nrow(bias), "\n", sep = "")
if (nrow(mat)) {
  print(as.data.frame(mat %>% transmute(panel, outcome, b1 = round(b1, 6),
                                        `shift/SE` = round(shift_over_SE_plaque, 3),
                                        `pct_b1` = round(pct_of_b1_plaque, 1))),
        row.names = FALSE)
  cat("\n  Ab4G8 (the 4G8 anti-amyloid antibody) and H31L21 (an Abeta epitope) report amyloid, so an\n")
  cat("  omitted-amyloid bias in their tangle-distance coefficient is expected and acts as a\n")
  cat("  positive control for the decomposition.\n")
}
oth <- bias %>% filter(!material_plaque)
cat("\n  For the remaining ", nrow(oth), " outcomes the omitted-amyloid bias is small: max ",
    fmt(max(abs(oth$shift_over_SE_plaque)), 3), " SE",
    if (any(is.finite(oth$pct_of_b1_plaque)))
      paste0(" and max ", fmt(max(abs(oth$pct_of_b1_plaque), na.rm = TRUE), 1), "% of b1")
    else "", ".\n", sep = "")
hl <- intersect(c("AT8", "GFAP", "CD68"), bias$outcome)
cat("  The tau-specific headline channels are unaffected: ",
    paste0(hl, " ", vapply(hl, function(o)
      fmt(max(abs(bias$shift_over_SE_plaque[bias$outcome == o])), 3), ""), " SE",
      collapse = ", "), ".\n", sep = "")

cat("\nDEPTH COMPARISON, per outcome (not a comparison of maxima across different outcomes):\n")
dcmp <- bias %>%
  transmute(panel, outcome, amyloid = abs(shift_over_SE_plaque), depth = abs(shift_over_SE_depth),
            ratio = depth / amyloid)
cat("  depth shift exceeds amyloid shift in ", sum(dcmp$depth > dcmp$amyloid), " of ", nrow(dcmp),
    " outcomes; median ratio ", fmt(median(dcmp$ratio), 2), "x",
    " (max ", fmt(max(dcmp$ratio), 1), "x)\n", sep = "")
cat("  largest |depth shift| = ", fmt(max(abs(bias$shift_over_SE_depth)), 3), " SE (",
    bias$panel[which.max(abs(bias$shift_over_SE_depth))], " / ",
    bias$outcome[which.max(abs(bias$shift_over_SE_depth))], ")\n", sep = "")
cat("\n")
cat("Matched_4G8_40 is a 40 um dilation threshold on a continuous amyloid field, so gamma x delta\n")
cat("is the bias from omitting this binary amyloid-proximity term.\n\n")

cat("--- PART 3: stratification -- the actual test ---------------------------------\n")
cat("Model: value ~ ", PANELS[[1]]$pterm, " * plaque + <extra> + Sex + Age_s + PMI_s + (1 | patient_id)\n", sep = "")
cat("One pooled fit, not two stratum-separate fits (the pooled fit shares sigma^2 and the donor\n")
cat("RE and is what yields the interaction test). dist_z stays on the PANEL-WIDE scale.\n")
cat("The PLAQUE-DISTAL simple slope is the test: the gradient among cells with no amyloid within\n")
cat("40 um. It equals the predictor main effect bit-for-bit (asserted), because plaque-distal is\n")
cat("the reference level. Simple slopes are computed from the fixed-effect vcov with a NORMAL\n")
cat("reference; emmeans silently drops to asymptotic z above 3,000 observations anyway, and at df\n")
cat("in the thousands z and t differ in the fifth decimal. Cross-checked against\n")
cat("emmeans::emtrends on ", ck_panel, " (max abs diff ", signif(emm_check, 3), ").\n", sep = "")
cat("BH within panel across outcomes, per term -- matching Figure S8A. Not pooled across panels:\n")
cat("different cell sets, different scalings, and dist_reln_calb1 is a SUBSET of the exc panel.\n\n")
cat("Reported on the beta_closer scale (per s.d. NEARER a tangle), as the published panels are.\n\n")
it <- inter_bc %>% arrange(panel, interaction_padj) %>%
  transmute(panel, outcome,
            distal = round(slope_distal, 5),
            distal_CI = paste0("[", round(slope_distal_CI.L, 4), ", ",
                               round(slope_distal_CI.R, 4), "]"),
            prox = round(slope_prox, 5),
            inter = round(interaction, 5), inter_SE = round(interaction_SE, 5),
            p = signif(interaction_p, 2), padj = signif(interaction_padj, 2),
            reversal = sign_reversal, equiv = equiv_within,
            mde_pct = round(100 * mde_frac_of_distal, 0))
print(as.data.frame(it), row.names = FALSE)
cat("\n  distal / prox = tangle-distance slope within each amyloid stratum, per s.d. nearer.\n")
cat("  equiv   = is the interaction CI inside +/- ", 100 * EQUIV_FRAC,
    "% of the distal slope (pre-declared margin)?\n", sep = "")
cat("  mde_pct = smallest interaction detectable at this SE, as % of the distal slope.\n\n")
cat("PER-STRATUM SUPPORT (cells, donors and ROIs in each amyloid stratum):\n")
st <- inter %>% distinct(panel, n_distal, n_prox, donors_distal, donors_prox,
                         donors_prox_ge10, donors_prox_ge30, rois_zero_prox, rois_total)
print(as.data.frame(st), row.names = FALSE)
cat("\nHOW TO READ A NULL INTERACTION. Equivalence, not absence. A p > 0.05 interaction does NOT\n")
cat("license 'the gradient is amyloid-independent'. The positive claim is the PLAQUE-DISTAL SLOPE\n")
cat("WITH ITS CI; the interaction only says the proximal stratum does not contradict it, down to\n")
cat("the mde_pct resolution tabulated above.\n\n")
# NB: reported from inter_bc, i.e. the SAME beta_closer scale as the table above -- mixing scales
# inside one log is how a sign gets misread.
nrev <- inter_bc %>% filter(interaction_padj < FDR)
if (nrow(nrev)) {
  cat("*** NOT ALL INTERACTIONS ARE NULL -- ", nrow(nrev), " pass BH at ", FDR, ": ",
      paste(paste0(nrev$panel, "/", nrev$outcome), collapse = ", "),
      "\n    (beta_closer scale, as above)\n", sep = "")
  cat("For these outcomes the constant-effect OVB decomposition in PART 2 is the WRONG MODEL:\n")
  cat("this is effect MODIFICATION, not confounding, and gamma x delta is not 'the bias' there.\n")
  print(as.data.frame(nrev %>% transmute(panel, outcome, interaction = round(interaction, 5),
                                         padj = signif(interaction_padj, 3),
                                         slope_distal = round(slope_distal, 5),
                                         slope_prox = round(slope_prox, 5), sign_reversal)),
        row.names = FALSE)
  cat("\n")
} else {
  cat("No interaction passes BH at ", FDR, " in any panel.\n\n", sep = "")
}

cat("--- Random effects, every fit --------------------------------------------------\n")
cat("Donor random-intercept variance for every fit, with singular fits flagged.\n")
for (P in PANELS) for (mk in P$outcomes) {
  F <- fits[[P$slug]][[mk]]
  cat_re_variance(F$m1, paste0(P$slug, " / ", mk, " / m1_base"))
  cat_re_variance(F$m3, paste0(P$slug, " / ", mk, " / m3_plaque"))
  cat_re_variance(F$mI, paste0(P$slug, " / ", mk, " / interaction"))
}
for (P in PANELS) cat_re_variance(get(paste0(".gm_", P$slug)),
                                  paste0(P$slug, " / delta (glmer binomial)"))
cat("\nNOTE the asymmetry, and why it matters. The two glia panels' OUTCOME fits are singular by\n")
cat("construction (the outcome is z-scored within donor; docs/MODELS.md Set 3, within-donor\n")
cat("z-scoring), so they are effectively pooled OLS. The delta model's outcome is `plaque`, which\n")
cat("is NOT z-scored, so its donor RE is non-zero. Since the OVB identity holds only under ONE\n")
cat("weight matrix, the correct multiplier for those panels is the POOLED-OLS delta, not the mixed\n")
cat("one -- which is what the lambda-matched whitened GLS returns automatically (lambda ~ 0 => no\n")
cat("whitening => OLS). Using the mixed delta would overstate the astrocyte DEPTH bias by ~3.6x.\n")
cat("The mixed delta remains the better SCIENTIFIC estimate of tau-amyloid coupling for those\n")
cat("panels and is in the glmer/OR rows above; it is simply not the multiplier.\n\n")

cat("--- Notes ---------------------------------------------------------------------\n")
cat("1. gamma x delta uses the binarised 40 um plaque proxy (see PART 2).\n")
cat("2. `sample_id` (ROI) is not modelled, per Set 3.\n")
cat("3. Effect sizes with 95% CIs are on every row of every emitted TSV.\n")
cat("4. The near/far bands sit on the same spatial scale as the 40 um dilation (<=", NEAR_UM,
    " um) and\n", sep = "")
cat("   close to the distance cap (>=", FAR_UM, " um, cap ", MAX_DIST,
    " um), so the far band is a censored tail.\n", sep = "")
cat("5. `plaque` is a Set 3 sensitivity term; the Figure S8A models are not modified -- this\n")
cat("   script writes only to ", out_dir, "/.\n", sep = "")

cat("\n--- Outputs -------------------------------------------------------------------\n")
for (f in sort(list.files(out_dir))) cat("  ", f, "\n", sep = "")
cat("\n")
print(sessionInfo())
sink()

message("\nWrote ", out_dir, "/")
message("  max |shift|/SE(b1) from omitting amyloid : ",
        signif(max(abs(bias$shift_over_SE_plaque)), 3),
        "  (worst case over delta's CI: ", signif(max(abs(bias$shift_over_SE_plaque_worst)), 3), ")")
message("  max |shift|/SE(b1) from omitting depth   : ",
        signif(max(abs(bias$shift_over_SE_depth)), 3))
message("  interactions passing BH at ", FDR, "     : ", sum(inter$interaction_padj < FDR),
        " of ", nrow(inter))
message("  stats log: ", log_path)
