#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# imc_covariate_sensitivity.R
#
# Figure panels: S8A, S8B, S8C
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# imc_covariate_sensitivity.R
#
# SUPPLEMENTARY covariate-sensitivity analysis for the IMC (EC project) panels. Asks one question:
# what happens to the reported effect when the two covariates that docs/MODELS.md classifies as
# "sensitivity terms, never in a headline" are added to the model?
#
#   plaque   = Matched_4G8_40, a logical in colData(spe): the nucleus lies inside the 40 um-DILATED
#              4G8+ amyloid mask. Built in IMCDataAnalysis-main/14-PlaqueNiche.Rmd:19-59 from
#              <EC_IMC_Project>/4G8_40_Match/*.csv. ~43.3k TRUE of 215,739 cells. Coded as in the
#              other IMC scripts, so the coefficient row is "plaqueplaque-prox".
#   depth_s  = relative cortical depth, i.e. distance from the PIAL SURFACE as a proportion. NOT
#              stored in the object; recomputed with the same per-ROI construction used in
#              imc_phf1_glia_distance.R
#                  scaled_Y = 100 * (y - min(y)) / (max(y) - min(y))     per sample_id (ROI)
#              0 = pia (the stored `layer` factor cuts scaled_Y at 0/10/55/58/100 with L1 = 0-10).
#              Here carried as depth_prop = scaled_Y/100 and modelled as
#              depth_s = scale(depth_prop) within each cell set.
#
# Model m1 below is the canonical headline model of each panel and is checked against the
# published estimates before the adjusted fits are used (see the "REGRESSION CHECK" section).
#
# Seven panels. Cohort filter, anchors, distance construction, cap, marker set, sign convention and
# figure style are copied from the ORIGIN script of each panel and must not drift from it:
#
#   slug                       origin                                predictor    extra
#   dist_reln_calb1            imc_phf1_distance_markers_subsets.R   dist_z       -
#   dist_exc_subcluster_adj    imc_phf1_distance_markers_subsets.R   dist_z       subcluster
#   state_ref_all_exc          imc_phf1_state_markers_subsets.R      grp          -
#   state_exc_subcluster_adj   imc_phf1_state_markers_subsets.R      grp          subcluster
#   state_reln_calb1           imc_phf1_state_markers_subsets.R      grp          -
#   cd68_microglia_distance    imc_phf1_glia_distance.R              dist_scaled  -
#   gfap_astrocytes_distance   imc_phf1_glia_distance.R              dist_scaled  -
#
# Model ladder, fitted per panel x outcome (4 fits; lmerTest::lmer, REML, random INTERCEPT only):
#   m1_base    value ~ P + extra +                     Sex + Age_s + PMI_s + (1 | patient_id)
#   m2_depth   value ~ P + extra + depth_s +           Sex + Age_s + PMI_s + (1 | patient_id)
#   m3_plaque  value ~ P + extra +           plaque +  Sex + Age_s + PMI_s + (1 | patient_id)
#   m4_full    value ~ P + extra + depth_s + plaque +  Sex + Age_s + PMI_s + (1 | patient_id)
#
# Outputs, all under plots/imc_covariate_sensitivity/ and all with the standard TRIPLE:
#   plot_<slug>_volcano_full.pdf / plot_<slug>_forest_full.pdf   the 5 marker panels redrawn from m4
#   plot_<slug>_full.pdf                                        the 2 glia curves, FDR from m4
#   plot_<slug>_ladder.pdf                                      m1..m4 side by side, per outcome
#   plot_beta3d_<slug>.pdf / plot_beta3d_combined.pdf            3-D: tau x plaque x depth
#   plot_beta_pairs_<slug>.pdf / plot_beta_pairs_combined.pdf    the readable 2-D companion
#   stats_covariate_ladder.tsv                                  master long table, every coefficient
#   stats_beta3d.tsv                                            the 3-D plot's exact coordinates
#
# Run locally, NOT on the HPC: spe.rds lives under <IMC_ROOT>. Optionally pass panel slugs as
# arguments to run a subset, e.g.
#   Rscript R/imc_covariate_sensitivity.R dist_reln_calb1

suppressPackageStartupMessages({
  library(SpatialExperiment); library(dplyr); library(tibble); library(tidyr)
  library(ggplot2); library(ggrepel); library(lmerTest); library(RANN)
  library(grid); library(patchwork); library(scatterplot3d)
})

# --- project root. The RDS share has been mounted at more than one path on this machine, so
# --- resolve against a candidate list rather than a single hard-coded fallback.
root_candidates <- c(
  "<PROJECT_ROOT>/phf1_v2",
  "<PROJECT_ROOT>/phf1_v2",
  "<PROJECT_ROOT>"
)
root <- root_candidates[dir.exists(root_candidates)][1]
if (is.na(root)) stop("Project root not found. Mount the RDS share.")
setwd(root)
source("R/palettes.R")
source("R/imc_utils.R")   # re_variance_summary(), cat_re_variance(), re_total_sd()

out_dir <- "plots/imc_covariate_sensitivity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Constants -- every one copied from the origin script that owns it
# -------------------------------------------------------------------
FDR          <- 0.05
LFC_MIN      <- 0.25       # state panels only (imc_phf1_state_markers_subsets.R)
PADJ_FLOOR   <- 1e-300     # state panels only
MAX_DIST     <- 300
EDGE_BUFFER  <- 0
MIN_ANCHORS  <- 1
WINDOW_UM    <- 75         # glia rolling mean
HALF_WINDOW  <- WINDOW_UM / 2
N_GRID       <- 200
MIN_WINDOW_N <- 10
RING_BREAKS  <- c(0, 25, 50, 100, 150, 200, 250, MAX_DIST)
SIG_BASIS    <- "padj<0.05"
SIG_ALPHA    <- 0.05

VULN_CLUSTERS <- c("Excitatory neuron cluster 4 (RELN)",
                   "Excitatory neuron cluster 2 (CALB1)")
MICRO <- c("Microglia cluster 1 (IBA1)", "Microglia cluster 2 (IBA1, CD68)")
ASTRO <- c("Astrocyte cluster 1 (GFAP, S100b)", "Astrocyte cluster 2 (S100b)")

# Ladder canvas. 9.2 cm so two ladders sit side by side within an 18.4 cm double column. The
# legend sits under the panel -- at this width a right-hand legend leaves only ~3.5 cm for the
# data. Heights carry the three legend rows this costs.
LADDER_W_CM      <- 9.2
LADDER_H_CM      <- 7.6    # 10 markers
LADDER_H_GLIA_CM <- 5.0    # 1 outcome

MODEL_LEVELS <- c("m1_base", "m2_depth", "m3_plaque", "m4_full")
# Legend text only (used by make_ladder). Kept short because the ladder legend sits UNDER a 9.2 cm
# panel, where "base (canonical headline)" alone would force the row to wrap.
MODEL_LABELS <- c(m1_base   = "base (headline)",
                  m2_depth  = "+ depth",
                  m3_plaque = "+ plaque",
                  m4_full   = "+ depth + plaque")
MODEL_EXTRA  <- list(m1_base   = character(0),
                     m2_depth  = "depth_s",
                     m3_plaque = "plaque",
                     m4_full   = c("depth_s", "plaque"))

# Okabe-Ito (colourblind-safe). One colour per
# COVARIATE AXIS, reused for the matching CI whisker in the 3-D view so a whisker's colour says
# which axis it belongs to.
COL_PROX   <- "#D55E00"   # vermillion  -- tau proximity (distance, or PHF1+ vs -)
COL_PLAQUE <- "#0072B2"   # blue        -- amyloid proximity
COL_DEPTH  <- "#009E73"   # bluish green-- cortical depth
MODEL_COLOURS <- c(m1_base = "#000000", m2_depth = COL_DEPTH,
                   m3_plaque = COL_PLAQUE, m4_full = "#CC79A7")
# 7 panels; Okabe-Ito minus the pale yellow, which is illegible as a small point on white
PANEL_COLOURS <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
                   "#D55E00", "#CC79A7", "#000000")

# Expected canonical values for the REGRESSION CHECK. The five marker panels are checked against
# their published source-data TSVs; the two glia panels have no per-marker TSV, so the numbers are
# taken from plots/imc_phf1_glia_distance/stats_*.txt, which quote 5 decimal places -- hence the
# looser tolerance for those two.
REF_TSV <- c(
  dist_reln_calb1          = "plots/imc_phf1_distance_markers_subsets/source_data_reln_calb1.tsv",
  dist_exc_subcluster_adj  = "plots/imc_phf1_distance_markers_subsets/source_data_exc_subcluster_adj.tsv",
  state_ref_all_exc        = "plots/imc_phf1_state_markers_subsets/source_data_ref_all_exc.tsv",
  state_exc_subcluster_adj = "plots/imc_phf1_state_markers_subsets/source_data_exc_subcluster_adj.tsv",
  state_reln_calb1         = "plots/imc_phf1_state_markers_subsets/source_data_reln_calb1.tsv"
)
REF_COL <- c(dist_reln_calb1 = "beta_closer", dist_exc_subcluster_adj = "beta_closer",
             state_ref_all_exc = "est", state_exc_subcluster_adj = "est",
             state_reln_calb1 = "est")
REF_GLIA <- tibble::tribble(
  ~panel,                     ~outcome, ~est,      ~n_cells, ~n_donors,
  "gfap_astrocytes_distance", "GFAP",    0.05446,   26832,    44,
  "cd68_microglia_distance",  "CD68",   -0.01743,   10970,    43
)
TOL_TSV  <- 1e-6
TOL_GLIA <- 1e-4

PANEL_SELECT <- commandArgs(trailingOnly = TRUE)
if (!length(PANEL_SELECT)) PANEL_SELECT <- NULL

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
# Cells. Two cohort frames, because the two origin scripts build theirs differently and m1 must
# reproduce both exactly.
#   base_all   -- all celltypes (the distance and glia scripts). Age_s/PMI_s scaled HERE.
#   base_state -- excitatory only (the state script). Age_s/PMI_s scaled on the SUBSET, as it does.
# Note Age_s/PMI_s are linear rescalings of Age/PMI, so which frame they are scaled on cannot move
# any other coefficient; the two frames are kept apart for exact reproduction, not for statistics.
# ===================================================================
cd <- as.data.frame(colData(spe))
cd$cell_id <- colnames(spe)
sc <- spatialCoords(spe); cd$x <- sc[, "Pos_X"]; cd$y <- sc[, "Pos_Y"]

min_rois_per_patient <- 3
cv_removal <- cd %>% distinct(patient_id, sample_id) %>% count(patient_id) %>%
  filter(n < min_rois_per_patient) %>% pull(patient_id)

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

# ROI geometry. scaled_Y / pia_dist_um / edge_dist are properties of the ROI, so they are computed
# over ALL cells in the ROI and then carried into every cell set -- computing them on a celltype
# subset would make min(y)/max(y) depend on which cells were selected.
base_all <- base_all %>%
  group_by(sample_id) %>%
  mutate(scaled_Y    = 100 * (y - min(y)) / (max(y) - min(y)),
         pia_dist_um = y - min(y),
         edge_dist   = pmin(x - min(x), max(x) - x, y - min(y), max(y) - y)) %>%
  ungroup() %>%
  mutate(depth_prop = scaled_Y / 100,     # 0 = pia, 1 = deep border
         # amyloid proximity, defined as in the other IMC scripts:
         # nucleus within the 40 um-dilated 4G8+ plaque mask
         plaque = factor(ifelse(Matched_4G8_40, "plaque-prox", "plaque-distal"),
                         levels = c("plaque-distal", "plaque-prox")),
         Age_s = as.numeric(scale(Age)), PMI_s = as.numeric(scale(PMI)))
# the bounds are exact by construction; the tolerance is for floating point on the max-y cell of an
# ROI, where 100*(y-min)/(max-min)/100 can land a few ulp above 1
stopifnot(all(is.finite(base_all$depth_prop)),
          min(base_all$depth_prop) >= -1e-9, max(base_all$depth_prop) <= 1 + 1e-9)

GEOM <- base_all %>% select(cell_id, scaled_Y, depth_prop, pia_dist_um, edge_dist, plaque)

base_state <- COHORT(cd) %>%
  filter(grepl("^Excitatory", celltype_clusters)) %>%
  mutate(grp = factor(ifelse(PHF1_Otsu == "PHF1_pos", "PHF1+", "PHF1-"),
                      levels = c("PHF1-", "PHF1+")),
         Age_s = as.numeric(scale(Age)), PMI_s = as.numeric(scale(PMI))) %>%
  left_join(GEOM, by = "cell_id")
stopifnot(!anyNA(base_state$depth_prop), !anyNA(base_state$plaque))
base_state$subcluster <- relevel(factor(base_state$celltype_clusters),
                                 ref = names(sort(table(base_state$celltype_clusters),
                                                  decreasing = TRUE))[1])
stopifnot(all(VULN_CLUSTERS %in% levels(base_state$subcluster)))

# --- distance to the nearest PHF1+ NEURON, within ROI (as in imc_phf1_distance_markers_subsets.R)
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
# --- all capped cells rather than per glial subset (as in imc_phf1_glia_distance.R), which is why
# --- the two glia panels share one dist_sd.
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

# depth_s: standardised within each modelled cell set, exactly as dist_z is. A linear rescaling, so
# it moves no other coefficient; the per-set SD is logged so the "per s.d." units are unambiguous.
add_depth_s <- function(md) { md$depth_s <- as.numeric(scale(md$depth_prop)); md }
md_dist_vuln  <- add_depth_s(md_dist_vuln)
md_dist_exc   <- add_depth_s(md_dist_exc)
md_glia_micro <- add_depth_s(md_glia_micro)
md_glia_astro <- add_depth_s(md_glia_astro)
md_state_all  <- add_depth_s(base_state)
md_state_vuln <- add_depth_s(base_state %>% filter(celltype_clusters %in% VULN_CLUSTERS))

# ===================================================================
# The model ladder
# ===================================================================
grab <- function(m, term) {
  if (is.null(m)) return(rep(NA_real_, 5))
  co <- summary(m)$coefficients
  if (!term %in% rownames(co)) return(rep(NA_real_, 5))
  c(co[term, "Estimate"], co[term, "Std. Error"], co[term, "df"],
    co[term, "t value"], co[term, "Pr(>|t|)"])
}

# The three terms of interest and the coefficient row each one occupies. `predictor` is filled in
# per panel because it is dist_z / grpPHF1+ / dist_scaled depending on the origin script.
TERM_ROWS <- c(plaque = "plaqueplaque-prox", depth = "depth_s")

#' Fit m1..m4 for one panel and return the coefficients in long form.
#'
#' @param md      modelled cell set; must carry the predictor, depth_s, plaque, Sex, Age_s, PMI_s,
#'                patient_id, and `extra` if used.
#' @param outcomes character vector of outcome names (10 state markers, or one glial marker).
#' @param prep    function(md, outcome) -> data frame with a finite `value` column.
#' @param pterm   the predictor as it appears in the FORMULA ("dist_z" / "grp" / "dist_scaled").
#' @param prow    the predictor's coefficient ROW name ("dist_z" / "grpPHF1+" / "dist_scaled").
#' @param extra   extra fixed effects present in every model of the ladder (e.g. "subcluster").
#' @param ref_fun function(d) -> the arcsinh reference mean for the fold-change back-transform,
#'                or NULL when the outcome is not on an arcsinh scale (the glia z-scores).
fit_ladder <- function(md, outcomes, prep, pterm, prow, extra = character(0), ref_fun = NULL) {
  forms <- lapply(MODEL_LEVELS, function(mn)
    as.formula(paste("value ~", paste(c(pterm, extra, MODEL_EXTRA[[mn]],
                                        "Sex", "Age_s", "PMI_s", "(1 | patient_id)"),
                                      collapse = " + "))))
  names(forms) <- MODEL_LEVELS

  bind_rows(lapply(outcomes, function(mk) {
    d <- prep(md, mk)
    refmean <- if (is.null(ref_fun)) NA_real_ else ref_fun(d)
    bind_rows(lapply(MODEL_LEVELS, function(mn) {
      m <- tryCatch(lmerTest::lmer(forms[[mn]], data = d, REML = TRUE), error = function(e) NULL)
      sd_tot <- re_total_sd(m)
      rv <- re_variance_summary(m); rv <- rv[!is.na(rv$group) & rv$group == "patient_id", ,
                                            drop = FALSE]
      want <- c(predictor = prow, TERM_ROWS)
      # a term absent from this model of the ladder simply yields no row
      keep <- names(want)[names(want) == "predictor" |
                          want %in% c(if ("depth_s" %in% MODEL_EXTRA[[mn]]) "depth_s",
                                      if ("plaque"  %in% MODEL_EXTRA[[mn]]) "plaqueplaque-prox")]
      bind_rows(lapply(keep, function(tn) {
        a <- grab(m, want[[tn]])
        tibble(outcome = mk, model = mn, term = tn, coef_row = want[[tn]],
               est = a[1], SE = a[2], df = a[3], t = a[4], pval = a[5],
               sd_tot = sd_tot, refmean = refmean,
               n_cells = nrow(d), n_donors = dplyr::n_distinct(d$patient_id),
               n_rois = dplyr::n_distinct(d$sample_id),
               re_var_donor = if (nrow(rv)) rv$vcov[1] else NA_real_,
               re_pct_donor = if (nrow(rv)) rv$pct_of_total[1] else NA_real_,
               singular = if (is.null(m)) NA else lme4::isSingular(m))
      }))
    }))
  }))
}

# prep functions -----------------------------------------------------
prep_marker <- function(md, mk) {
  d <- md
  d$value <- as.numeric(EM[mk, md$cell_id])
  d[is.finite(d$value), , drop = FALSE]
}
# glia: within-donor z of the arcsinh intensity, computed on the analysed population
# (as in imc_phf1_glia_distance.R). This is why the donor random-intercept variance is zero.
prep_glia <- function(md, mk) {
  d <- md
  d$value_raw <- d[[mk]]
  d <- d[is.finite(d$value_raw), , drop = FALSE]
  d <- d %>% group_by(patient_id) %>% mutate(value = as.numeric(scale(value_raw))) %>% ungroup()
  as.data.frame(d[is.finite(d$value), , drop = FALSE])
}

# ===================================================================
# Panel table
# ===================================================================
# `flip` orients the predictor toward MORE tau: the marker panels already report per s.d. NEARER a
# tangle (beta_closer = -beta(dist_z)), the state panels' grp contrast is already PHF1+ minus PHF1-,
# and the glia panels natively report per s.d. FURTHER, so only they are flipped. Every panel's
# `beta_prox` therefore points the same way, so all panels share one 3-D axis.
PANELS <- list(
  list(slug = "dist_reln_calb1", kind = "marker", style = "dist", md = md_dist_vuln,
       outcomes = markers, prep = prep_marker, pterm = "dist_z", prow = "dist_z",
       extra = character(0), flip = -1,
       ref_fun = function(d) mean(d$value),
       origin = "R/imc_phf1_distance_markers_subsets.R",
       title = "Covariate sensitivity: distance gradient in RELN+ and CALB1+ excitatory neurons (IMC)"),
  list(slug = "dist_exc_subcluster_adj", kind = "marker", style = "dist", md = md_dist_exc,
       outcomes = markers, prep = prep_marker, pterm = "dist_z", prow = "dist_z",
       extra = "subcluster", flip = -1,
       ref_fun = function(d) mean(d$value),
       origin = "R/imc_phf1_distance_markers_subsets.R",
       title = "Covariate sensitivity: distance gradient in excitatory neurons, subcluster-adjusted (IMC)"),
  list(slug = "state_ref_all_exc", kind = "marker", style = "state", md = md_state_all,
       outcomes = markers, prep = prep_marker, pterm = "grp", prow = "grpPHF1+",
       extra = character(0), flip = 1,
       ref_fun = function(d) mean(d$value[d$grp == "PHF1-"]),
       origin = "R/imc_phf1_state_markers_subsets.R",
       title = "Covariate sensitivity: PHF1+ vs PHF1- state markers, all excitatory neurons (IMC)"),
  list(slug = "state_exc_subcluster_adj", kind = "marker", style = "state", md = md_state_all,
       outcomes = markers, prep = prep_marker, pterm = "grp", prow = "grpPHF1+",
       extra = "subcluster", flip = 1,
       ref_fun = function(d) mean(d$value[d$grp == "PHF1-"]),
       origin = "R/imc_phf1_state_markers_subsets.R",
       title = "Covariate sensitivity: PHF1+ vs PHF1- state markers, subcluster-adjusted (IMC)"),
  list(slug = "state_reln_calb1", kind = "marker", style = "state", md = md_state_vuln,
       outcomes = markers, prep = prep_marker, pterm = "grp", prow = "grpPHF1+",
       extra = character(0), flip = 1,
       ref_fun = function(d) mean(d$value[d$grp == "PHF1-"]),
       origin = "R/imc_phf1_state_markers_subsets.R",
       title = "Covariate sensitivity: PHF1+ vs PHF1- state markers, RELN+ and CALB1+ neurons (IMC)"),
  list(slug = "cd68_microglia_distance", kind = "glia", style = "glia", md = md_glia_micro,
       outcomes = "CD68", prep = prep_glia, pterm = "dist_scaled", prow = "dist_scaled",
       extra = character(0), flip = -1, ref_fun = NULL,
       series = "Microglia", ylab = "CD68 (within-donor z)",
       origin = "R/imc_phf1_glia_distance.R",
       title = "Covariate sensitivity: CD68 in pooled microglia vs distance to nearest PHF1+ neuron (IMC)"),
  list(slug = "gfap_astrocytes_distance", kind = "glia", style = "glia", md = md_glia_astro,
       outcomes = "GFAP", prep = prep_glia, pterm = "dist_scaled", prow = "dist_scaled",
       extra = character(0), flip = -1, ref_fun = NULL,
       series = "Astrocytes", ylab = "GFAP (within-donor z)",
       origin = "R/imc_phf1_glia_distance.R",
       title = "Covariate sensitivity: GFAP in pooled astrocytes vs distance to nearest PHF1+ neuron (IMC)")
)
names(PANELS) <- vapply(PANELS, `[[`, "", "slug")
if (!is.null(PANEL_SELECT)) {
  unknown <- setdiff(PANEL_SELECT, names(PANELS))
  if (length(unknown)) stop("unknown panel slug(s): ", paste(unknown, collapse = ", "))
  PANELS <- PANELS[PANEL_SELECT]
}

# ===================================================================
# Fit everything
# ===================================================================
# BH is applied WITHIN a panel, across the outcomes, SEPARATELY for each (model, term) -- so the
# base model's predictor p-values are adjusted among themselves, the full model's among themselves,
# and the plaque and depth coefficients each get their own family. The glia panels have one
# outcome, so BH is a no-op there and is kept explicit so the basis is unambiguous.
run_panel <- function(P) {
  message("== ", P$slug, "  (", nrow(P$md), " cells, ", length(P$outcomes), " outcome(s))")
  raw <- fit_ladder(P$md, P$outcomes, P$prep, P$pterm, P$prow, P$extra, P$ref_fun)

  res <- raw %>%
    mutate(panel = P$slug,
           # orientation: predictor toward MORE tau; plaque and depth are never flipped
           orient_sign = ifelse(term == "predictor", P$flip, 1),
           est_or = est * orient_sign,
           crit = stats::qt(0.975, ifelse(is.na(df), Inf, df)),
           CI.L = est_or - crit * SE,
           CI.R = est_or + crit * SE,
           # aliases so make_ladder() can address either currency with one column-stem argument
           est_or.L = CI.L, est_or.R = CI.R,
           cohens_d = est_or / sd_tot,
           d.L = CI.L / sd_tot, d.R = CI.R / sd_tot,
           # exact inverse of asinh, so a ratio on the intensity scale is a genuine fold change.
           # NA for the glia panels: their outcome is a within-donor z, not an arcsinh intensity.
           logFC   = log2(sinh(refmean + est_or) / sinh(refmean)),
           logFC.L = log2(sinh(refmean + CI.L)   / sinh(refmean)),
           logFC.R = log2(sinh(refmean + CI.R)   / sinh(refmean))) %>%
    group_by(model, term) %>%
    mutate(padj = if (P$style == "state") pmax(p.adjust(pval, method = "BH"), PADJ_FLOOR)
                  else p.adjust(pval, method = "BH")) %>%
    ungroup() %>%
    mutate(model = factor(model, levels = MODEL_LEVELS),
           term_label = recode(term, predictor = "tau proximity",
                               plaque = "plaque-prox vs distal", depth = "deeper (away from pia)")) %>%
    arrange(term, model, pval)

  # the per-cell frame is needed again for the glia rolling curve
  d_glia <- if (P$kind == "glia") P$prep(P$md, P$outcomes[1]) else NULL
  list(panel = P, res = res, d = d_glia)
}

RUNS <- lapply(PANELS, run_panel)
ALL  <- bind_rows(lapply(RUNS, `[[`, "res"))

# ===================================================================
# REGRESSION CHECK -- m1 must reproduce the published headline
# ===================================================================
check_rows <- bind_rows(lapply(RUNS, function(R) {
  P <- R$panel
  mine <- R$res %>% filter(term == "predictor", model == "m1_base") %>%
    select(panel, outcome, mine = est_or, n_cells, n_donors)
  if (P$kind == "glia") {
    ref <- REF_GLIA %>% filter(panel == P$slug) %>%
      select(outcome, ref = est, ref_n_cells = n_cells, ref_n_donors = n_donors)
    # the published glia estimate is per s.d. FURTHER; est_or is per s.d. NEARER, so un-flip
    mine %>% left_join(ref, by = "outcome") %>%
      mutate(mine = mine * P$flip, source = "stats_*.txt (5 d.p.)", tol = TOL_GLIA)
  } else {
    f <- REF_TSV[[P$slug]]
    if (is.null(f) || !file.exists(f)) return(tibble())
    rr <- read.delim(f, check.names = FALSE)
    ref <- tibble(outcome = as.character(rr$gene), ref = as.numeric(rr[[REF_COL[[P$slug]]]]))
    mine %>% left_join(ref, by = "outcome") %>%
      mutate(ref_n_cells = NA_integer_, ref_n_donors = NA_integer_,
             source = basename(f), tol = TOL_TSV)
  }
}))
if (nrow(check_rows)) {
  check_rows <- check_rows %>%
    mutate(delta = mine - ref, ok = is.finite(delta) & abs(delta) <= tol)
} else {
  check_rows <- tibble(panel = character(0), outcome = character(0), mine = numeric(0),
                       ref = numeric(0), delta = numeric(0), tol = numeric(0), ok = logical(0),
                       n_cells = integer(0), n_donors = integer(0),
                       ref_n_cells = integer(0), ref_n_donors = integer(0), source = character(0))
}

if (nrow(check_rows) && !all(check_rows$ok, na.rm = TRUE)) {
  bad <- check_rows %>% filter(!ok)
  warning("REGRESSION CHECK FAILED for ", nrow(bad), " outcome(s): m1 does not reproduce the ",
          "published headline. The cell set was rebuilt wrongly and every adjusted fit below is ",
          "suspect. See stats_covariate_sensitivity_overview.txt.", call. = FALSE)
  print(as.data.frame(bad %>% select(panel, outcome, mine, ref, delta, source)), row.names = FALSE)
}

# ===================================================================
# Figure helpers
# ===================================================================
wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE)

XLAB_DIST  <- expression(Log[2]*" FC per s.d. nearer tangle")
XLAB_STATE <- expression(Log[2]*" (fold-change), PHF1+ / PHF1-")

# --- forest / volcano. Copied from imc_phf1_distance_markers_subsets.R and
# --- imc_phf1_state_markers_subsets.R, parameterised on the column quadruple.
# SIG_LEVELS as a fixed factor level set, used by every legend that can lose a level. A guide key is
# built from the DATA, so if every marker is significant the "n.s." entry silently vanishes; making
# the variable a factor with both levels and passing drop = FALSE keeps the legend complete and
# comparable across panels.
SIG_LEVELS <- c("padj < 0.05", "n.s.")
sig_lab_factor <- function(ok) factor(ifelse(ok, SIG_LEVELS[1], SIG_LEVELS[2]), levels = SIG_LEVELS)
SIG_COLOURS <- setNames(c("#BD0026", "grey65"), SIG_LEVELS)

#' Carry one undrawn row per absent factor level, so its legend key gets a GLYPH.
#'
#' `drop = FALSE` plus `limits` is NOT enough on its own: it keeps the absent level's LABEL but the
#' key is built from LAYER DATA, so with nothing to draw from the glyph is simply missing and the
#' legend reads as a label with blank space beside it. Verified here on the shape guide when every
#' marker clears the effect-size floor. Same problem and same fix as `pad_levels()` below, which the
#' glia curve inherited from R/phf1_channel_intensity_exc.R -- generalised to any factor column so
#' the forest, ladder and pairs guides can all use it. The padded row plots nothing because its
#' positional values are NA, at the cost of ggplot's "Removed N rows containing missing values"
#' message, which is expected rather than a fault.
pad_factor_levels <- function(d, col, value_cols) {
  f <- d[[col]]
  if (!is.factor(f)) return(d)
  miss <- setdiff(levels(f), as.character(f))
  if (!length(miss)) return(d)
  filler <- d[rep(1L, length(miss)), , drop = FALSE]
  filler[[col]] <- factor(miss, levels = levels(f))
  for (v in intersect(value_cols, names(filler))) filler[[v]] <- NA_real_
  dplyr::bind_rows(as.data.frame(d), as.data.frame(filler))
}

make_forest <- function(d, xlab, lfc_min = 0, at8_label = FALSE) {
  fpd <- d %>%
    transmute(gene = outcome, x = logFC, lo = logFC.L, hi = logFC.R, padj_use = padj) %>%
    mutate(gene_lab = if (at8_label) ifelse(gene == "AT8", "AT8 (pos. control)", gene) else gene) %>%
    mutate(gene_lab = factor(gene_lab, levels = rev(gene_lab[order(x)])),
           sig = sig_lab_factor(!is.na(padj_use) & padj_use < FDR & abs(x) >= lfc_min)) %>%
    pad_factor_levels("sig", c("x", "lo", "hi"))
  ggplot(fpd, aes(x = x, y = gene_lab, colour = sig)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.4) +
    geom_point(size = 1.7) +
    scale_colour_manual(values = SIG_COLOURS, name = NULL, drop = FALSE, limits = SIG_LEVELS) +
    labs(x = xlab, y = NULL) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "top", plot.margin = margin(4, 8, 4, 4, unit = "pt"))
}

make_volcano <- function(dt, xlab, lfc_min = 0) {
  dt <- as.data.frame(dt)
  dt <- dt[!is.na(dt$padj) & !is.nan(dt$logFC), ]
  dt <- dt[order(dt$padj), ]
  dt$de <- "Not sig"
  dt$de[dt$padj <= FDR & dt$logFC >=  lfc_min & dt$logFC > 0] <- "Up"
  dt$de[dt$padj <= FDR & dt$logFC <= -lfc_min & dt$logFC <= 0] <- "Down"
  dt$de <- factor(dt$de, levels = c("Up", "Down", "Not sig"))
  dt$lab <- ifelse(dt$de != "Not sig", dt$outcome, "")
  if (any(dt$padj == 0)) dt$padj[dt$padj == 0] <- min(dt$padj[dt$padj > 0])
  max_x <- max(abs(dt$logFC), na.rm = TRUE) * 1.15
  max_y <- max(-log10(dt$padj), na.rm = TRUE) * 1.12
  p <- ggplot(dt, aes(x = logFC, y = -log10(padj)))
  if (lfc_min > 0) {
    p <- p +
      annotate("rect", xmin = -lfc_min, xmax = lfc_min, ymin = -Inf, ymax = Inf,
               fill = "grey55", alpha = 0.12) +
      geom_vline(xintercept = c(-lfc_min, lfc_min), linetype = 2, linewidth = 0.3, alpha = 0.5)
  } else {
    p <- p + geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, alpha = 0.5)
  }
  p +
    geom_hline(yintercept = -log10(FDR), linetype = 2, linewidth = 0.3, alpha = 0.5) +
    geom_point(aes(colour = de), size = 1.4, alpha = 0.9, show.legend = FALSE) +
    ggrepel::geom_text_repel(
      aes(label = lab, colour = de), size = 2, max.overlaps = Inf,
      min.segment.length = 0, segment.size = 0.2, segment.colour = "grey45",
      box.padding = 0.45, point.padding = 0.15, force = 4, force_pull = 0.5,
      max.iter = 20000, max.time = 2, seed = 42, na.rm = TRUE, show.legend = FALSE) +
    scale_colour_manual(values = c("Up" = "#DC0000FF", "Down" = "#3C5488FF",
                                   "Not sig" = "grey70"), guide = "none") +
    scale_x_continuous(limits = c(-max_x, max_x), oob = scales::squish) +
    scale_y_continuous(limits = c(0, max_y), oob = scales::squish) +
    labs(x = xlab, y = expression("-" * log[10] * " (adjusted p-value)")) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.margin = margin(6, 8, 6, 8, unit = "pt"))
}

# --- the ladder: m1..m4 side by side for every outcome. This is the panel that actually answers
# --- "what do the covariates do", because a sign flip between models is visible at a glance.
make_ladder <- function(d, xcol, xlab, lfc_min = 0) {
  fpd <- d %>%
    transmute(gene = outcome, model,
              x = .data[[xcol]], lo = .data[[paste0(xcol, ".L")]],
              hi = .data[[paste0(xcol, ".R")]], padj_use = padj) %>%
    mutate(gene_lab = factor(gene, levels = rev(sort(unique(gene)))),
           model_lab = factor(MODEL_LABELS[as.character(model)],
                              levels = unname(MODEL_LABELS)),
           # lfc_min > 0 adds the effect-size floor to the significance call, so a filled dot means
           # padj < FDR AND |logFC| >= lfc_min -- the same rule the state volcano and forest use.
           sig = sig_lab_factor(!is.na(padj_use) & padj_use < FDR & abs(x) >= lfc_min)) %>%
    pad_factor_levels("sig", c("x", "lo", "hi")) %>%
    pad_factor_levels("model_lab", c("x", "lo", "hi"))
  # group = model_lab EXPLICITLY on both layers, and one shared position_dodge object. Without it
  # the point layer's default group is interaction(model_lab, shape) -- 8 dodge slots -- while the
  # errorbar layer has 4, so ggplot dodges them to different offsets and a dot ends up sitting next
  # to a CI line of a different colour.
  pd <- position_dodge(width = 0.7)
  p <- ggplot(fpd, aes(x = x, y = gene_lab, colour = model_lab, group = model_lab))
  if (lfc_min > 0) {
    # shaded dead-band, drawn FIRST so it sits behind the estimates
    p <- p +
      annotate("rect", xmin = -lfc_min, xmax = lfc_min, ymin = -Inf, ymax = Inf,
               fill = "grey55", alpha = 0.12) +
      geom_vline(xintercept = c(-lfc_min, lfc_min), linetype = 2, linewidth = 0.3,
                 colour = "grey50")
  }
  p +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.35,
                  position = pd) +
    geom_point(aes(shape = sig), size = 1.3, position = pd) +
    scale_colour_manual(values = setNames(unname(MODEL_COLOURS[MODEL_LEVELS]),
                                          unname(MODEL_LABELS[MODEL_LEVELS])),
                        name = "Model", drop = FALSE,
                        limits = unname(MODEL_LABELS[MODEL_LEVELS])) +
    scale_shape_manual(values = setNames(c(16, 1), SIG_LEVELS), name = NULL,
                       drop = FALSE, limits = SIG_LEVELS) +
    # Legend BELOW the panel, stacked: at LADDER_W_CM the panel would otherwise be left ~3.5 cm by a
    # right-hand legend. Model on two rows of two, padj on one, so the whole block is three compact
    # rows -- side by side at the bottom the two guides do not fit in 9.2 cm.
    guides(colour = guide_legend(order = 1, nrow = 2, byrow = TRUE),
           shape  = guide_legend(order = 2, nrow = 1)) +
    labs(x = xlab, y = NULL) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "bottom", legend.box = "vertical",
          legend.key.height = grid::unit(7, "pt"),
          legend.key.width = grid::unit(7, "pt"),
          legend.box.spacing = grid::unit(2, "pt"),
          legend.spacing.y = grid::unit(1, "pt"),
          legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(4, 8, 4, 4, unit = "pt"))
}

# --- glia rolling-mean curve. roll_mean / blank_sparse / sig_factor / pad_levels / make_panel /
# --- save_panel are copied VERBATIM from imc_phf1_glia_distance.R so the redrawn panel is
# --- pixel-identical to the published one except for the FDR linetype, which here comes from m4.
roll_mean <- function(x, y, grid, half_window) {
  out <- lapply(grid, function(g) {
    idx <- which(x >= g - half_window & x <= g + half_window); n <- length(idx)
    if (n < 1) return(c(roll_mean = NA_real_, sem = NA_real_, n_window = 0))
    m <- mean(y[idx]); s <- if (n > 1) stats::sd(y[idx]) / sqrt(n) else NA_real_
    c(roll_mean = m, sem = s, n_window = n)
  })
  as.data.frame(do.call(rbind, out))
}
blank_sparse <- function(r) {
  r$roll_mean[r$n_window < MIN_WINDOW_N] <- NA
  r$sem[r$n_window < MIN_WINDOW_N]       <- NA
  r
}
sig_factor <- function(x) factor(ifelse(x, SIG_BASIS, "ns"), levels = c(SIG_BASIS, "ns"))
pad_levels <- function(d, col = "sig") {
  miss <- setdiff(c(SIG_BASIS, "ns"), as.character(d[[col]]))
  if (!length(miss)) return(d)
  filler <- d[rep(1L, length(miss)), , drop = FALSE]
  filler[[col]] <- factor(miss, levels = c(SIG_BASIS, "ns"))
  for (v in intersect(c("roll_mean", "sem"), names(filler))) filler[[v]] <- NA_real_
  dplyr::bind_rows(d, filler)
}
XLAB_UM <- expression("Distance to PHF1+ neuron (" * mu * "m)")
make_panel <- function(roll_df, ylab, legend_title, cols) {
  roll_df <- pad_levels(roll_df, "sig")
  p <- ggplot(roll_df, aes(dist_um, roll_mean, colour = series, fill = series,
                           group = interaction(series, sig))) +
    geom_ribbon(aes(ymin = roll_mean - 1.96 * sem, ymax = roll_mean + 1.96 * sem),
                alpha = 0.15, colour = NA) +
    geom_line(aes(linetype = sig, linewidth = sig)) +
    scale_colour_manual(values = cols, name = legend_title, drop = FALSE) +
    scale_fill_manual(values = cols, guide = "none", drop = FALSE) +
    scale_linetype_manual(values = setNames(c("solid", "dashed"), c(SIG_BASIS, "ns")),
                          name = "Model FDR", drop = FALSE, limits = c(SIG_BASIS, "ns")) +
    scale_linewidth_manual(values = setNames(c(0.8, 0.4), c(SIG_BASIS, "ns")),
                           name = "Model FDR", drop = FALSE, limits = c(SIG_BASIS, "ns")) +
    guides(colour = guide_legend(order = 1),
           linetype = guide_legend(order = 2), linewidth = guide_legend(order = 2)) +
    labs(x = XLAB_UM, y = ylab) +
    coord_cartesian(xlim = c(0, MAX_DIST)) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.text = element_text(size = 6),
          legend.key.width = grid::unit(18, "pt"), legend.key.height = grid::unit(8, "pt"),
          legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
  g <- ggplot2::ggplotGrob(p)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  g$widths[pcol] <- grid::unit(1.6, "in")
  g
}
save_panel <- function(g, file, total_w_in = NULL, total_h_in = 2.1) {
  if (!is.null(total_w_in)) {
    pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
    furniture <- sum(grid::convertWidth(g$widths[-pcol], "in", valueOnly = TRUE))
    g$widths[pcol] <- grid::unit(max(total_w_in - furniture, 0.6), "in")
  }
  ggsave(file.path(out_dir, file), g,
         width = if (is.null(total_w_in))
                   grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE) else total_w_in,
         height = total_h_in, units = "in", device = "pdf")
}
PANEL_TOTAL_W_IN <- (12.7 / 2) / 2.54     # 2.500 in, matches the published glia panels
PANEL_TOTAL_H_IN <- 4.5 / 2.54            # 1.772 in
MODULE_COLOURS <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                    "#A65628", "#F781BF", "#666666", "#66C2A5", "#FC8D62")

ring_table <- function(d, value_col, series_label) {
  d %>%
    mutate(dist_ring = cut(dist_um, breaks = RING_BREAKS, include.lowest = TRUE, right = TRUE),
           .value = .data[[value_col]]) %>%
    filter(!is.na(dist_ring)) %>%
    group_by(patient_id, dist_ring) %>%
    summarise(mean_value = mean(.value), sd_value = sd(.value), n_cells = n(), .groups = "drop") %>%
    mutate(series = series_label,
           ring_mid_um = (RING_BREAKS[as.integer(dist_ring)] +
                            RING_BREAKS[as.integer(dist_ring) + 1]) / 2) %>%
    select(series, patient_id, dist_ring, ring_mid_um, mean_value, sd_value, n_cells)
}

# ===================================================================
# The 3-D coefficient view
# ===================================================================
# One point per outcome. Coordinates are Cohen's d for the three coefficients of the SAME model
# (m4_full), d rather than the raw beta so the arcsinh-scale marker panels and the within-donor-z
# glia panels can share one cube. Whisker colour identifies its axis.
beta3d_table <- function(res) {
  res %>% filter(model == "m4_full") %>%
    select(panel, outcome, n_cells, n_donors, term,
           est_or, CI.L, CI.R, cohens_d, d.L, d.R, pval, padj) %>%
    pivot_wider(id_cols = c(panel, outcome, n_cells, n_donors),
                names_from = term,
                values_from = c(est_or, CI.L, CI.R, cohens_d, d.L, d.R, pval, padj)) %>%
    # label only the outcomes whose tau coefficient survives; a 52-point cube labelled in full is
    # unreadable, and the pairs companion carries every label anyway
    mutate(lab3d = ifelse(!is.na(padj_predictor) & padj_predictor < FDR, outcome, ""))
}

pad_lim <- function(...) {
  v <- c(...); v <- v[is.finite(v)]
  if (!length(v)) return(c(-1, 1))
  m <- max(abs(v)) * 1.12
  if (!is.finite(m) || m <= 0) m <- 1
  c(-m, m)
}

beta3d_plot <- function(b3, file, point_col, label_col = "outcome",
                        legend_lab = NULL, legend_cols = NULL,
                        xlab = "d: nearer tangle") {
  b3 <- as.data.frame(b3)
  xl <- pad_lim(b3$cohens_d_predictor, b3$d.L_predictor, b3$d.R_predictor)
  yl <- pad_lim(b3$cohens_d_plaque,    b3$d.L_plaque,    b3$d.R_plaque)
  zl <- pad_lim(b3$cohens_d_depth,     b3$d.L_depth,     b3$d.R_depth)

  # A cube needs room that a 3.5 in subfigure does not have: three axis titles, three sets of tick
  # labels (the y ones sit OUTSIDE the box on the right, so the right margin is not optional), point
  # labels, and two legends. 4.6 x 4.4 in with a reserved bottom strip for the legends is the
  # smallest canvas at which nothing is clipped. This is the one figure here that is not a
  # Nature-compact subfigure, by necessity.
  pdf(file.path(out_dir, file), width = 4.6, height = 4.8, useDingbats = FALSE)
  s3d <- scatterplot3d(
    x = b3$cohens_d_predictor, y = b3$cohens_d_plaque, z = b3$cohens_d_depth,
    xlim = xl, ylim = yl, zlim = zl,
    type = "p", pch = 16, cex.symbols = 0.7,
    color = point_col, angle = 52, scale.y = 0.9, grid = TRUE, box = TRUE,
    xlab = xlab, ylab = "d: plaque-prox", zlab = "d: deeper (from pia)",
    cex.axis = 0.5, cex.lab = 0.65, tick.marks = TRUE,
    mar = c(4.2, 3.0, 3.2, 2.6))

  seg3 <- function(x0, y0, z0, x1, y1, z1, ...) {
    p <- s3d$xyz.convert(c(x0, x1), c(y0, y1), c(z0, z1))
    segments(p$x[1], p$y[1], p$x[2], p$y[2], ...)
  }
  # zero reference lines, one per axis, drawn on the floor / back wall. Dotted grey: the null.
  seg3(0, yl[1], zl[1], 0, yl[2], zl[1], col = "grey45", lty = 3, lwd = 0.6)
  seg3(xl[1], 0, zl[1], xl[2], 0, zl[1], col = "grey45", lty = 3, lwd = 0.6)
  seg3(xl[1], yl[1], 0, xl[2], yl[1], 0, col = "grey45", lty = 3, lwd = 0.6)

  # faint drop line to the floor, for depth perception. Drawn by hand rather than via type = "h"
  # because scatterplot3d colours the h-lines like the points, and a thick coloured vertical then
  # reads as a fourth CI whisker.
  for (i in seq_len(nrow(b3))) {
    x <- b3$cohens_d_predictor[i]; y <- b3$cohens_d_plaque[i]; z <- b3$cohens_d_depth[i]
    if (!all(is.finite(c(x, y, z)))) next
    seg3(x, y, zl[1], x, y, z, col = "grey80", lty = 1, lwd = 0.4)
  }
  # 95% CI whiskers on all three axes, coloured by axis
  for (i in seq_len(nrow(b3))) {
    x <- b3$cohens_d_predictor[i]; y <- b3$cohens_d_plaque[i]; z <- b3$cohens_d_depth[i]
    if (!all(is.finite(c(x, y, z)))) next
    seg3(b3$d.L_predictor[i], y, z, b3$d.R_predictor[i], y, z, col = COL_PROX,   lwd = 0.9)
    seg3(x, b3$d.L_plaque[i], z, x, b3$d.R_plaque[i], z,       col = COL_PLAQUE, lwd = 0.9)
    seg3(x, y, b3$d.L_depth[i], x, y, b3$d.R_depth[i],         col = COL_DEPTH,  lwd = 0.9)
  }
  # points last so the whiskers do not cover them
  pts <- s3d$xyz.convert(b3$cohens_d_predictor, b3$cohens_d_plaque, b3$cohens_d_depth)
  points(pts$x, pts$y, pch = 16, cex = 0.7, col = point_col)
  lab <- b3[[label_col]]
  if (any(nzchar(lab))) {
    # alternate the label side so a tight cluster of points does not stack every label on the same
    # spot. Base graphics has no repel; the pairs companion is the figure to read when this is dense.
    side <- ifelse(seq_along(lab) %% 2 == 0, 4, 2)
    text(pts$x, pts$y, labels = lab, cex = 0.45, pos = side, offset = 0.3,
         col = "grey15", xpd = NA)
  }

  # Legends go ABOVE the cube, in the reserved top margin. The bottom strip is not usable:
  # scatterplot3d puts the x-axis title there and a negative "bottom" inset large enough to clear it
  # lands the legend off-canvas. Negative inset from "top" moves upward, so -0.105 is the upper row.
  graphics::legend("top", bty = "n", cex = 0.45, horiz = TRUE, inset = c(0, -0.105),
                   legend = c("95% CI: tau", "95% CI: plaque", "95% CI: depth"),
                   lty = 1, lwd = 1.1, col = c(COL_PROX, COL_PLAQUE, COL_DEPTH), xpd = NA)
  if (!is.null(legend_lab))
    graphics::legend("top", bty = "n", cex = 0.45, inset = c(0, -0.055), pch = 16,
                     ncol = min(4, length(legend_lab)),
                     legend = legend_lab, col = legend_cols, xpd = NA)
  dev.off()
}

# --- the 2-D companion. Three pairwise views of the same three coefficients, with CI crossbars on
# --- both axes.
pairs_panel <- function(b3, xn, yn, xlab, ylab, colour_by, colour_vals, colour_name,
                        label_col = "outcome") {
  d <- b3 %>% transmute(lab = .data[[label_col]], grp = .data[[colour_by]],
                        x  = .data[[paste0("cohens_d_", xn)]],
                        xl = .data[[paste0("d.L_", xn)]], xr = .data[[paste0("d.R_", xn)]],
                        y  = .data[[paste0("cohens_d_", yn)]],
                        yl = .data[[paste0("d.L_", yn)]], yr = .data[[paste0("d.R_", yn)]]) %>%
    pad_factor_levels("grp", c("x", "xl", "xr", "y", "yl", "yr"))
  ggplot(d, aes(x, y, colour = grp)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
    geom_linerange(aes(ymin = yl, ymax = yr), linewidth = 0.3, alpha = 0.7) +
    geom_linerange(aes(xmin = xl, xmax = xr), orientation = "y",
                   linewidth = 0.3, alpha = 0.7) +
    geom_point(size = 1.2) +
    ggrepel::geom_text_repel(aes(label = lab), size = 1.7, max.overlaps = Inf,
                             min.segment.length = 0, segment.size = 0.2,
                             segment.colour = "grey55", box.padding = 0.35,
                             seed = 42, show.legend = FALSE) +
    scale_colour_manual(values = colour_vals, name = colour_name, drop = FALSE,
                        limits = if (is.factor(d$grp)) levels(d$grp) else names(colour_vals)) +
    labs(x = xlab, y = ylab) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.margin = margin(4, 6, 4, 4, unit = "pt"))
}

make_pairs <- function(b3, colour_by, colour_vals, colour_name, xlab_prox,
                       label_col = "outcome") {
  p1 <- pairs_panel(b3, "predictor", "plaque", xlab_prox, "d, plaque-prox vs distal",
                    colour_by, colour_vals, colour_name, label_col)
  p2 <- pairs_panel(b3, "predictor", "depth", xlab_prox, "d, per s.d. deeper",
                    colour_by, colour_vals, colour_name, label_col)
  p3 <- pairs_panel(b3, "plaque", "depth", "d, plaque-prox vs distal", "d, per s.d. deeper",
                    colour_by, colour_vals, colour_name, label_col)
  (p1 | p2 | p3) + patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

# ===================================================================
# Emit, per panel
# ===================================================================
emit_panel <- function(R) {
  P <- R$panel; res <- R$res; slug <- P$slug
  is_state <- P$style == "state"
  xlab_fc  <- if (is_state) XLAB_STATE else XLAB_DIST
  lfc_min  <- if (is_state) LFC_MIN else 0
  full <- res %>% filter(term == "predictor", model == "m4_full")
  base <- res %>% filter(term == "predictor", model == "m1_base")
  pred <- res %>% filter(term == "predictor")

  # ---- figures
  if (P$kind == "marker") {
    ggsave(file.path(out_dir, sprintf("plot_%s_forest_full.pdf", slug)),
           make_forest(full, xlab_fc, lfc_min, at8_label = is_state),
           width = 8.5, height = 6.4, units = "cm", device = "pdf")
    ggsave(file.path(out_dir, sprintf("plot_%s_volcano_full.pdf", slug)),
           make_volcano(full, xlab_fc, lfc_min),
           width = 6.3, height = 5.4, units = "cm", device = "pdf")
    ggsave(file.path(out_dir, sprintf("plot_%s_ladder.pdf", slug)),
           make_ladder(pred, "logFC", xlab_fc, lfc_min),
           width = LADDER_W_CM, height = LADDER_H_CM, units = "cm", device = "pdf")
  } else {
    d <- R$d
    sig_full <- isTRUE(full$padj[1] < SIG_ALPHA)
    roll <- blank_sparse(roll_mean(d$dist_um, d$value, seq(0, MAX_DIST, length.out = N_GRID),
                                  HALF_WINDOW)) %>%
      mutate(dist_um = seq(0, MAX_DIST, length.out = N_GRID),
             series = P$series, sig = sig_factor(sig_full))
    save_panel(make_panel(roll, P$ylab, "Celltype", setNames(MODULE_COLOURS[1], P$series)),
               sprintf("plot_%s_full.pdf", slug), PANEL_TOTAL_W_IN, PANEL_TOTAL_H_IN)
    wr(roll %>% mutate(window_um = WINDOW_UM, significant = sig) %>%
         select(series, window_um, significant, dist_um, roll_mean, sem, n_window),
       sprintf("source_data_%s_rollmean.tsv", slug))
    wr(ring_table(d, "value_raw", P$series), sprintf("source_data_%s_persample.tsv", slug))
    ggsave(file.path(out_dir, sprintf("plot_%s_ladder.pdf", slug)),
           make_ladder(pred, "est_or", "Beta per s.d. nearer tangle (within-donor z)"),
           width = LADDER_W_CM, height = LADDER_H_GLIA_CM, units = "cm", device = "pdf")
  }

  # 3-D + pairs, this panel only. Point colour = whether the tau coefficient survives in m4.
  b3 <- beta3d_table(res)
  if (nrow(b3) >= 2) {
    sig <- !is.na(b3$padj_predictor) & b3$padj_predictor < FDR
    beta3d_plot(b3, sprintf("plot_beta3d_%s.pdf", slug),
                point_col = ifelse(sig, "#BD0026", "grey55"),
                legend_lab = c("tau padj < 0.05", "n.s."), legend_cols = c("#BD0026", "grey55"),
                xlab = if (is_state) "d: PHF1+ vs PHF1-" else "d: nearer tangle")
    b3$sig_lab <- sig_lab_factor(sig)
    ggsave(file.path(out_dir, sprintf("plot_beta_pairs_%s.pdf", slug)),
           make_pairs(b3, "sig_lab", SIG_COLOURS, "tau term (m4)",
                      if (is_state) "d, PHF1+ vs PHF1-" else "d, per s.d. nearer tangle"),
           width = 16.0, height = 6.2, units = "cm", device = "pdf")
  }

  # ---- tables
  wr(res, sprintf("source_data_%s.tsv", slug))
  wr(res %>% filter(term == "predictor") %>%
       select(outcome, model, est_or, CI.L, CI.R, cohens_d, d.L, d.R,
              logFC, logFC.L, logFC.R, pval, padj, n_cells, n_donors),
     sprintf("stats_%s_effectsize.tsv", slug))

  # ---- log
  md <- P$md
  depth_sd <- sd(md$depth_prop)
  sink(file.path(out_dir, sprintf("stats_%s.txt", slug)))
  cat(P$title, "\n"); cat(strrep("=", nchar(P$title)), "\n", sep = "")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Origin script (cell set, distance, sign convention, figure style):", P$origin, "\n\n")

  cat("PURPOSE\n-------\n")
  cat("SUPPLEMENTARY covariate sensitivity. m1_base IS the canonical headline of the origin script,\n")
  cat("reproduced here so the adjusted fits are anchored to it. docs/MODELS.md classifies BOTH added\n")
  cat("covariates as sensitivity terms that never enter a headline; nothing here changes that.\n\n")

  cat("MODEL LADDER (lmerTest::lmer, REML, random INTERCEPT only)\n")
  for (mn in MODEL_LEVELS)
    cat(sprintf("  %-10s value ~ %s\n", mn,
                paste(c(P$pterm, P$extra, MODEL_EXTRA[[mn]],
                        "Sex", "Age_s", "PMI_s", "(1 | patient_id)"), collapse = " + ")))
  cat("\nCOVARIATES ADDED\n----------------\n")
  cat("  plaque  = Matched_4G8_40, a logical in colData(spe): nucleus inside the 40 um-DILATED 4G8+\n")
  cat("            amyloid mask (IMCDataAnalysis-main/14-PlaqueNiche.Rmd:19-59). Coded\n")
  cat("            factor(c('plaque-distal','plaque-prox')), so the coefficient row is\n")
  cat("            'plaqueplaque-prox' and the estimate is prox MINUS distal.\n")
  cat(sprintf("            In this cell set: %d plaque-prox, %d plaque-distal (%.1f%% prox).\n",
              sum(md$plaque == "plaque-prox"), sum(md$plaque == "plaque-distal"),
              100 * mean(md$plaque == "plaque-prox")))
  cat("  depth_s = relative cortical depth, i.e. DISTANCE FROM THE PIAL SURFACE AS A PROPORTION.\n")
  cat("            Not stored in the object; recomputed per ROI as\n")
  cat("              scaled_Y = 100 * (y - min(y)) / (max(y) - min(y))   per sample_id (ROI)\n")
  cat("            then depth_prop = scaled_Y/100 and depth_s = scale(depth_prop) within this set.\n")
  cat("            0 = PIA. Positive coefficient = higher DEEPER, away from the pia. Never flipped.\n")
  cat(sprintf("            depth_prop: median %.3f, IQR %.3f-%.3f; sd = %.4f, so one s.d. of\n",
              median(md$depth_prop), quantile(md$depth_prop, .25), quantile(md$depth_prop, .75),
              depth_sd))
  cat(sprintf("            depth_s = %.1f%% of the ROI's pia-to-deep extent. Divide a depth_s\n",
              100 * depth_sd))
  cat(sprintf("            coefficient by %.4f to read it per unit proportion, or by %.4f for\n",
              1 / depth_sd, 1 / (10 * depth_sd)))
  cat("            per 10 percentage points of depth.\n")
  cat("            NOTE: scaled_Y is a proportion of the ROI's OWN y-extent, not of true\n")
  cat("            pia-to-white-matter thickness. An ROI that does not span the full cortical\n")
  cat("            ribbon has its depths stretched onto 0-1 regardless. The raw pial distance in\n")
  cat("            microns (pia_dist_um = y - min(y)) is reported below for comparison but is NOT\n")
  cat("            modelled, because it is not comparable across ROIs of different heights.\n")
  cat(sprintf("            pia_dist_um: median %.0f um, range %.0f-%.0f um.\n",
              median(md$pia_dist_um), min(md$pia_dist_um), max(md$pia_dist_um)))
  cat("            depth_s is standardised WITHIN this cell set, exactly as dist_z is, so its\n")
  cat("            magnitude is in this set's own s.d. units and is not comparable across panels.\n")
  cat("            Standardising is a linear rescaling and moves no other coefficient; it is done\n")
  cat("            so the covariate is on a unit scale alongside dist_z (the reason\n")
  cat("            imc_phf1_glia_distance.R z-scores it: glmer otherwise warns 'Rescale\n")
  cat("            variables?' when scaled_Y is left on 0-100).\n\n")

  cat("CELL SET\n--------\n")
  cat(sprintf("  Cells %d | donors %d | ROIs %d\n", nrow(md),
              dplyr::n_distinct(md$patient_id), dplyr::n_distinct(md$sample_id)))
  cat("  Subclusters:\n"); print(sort(table(md$celltype_clusters), decreasing = TRUE))
  if (!is.null(P$extra) && "subcluster" %in% P$extra)
    cat("  Reference subcluster:", levels(md$subcluster)[1], "\n")
  if (P$pterm %in% c("dist_z", "dist_scaled")) {
    cat(sprintf("\n  Distance to nearest PHF1+ NEURON, within ROI, cap %d um.\n", MAX_DIST))
    cat("  Transform log(dist)/sd(log(dist)) -- plain log, /SD, NOT centred (canonical).\n")
    if (P$pterm == "dist_z") {
      cat(sprintf("    dist_sd = %.5f, computed WITHIN this cell set.\n", sd(log(md$dist_um))))
    } else {
      cat(sprintf("    dist_sd = %.5f, computed ONCE over all capped cells (not per glial subset),\n",
                  DIST_SD_GLIA))
      cat("    as in imc_phf1_glia_distance.R, so both glia panels share one dist_sd.\n")
    }
  }
  if (P$pterm == "grp")
    cat(sprintf("\n  PHF1+ %d | PHF1- %d cells.\n", sum(md$grp == "PHF1+"), sum(md$grp == "PHF1-")))

  cat("\nCOLLINEARITY OF THE THREE PREDICTORS\n------------------------------------\n")
  cat("If the added covariates are strongly related to the predictor of interest, an attenuated\n")
  cat("coefficient in m4 is shared variance being reassigned, not the effect being explained away.\n")
  pvec <- if (P$pterm == "grp") as.numeric(md$grp == "PHF1+") else md[[P$pterm]]
  plq  <- as.numeric(md$plaque == "plaque-prox")
  cat(sprintf("  Spearman rho(predictor, depth_prop) = %+.4f\n",
              suppressWarnings(cor(pvec, md$depth_prop, method = "spearman"))))
  cat(sprintf("  Spearman rho(predictor, plaque)     = %+.4f\n",
              suppressWarnings(cor(pvec, plq, method = "spearman"))))
  cat(sprintf("  Spearman rho(depth_prop, plaque)    = %+.4f\n",
              suppressWarnings(cor(md$depth_prop, plq, method = "spearman"))))
  cat(sprintf("  mean depth_prop: plaque-prox %.3f vs plaque-distal %.3f\n",
              mean(md$depth_prop[plq == 1]), mean(md$depth_prop[plq == 0])))

  cat("\nSIGN CONVENTIONS\n----------------\n")
  if (P$flip == -1 && P$pterm != "grp") {
    cat("  Predictor reported per s.d. NEARER a tangle (est_or = -beta), so POSITIVE = higher near\n")
    cat("  tangles.")
    if (P$kind == "glia")
      cat(" NOTE the origin script reports the OPPOSITE orientation (positive = higher\n  FURTHER); it is flipped here so all seven panels point the same way and can share one\n  3-D axis. The origin-orientation value is -est_or.\n")
    else cat(" Same orientation as the origin script.\n")
  } else {
    cat("  Predictor is the grp contrast, PHF1+ minus PHF1-, read straight off the coefficient\n")
    cat("  table (deliberately NOT emmeans::pairs(), which returns PHF1- minus PHF1+).\n")
  }
  cat("  plaque and depth are NEVER flipped: plaque = prox minus distal, depth = per s.d. deeper.\n")
  if (P$kind == "marker") {
    cat("  logFC = log2( sinh(ref + beta) / sinh(ref) ); the model is fitted on asinh intensity so\n")
    cat("    sinh() is the exact inverse and the ratio is a genuine fold change. ref =")
    cat(if (is_state) " observed mean\n    arcsinh in the PHF1- cells.\n" else " mean arcsinh over the modelled cells.\n")
  } else {
    cat("  logFC is NA: the outcome is a within-donor z-score, not an arcsinh intensity, so there\n")
    cat("    is no intensity scale to express a fold change on. Read est_or and cohens_d.\n")
  }
  cat("  cohens_d = est_or / sqrt(sum of all VarCorr variances) FROM THE SAME MODEL as est_or.\n")

  cat("\nMULTIPLE TESTING\n----------------\n")
  cat(sprintf("  BH within this panel across the %d outcome(s), SEPARATELY for each (model, term):\n",
              length(P$outcomes)))
  cat("  the base model's predictor p-values are adjusted among themselves, the full model's among\n")
  cat("  themselves, and the plaque and depth coefficients each form their own family.\n")
  if (length(P$outcomes) == 1) cat("  One outcome here, so BH is a no-op; kept explicit.\n")
  if (is_state) cat(sprintf("  padj floored at %s; effect-size floor |log2FC| >= %.2f.\n",
                            format(PADJ_FLOOR, scientific = TRUE), LFC_MIN))

  cat("\n=== PREDICTOR OF INTEREST ACROSS THE LADDER ===\n")
  cat("The question this script exists to answer. Compare m1_base with m4_full.\n")
  tab <- pred %>%
    transmute(outcome, model, est = est_or, CI.L, CI.R, cohens_d, logFC, pval, padj) %>%
    arrange(outcome, model)
  print(as.data.frame(tab), row.names = FALSE, digits = 3)

  cat("\n=== SHIFT FROM BASE TO FULLY ADJUSTED ===\n")
  shift <- base %>% select(outcome, base_est = est_or, base_padj = padj) %>%
    left_join(full %>% select(outcome, full_est = est_or, full_padj = padj), by = "outcome") %>%
    mutate(delta = full_est - base_est,
           pct_of_base = ifelse(base_est != 0, 100 * delta / abs(base_est), NA_real_),
           sign_flip = sign(full_est) != sign(base_est),
           sig_base = base_padj < FDR, sig_full = full_padj < FDR,
           lost_sig = sig_base & !sig_full, gained_sig = !sig_base & sig_full) %>%
    arrange(desc(abs(pct_of_base)))
  print(as.data.frame(shift), row.names = FALSE, digits = 3)
  cat(sprintf("\n  Sign flips: %d | lost significance: %d | gained significance: %d | of %d\n",
              sum(shift$sign_flip, na.rm = TRUE), sum(shift$lost_sig, na.rm = TRUE),
              sum(shift$gained_sig, na.rm = TRUE), nrow(shift)))
  cat("  A large shift means the headline coefficient was partly carrying amyloid proximity or\n")
  cat("  laminar position. A small shift means the two covariates are close to orthogonal to the\n")
  cat("  predictor here -- read it together with the collinearity block above.\n")

  cat("\n=== THE ADDED COVARIATES' OWN EFFECTS (from m4_full) ===\n")
  print(as.data.frame(res %>% filter(term != "predictor", model == "m4_full") %>%
                        transmute(outcome, term = term_label, est = est_or, CI.L, CI.R,
                                  cohens_d, pval, padj) %>%
                        arrange(term, pval)), row.names = FALSE, digits = 3)

  cat("\n=== RANDOM EFFECTS (donor) ===\n")
  cat("A variance at ~0 with singular TRUE means the donor random-intercept variance is estimated\n")
  cat("at zero.\n")
  print(as.data.frame(pred %>% select(outcome, model, re_var_donor, re_pct_donor, singular) %>%
                        arrange(outcome, model)), row.names = FALSE, digits = 4)
  if (P$kind == "glia") {
    cat("\nGlia panels: the outcome is z-scored within donor, so it has no between-donor variance;\n")
    cat("the donor random-intercept variance is estimated at zero and the fit is equivalent to\n")
    cat("pooled least squares on the standardised outcome (docs/MODELS.md, Set 3, within-donor\n")
    cat("z-scoring). The z-score standardises the outcome, not the predictor, so the donor-level\n")
    cat("covariates are retained. Full block for m1 and m4 below.\n")
    cat("lme4's 'Model failed to converge with 1 negative eigenvalue' warning on these fits\n")
    cat("reflects the same zero donor variance seen from the Hessian.\n")
  }
  # fit_ladder deliberately does not retain the merMod objects (10 markers x 4 models x 7 panels is
  # a lot of memory to hold for a printout); the per-row variance, % of total and singularity flag
  # it DOES retain are in the table above. cat_re_variance() needs a real object, so m1 and m4 are
  # refitted here for the FIRST outcome only, for the donor random-effect block reported in every
  # IMC log (R/imc_utils.R).
  for (mn in c("m1_base", "m4_full")) {
    f_re <- as.formula(paste("value ~", paste(c(P$pterm, P$extra, MODEL_EXTRA[[mn]],
                                                "Sex", "Age_s", "PMI_s", "(1 | patient_id)"),
                                              collapse = " + ")))
    m_re <- tryCatch(lmerTest::lmer(f_re, data = P$prep(P$md, P$outcomes[1]), REML = TRUE),
                     error = function(e) NULL)
    cat_re_variance(m_re, sprintf("Random effects (donor), %s, outcome %s",
                                  mn, P$outcomes[1]))
  }

  cat("\nNOTES\n-----\n")
  cat("- These are SENSITIVITY fits. The headline model of every panel is m1_base.\n")
  cat("- AT8 is a POSITIVE CONTROL: it stains the same tau species as PHF1.\n")
  cat("- IMC ROIs are cortical strips (median 478 x 2525 um), so distance to a tangle and laminar\n")
  cat("  position can covary; the depth term is added to test this.\n")
  cat("- The 40 um dilation means 'plaque-prox' is a niche label, not contact: a nucleus up to\n")
  cat("  40 um from 4G8+ signal counts as prox.\n")
  if (P$kind == "glia") {
    cat("- This panel pools two clusters and every model above leaves them pooled, matching\n")
    cat("  imc_phf1_glia_distance.R; no subcluster term is included.\n")
  }
  cat("\n=== sessionInfo() ===\n"); print(sessionInfo())
  sink()
  invisible(shift)
}

for (R in RUNS) emit_panel(R)

# ===================================================================
# Cross-panel outputs
# ===================================================================
# A SUBSET run must not clobber the full run's cross-panel artefacts. These four files and the two
# combined figures are only meaningful over all seven panels, so a subset run writes them under a
# _subset suffix instead -- otherwise `Rscript ... dist_reln_calb1` would silently leave
# stats_beta3d.tsv and plot_beta3d_combined.pdf holding one panel while looking complete.
SUB <- !is.null(PANEL_SELECT)
xf <- function(stem, ext) sprintf("%s%s.%s", stem, if (SUB) "_subset" else "", ext)
if (SUB) message("SUBSET run: cross-panel outputs written with a _subset suffix, not overwriting ",
                 "the full-run files.")

wr(ALL %>% select(panel, outcome, model, term, term_label, coef_row, est, est_or, SE, df, t,
                  pval, padj, CI.L, CI.R, cohens_d, d.L, d.R, logFC, logFC.L, logFC.R,
                  sd_tot, refmean, n_cells, n_donors, n_rois,
                  re_var_donor, re_pct_donor, singular),
   xf("stats_covariate_ladder", "tsv"))

B3_ALL <- beta3d_table(ALL)
wr(B3_ALL, xf("stats_beta3d", "tsv"))

if (nrow(B3_ALL) >= 2) {
  pan <- unique(B3_ALL$panel)
  cols <- setNames(PANEL_COLOURS[seq_along(pan)], pan)
  B3_ALL$panel_lab <- factor(B3_ALL$panel, levels = pan)
  beta3d_plot(B3_ALL, xf("plot_beta3d_combined", "pdf"),
              point_col = unname(cols[B3_ALL$panel]), label_col = "lab3d",
              legend_lab = pan, legend_cols = unname(cols[pan]),
              xlab = "d: tau proximity")
  ggsave(file.path(out_dir, xf("plot_beta_pairs_combined", "pdf")),
         make_pairs(B3_ALL, "panel_lab", cols, "Panel", "d, tau proximity",
                    label_col = "lab3d"),
         width = 16.0, height = 7.4, units = "cm", device = "pdf")
}

# ---- overview log, including the regression check
sink(file.path(out_dir, xf("stats_covariate_sensitivity_overview", "txt")))
cat("IMC covariate sensitivity -- overview and regression check\n")
cat("=========================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Script: R/imc_covariate_sensitivity.R\n")
cat("Panels run:", paste(names(PANELS), collapse = ", "), "\n\n")

cat("REGRESSION CHECK\n----------------\n")
cat("m1_base is the canonical headline of each origin script and is checked against the\n")
cat("published estimate.\n")
cat("Marker panels are checked against their published source-data TSVs (tolerance 1e-6); the two\n")
cat("glia panels have no per-marker TSV, so they are checked against the estimates quoted to 5\n")
cat("d.p. in plots/imc_phf1_glia_distance/stats_*.txt (tolerance 1e-4).\n")
cat("Glia rows are compared in the ORIGIN orientation (per s.d. FURTHER), which is why `mine` is\n")
cat("un-flipped there.\n\n")
if (nrow(check_rows)) {
  print(as.data.frame(check_rows %>%
    select(panel, outcome, mine, ref, delta, tol, ok, n_cells, n_donors,
           ref_n_cells, ref_n_donors, source)), row.names = FALSE, digits = 8)
  cat(sprintf("\n  PASS %d / %d\n", sum(check_rows$ok, na.rm = TRUE), nrow(check_rows)))
  if (!all(check_rows$ok, na.rm = TRUE))
    cat("  *** FAILED -- do not use the adjusted fits until this is resolved. ***\n")
} else cat("  (no reference available for the selected panels)\n")

cat("\nCELL-SET SIZES\n--------------\n")
for (R in RUNS) {
  md <- R$panel$md
  cat(sprintf("  %-26s cells %6d  donors %3d  ROIs %4d  depth_prop sd %.4f", R$panel$slug,
              nrow(md), dplyr::n_distinct(md$patient_id), dplyr::n_distinct(md$sample_id),
              sd(md$depth_prop)))
  if (R$panel$pterm == "dist_z") cat(sprintf("  dist_sd %.5f", sd(log(md$dist_um))))
  if (R$panel$pterm == "dist_scaled") cat(sprintf("  dist_sd %.5f (shared)", DIST_SD_GLIA))
  cat("\n")
}

cat("\nBASE -> FULLY ADJUSTED, ALL PANELS\n----------------------------------\n")
cat("Counts of outcomes whose predictor coefficient changes materially once plaque and depth are\n")
cat("both in the model.\n")
summ <- ALL %>% filter(term == "predictor", model %in% c("m1_base", "m4_full")) %>%
  select(panel, outcome, model, est_or, padj) %>%
  pivot_wider(names_from = model, values_from = c(est_or, padj)) %>%
  mutate(delta = est_or_m4_full - est_or_m1_base,
         pct = ifelse(est_or_m1_base != 0, 100 * delta / abs(est_or_m1_base), NA_real_),
         sign_flip = sign(est_or_m4_full) != sign(est_or_m1_base),
         lost_sig = padj_m1_base < FDR & !(padj_m4_full < FDR),
         gained_sig = !(padj_m1_base < FDR) & padj_m4_full < FDR)
print(as.data.frame(summ %>% group_by(panel) %>%
  summarise(n = n(), sig_base = sum(padj_m1_base < FDR), sig_full = sum(padj_m4_full < FDR),
            sign_flips = sum(sign_flip, na.rm = TRUE),
            lost = sum(lost_sig, na.rm = TRUE), gained = sum(gained_sig, na.rm = TRUE),
            median_abs_pct_shift = median(abs(pct), na.rm = TRUE), .groups = "drop")),
  row.names = FALSE, digits = 3)

cat("\nOutcomes with a sign flip or a change in significance:\n")
flagged <- summ %>% filter(sign_flip | lost_sig | gained_sig)
if (nrow(flagged)) print(as.data.frame(flagged), row.names = FALSE, digits = 3) else
  cat("  (none)\n")

cat("\nTHE ADDED COVARIATES' OWN EFFECTS (m4_full, padj < 0.05 only)\n")
cat("------------------------------------------------------------\n")
covs <- ALL %>% filter(term != "predictor", model == "m4_full", padj < FDR) %>%
  transmute(panel, outcome, term = term_label, est = est_or, CI.L, CI.R, cohens_d, padj) %>%
  arrange(term, panel, padj)
if (nrow(covs)) print(as.data.frame(covs), row.names = FALSE, digits = 3) else
  cat("  (none reach padj < 0.05)\n")

cat("\nFIGURE INVENTORY\n----------------\n")
cat("  plot_<slug>_forest_full.pdf / _volcano_full.pdf   marker panels redrawn from m4_full\n")
cat("  plot_<slug>_full.pdf                             glia curves; the CURVE IS UNCHANGED (it is\n")
cat("                                                   model-free), only the FDR linetype comes\n")
cat("                                                   from m4_full\n")
cat("  plot_<slug>_ladder.pdf                           m1..m4 side by side -- the panel that\n")
cat("                                                   actually shows what the covariates do\n")
cat("  plot_beta3d_*.pdf                                Cohen's d for the three coefficients of\n")
cat("                                                   m4_full, with 95% CI whiskers on all three\n")
cat("                                                   axes (whisker colour = axis). d rather than\n")
cat("                                                   raw beta so arcsinh-scale marker panels and\n")
cat("                                                   within-donor-z glia panels share one cube.\n")
cat("  plot_beta_pairs_*.pdf                            the same three coefficients as three 2-D\n")
cat("                                                   scatters -- readable at print size, which a\n")
cat("                                                   cube with three sets of whiskers is not\n")
cat("\nNOTE ON THE 3-D AXES: all three axes are Cohen's d, so a point far from the origin on the tau\n")
cat("axis and at zero on the other two is an effect specific to tangle proximity; a point out along\n")
cat("the plaque or depth axis is a marker tracking amyloid or lamina instead. Because the axes are\n")
cat("standardised by the OUTCOME's total SD from the same model, they are comparable in kind but\n")
cat("each predictor's 'one unit' still differs (one s.d. of log-distance, one plaque category step,\n")
cat("one s.d. of relative depth) -- and dist_sd/depth_sd are per-cell-set, so magnitudes are not\n")
cat("comparable ACROSS panels. See the per-panel logs.\n")
cat("\n=== sessionInfo() ===\n"); print(sessionInfo())
sink()

cat("\n== done ==\n")
cat("Output:", out_dir, "\n")
if (nrow(check_rows))
  cat(sprintf("Regression check: %d / %d pass\n", sum(check_rows$ok, na.rm = TRUE),
              nrow(check_rows)))
print(as.data.frame(ALL %>% filter(term == "predictor", model %in% c("m1_base", "m4_full")) %>%
                      select(panel, outcome, model, est_or, padj) %>%
                      arrange(panel, outcome, model)), row.names = FALSE, digits = 3)
