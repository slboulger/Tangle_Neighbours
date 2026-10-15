#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# imc_cluster_loss_vs_phf1_intensity.R
#
# Figure panels: 4B
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# imc_cluster_loss_vs_phf1_intensity.R
#
# Are the EC neuronal clusters that are LOST most in AD the ones that carry the most PHF1 in AD?
# Sister analysis to R/imc_cluster_loss_vs_tangle_distance.R (Fig. 4A) -- same x, same test, same
# plot style; only y changes.
#
#   x = % difference in columnar density, AD vs control, per cluster -- the cluster percentage-change
#       bar of Fig. 2a of the previously published IMC study of this cohort (see README), computed
#       as in its analysis code (IMCDataAnalysis-main/12b-Cluster_Percentage_Change.Rmd:39-103). The
#       x code below is identical to the distance script and gated against the IMC study's values.
#   y = mean within-donor z-scored PHF1 intensity per cluster, AD donors only.
#
# ---------------------------------------------------------------------------
# PHF1 NORMALISATION
# ---------------------------------------------------------------------------
# Intensity = assay "exprs" = asinh(counts / 1) (04-read_data.Rmd:311), the assay every IMC statistical
# analysis models (NOT exprs_masked, 79% NA). It is then Z-SCORED WITHIN DONOR over the analysed
# population, as in R/imc_phf1_glia_distance.R (`group_by(patient_id) %>% mutate(scale(...))`).
# Without it, donors with heavier tangle burden also stain globally brighter and between-donor
# variation swamps within-donor contrasts.
#   Analysed population = every cell of the 10 neuronal clusters in each AD donor, PHF1+ and PHF1-
#   alike, so a cluster's y says where its cells sit in that donor's neuronal PHF1 distribution.
#   Glia are excluded from the reference because they are excluded from the analysis.
# A cluster mean of z combines TWO things: how many of its cells bear tangles (the PHF1_Otsu+
# fraction is reported alongside) and how bright its PHF1- cells are. The `neg` row separates the
# two: PHF1- cells only, same z scale.
#
# ---------------------------------------------------------------------------
# THE TEST -- identical to the distance script
# ---------------------------------------------------------------------------
# Spearman rho across the 10 neuronal clusters; donor-level STRATIFIED bootstrap (control and AD
# donors resampled separately) recomputing x and y from the same resampled donors; percentile 95% CI;
# CI-inversion bootstrap p; cluster-label permutation p as the fixed-point reference. Why not
# Deming/York: see the distance script header (correlated x/y errors via shared AD donors, clusters
# sharing donors, the control-mean denominator).
#
# SIGN: rho < 0 = the most-lost clusters carry the MOST PHF1. This is the opposite sign convention
# to the distance script, where nearer = smaller y.
#
# Pre-specified sensitivities (rows of stats_*_sensitivity.tsv):
#   neg      -- PHF1_Otsu-negative cells only (no tangle-bearers).
#   layer    -- z minus the mean z of all neuronal cells in the same ROI x layer.
#   depth    -- z minus a per-ROI natural spline (df = DEPTH_DF) of z on relative cortical depth;
#               df sweep in depth_df* rows. Same construction as the distance script.
#   age, zerofill, exc_only(_layer/_depth), leave-one-cluster-out -- as the distance script.
#
# ---------------------------------------------------------------------------
# Alternative outcomes (not drawn in the paper):
# --outcome=pct_plaque: y = % of the cluster's cells with Matched_4G8_40 (nucleus within the 40 um-dilated
# 4G8+ plaque mask). Same construction as pct_pos; outputs plots/imc_cluster_loss_vs_plaque_pct/.
# --outcome=pct_pos: y = % of the cluster's cells that are PHF1_Otsu+ (the per-donor Otsu call on this
# same channel, floored at 0.88), per AD donor, then averaged over donors.
# No z-score: a proportion is already on a common scale across donors. layer/depth centre the 0/100
# indicator exactly as they centre z, so those rows read as "percentage points above the PHF1+ rate of
# neuronal cells at the same layer/depth of the same ROI". There is no `neg` row (it would be 0).
# The denominator of % PHF1+ is the surviving cells: preferential death of tangle-free cells of a
# cluster raises it, death of tangle-bearing cells (ghost tangles are not segmented) lowers it.
# Outputs go to plots/imc_cluster_loss_vs_phf1_pct/.
#
# Input: <IMC_ROOT>/spe.rds. Reproduction check (GATE 1): if
# <IMC_ROOT>/Figures/SupplementaryTable4_patient_percent_change.csv -- Supplementary Table 4 of the
# IMC study, the per-donor % change per cluster -- is present, x must match it.
#
# Run: Rscript R/imc_cluster_loss_vs_phf1_intensity.R
#      Rscript R/imc_cluster_loss_vs_phf1_intensity.R --outcome=pct_pos
suppressPackageStartupMessages({
  library(SpatialExperiment); library(dplyr); library(tibble); library(tidyr)
  library(ggplot2); library(splines)
})

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")

args <- commandArgs(trailingOnly = TRUE)
OUTCOME <- sub("^--outcome=", "", grep("^--outcome=", args, value = TRUE))
OUTCOME <- if (length(OUTCOME)) OUTCOME else "intensity"
stopifnot(OUTCOME %in% c("intensity", "pct_pos", "pct_plaque"))
IS_INT <- OUTCOME == "intensity"
IS_PLQ <- OUTCOME == "pct_plaque"
stem <- switch(OUTCOME, intensity = "cluster_loss_vs_phf1_intensity", pct_pos = "cluster_loss_vs_phf1_pct",
               pct_plaque = "cluster_loss_vs_plaque_pct")
out_dir <- paste0("plots/imc_", stem)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE)

B           <- 9999     # donor bootstrap replicates
B_PERM      <- 199999   # cluster-label permutations (k = 10; Monte Carlo)
DEPTH_DF    <- 6        # spline df for the per-ROI depth expectation of z
DEPTH_DF_SWEEP <- c(3, 10, 15)
MIN_CELLS   <- 5        # cells a donor x cluster needs to contribute a y
SEED        <- 1
EXCLUDE     <- c("Artefact cluster", "Unassigned cluster")
MARKER      <- "PHF1"
set.seed(SEED)

localwd <- "<IMC_ROOT>/"
ref_pct <- paste0(localwd, "Figures/SupplementaryTable4_patient_percent_change.csv")
spe <- readRDS(paste0(localwd, "spe.rds"))
stopifnot(MARKER %in% rownames(spe), "exprs" %in% assayNames(spe))

cd <- as.data.frame(colData(spe))
cd$cell_id <- colnames(spe)
sc <- spatialCoords(spe)
cd$x <- sc[, "Pos_X"]; cd$y <- sc[, "Pos_Y"]
cd <- cd %>% mutate(patient_id = as.character(patient_id), sample_id = as.character(sample_id),
                    celltype_clusters = as.character(celltype_clusters))
# relative cortical depth: a property of the ROI, so computed over ALL its cells before any filter
cd <- cd %>% group_by(sample_id) %>% mutate(scaled_Y = 100 * (y - min(y)) / (max(y) - min(y))) %>% ungroup()

cv_removal <- cd %>% distinct(patient_id, sample_id) %>% count(patient_id) %>%
  filter(n < 3) %>% pull(patient_id)

short_lab <- function(cl) sub("^.*\\((.*)\\)$", "\\1", cl)
cl_class  <- function(cl) ifelse(grepl("^Excitatory", cl), "Excitatory",
                          ifelse(grepl("^Inhibitory", cl), "Inhibitory", "Glia"))

# ===================================================================
# x: columnar density, as in the IMC study's 12b-Cluster_Percentage_Change.Rmd:39-71
# (identical to imc_cluster_loss_vs_tangle_distance.R; GATE 1 holds both to the IMC study's values)
# ===================================================================
sample_meta <- cd %>% select(sample_id, patient_id, Age, Braak, BraakGroup, width_px) %>% distinct()
stopifnot(!anyDuplicated(sample_meta$sample_id))

roi_dens <- cd %>% count(sample_id, celltype_clusters, name = "n_cells") %>%
  left_join(sample_meta, by = "sample_id") %>%
  mutate(counts_per_mm = n_cells / width_px * 1e3) %>%
  filter(!patient_id %in% cv_removal, !celltype_clusters %in% EXCLUDE)

roi_dens_zf <- tidyr::expand_grid(sample_id = unique(roi_dens$sample_id),
                                  celltype_clusters = sort(unique(roi_dens$celltype_clusters))) %>%
  left_join(roi_dens %>% select(sample_id, celltype_clusters, n_cells),
            by = c("sample_id", "celltype_clusters")) %>%
  mutate(n_cells = coalesce(n_cells, 0L)) %>%
  left_join(sample_meta, by = "sample_id") %>%
  mutate(counts_per_mm = n_cells / width_px * 1e3)

donor_mean <- function(rd) rd %>% group_by(patient_id, celltype_clusters) %>%
  summarise(dens = mean(counts_per_mm), .groups = "drop")

donor_meta <- sample_meta %>% filter(!patient_id %in% cv_removal) %>%
  distinct(patient_id, Age, Braak, BraakGroup) %>%
  mutate(AD = Braak >= 2)                          # Rmd: braak_cutoff <- 2
stopifnot(!anyDuplicated(donor_meta$patient_id))
ctrl_ids <- sort(donor_meta$patient_id[!donor_meta$AD])
ad_ids   <- sort(donor_meta$patient_id[donor_meta$AD])
stopifnot(length(ctrl_ids) == 18, length(ad_ids) == 44,
          all(donor_meta$BraakGroup[donor_meta$AD] != "Braak_0_1"))

clusters_all <- sort(unique(roi_dens$celltype_clusters))
stopifnot(length(clusters_all) == 14)
neur  <- clusters_all[cl_class(clusters_all) != "Glia"]
exc   <- clusters_all[cl_class(clusters_all) == "Excitatory"]
stopifnot(length(neur) == 10, length(exc) == 4)

to_mat <- function(dd, ids, val = "dens") {
  m <- matrix(NA_real_, length(ids), length(clusters_all), dimnames = list(ids, clusters_all))
  dd <- dd[dd$patient_id %in% ids, ]
  m[cbind(match(dd$patient_id, ids), match(dd$celltype_clusters, clusters_all))] <- dd[[val]]
  m
}
dm      <- donor_mean(roi_dens);    dm_zf <- donor_mean(roi_dens_zf)
D_ctrl  <- to_mat(dm, ctrl_ids);    D_ad  <- to_mat(dm, ad_ids)
Dz_ctrl <- to_mat(dm_zf, ctrl_ids); Dz_ad <- to_mat(dm_zf, ad_ids)

pct_change <- function(Dc, Da) 100 * (colMeans(Da, na.rm = TRUE) / colMeans(Dc, na.rm = TRUE) - 1)

# ---- GATE 1: reproduce the IMC study's Fig. 2a bars (mean + SEM, donor counts) ----
# and, when its Supplementary Table 4 CSV is present, the per-cluster mean of its per-donor values
x_hat <- pct_change(D_ctrl, D_ad)
ctrl_mean <- colMeans(D_ctrl, na.rm = TRUE)
per_donor_pct <- sweep(sweep(D_ad, 2, ctrl_mean, "-"), 2, ctrl_mean, "/") * 100
sem_pub <- apply(per_donor_pct, 2, function(v) sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v))))
n_ad_dens <- colSums(!is.na(D_ad)); n_ctrl_dens <- colSums(!is.na(D_ctrl))
stopifnot(abs(x_hat - colMeans(per_donor_pct, na.rm = TRUE)) < 1e-9)
exp_pub <- c("Excitatory neuron cluster 2 (CALB1)" = -49.7, "Inhibitory neuron cluster 2 (PVALB)" = -20.0,
             "Excitatory neuron cluster 3 (RORB)" = -18.0)
exp_sem <- c(3.9, 6.9, 4.3)
stopifnot(abs(x_hat[names(exp_pub)] - exp_pub) < 0.05, abs(sem_pub[names(exp_pub)] - exp_sem) < 0.05,
          n_ad_dens["Excitatory neuron cluster 2 (CALB1)"] == 43, all(n_ctrl_dens >= 17))
gate1_csv <- NA_real_
if (file.exists(ref_pct)) {
  ref <- read.csv(ref_pct, check.names = FALSE, colClasses = c(ID = "character")) %>%
    group_by(Cluster) %>% summarise(pct = mean(`Percentage Change`, na.rm = TRUE), n = n())
  gate1_csv <- max(abs(ref$pct - x_hat[ref$Cluster]))
  stopifnot(gate1_csv < 1e-6, all(ref$n == n_ad_dens[ref$Cluster]))
}

# ===================================================================
# y: within-donor z-scored PHF1 (exprs = asinh(counts/1)), AD donors, neuronal clusters
# ===================================================================
# same cohort frame as the distance script and imc_phf1_distance_markers_subsets.R
base <- cd %>%
  filter(!celltype_clusters %in% EXCLUDE, !patient_id %in% cv_removal,
         BraakGroup != "Braak_0_1", !is.na(PHF1_Otsu))
cells <- base %>% filter(celltype_clusters %in% neur)
cells$phf1_raw <- as.numeric(assay(spe, "exprs")[MARKER, cells$cell_id])
stopifnot(all(is.finite(cells$phf1_raw)))
# `z` is the per-cell value every downstream step averages: the within-donor z-scored intensity, or
# for --outcome=pct_pos the 0/100 PHF1_Otsu+ indicator (its donor x cluster mean is the % PHF1+).
# for --outcome=pct_plaque the 0/100 Matched_4G8_40 indicator (nucleus within the 40 um-dilated 4G8+ plaque
# mask, matched in Steinbock -- the Set 3 plaque sensitivity term), so the mean is the % plaque-adjacent.
stopifnot(!IS_PLQ || !anyNA(cells$Matched_4G8_40))
cells <- cells %>% group_by(patient_id) %>%
  mutate(z = if (IS_INT) as.numeric(scale(phf1_raw)) else
               if (IS_PLQ) 100 * as.numeric(Matched_4G8_40) else 100 * (PHF1_Otsu == "PHF1_pos")) %>% ungroup()

# ---- GATE 2: cohort and z-score sanity ----
zchk <- cells %>% group_by(patient_id) %>% summarise(m = mean(z), s = sd(z), .groups = "drop")
pos_vs_neg <- cells %>% group_by(patient_id) %>%
  summarise(dz = mean(z[PHF1_Otsu == "PHF1_pos"]) - mean(z[PHF1_Otsu == "PHF1_neg"]),
            n_pos = sum(PHF1_Otsu == "PHF1_pos"), .groups = "drop")
g2 <- c(donors = n_distinct(cells$patient_id), rois = n_distinct(cells$sample_id),
        cells = nrow(cells), phf1_pos = sum(cells$PHF1_Otsu == "PHF1_pos"))
stopifnot(g2["donors"] == 44, g2["rois"] == 132,
          !IS_INT || (max(abs(zchk$m)) < 1e-9 && max(abs(zchk$s - 1)) < 1e-9),
          # the Otsu call is made on this same channel, so PHF1+ cells must be brighter in every donor
          # that has any
          IS_PLQ || all(pos_vs_neg$dz[pos_vs_neg$n_pos > 0] > 0))   # PHF1 outcomes only

# layer: centre z on all neuronal cells of the same ROI x layer
layer_ref <- cells %>% filter(!is.na(layer)) %>% group_by(sample_id, layer) %>%
  summarise(z_ref = mean(z), .groups = "drop")
cells <- cells %>% left_join(layer_ref, by = c("sample_id", "layer")) %>% mutate(z_layer = z - z_ref)

# depth: residual of z from a per-ROI natural spline of z on relative depth (all neuronal cells)
depth_resid <- function(k) {
  out <- rep(NA_real_, nrow(cells))
  for (s_ in unique(cells$sample_id)) {
    i <- which(cells$sample_id == s_)
    if (length(i) < 10 * k) next
    out[i] <- residuals(lm(cells$z[i] ~ ns(cells$scaled_Y[i], df = k)))
  }
  out
}
cells$z_depth <- depth_resid(DEPTH_DF)
n_roi_depth <- n_distinct(cells$sample_id[!is.na(cells$z_depth)])
sweep_nm <- paste0("depth_df", DEPTH_DF_SWEEP)
for (k in DEPTH_DF_SWEEP) cells[[paste0("z_depth_df", k)]] <- depth_resid(k)

donor_y <- cells %>% group_by(patient_id, celltype_clusters) %>%
  summarise(n_cells = n(), n_neg = sum(PHF1_Otsu == "PHF1_neg"),
            raw   = mean(z),
            neg   = if (any(PHF1_Otsu == "PHF1_neg")) mean(z[PHF1_Otsu == "PHF1_neg"]) else NA_real_,
            layer = mean(z_layer, na.rm = TRUE),
            depth = mean(z_depth, na.rm = TRUE),
            across(all_of(paste0("z_", sweep_nm)), ~ mean(.x, na.rm = TRUE), .names = "{sub('^z_', '', .col)}"),
            .groups = "drop") %>%
  mutate(across(c(raw, layer, depth, all_of(sweep_nm)), ~ ifelse(n_cells >= MIN_CELLS, .x, NA_real_)),
         neg = ifelse(n_neg >= MIN_CELLS, neg, NA_real_))
y_vars <- c("raw", if (IS_INT) "neg", "layer", "depth", sweep_nm)
Y <- lapply(setNames(nm = y_vars),
            function(v) to_mat(donor_y, ad_ids, v))
n_ad_y <- colSums(!is.na(Y$raw))[neur]
low_cov <- names(n_ad_y)[n_ad_y < 30]

phf1_frac <- cells %>% group_by(patient_id, celltype_clusters) %>%
  summarise(f = mean(PHF1_Otsu == "PHF1_pos"), .groups = "drop") %>%
  group_by(celltype_clusters) %>% summarise(phf1_pos_frac = mean(f), .groups = "drop")

# ===================================================================
# age-adjusted x: OLS on donor means, dens ~ AD + Age, per cluster; % = beta / control EMM
# ===================================================================
D_all <- rbind(D_ctrl, D_ad)
age_v <- donor_meta$Age[match(rownames(D_all), donor_meta$patient_id)]
ad_v  <- as.numeric(rownames(D_all) %in% ad_ids)
age_bar <- mean(age_v)
pct_age <- function(rows) {
  sapply(clusters_all, function(cl) {
    yv <- D_all[rows, cl]; ok <- !is.na(yv)
    X  <- cbind(1, ad_v[rows][ok], age_v[rows][ok])
    b  <- lm.fit(X, yv[ok])$coefficients
    unname(100 * b[2] / (b[1] + b[3] * age_bar))
  })
}

# ===================================================================
# Bootstrap
# ===================================================================
spear <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))
ic <- seq_along(ctrl_ids); ia <- seq_along(ad_ids)
n_c <- length(ic); n_a <- length(ia)

stat_once <- function(sc_, sa_) {
  x  <- pct_change(D_ctrl[sc_, , drop = FALSE], D_ad[sa_, , drop = FALSE])
  xz <- pct_change(Dz_ctrl[sc_, , drop = FALSE], Dz_ad[sa_, , drop = FALSE])
  xa <- pct_age(c(sc_, n_c + sa_))
  ys <- lapply(Y, function(m) colMeans(m[sa_, , drop = FALSE], na.rm = TRUE))
  list(x = x, xz = xz, xa = xa, y = ys)
}
rho_set <- function(s) c(
  headline = spear(s$x[neur],  s$y$raw[neur]),
  if (IS_INT) c(neg = spear(s$x[neur], s$y$neg[neur])),
  layer    = spear(s$x[neur],  s$y$layer[neur]),
  depth    = spear(s$x[neur],  s$y$depth[neur]),
  age      = spear(s$xa[neur], s$y$raw[neur]),
  zerofill = spear(s$xz[neur], s$y$raw[neur]),
  exc_only = spear(s$x[exc],   s$y$raw[exc]),
  exc_only_layer = spear(s$x[exc], s$y$layer[exc]),
  exc_only_depth = spear(s$x[exc], s$y$depth[exc]),
  setNames(sapply(sweep_nm, function(v) spear(s$x[neur], s$y[[v]][neur])), sweep_nm),
  setNames(sapply(neur, function(cl) spear(s$x[setdiff(neur, cl)], s$y$raw[setdiff(neur, cl)])),
           paste0("loco:", neur)))

obs <- stat_once(ic, ia)
rho_obs <- rho_set(obs)

boot_x <- matrix(NA_real_, B, length(neur), dimnames = list(NULL, neur))
boot_y <- lapply(Y, function(m) boot_x)
boot_rho <- matrix(NA_real_, B, length(rho_obs), dimnames = list(NULL, names(rho_obs)))
null_rho <- numeric(B)
for (b in seq_len(B)) {
  s <- stat_once(sample(ic, n_c, TRUE), sample(ia, n_a, TRUE))
  boot_rho[b, ] <- rho_set(s)
  boot_x[b, ] <- s$x[neur]
  for (v in names(Y)) boot_y[[v]][b, ] <- s$y[[v]][neur]
  null_rho[b] <- spear(s$x[neur], sample(s$y$raw[neur]))   # GATE 3a: shuffled y per replicate
}

boot_p <- function(r) { r <- r[!is.na(r)]; min(1, 2 * min((sum(r <= 0) + 1), (sum(r >= 0) + 1)) / (length(r) + 1)) }
perm_p <- function(x, y, n_perm) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; r0 <- spear(x, y); k <- length(x)
  if (k <= 8) {                                    # exact enumeration
    perms <- function(v) if (length(v) <= 1) list(v) else
      do.call(c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(p) c(v[i], p))))
    r <- sapply(perms(seq_len(k)), function(p) spear(x, y[p]))
    return(c(p = mean(abs(r) >= abs(r0) - 1e-12), n = length(r)))
  }
  rx <- rank(x); ry <- rank(y)
  r <- vapply(seq_len(n_perm), function(i) cor(rx, sample(ry)), 0)
  c(p = (sum(abs(r) >= abs(r0) - 1e-12) + 1) / (n_perm + 1), n = n_perm)
}
y_for <- c(headline = "raw", neg = "neg", layer = "layer", depth = "depth", age = "raw", zerofill = "raw",
           exc_only = "raw", exc_only_layer = "layer", exc_only_depth = "depth", setNames(sweep_nm, sweep_nm))
x_for <- c(headline = "x", neg = "x", layer = "x", depth = "x", age = "xa", zerofill = "xz",
           exc_only = "x", exc_only_layer = "x", exc_only_depth = "x",
           setNames(rep("x", length(sweep_nm)), sweep_nm))
SIGN_TXT <- if (IS_INT) "carry the most PHF1" else if (IS_PLQ) "have the most cells within 40 um of a plaque" else
  "carry the most PHF1 (highest % PHF1+)"
YW <- if (IS_INT) "z" else if (IS_PLQ) "% plaque-adjacent" else "% PHF1+"
analysis_desc <- c(
  headline = if (IS_INT) "10 neuronal clusters; y = mean within-donor z PHF1, all cells" else
             if (IS_PLQ) "10 neuronal clusters; y = % of the cluster's cells within 40 um of a 4G8+ plaque" else
             "10 neuronal clusters; y = % of the cluster's cells that are PHF1_Otsu+",
  if (IS_INT) c(neg = "y = PHF1_Otsu-negative cells only (same z scale)"),
  layer    = sprintf("y = %s centred on all neuronal cells in the same ROI x layer", YW),
  depth    = sprintf("y = %s minus per-ROI spline (df %d) of it on relative cortical depth", YW, DEPTH_DF),
  age      = "x = AD + Age fit on donor means, % of control EMM at mean age",
  zerofill = "x with absent ROI x cluster counted as 0",
  exc_only = "4 excitatory clusters; y raw (CANNOT reach p<0.05: min p = 2/24)",
  exc_only_layer = "4 excitatory clusters; y layer-centred (direction check only)",
  exc_only_depth = "4 excitatory clusters; y depth-adjusted (direction check only)",
  setNames(sprintf("as depth, spline df = %d (df sensitivity)", DEPTH_DF_SWEEP), sweep_nm))

y_for <- y_for[names(analysis_desc)]; x_for <- x_for[names(analysis_desc)]
sens <- bind_rows(lapply(names(analysis_desc), function(a) {
  cls <- if (grepl("^exc_only", a)) exc else neur
  pp  <- perm_p(obs[[x_for[a]]][cls], obs$y[[y_for[a]]][cls], B_PERM)
  r   <- boot_rho[, a]
  tibble(analysis = a, description = analysis_desc[[a]], k = length(cls), rho = rho_obs[[a]],
         ci_lo = unname(quantile(r, 0.025, na.rm = TRUE)), ci_hi = unname(quantile(r, 0.975, na.rm = TRUE)),
         p_boot = boot_p(r), p_perm = unname(pp["p"]), n_perm = unname(pp["n"]),
         frac_boot_neg = mean(r < 0, na.rm = TRUE), n_boot_na = sum(is.na(r)))
}))

loco <- tibble(dropped = neur, label = short_lab(neur), class = cl_class(neur),
               rho = rho_obs[paste0("loco:", neur)],
               ci_lo = apply(boot_rho[, paste0("loco:", neur)], 2, quantile, 0.025, na.rm = TRUE),
               ci_hi = apply(boot_rho[, paste0("loco:", neur)], 2, quantile, 0.975, na.rm = TRUE))

# ---- GATE 3b: planted relation recovered ----
# Synthetic donor-level z with the observed per-cluster donor SD, cluster means a monotone
# (more loss -> more PHF1) function of the observed x. The CI must sit entirely below 0.
sd_y <- apply(Y$raw[, neur], 2, sd, na.rm = TRUE)
# 0.5 within-donor per-cell SDs per SD of x (= 0.5 z for intensity)
plant_unit <- if (IS_INT) 1 else mean(tapply(cells$z, cells$patient_id, sd))
mu_plant <- mean(obs$y$raw[neur]) - 0.5 * plant_unit * as.numeric(scale(obs$x[neur]))
Yp <- Y$raw; Yp[, neur] <- sapply(seq_along(neur), function(j)
  ifelse(is.na(Y$raw[, neur[j]]), NA, rnorm(nrow(Yp), mu_plant[j], sd_y[j])))
plant_r <- replicate(1999, {
  sc_ <- sample(ic, n_c, TRUE); sa_ <- sample(ia, n_a, TRUE)
  spear(pct_change(D_ctrl[sc_, ], D_ad[sa_, ])[neur], colMeans(Yp[sa_, neur], na.rm = TRUE))
})
plant_ci <- quantile(plant_r, c(0.025, 0.975), na.rm = TRUE)
if (plant_ci[2] >= 0) stop("GATE 3b failed: a planted monotone relation was not recovered")
if (abs(mean(null_rho, na.rm = TRUE)) > 0.02) stop("GATE 3a failed: shuffled-y rho not centred on 0")

# ===================================================================
# Source data (per cluster) + per-donor audit table
# ===================================================================
qlo <- function(m) apply(m, 2, quantile, 0.025, na.rm = TRUE)
qhi <- function(m) apply(m, 2, quantile, 0.975, na.rm = TRUE)
panel_of <- if (IS_INT) c(raw = "phf1_z", neg = "phf1_z_negative_cells", layer = "phf1_z_layer_relative",
                          depth = "phf1_z_depth_relative") else
  if (IS_PLQ) c(raw = "pct_plaque", layer = "pct_plaque_layer_excess", depth = "pct_plaque_depth_excess") else
  c(raw = "pct_phf1_pos", layer = "pct_phf1_pos_layer_excess", depth = "pct_phf1_pos_depth_excess")
pts <- bind_rows(lapply(names(panel_of), function(v) {
  tibble(panel = panel_of[[v]],
         cluster = neur, label = short_lab(neur), class = cl_class(neur),
         pct_change = obs$x[neur], pct_lo = qlo(boot_x), pct_hi = qhi(boot_x),
         pct_sem_published = sem_pub[neur],
         y = obs$y[[v]][neur], y_lo = qlo(boot_y[[v]]), y_hi = qhi(boot_y[[v]]),
         n_ctrl_donors = n_ctrl_dens[neur], n_ad_donors_density = n_ad_dens[neur],
         n_ad_donors_phf1 = colSums(!is.na(Y[[v]]))[neur],
         n_cells = as.integer(table(factor(cells$celltype_clusters, levels = neur))))
})) %>% left_join(phf1_frac, by = c("cluster" = "celltype_clusters"))

audit <- full_join(
  dm %>% mutate(group = ifelse(patient_id %in% ad_ids, "AD", "Control")),
  donor_y %>% select(-any_of(if (IS_INT) character(0) else "neg")) %>%
    rename_with(~ paste0(if (IS_INT) "phf1_z_" else if (IS_PLQ) "pct_plaque_" else "pct_phf1_pos_", sub("^raw$", "all", .x)),
                any_of(c("raw", "neg", "layer", "depth"))),
  by = c("patient_id", "celltype_clusters")) %>%
  filter(celltype_clusters %in% neur) %>%
  rename(columnar_density_per_mm = dens) %>% arrange(celltype_clusters, group, patient_id)
wr(audit, paste0("source_data_", stem, "_donor.tsv"))
wr(sens,  paste0("stats_", stem, "_sensitivity.tsv"))
wr(loco,  paste0("stats_", stem, "_loco.tsv"))

# ===================================================================
# Plot -- identical style to the distance script
# ===================================================================
# IMC study cluster colours: IMCDataAnalysis-main/19-ColoursForPublication.Rmd:46-69.
cluster_cols <- c(
  "Excitatory neuron cluster 1 (MAP2)"      = "#00796B",
  "Excitatory neuron cluster 2 (CALB1)"     = "#26A69A",
  "Excitatory neuron cluster 3 (RORB)"      = "#66BB6A",
  "Excitatory neuron cluster 4 (RELN)"      = "#A5D6A7",
  "Inhibitory neuron cluster 1 (GAD1)"      = "#E65100",
  "Inhibitory neuron cluster 2 (PVALB)"     = "#EF6C00",
  "Inhibitory neuron cluster 3 (VIP)"       = "#FB8C00",
  "Inhibitory neuron cluster 4 (SST)"       = "#FFA726",
  "Inhibitory neuron cluster 5 (SST, RELN)" = "#FDD835",
  "Inhibitory neuron cluster 6 (CALB2)"     = "#FFFB99")
pal_meta <- tryCatch(metadata(spe)$color_vectors$celltype_clusters, error = function(e) NULL)
if (!is.null(pal_meta)) stopifnot(identical(toupper(unname(pal_meta[neur])), toupper(unname(cluster_cols[neur]))))
stopifnot(setequal(names(cluster_cols), neur))

mk <- function(d, ylab, hline = 0) {
  d$cluster <- factor(d$cluster, levels = names(cluster_cols))
  g <- ggplot(d, aes(pct_change, y)) +
    geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey60")
  if (!is.null(hline)) g <- g + geom_hline(yintercept = hline, linewidth = 0.25, linetype = "dashed", colour = "grey60")
  g +
    # outlined crosshairs, square caps so the outline wraps the bar ends (see distance script)
    geom_segment(aes(x = pct_change, xend = pct_change, y = y_lo, yend = y_hi),
                 linewidth = 0.75, colour = "grey25", lineend = "square") +
    geom_segment(aes(x = pct_lo, xend = pct_hi, y = y, yend = y),
                 linewidth = 0.75, colour = "grey25", lineend = "square") +
    geom_segment(aes(x = pct_change, xend = pct_change, y = y_lo, yend = y_hi, colour = cluster),
                 linewidth = 0.45, lineend = "square") +
    geom_segment(aes(x = pct_lo, xend = pct_hi, y = y, yend = y, colour = cluster),
                 linewidth = 0.45, lineend = "square") +
    geom_point(aes(fill = cluster), shape = 21, size = 2, stroke = 0.3, colour = "grey15") +
    scale_fill_manual(values = cluster_cols, name = NULL, drop = FALSE,
                      labels = function(x) sub("^Inhibitory neuron cluster ", "Inh ", sub("^Excitatory neuron cluster ", "Exc ", x))) +
    scale_colour_manual(values = cluster_cols, guide = "none", drop = FALSE) +
    labs(x = "Difference in columnar density,\nAD vs control (%)", y = ylab) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7),
          plot.margin = margin(4, 4, 4, 4), legend.position = "right",
          legend.key.size = grid::unit(8, "pt"), legend.text = element_text(size = 6))
}
panels <- if (IS_INT) list(
  list(panel = "phf1_z",                name = "raw",   analysis = "headline", hline = 0,
       ylab = "PHF1 intensity (within-donor z)"),
  list(panel = "phf1_z_negative_cells", name = "neg",   analysis = "neg", hline = 0,
       ylab = "PHF1 intensity, PHF1- cells (z)"),
  list(panel = "phf1_z_layer_relative", name = "layer", analysis = "layer", hline = 0,
       ylab = "PHF1 intensity vs same-layer cells (z)"),
  list(panel = "phf1_z_depth_relative", name = "depth", analysis = "depth", hline = 0,
       ylab = "PHF1 intensity vs same-depth cells (z)")) else if (IS_PLQ) list(
  list(panel = "pct_plaque",              name = "raw",   analysis = "headline", hline = NULL,
       ylab = expression(atop("Cells within 40 " * mu * "m", "of a plaque (%)"))),
  list(panel = "pct_plaque_layer_excess", name = "layer", analysis = "layer", hline = 0,
       ylab = expression(atop("Plaque-adjacent vs", "same-layer cells (% points)"))),
  list(panel = "pct_plaque_depth_excess", name = "depth", analysis = "depth", hline = 0,
       ylab = expression(atop("Plaque-adjacent vs", "same-depth cells (% points)")))) else list(
  list(panel = "pct_phf1_pos",              name = "raw",   analysis = "headline", hline = NULL,
       ylab = "PHF1+ cells (% of cluster)"),
  list(panel = "pct_phf1_pos_layer_excess", name = "layer", analysis = "layer", hline = 0,
       ylab = "PHF1+ vs same-layer cells (% points)"),
  list(panel = "pct_phf1_pos_depth_excess", name = "depth", analysis = "depth", hline = 0,
       ylab = "PHF1+ vs same-depth cells (% points)"))
for (pn in panels) {
  nm <- paste0(stem, "_", pn$name)
  d  <- filter(pts, panel == pn$panel)
  ggsave(file.path(out_dir, paste0("plot_", nm, ".pdf")), mk(d, pn$ylab, pn$hline),
         width = 9.5, height = 6, units = "cm", device = "pdf")
  wr(d, paste0("source_data_", nm, ".tsv"))
  r <- sens[sens$analysis == pn$analysis, ]
  writeLines(c(
    sprintf("plot_%s.pdf -- %s", nm, r$description),
    "Spearman rho across the 10 neuronal clusters; donor-level stratified bootstrap (B = 9999).",
    sprintf("rho = %.3f, 95%% CI [%.3f, %.3f]; p(boot, CI inversion) = %s; p(perm, cluster labels, %d draws) = %s",
            r$rho, r$ci_lo, r$ci_hi, signif(r$p_boot, 3), r$n_perm, signif(r$p_perm, 3)),
    sprintf("rho < 0 = the most-lost clusters %s.", SIGN_TXT),
    "Points = cluster estimate; crosshairs = 95% donor-bootstrap intervals on both axes.",
    sprintf("Full log with design, gates, sensitivities and sessionInfo: stats_%s.txt", stem)),
    file.path(out_dir, paste0("stats_", nm, ".txt")))
}

# ===================================================================
# Stats log
# ===================================================================
sink(file.path(out_dir, paste0("stats_", stem, ".txt")))
cat(if (IS_INT) "IMC (EC): cluster % loss in AD vs cluster PHF1 intensity in AD\n" else
    if (IS_PLQ) "IMC (EC): cluster % loss in AD vs % of the cluster's cells within 40 um of an amyloid plaque in AD\n" else
    "IMC (EC): cluster % loss in AD vs % of the cluster's cells that are PHF1+ in AD\n")
cat("Script: R/imc_cluster_loss_vs_phf1_intensity.R --outcome=", OUTCOME, " | run ", format(Sys.time()), "\n", sep = "")
cat("Source: ", localwd, "spe.rds\n\n", sep = "")
cat("=== Design ===\n")
cat("Unit of inference: cluster (k = 10 neuronal: Excitatory 1-4 + Inhibitory 1-6; glia excluded).\n")
cat("x = 100*(mean AD donor density / mean control donor density - 1), columnar density = cells/mm\n")
cat("    ROI width, per ROI then per donor (IMC study Fig. 2a, 12b-Cluster_Percentage_Change.Rmd).\n")
if (IS_INT) {
  cat("y = mean over AD donors of the donor-mean PHF1 z, where PHF1 = assay exprs = asinh(counts/1)\n")
  cat("    Z-SCORED WITHIN DONOR over all cells of the 10 neuronal clusters (PHF1+ and PHF1- alike),\n")
  cat("    as in imc_phf1_glia_distance.R. Donor x cluster needs >=", MIN_CELLS, "cells.\n")
} else if (IS_PLQ) {
  cat("y = mean over AD donors of the donor's % of the cluster's cells with Matched_4G8_40 = TRUE: nucleus\n")
  cat("    within the 40 um-dilated Fiji 4G8+ plaque mask (matched in Steinbock), the Set 3 plaque\n")
  cat("    sensitivity term. Donor x cluster needs >=", MIN_CELLS, "cells. layer/depth rows: percentage points above\n")
  cat("    the plaque-adjacent rate of all neuronal cells at the same layer / relative depth of the same ROI.\n")
} else {
  cat("y = mean over AD donors of the donor's % of the cluster's cells that are PHF1_Otsu+ (per-donor\n")
  cat("    Otsu on asinh PHF1, floored at 0.88 -- the call that defines the tangle anchors elsewhere).\n")
  cat("    Donor x cluster needs >=", MIN_CELLS, "cells. layer/depth rows: percentage points above the PHF1+\n")
  cat("    rate of all neuronal cells at the same layer / relative depth of the same ROI.\n")
}
cat("Test: Spearman rho across clusters; donor-level STRATIFIED bootstrap (control and AD donors\n")
cat("  resampled separately, B =", B, ") recomputing x and y from the same resampled donors.\n")
cat("  Effect size = rho with percentile 95% CI. p(boot) = CI-inversion; p(perm) = cluster-label\n")
cat("  permutation of y (", B_PERM, " draws; exact enumeration for k <= 8), which treats the point\n", sep = "")
cat("  estimates as fixed and ignores their uncertainty.\n")
cat("Power: at k = 10, two-sided p < 0.05 needs |rho| >~ 0.65. At k = 4 it is unreachable.\n")
cat("SIGN: rho < 0 = the most-lost clusters (most negative x) ", SIGN_TXT, ".\n\n", sep = "")

cat("=== Cohort ===\n")
cat(sprintf("Donors: %d control (Braak 0-1), %d AD (Braak 3-6); 3-ROI filter removed %d donors.\n",
            length(ctrl_ids), length(ad_ids), length(cv_removal)))
cat(sprintf("PHF1 frame: %d AD donors, %d ROIs, %d neuronal cells (%d PHF1_Otsu+).\n",
            g2["donors"], g2["rois"], g2["cells"], g2["phf1_pos"]))
cat(sprintf("GATE 1 (IMC study Fig. 2a reproduced): CALB1 %.1f +/- %.1f, PVALB %.1f +/- %.1f, RORB %.1f +/- %.1f;\n",
            x_hat[names(exp_pub)[1]], sem_pub[names(exp_pub)[1]], x_hat[names(exp_pub)[2]],
            sem_pub[names(exp_pub)[2]], x_hat[names(exp_pub)[3]], sem_pub[names(exp_pub)[3]]))
cat(sprintf("  max |diff| vs SupplementaryTable4_patient_percent_change.csv = %s. PASS\n", format(gate1_csv, digits = 3)))
if (IS_INT) {
  cat(sprintf("GATE 2 (z = mean 0, SD 1 within every donor; PHF1_Otsu+ brighter than PHF1- in all %d donors\n",
              sum(pos_vs_neg$n_pos > 0)))
  cat(sprintf("  with any; median difference %.2f z). PASS\n", median(pos_vs_neg$dz[pos_vs_neg$n_pos > 0])))
} else cat("GATE 2 (cohort: 44 AD donors, 132 ROIs). PASS\n")
cat(sprintf("GATE 3a (shuffled y per replicate): mean rho = %.4f (should be ~0). PASS\n", mean(null_rho, na.rm = TRUE)))
cat(sprintf("GATE 3b (planted monotone relation, 0.5 within-donor cell SD (%.2f) per SD of x): 95%% CI [%.2f, %.2f] excludes 0. PASS\n",
            plant_unit * 0.5, plant_ci[1], plant_ci[2]))
cat(sprintf("Clusters with < 30 AD donors contributing a y: %s\n\n",
            if (length(low_cov)) paste(short_lab(low_cov), collapse = ", ") else "none"))

for (v in names(panel_of)) {
  cat(sprintf("=== Per-cluster values: %s (95%% donor-bootstrap intervals) ===\n", panel_of[[v]]))
  print(as.data.frame(pts %>% filter(panel == panel_of[[v]]) %>%
    transmute(label, class, pct = round(pct_change, 1), pct_lo = round(pct_lo, 1), pct_hi = round(pct_hi, 1),
              y = round(y, 3), y_lo = round(y_lo, 3), y_hi = round(y_hi, 3),
              n_ad = n_ad_donors_phf1, phf1_frac = round(phf1_pos_frac, 3))), row.names = FALSE)
  cat("\n")
}

cat("=== HEADLINE and pre-specified sensitivities (effect size = Spearman rho, 95% CI) ===\n")
print(as.data.frame(sens %>% mutate(across(c(rho, ci_lo, ci_hi, frac_boot_neg), ~ round(.x, 3)),
                                    across(c(p_boot, p_perm), ~ signif(.x, 3)))), row.names = FALSE)
cat("\nSIGN: rho < 0 = the most-lost clusters ", SIGN_TXT, ".\n", sep = "")
cat(sprintf("depth_df* rows: the depth adjustment refitted at df = %s.\n",
            paste(DEPTH_DF_SWEEP, collapse = ", ")))
cat("exc_only rows: k = 4, so p cannot fall below 2/24 = 0.083. Read rho and its CI as a direction\n")
cat("check; an n.s. p there is NOT evidence of no association.\n")

cat("\n=== Leave-one-cluster-out (headline) ===\n")
print(as.data.frame(loco %>% mutate(across(c(rho, ci_lo, ci_hi), ~ round(.x, 3)))), row.names = FALSE)

cat("\n=== Notes ===\n")
cat(" - Ecological association across 10 clusters: it says which POPULATIONS are both lost and\n")
cat(if (IS_PLQ) "   plaque-adjacent, not that plaques cause loss in any cell.\n" else
    "   PHF1-rich, not that PHF1 causes loss in any cell.\n")
if (IS_PLQ) {
  cat(" - Plaque adjacency is a property of WHERE a cell sits, not of the cell, so it carries no self-anchor\n")
  cat("   coupling (plaques are not cells of the cluster). Under death independent of plaques the % is\n")
  cat("   unchanged by loss; preferential death NEAR plaques lowers it (survivors are the far cells).\n")
  cat(" - Amyloid is laminar (and ROI coverage of the 4G8 mask varies); read the layer/depth rows.\n")
  cat(" - 4G8 detects APP/Abeta 17-24 and labels diffuse as well as dense-core plaques; the 40 um niche\n")
  cat("   is a fixed dilation and was not tuned here.\n")
} else if (IS_INT) {
  cat(" - A cluster mean of z combines the tangle-bearing fraction (phf1_frac) with the brightness of\n")
  cat("   PHF1- cells; the `neg` row isolates the latter.\n")
} else {
  cat(" - The denominator of % PHF1+ is the surviving cells. Preferential\n")
  cat("   death of tangle-FREE cells raises it with loss; death of tangle-BEARING cells lowers it. The\n")
  cat("   correlation therefore cannot separate 'tangles cause loss' from 'loss concentrates tangles'.\n")
}
if (!IS_PLQ) {
  cat(" - Survivorship: y is measured on SURVIVING cells. Tangle-bearing neurons that died (ghost\n")
  cat("   tangles excluded by segmentation) are not counted.\n")
}
if (IS_INT) {
  cat(" - The z-score is within donor over the neuronal population, so y is a cluster's position in its\n")
  cat("   donor's neuronal PHF1 distribution; it carries no between-donor burden information.\n")
}
cat(" - x is age-unadjusted, as in IMC study Fig. 2a; see the 'age' row.\n")
cat("\n")
print(sessionInfo())
sink()
cat("Done ->", out_dir, "\n")
