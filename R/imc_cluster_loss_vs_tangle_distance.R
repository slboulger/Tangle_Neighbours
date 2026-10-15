#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# imc_cluster_loss_vs_tangle_distance.R
#
# Figure panels: 4A
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# imc_cluster_loss_vs_tangle_distance.R
#
# Are the EC neuronal clusters that are LOST most in AD the ones that sit CLOSEST to
# tangle-bearing neurons in AD tissue?
#
#   x = % change in columnar density, AD vs control, per cluster -- the cluster percentage-change
#       bar of Fig. 2a of the previously published IMC study of this cohort (see README),
#       computed as in its analysis code (IMCDataAnalysis-main/12b-Cluster_Percentage_Change.Rmd:39-103)
#   y = geometric-mean distance to the nearest PHF1+ neuron, per cluster, AD donors only
#
# ---------------------------------------------------------------------------
# THE TEST: donor-level stratified bootstrap of Spearman rho across clusters
# ---------------------------------------------------------------------------
# The unit of inference is the CLUSTER (k = 10 neuronal clusters). The uncertainty in each
# cluster's x and y is DONOR sampling, and it is propagated by resampling donors -- control and AD
# separately -- and recomputing BOTH axes from the same resampled donors in every replicate.
# Why not Deming / York / BCES:
#   * x and y errors are correlated: the same 44 AD donors feed both axes.
#   * the clusters are not independent points: they share donors, ROIs and the PHF1+ anchors.
#     Resampling whole donors keeps that cross-cluster correlation; per-point SEs cannot.
#   * x's uncertainty must include the control mean (its denominator), which a per-cluster SEM
#     over AD donors ignores. Resampling control donors carries it.
#   * York/Deming treat all scatter as measurement error (no intrinsic scatter).
# Spearman because CALB1 (-50%) is an extreme point on x and nothing says x and y are linear.
# Reported: rho + percentile 95% CI (the effect size), CI-inversion bootstrap p, and a cluster-label
# permutation p as the fixed-point reference (ignores measurement error; answers "is this ordering
# unusual among these 10 clusters").
#
# Scope: headline = 10 neuronal clusters (Exc 1-4 + Inh 1-6), glia excluded. Sensitivity = 4
# excitatory clusters only; with k = 4 the minimum two-sided permutation p is 2/24, so those rows
# are read as a direction check.
#
# Pre-specified sensitivities (each is a row of stats_*_sensitivity.tsv):
#   layer   -- y = cell log d minus the mean log d of all PHF1- cells in the same ROI x layer.
#              Removes laminar position AND donor/ROI tangle burden. Tangles are layer-II enriched,
#              CALB1 is a layer-II population and ROIs are 478 x 2525 um strips, so this row asks
#              whether the association is carried by laminar position.
#   depth   -- y = cell log d minus its ROI's expected log d at that relative cortical depth, from a
#              per-ROI natural spline (df = DEPTH_DF) of log d on scaled_Y over all PHF1- cells.
#              scaled_Y = 100*(y - min y)/(max y - min y) over ALL cells in the ROI, 0 = pia, as in
#              imc_covariate_sensitivity.R. A continuous version of `layer`, which cannot separate
#              layer II from III because the stored layer factor merges them as 2/3.
#   other   -- y from anchors in OTHER clusters only. Same-type cells aggregate, so a cluster that
#              carries many tangles itself sits near PHF1+ cells trivially.
#   cap300  -- cells > 300 um dropped, the IMC Set 3 gradient-model cap (docs/MODELS.md). The
#              headline is UNCAPPED because a cap censors a mean distance.
#   age     -- x from an AD + Age fit (OLS on donor means == the ROI-level lmer fixed effect for a
#              donor-level predictor under 3-ROI balance) / control EMM at the cohort mean age.
#   zerofill-- ROI x cluster combinations with no cells counted as 0 (the Rmd drops them, which is
#              why CALB1 has 43 AD donors, not 44).
#   exc_only-- 4 excitatory clusters.
#
# Distance construction is the one used by the IMC distance analyses
# (imc_phf1_distance_markers_subsets.R): 3-ROI filter, Artefact/Unassigned dropped,
# AD = BraakGroup != Braak_0_1, anchors = PHF1_pos in any neuron cluster, targets PHF1_neg,
# nearest anchor WITHIN sample_id (Pos_X/Pos_Y are per-ROI local).
# Natural log, never log1p. Geometric means reported back in um.
#
# Input: <IMC_ROOT>/spe.rds. Reproduction check (GATE 1): if
# <IMC_ROOT>/Figures/SupplementaryTable4_patient_percent_change.csv -- Supplementary Table 4 of the
# IMC study, the per-donor % change per cluster -- is present, x must match it.
#
# Run: Rscript R/imc_cluster_loss_vs_tangle_distance.R
suppressPackageStartupMessages({
  library(SpatialExperiment); library(dplyr); library(tibble); library(tidyr)
  library(ggplot2); library(RANN); library(splines)
})

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")
source("R/imc_nn3_utils.R")   # nearest_dist_imc()

out_dir <- "plots/imc_cluster_loss_vs_tangle_distance"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stem <- "cluster_loss_vs_tangle_distance"
wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE)

B           <- 9999     # donor bootstrap replicates
B_PERM      <- 199999   # cluster-label permutations (k = 10; Monte Carlo, SE ~0.0005 at p = 0.05)
MAX_DIST    <- 300      # cap300 sensitivity only
DEPTH_DF    <- 6        # spline df for the per-ROI depth expectation of log d
DEPTH_DF_SWEEP <- c(3, 10, 15)   # df sensitivity for the depth adjustment
MIN_CELLS   <- 5        # PHF1- cells a donor x cluster needs to contribute a distance
SEED        <- 1
EXCLUDE     <- c("Artefact cluster", "Unassigned cluster")
set.seed(SEED)

localwd <- "<IMC_ROOT>/"
ref_pct <- paste0(localwd, "Figures/SupplementaryTable4_patient_percent_change.csv")
spe <- readRDS(paste0(localwd, "spe.rds"))

cd <- as.data.frame(colData(spe))
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
# ===================================================================
sample_meta <- cd %>% select(sample_id, patient_id, Age, Braak, BraakGroup, width_px) %>% distinct()
stopifnot(!anyDuplicated(sample_meta$sample_id))

roi_dens <- cd %>% count(sample_id, celltype_clusters, name = "n_cells") %>%
  left_join(sample_meta, by = "sample_id") %>%
  mutate(counts_per_mm = n_cells / width_px * 1e3) %>%
  filter(!patient_id %in% cv_removal, !celltype_clusters %in% EXCLUDE)

# zerofill sensitivity: every kept ROI x every cluster, absent = 0
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
stopifnot(abs(x_hat - colMeans(per_donor_pct, na.rm = TRUE)) < 1e-9)   # ratio of means == mean of %
gate1 <- tibble(cluster = clusters_all, pct = x_hat, sem = sem_pub)
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
# y: distance to nearest PHF1+ neuron, AD donors, within ROI
# ===================================================================
base <- cd %>%
  filter(!celltype_clusters %in% EXCLUDE, !patient_id %in% cv_removal,
         BraakGroup != "Braak_0_1", !is.na(PHF1_Otsu))
anchors <- base %>% filter(PHF1_Otsu == "PHF1_pos", grepl("neuron", celltype_clusters, ignore.case = TRUE))
negs    <- base %>% filter(PHF1_Otsu == "PHF1_neg")
negs$d  <- nearest_dist_imc(negs$x, negs$y, negs$sample_id, anchors$x, anchors$y, anchors$sample_id)

# other-cluster anchors: a cell's nearest PHF1+ neuron NOT of its own cluster
negs$d_other <- NA_real_
for (cl in neur) {
  i <- which(negs$celltype_clusters == cl); a <- anchors[anchors$celltype_clusters != cl, ]
  negs$d_other[i] <- nearest_dist_imc(negs$x[i], negs$y[i], negs$sample_id[i], a$x, a$y, a$sample_id)
}

# ---- GATE 2: expected cohort/anchor/target counts of the IMC distance cohort ----
exc_t <- negs %>% filter(grepl("^Excitatory", celltype_clusters))
# ROIs are counted WITH an anchor (131 of 132; BBN_9110_4 has none); exc_targets counts ALL PHF1- Exc
# targets including that ROI's 64 cells with no distance (34563); exc_le300 is the <= 300 um count.
g2 <- c(anchors = nrow(anchors), donors = n_distinct(base$patient_id),
        rois = n_distinct(exc_t$sample_id[!is.na(exc_t$d)]),
        rois_no_anchor = length(setdiff(unique(exc_t$sample_id), unique(anchors$sample_id))),
        exc_targets = nrow(exc_t), exc_le300 = sum(exc_t$d <= MAX_DIST, na.rm = TRUE))
g2_exp <- c(anchors = 2467, donors = 44, rois = 131, rois_no_anchor = 1, exc_targets = 34563,
            exc_le300 = 26653)
if (!all(unname(g2) == unname(g2_exp))) {
  print(rbind(observed = g2, expected = g2_exp)); stop("GATE 2 failed: cohort differs from the IMC distance cohort")
}
stopifnot(all(negs$d > 0, na.rm = TRUE), all(negs$d_other > 0, na.rm = TRUE))   # log(), never log1p

negs <- negs %>% mutate(ld = log(d), ld_other = log(d_other))
layer_ref <- negs %>% filter(!is.na(ld), !is.na(layer)) %>%
  group_by(sample_id, layer) %>% summarise(ld_ref = mean(ld), n_ref = n(), .groups = "drop")
negs <- negs %>% left_join(layer_ref, by = c("sample_id", "layer")) %>% mutate(ld_layer = ld - ld_ref)

# depth: residual of log d from a per-ROI smooth of log d on relative depth (all PHF1- cells)
negs$ld_depth <- NA_real_
for (s_ in unique(negs$sample_id)) {
  i <- which(negs$sample_id == s_ & !is.na(negs$ld))
  if (length(i) < 10 * DEPTH_DF) next
  negs$ld_depth[i] <- residuals(lm(negs$ld[i] ~ ns(negs$scaled_Y[i], df = DEPTH_DF)))
}
n_roi_depth <- n_distinct(negs$sample_id[!is.na(negs$ld_depth)])
for (k in DEPTH_DF_SWEEP) {
  col <- paste0("ld_depth_df", k); negs[[col]] <- NA_real_
  for (s_ in unique(negs$sample_id)) {
    i <- which(negs$sample_id == s_ & !is.na(negs$ld))
    if (length(i) < 10 * k) next
    negs[[col]][i] <- residuals(lm(negs$ld[i] ~ ns(negs$scaled_Y[i], df = k)))
  }
}
sweep_nm <- paste0("depth_df", DEPTH_DF_SWEEP)

tg <- negs %>% filter(celltype_clusters %in% neur, !is.na(ld))
donor_y <- tg %>% group_by(patient_id, celltype_clusters) %>%
  summarise(n_cells     = n(),
            raw         = mean(ld),
            layer       = mean(ld_layer, na.rm = TRUE),
            depth       = mean(ld_depth, na.rm = TRUE),
            other       = if (all(is.na(ld_other))) NA_real_ else mean(ld_other, na.rm = TRUE),
            cap300      = if (any(d <= MAX_DIST)) mean(ld[d <= MAX_DIST]) else NA_real_,
            n_cap300    = sum(d <= MAX_DIST),
            .groups = "drop") %>%
  mutate(across(c(raw, layer, depth, other), ~ ifelse(n_cells >= MIN_CELLS, .x, NA_real_)),
         cap300 = ifelse(n_cap300 >= MIN_CELLS, cap300, NA_real_))
donor_y <- donor_y %>% left_join(
  tg %>% group_by(patient_id, celltype_clusters) %>%
    summarise(across(all_of(paste0("ld_", sweep_nm)), ~ mean(.x, na.rm = TRUE), .names = "{sub('^ld_', '', .col)}"),
              .groups = "drop"),
  by = c("patient_id", "celltype_clusters")) %>%
  mutate(across(all_of(sweep_nm), ~ ifelse(n_cells >= MIN_CELLS, .x, NA_real_)))
Y <- lapply(setNames(nm = c("raw", "layer", "depth", "other", "cap300", sweep_nm)),
            function(v) to_mat(donor_y, ad_ids, v))
n_ad_dist <- colSums(!is.na(Y$raw))[neur]
low_cov <- names(n_ad_dist)[n_ad_dist < 30]

# PHF1+ fraction per cluster (descriptive; donor-averaged)
phf1_frac <- base %>% filter(celltype_clusters %in% neur) %>%
  group_by(patient_id, celltype_clusters) %>%
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
  layer    = spear(s$x[neur],  s$y$layer[neur]),
  depth    = spear(s$x[neur],  s$y$depth[neur]),
  other    = spear(s$x[neur],  s$y$other[neur]),
  cap300   = spear(s$x[neur],  s$y$cap300[neur]),
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
y_for <- c(headline = "raw", layer = "layer", depth = "depth", other = "other", cap300 = "cap300", age = "raw",
           zerofill = "raw", exc_only = "raw", exc_only_layer = "layer",
           exc_only_depth = "depth", setNames(sweep_nm, sweep_nm))
x_for <- c(headline = "x", layer = "x", depth = "x", other = "x", cap300 = "x", age = "xa", zerofill = "xz",
           exc_only = "x", exc_only_layer = "x", exc_only_depth = "x",
           setNames(rep("x", length(sweep_nm)), sweep_nm))
analysis_desc <- c(
  headline = "10 neuronal clusters; y = raw distance (uncapped)",
  layer    = "y = log d centred on all PHF1- cells in the same ROI x layer",
  depth    = sprintf("y = log d minus per-ROI spline (df %d) of log d on relative cortical depth", DEPTH_DF),
  other    = "y = distance to nearest PHF1+ neuron of a DIFFERENT cluster",
  cap300   = "y from cells <= 300 um only (Set 3 cap)",
  age      = "x = AD + Age fit on donor means, % of control EMM at mean age",
  zerofill = "x with absent ROI x cluster counted as 0",
  exc_only = "4 excitatory clusters; y raw (CANNOT reach p<0.05: min p = 2/24)",
  exc_only_layer = "4 excitatory clusters; y layer-centred (direction check only)",
  exc_only_depth = "4 excitatory clusters; y depth-adjusted (direction check only)",
  setNames(sprintf("as depth, spline df = %d (df sensitivity)", DEPTH_DF_SWEEP), sweep_nm))

sens <- bind_rows(lapply(names(analysis_desc), function(a) {
  cls <- if (grepl("^exc_only", a)) exc else neur
  pp  <- perm_p(obs[[x_for[a]]][cls], if (y_for[a] == "raw") obs$y$raw[cls] else obs$y[[y_for[a]]][cls], B_PERM)
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
# Synthetic donor-level log distances with the observed per-cluster donor SD, whose cluster means
# are a monotone (decreasing-loss -> nearer) function of the observed x. CI must exclude 0.
sd_y <- apply(Y$raw[, neur], 2, sd, na.rm = TRUE)
mu_plant <- mean(obs$y$raw[neur]) + 0.5 * as.numeric(scale(obs$x[neur]))   # 0.5 log-units per SD of x
Yp <- Y$raw; Yp[, neur] <- sapply(seq_along(neur), function(j)
  ifelse(is.na(Y$raw[, neur[j]]), NA, rnorm(nrow(Yp), mu_plant[j], sd_y[j])))
plant_r <- replicate(1999, {
  sc_ <- sample(ic, n_c, TRUE); sa_ <- sample(ia, n_a, TRUE)
  spear(pct_change(D_ctrl[sc_, ], D_ad[sa_, ])[neur], colMeans(Yp[sa_, neur], na.rm = TRUE))
})
plant_ci <- quantile(plant_r, c(0.025, 0.975), na.rm = TRUE)
if (plant_ci[1] <= 0) stop("GATE 3b failed: a planted monotone relation was not recovered")
if (abs(mean(null_rho, na.rm = TRUE)) > 0.02) stop("GATE 3a failed: shuffled-y rho not centred on 0")

# ===================================================================
# Source data (per cluster) + per-donor audit table
# ===================================================================
qlo <- function(m) apply(m, 2, quantile, 0.025, na.rm = TRUE)
qhi <- function(m) apply(m, 2, quantile, 0.975, na.rm = TRUE)
panel_of <- c(raw = "distance_um", layer = "layer_relative", depth = "depth_relative")
pts <- bind_rows(lapply(names(panel_of), function(v) {
  tibble(panel = panel_of[[v]],
         cluster = neur, label = short_lab(neur), class = cl_class(neur),
         pct_change = obs$x[neur], pct_lo = qlo(boot_x), pct_hi = qhi(boot_x),
         pct_sem_published = sem_pub[neur],
         # raw: geometric mean distance (um); layer/depth: ratio to same-ROI PHF1- cells at that layer/depth
         y = exp(obs$y[[v]][neur]), y_lo = exp(qlo(boot_y[[v]])), y_hi = exp(qhi(boot_y[[v]])),
         y_log_mean = obs$y[[v]][neur],
         n_ctrl_donors = n_ctrl_dens[neur], n_ad_donors_density = n_ad_dens[neur],
         n_ad_donors_distance = colSums(!is.na(Y[[v]]))[neur],
         n_phf1neg_cells = as.integer(table(factor(tg$celltype_clusters, levels = neur))))
})) %>% left_join(phf1_frac, by = c("cluster" = "celltype_clusters"))

audit <- full_join(
  bind_rows(dm %>% mutate(group = ifelse(patient_id %in% ad_ids, "AD", "Control"))),
  donor_y %>% rename(log_dist_raw = raw, log_dist_layer = layer, log_dist_depth = depth, log_dist_other = other,
                     log_dist_cap300 = cap300, n_phf1neg_cells = n_cells),
  by = c("patient_id", "celltype_clusters")) %>%
  filter(celltype_clusters %in% neur) %>%
  rename(columnar_density_per_mm = dens) %>% arrange(celltype_clusters, group, patient_id)
wr(audit, paste0("source_data_", stem, "_donor.tsv"))
wr(sens,  paste0("stats_", stem, "_sensitivity.tsv"))
wr(loco,  paste0("stats_", stem, "_loco.tsv"))

# ===================================================================
# Plot
# ===================================================================
# IMC study cluster colours: IMCDataAnalysis-main/19-ColoursForPublication.Rmd:46-69
# (metadata(spe)$color_vectors$celltype_clusters). Points and crosshairs take the cluster colour with
# a thin dark outline, because the Inh 5/6 yellows are illegible on white otherwise.
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

mk <- function(d, ylab, hline = NULL) {
  d$cluster <- factor(d$cluster, levels = names(cluster_cols))
  g <- ggplot(d, aes(pct_change, y)) +
    geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey60")
  if (!is.null(hline)) g <- g + geom_hline(yintercept = hline, linewidth = 0.25, linetype = "dashed", colour = "grey60")
  # outlined crosshairs: a wider dark segment underneath, the cluster colour on top. Both use SQUARE
  # caps, which extend each end by half the line width, so the outline also wraps the bar ends with
  # the same thickness as along the sides ((0.75 - 0.45) / 2 both ways). All outlines are drawn before
  # any colour so an outline never cuts across another cluster's bar.
  g + geom_segment(aes(x = pct_change, xend = pct_change, y = y_lo, yend = y_hi),
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
    scale_y_log10() +
    labs(x = "Difference in columnar density,\nAD vs control (%)", y = ylab) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7),
          plot.margin = margin(4, 4, 4, 4), legend.position = "right",
          legend.key.size = grid::unit(8, "pt"), legend.text = element_text(size = 6))
}
panels <- list(
  list(panel = "distance_um",    name = "raw",   analysis = "headline", hline = NULL,
       ylab = expression("Distance to nearest PHF1+ neuron (" * mu * "m)")),
  list(panel = "layer_relative", name = "layer", analysis = "layer", hline = 1,
       ylab = "Distance relative to same-layer cells"),
  list(panel = "depth_relative", name = "depth", analysis = "depth", hline = 1,
       ylab = "Distance relative to same-depth cells"))
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
    "Points = cluster estimate; crosshairs = 95% donor-bootstrap intervals on both axes.",
    sprintf("Full log with design, gates, sensitivities and sessionInfo: stats_%s.txt", stem)),
    file.path(out_dir, paste0("stats_", nm, ".txt")))
}

# ===================================================================
# Stats log
# ===================================================================
sink(file.path(out_dir, paste0("stats_", stem, ".txt")))
cat("IMC (EC): cluster % loss in AD vs cluster distance to nearest PHF1+ neuron in AD\n")
cat("Script: R/imc_cluster_loss_vs_tangle_distance.R | run", format(Sys.time()), "\n")
cat("Source: ", localwd, "spe.rds\n\n", sep = "")
cat("=== Design ===\n")
cat("Unit of inference: cluster (k = 10 neuronal: Excitatory 1-4 + Inhibitory 1-6; glia excluded).\n")
cat("x = 100*(mean AD donor density / mean control donor density - 1), columnar density = cells/mm\n")
cat("    ROI width, per ROI then per donor (IMC study Fig. 2a, 12b-Cluster_Percentage_Change.Rmd).\n")
cat("y = exp(mean over AD donors of the donor-mean log distance) to the nearest PHF1+ neuron,\n")
cat("    PHF1- cells of the cluster, within ROI, uncapped. Donor x cluster needs >=", MIN_CELLS, "cells.\n")
cat("Test: Spearman rho across clusters; donor-level STRATIFIED bootstrap (control and AD donors\n")
cat("  resampled separately, B =", B, ") recomputing x and y from the same resampled donors, so the\n")
cat("  uncertainty in both axes, their correlation (shared AD donors) and the control-mean\n")
cat("  denominator are all propagated. Effect size = rho with percentile 95% CI. p(boot) = CI-inversion;\n")
cat("  p(perm) = cluster-label permutation of y (", B_PERM, " draws; exact enumeration for k <= 8),\n", sep = "")
cat("  which treats the point estimates as fixed and ignores their uncertainty.\n")
cat(sprintf("Power: at k = 10, two-sided p < 0.05 needs |rho| >~ 0.65. At k = 4 it is unreachable.\n\n"))

cat("=== Cohort ===\n")
cat(sprintf("Donors: %d control (Braak 0-1), %d AD (Braak 3-6); 3-ROI filter removed %d donors.\n",
            length(ctrl_ids), length(ad_ids), length(cv_removal)))
cat(sprintf("Distance frame: %d anchors (PHF1+ neurons), %d ROIs, %d ROI(s) without an anchor.\n",
            g2["anchors"], g2["rois"], g2["rois_no_anchor"]))
cat(sprintf("GATE 1 (IMC study Fig. 2a reproduced): CALB1 %.1f +/- %.1f, PVALB %.1f +/- %.1f, RORB %.1f +/- %.1f;\n",
            x_hat[names(exp_pub)[1]], sem_pub[names(exp_pub)[1]], x_hat[names(exp_pub)[2]],
            sem_pub[names(exp_pub)[2]], x_hat[names(exp_pub)[3]], sem_pub[names(exp_pub)[3]]))
cat(sprintf("  max |diff| vs SupplementaryTable4_patient_percent_change.csv = %s. PASS\n", format(gate1_csv, digits = 3)))
cat("GATE 2 (cohort counts match the IMC distance cohort): PASS\n")
print(rbind(observed = g2, expected = g2_exp))
cat(sprintf("GATE 3a (shuffled y per replicate): mean rho = %.4f (should be ~0). PASS\n", mean(null_rho, na.rm = TRUE)))
cat(sprintf("GATE 3b (planted monotone relation, 0.5 log-units per SD of x): 95%% CI [%.2f, %.2f] excludes 0. PASS\n",
            plant_ci[1], plant_ci[2]))
cat(sprintf("Clusters with < 30 AD donors contributing a distance: %s\n\n",
            if (length(low_cov)) paste(short_lab(low_cov), collapse = ", ") else "none"))

cat("=== Per-cluster values (95% donor-bootstrap intervals) ===\n")
print(as.data.frame(pts %>% filter(panel == "distance_um") %>%
  transmute(label, class, pct = round(pct_change, 1), pct_lo = round(pct_lo, 1), pct_hi = round(pct_hi, 1),
            dist_um = round(y, 1), d_lo = round(y_lo, 1), d_hi = round(y_hi, 1),
            n_ctrl = n_ctrl_donors, n_ad = n_ad_donors_density, n_ad_dist = n_ad_donors_distance,
            phf1_frac = round(phf1_pos_frac, 3))), row.names = FALSE)
cat("\nLayer-relative y (ratio of the cluster's geometric-mean distance to that of all PHF1- cells in\n")
cat("the same ROI x layer; < 1 = nearer to tangles than its laminar neighbours):\n")
print(as.data.frame(pts %>% filter(panel == "layer_relative") %>%
  transmute(label, ratio = round(y, 3), lo = round(y_lo, 3), hi = round(y_hi, 3))), row.names = FALSE)

cat(sprintf("\nDepth-relative y (ratio to the ROI's spline-expected distance at that relative depth,\nper-ROI ns(scaled_Y, df = %d) over all PHF1- cells; %d ROIs):\n", DEPTH_DF, n_roi_depth))
print(as.data.frame(pts %>% filter(panel == "depth_relative") %>%
  transmute(label, ratio = round(y, 3), lo = round(y_lo, 3), hi = round(y_hi, 3))), row.names = FALSE)

cat("\n=== HEADLINE and pre-specified sensitivities (effect size = Spearman rho, 95% CI) ===\n")
print(as.data.frame(sens %>% mutate(across(c(rho, ci_lo, ci_hi, frac_boot_neg), ~ round(.x, 3)),
                                    across(c(p_boot, p_perm), ~ signif(.x, 3)))), row.names = FALSE)
cat("\nSign: rho > 0 = clusters with a more NEGATIVE % change (more loss) have SHORTER distances,\n")
cat("i.e. the most-lost clusters sit nearest to tangles. rho < 0 = the opposite.\n")
cat(sprintf("depth_df* rows: the depth adjustment refitted at df = %s.\n",
            paste(DEPTH_DF_SWEEP, collapse = ", ")))
cat("exc_only rows: k = 4, so p cannot fall below 2/24 = 0.083. Read rho and its CI as a direction\n")
cat("check; an n.s. p there is NOT evidence of no association.\n")

cat("\n=== Leave-one-cluster-out (headline) ===\n")
print(as.data.frame(loco %>% mutate(across(c(rho, ci_lo, ci_hi), ~ round(.x, 3)))), row.names = FALSE)

cat("\n=== Notes ===\n")
cat(" - Ecological association across 10 clusters: it says which POPULATIONS are both lost and\n")
cat("   tangle-proximal, not that proximity causes loss in any cell.\n")
cat(" - Survivorship: y is measured on SURVIVING PHF1- cells. If the cells nearest tangles died, a\n")
cat("   lost cluster's survivors are biased AWAY from tangles, which works against rho > 0.\n")
cat(" - Tangle-bearing cells are excluded from the targets but remain anchors; the 'other' row uses\n")
cat("   anchors of other clusters only.\n")
cat(" - The headline is uncapped (a cap censors a mean); the cap300 row applies the Set 3 300 um cap.\n")
cat(" - x is age-unadjusted, as in IMC study Fig. 2a; the 'age' row uses OLS on donor means, which equals the\n")
cat("   Rmd's ROI-level lmer fixed effect only under exact 3-ROI balance (approximately here).\n")
cat("\n")
print(sessionInfo())
sink()
cat("Done ->", out_dir, "\n")
