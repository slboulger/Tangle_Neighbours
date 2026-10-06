#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# phf1_intensity_gradient_neurons.R
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
# phf1_intensity_gradient_neurons.R
#
# Is sub-threshold PHF1 (ptau) signal GRADED in non-PHF1+ neurons around PHF1+
# neurons? Tests whether the per-cell PHF1 intensity of neurons that were NOT
# manually called PHF1+ decays with distance to the nearest PHF1+ (NFT) neuron --
# i.e. a spatial "halo" of tau, vs a purely cell-autonomous on/off state.
#
# DESIGN:
#   * Population = ALL neurons (celltype %in% neuron_order), manual PHF1 == FALSE,
#     within MAX_DIST um of a PHF1+ neuron.
#   * Outcome  = per-cell 95th-pct masked PHF1 intensity (phf1_intensity_p95),
#     robust-z normalised WITHIN each sample ((x - median)/MAD) -- confocal
#     exposure/gain differ per acquisition, so raw DN is only within-sample
#     comparable. This reproduces the phf1_intensity_p95_z column of
#     add_phf1_intensity.R, recomputed here so the script is self-contained.
#   * Models (mirroring plot_modulescore_vs_phf1_distance_modelp.R machinery):
#       pooled linear LMM (headline):
#         z_p95 ~ dist_scaled + celltype + Sex + Age_s + PMI_s + (1|sample_id)
#       pooled log-distance LMM: as above on log(distance), SD-scaled
#       pooled GAM: z_p95 ~ s(dist,k=5) + celltype + Sex + Age_s + PMI_s + s(sample_id,bs='re')
#       per-subtype linear LMM (sensitivity; drops the celltype term), BH across subtypes.
#     A negative, significant dist slope = graded halo. celltype fixed effect makes
#     the pooled slope composition-robust (within-subtype gradient).
#
# INFERENCE: donor random intercept (1|sample_id); n=9 donors (3 per Braak stage).
#
# OUTPUT (the standard triple), under OUTPUT_DIR:
#   plot_phf1_p95_gradient_neurons[_log|_gam].pdf   -- pooled fitted curve + 95% CI
#   plot_phf1_p95_gradient_neurons_rollmean.pdf     -- points + model-free rolling mean + 95% CI
#                                                      (line solid/thick if LOG model significant)
#   plot_phf1_p95_gradient_by_subtype.pdf           -- per-subtype fitted lines
#   plot_phf1_p95_gradient_by_subtype_rollmean.pdf  -- per-subtype 75um rolling means + 95% CI (no fit)
#   plot_phf1_p95_gradient_<FRONT_CT>_rollmean.pdf   -- single-cluster (Exc-IT-L2-3-CBLN2-HOPX) rolling mean
#   source_data_phf1_p95_gradient_neurons[...]_fit.tsv        -- exact drawn rows (fitted curves)
#   source_data_phf1_p95_gradient_neurons_rollmean.tsv        -- rolling-mean grid (mean, SEM, CI, n)
#   source_data_phf1_p95_gradient_by_subtype_rollmean.tsv     -- per-subtype rolling-mean grid
#   source_data_phf1_p95_gradient_neurons_persample.tsv       -- per-donor x ring medians
#   stats_phf1_p95_gradient_neurons.txt             -- formulae, counts, refs, sessionInfo
#   stats_phf1_p95_gradient_coeffs.tsv              -- structured coefficient table
#
# Run interactively (on HPC, where seu_PHF1.rds lives): edit the Settings block
# below (SEU_PATH, OUTPUT_DIR, MAX_DIST, ...), then source() or step through the file.

##  ............................................................................
##  Packages + setup                                                        ####
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(mgcv)
})

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")  # fig_theme, celltype_palette, neuron_order, ...

##  ............................................................................
##  Settings -- edit these, then run the script interactively (source/step)  ####
SEU_PATH          <- "seu_PHF1.rds"                                   # path to the Seurat object
OUTPUT_DIR        <- "plots/phf1_intensity_gradient_neurons_1000um"   # output dir (cap-specific)
MAX_DIST          <- 1000    # cap: model non-PHF1 neurons with dist_to_phf1_um <= this (um)
MIN_CELLS_SUBTYPE <- 100     # skip a per-subtype LMM below this many modelled cells
MIN_DONORS        <- 3       # skip a model below this many distinct donors
WINDOW_UM         <- 75      # rolling-mean window WIDTH (um) for the actual-data plot (matches module scripts)
POINT_MAX         <- 20000   # subsample this many cells for the scatter overlay (0 = no points)
SEED              <- 42      # provenance / reproducibility

FRAC_IN_BOUNDS_MIN <- 0.5    # mirror add_phf1_intensity.R: poorly-registered masks -> NA
SIG_ALPHA          <- 0.05
GAM_K              <- 5
P_FLOOR            <- .Machine$double.xmin   # clamp underflowing p (e.g. GAM p=0); record p_at_floor
RING_BREAKS        <- c(0, 50, 100, 200, 300, 500, 700, 1000)
HALF_WINDOW        <- WINDOW_UM / 2
set.seed(SEED)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

.ring_lower <- head(RING_BREAKS, -1)
.ring_upper <- tail(RING_BREAKS, -1)
RING_LAB    <- sprintf("%g-%g", .ring_lower, .ring_upper)
RING_MID    <- setNames((.ring_lower + .ring_upper) / 2, RING_LAB)

POOLED_COL   <- "#377EB8"   # single-series colour for the pooled fitted / rolling-mean curve
PHF1_REF_COL <- "#7F7F7F"   # grey for the PHF1+ neuron reference band
# Concise axis labels (x = distance to the nearest PHF1+ neuron; y = per-cell p95
# PHF1 intensity, within-sample robust z).
XLAB       <- expression("Distance to PHF1+ neuron (" * mu * "m)")
YLAB       <- "PHF1 intensity (z)"

##  ............................................................................
##  Helpers                                                                 ####

# First present column name from candidates (defensive against naming drift).
pick_col <- function(meta, candidates, what) {
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0)
    stop(sprintf("None of the %s columns (%s) found in seu meta.data",
                 what, paste(candidates, collapse = ", ")))
  hit[1]
}

# Robust z (median/MAD), identical to add_phf1_intensity.R::robust_z.
robust_z <- function(x) {
  m <- median(x, na.rm = TRUE); s <- mad(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) s <- sd(x, na.rm = TRUE)
  (x - m) / s
}

# Fit one LMM; return dist-slope stats + the model. coef_name is the distance term.
fit_lmer <- function(form, data, coef_name = "dist_scaled") {
  m <- tryCatch(lmerTest::lmer(form, data = data, REML = TRUE),
                error = function(e) { message("    lmer failed: ", conditionMessage(e)); NULL })
  if (is.null(m)) return(NULL)
  cf <- coef(summary(m))
  if (!coef_name %in% rownames(cf)) return(NULL)
  list(model = m,
       estimate = cf[coef_name, "Estimate"],
       se       = cf[coef_name, "Std. Error"],
       t        = cf[coef_name, "t value"],
       df       = cf[coef_name, "df"],
       p        = cf[coef_name, "Pr(>|t|)"],
       singular = isSingular(m))
}

# Fitted trend across a raw-um distance grid, Wald 95% CI from the fixed-effect
# vcov (RE variance ignored). Covariates held at means / reference factor levels.
# fe_form is the fixed-effect formula (must match the model's fixed effects).
predict_grid_lmer <- function(model, df, dist_sd, log_distance, x_cap, fe_form,
                              include_celltype) {
  lo      <- if (log_distance) min(df$dist_to_phf1_um) else 0
  grid_um <- seq(lo, x_cap, length.out = 200)
  dt      <- if (log_distance) log(grid_um) else grid_um
  nd <- data.frame(
    dist_scaled = dt / dist_sd,
    Sex         = factor(levels(df$Sex)[1], levels = levels(df$Sex)),
    Age_s       = mean(df$Age_s),
    PMI_s       = mean(df$PMI_s)
  )
  if (include_celltype)
    nd$celltype <- factor(levels(df$celltype)[1], levels = levels(df$celltype))
  X    <- model.matrix(fe_form, nd)
  beta <- lme4::fixef(model)
  X    <- X[, names(beta), drop = FALSE]
  fit  <- as.numeric(X %*% beta)
  V    <- as.matrix(vcov(model))
  se   <- sqrt(rowSums((X %*% V) * X))
  data.frame(dist_to_phf1_um = grid_um, fitted = fit,
             ci_lo = fit - 1.96 * se, ci_hi = fit + 1.96 * se)
}

# GAM smooth (mgcv::bam) on RAW distance + donor RE smooth. Returns edf/F/p + model.
.gam_dist_stat <- function(m) {
  st <- summary(m)$s.table
  rn <- grep("s\\(dist", rownames(st))
  if (length(rn) == 0) return(NULL)
  fcol <- if ("F" %in% colnames(st)) "F" else "Chi.sq"
  list(edf = st[rn[1], "edf"], stat = st[rn[1], fcol], p = st[rn[1], "p-value"])
}
fit_gam <- function(form, data) {
  m <- tryCatch(mgcv::bam(form, data = data, method = "fREML", discrete = TRUE),
                error = function(e) { message("    gam failed: ", conditionMessage(e)); NULL })
  if (is.null(m)) return(NULL)
  s <- .gam_dist_stat(m); if (is.null(s)) return(NULL)
  list(model = m, edf = s$edf, stat = s$stat, p = s$p)
}
predict_grid_gam <- function(model, df, x_cap) {
  grid_um <- seq(0, x_cap, length.out = 200)
  nd <- data.frame(
    dist      = grid_um,
    Sex       = factor(levels(df$Sex)[1], levels = levels(df$Sex)),
    Age_s     = mean(df$Age_s),
    PMI_s     = mean(df$PMI_s),
    celltype  = factor(levels(df$celltype)[1], levels = levels(df$celltype)),
    sample_id = factor(levels(df$sample_id)[1], levels = levels(df$sample_id))
  )
  pr  <- predict(model, newdata = nd, se.fit = TRUE, exclude = "s(sample_id)")
  fit <- as.numeric(pr$fit); se <- as.numeric(pr$se.fit)
  data.frame(dist_to_phf1_um = grid_um, fitted = fit,
             ci_lo = fit - 1.96 * se, ci_hi = fit + 1.96 * se)
}

# Save the pooled fitted-curve figure (single series; solid/bold if significant).
save_pooled_plot <- function(fit_df, significant, x_cap, out_pdf) {
  lw <- if (isTRUE(significant)) 0.8 else 0.4
  lt <- if (isTRUE(significant)) "solid" else "dashed"
  p <- ggplot(fit_df, aes(dist_to_phf1_um, fitted)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = POOLED_COL) +
    geom_line(colour = POOLED_COL, linewidth = lw, linetype = lt) +
    labs(x = XLAB, y = YLAB) +
    coord_cartesian(xlim = c(0, x_cap)) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.margin  = margin(t = 4, r = 6, b = 4, l = 5))
  g <- ggplot2::ggplotGrob(p)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  g$widths[pcol] <- grid::unit(1.4, "in")
  ggsave(out_pdf, g, width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
         height = 2.1, units = "in", device = "pdf")
}

# Model-free rolling (sliding-window) mean of y over x, evaluated on a grid
# (mean, +/- 1 SEM; 95% CI = +/- 1.96 SEM). Verbatim from
# plot_phf1_module_vs_phf1_distance_modelp.R::roll_mean.
roll_mean <- function(x, y, grid, half_window) {
  out <- lapply(grid, function(g) {
    idx <- which(x >= g - half_window & x <= g + half_window); n <- length(idx)
    if (n < 1) return(c(roll_mean = NA_real_, sem = NA_real_, n_window = 0))
    m <- mean(y[idx]); s <- if (n > 1) stats::sd(y[idx]) / sqrt(n) else NA_real_
    c(roll_mean = m, sem = s, n_window = n)
  })
  as.data.frame(do.call(rbind, out))
}

# Mean + 95% CI of a value over a subset (the PHF1+ neuron reference band).
mean_ci <- function(v) {
  v <- v[!is.na(v)]; n <- length(v)
  if (n < 1) return(NULL)
  m <- mean(v); sem <- if (n > 1) stats::sd(v) / sqrt(n) else NA_real_
  list(mean = m, sem = sem, ci_lo = m - 1.96 * sem, ci_hi = m + 1.96 * sem, n = n)
}

# Nakagawa marginal R2 (fixed-effects variance explained) for a fitted lmer.
r2m <- function(m) {
  if (is.null(m)) return(NA_real_)
  vf <- stats::var(as.numeric(stats::model.matrix(m) %*% lme4::fixef(m)))
  vc <- as.data.frame(lme4::VarCorr(m))
  vr <- sum(vc$vcov[vc$grp != "Residual"]); ve <- vc$vcov[vc$grp == "Residual"]
  as.numeric(vf / (vf + vr + ve))
}
# Clamp an underflowing p to P_FLOOR (returns list: clamped p + logical flag).
clamp_p <- function(p) {
  floored <- is.finite(p) && p < P_FLOOR
  list(p = if (floored) P_FLOOR else p, floored = floored)
}

##  ............................................................................
##  Load object, compute within-sample robust-z of p95 (over ALL cells)     ####
cat("Loading Seurat object (this is large)...\n")
seu <- readRDS(SEU_PATH)
md_all <- seu@meta.data
rm(seu); invisible(gc())

req <- c("celltype", "PHF1", "dist_to_phf1_um", "phf1_intensity_p95", "sample_id")
miss <- setdiff(req, colnames(md_all))
if (length(miss))
  stop(sprintf("seu@meta.data missing required columns: %s.\nRun add_phf1_intensity.R first.",
               paste(miss, collapse = ", ")))

sex_col <- pick_col(md_all, c("Sex"), "Sex")
age_col <- pick_col(md_all, c("Age"), "Age")
pmi_col <- pick_col(md_all, c("PMI", "PostMortemInterval", "PostMortem_Interval"), "PMI")
cat(sprintf("Donor covariate columns: Sex=%s Age=%s PMI=%s\n", sex_col, age_col, pmi_col))

# NA-out poorly-registered masks exactly as add_phf1_intensity.R (belt & braces;
# they are usually already NA on the object).
if ("phf1_frac_in_bounds" %in% colnames(md_all)) {
  poor <- !is.na(md_all$phf1_frac_in_bounds) & md_all$phf1_frac_in_bounds < FRAC_IN_BOUNDS_MIN
  md_all$phf1_intensity_p95[poor] <- NA
  cat(sprintf("Set %d cells with frac_in_bounds < %.2f to NA p95.\n", sum(poor), FRAC_IN_BOUNDS_MIN))
}

# Within-sample robust z over ALL cells (matches phf1_intensity_p95_z provenance).
md_all <- md_all %>%
  group_by(sample_id) %>%
  mutate(p95_z = robust_z(phf1_intensity_p95)) %>%
  ungroup() %>%
  as.data.frame()

##  ............................................................................
##  PHF1+ reference (descriptive; the "distance 0" ceiling)                 ####
neuron_cts <- intersect(neuron_order, unique(as.character(md_all$celltype)))
cat("Neuron subtypes present:", paste(neuron_cts, collapse = ", "), "\n")

ref_tab <- md_all %>%
  filter(celltype %in% neuron_cts, !is.na(p95_z)) %>%
  mutate(group = ifelse(PHF1, "PHF1_pos_neuron", "PHF1_neg_neuron")) %>%
  group_by(group) %>%
  summarise(n = dplyr::n(), median_p95_z = median(p95_z),
            median_p95_raw = median(phf1_intensity_p95), .groups = "drop")
phf1_ref_z <- ref_tab$median_p95_z[ref_tab$group == "PHF1_pos_neuron"]
if (length(phf1_ref_z) == 0) phf1_ref_z <- NA_real_

##  ............................................................................
##  Build the modelled set: non-PHF1 neurons, within cap, complete cases    ####
dat <- md_all %>%
  transmute(
    cell_id         = rownames(md_all),
    celltype        = as.character(celltype),
    PHF1            = PHF1,
    dist_to_phf1_um = as.numeric(dist_to_phf1_um),
    p95_z           = p95_z,
    sample_id       = as.character(sample_id),
    Sex             = as.character(.data[[sex_col]]),
    Age             = as.numeric(.data[[age_col]]),
    PMI             = as.numeric(.data[[pmi_col]])
  ) %>%
  filter(celltype %in% neuron_cts,
         PHF1 == FALSE,                       # manual-negative neurons only
         !is.na(dist_to_phf1_um),             # also drops PHF1+ (dist is NA there)
         dist_to_phf1_um <= MAX_DIST,
         dist_to_phf1_um > 0,                 # log() defined; PHF1-neg dist is > 0
         !is.na(p95_z), !is.na(Age), !is.na(PMI), !is.na(Sex))

dat$Sex       <- droplevels(factor(dat$Sex))
dat$sample_id <- factor(dat$sample_id)
dat$celltype  <- factor(dat$celltype, levels = neuron_cts)
dat$celltype  <- droplevels(dat$celltype)
dat$Age_s     <- as.numeric(scale(dat$Age))
dat$PMI_s     <- as.numeric(scale(dat$PMI))

n_cells  <- nrow(dat)
n_donors <- dplyr::n_distinct(dat$sample_id)
cat(sprintf("Modelled non-PHF1 neurons (<= %g um): %d cells, %d donors.\n",
            MAX_DIST, n_cells, n_donors))
if (n_cells < 50)  stop("Too few modelled cells (<50).")
if (n_donors < MIN_DONORS) stop("Too few donors.")

coef_rows <- list()          # structured coefficient accumulator
pooled_log_sig <- NA         # set by run_pooled("log"); styles the rolling-mean line

##  ............................................................................
##  POOLED models: linear (headline), log, GAM                              ####
run_pooled <- function(model_type) {
  is_gam       <- model_type == "gam"
  log_distance <- model_type == "log"
  suffix       <- switch(model_type, linear = "", log = "_log", gam = "_gam")
  scale_label  <- switch(model_type, linear = "um", log = "log_um", gam = "gam")

  if (is_gam) {
    dat$dist <- dat$dist_to_phf1_um
    obs  <- fit_gam(p95_z ~ s(dist, k = GAM_K) + celltype + Sex + Age_s + PMI_s +
                      s(sample_id, bs = "re"), dat)
    if (is.null(obs)) { cat("  pooled GAM failed\n"); return(invisible(NULL)) }
    grid <- predict_grid_gam(obs$model, dat, MAX_DIST)
    dist_sd <- NA_real_
    pc <- clamp_p(obs$p)
    row <- tibble(model = "pooled", distance_scale = scale_label,
                  slope_per_sd = NA_real_, se = NA_real_, ci_lo_per_sd = NA_real_,
                  ci_hi_per_sd = NA_real_, slope_per_unit = NA_real_,
                  stat = obs$stat, df = obs$edf, p = pc$p, p_at_floor = pc$floored,
                  singular = NA, dist_sd = dist_sd,
                  r2_marginal = NA_real_, sigma_resid = NA_real_,
                  outcome_sd = stats::sd(dat$p95_z), n_cells = n_cells, n_donors = n_donors)
  } else {
    dt <- if (log_distance) log(dat$dist_to_phf1_um) else dat$dist_to_phf1_um
    dist_sd <- sd(dt); dat$dist_scaled <- dt / dist_sd
    fe_form <- ~ dist_scaled + celltype + Sex + Age_s + PMI_s
    obs <- fit_lmer(p95_z ~ dist_scaled + celltype + Sex + Age_s + PMI_s + (1 | sample_id), dat)
    if (is.null(obs)) { cat("  pooled LMM (", scale_label, ") failed\n"); return(invisible(NULL)) }
    grid <- predict_grid_lmer(obs$model, dat, dist_sd, log_distance, MAX_DIST,
                              fe_form, include_celltype = TRUE)
    pc <- clamp_p(obs$p)
    row <- tibble(model = "pooled", distance_scale = scale_label,
                  slope_per_sd = obs$estimate, se = obs$se,
                  ci_lo_per_sd = obs$estimate - 1.96 * obs$se,
                  ci_hi_per_sd = obs$estimate + 1.96 * obs$se,
                  slope_per_unit = obs$estimate / dist_sd,
                  stat = obs$t, df = obs$df, p = pc$p, p_at_floor = pc$floored,
                  singular = obs$singular, dist_sd = dist_sd,
                  r2_marginal = r2m(obs$model), sigma_resid = stats::sigma(obs$model),
                  outcome_sd = stats::sd(dat$p95_z), n_cells = n_cells, n_donors = n_donors)
  }
  # Pooled is one hypothesis viewed 3 ways -> no cross-model multiplicity; each
  # judged by its own p (headline = linear).
  significant <- is.finite(row$p) && row$p < SIG_ALPHA
  coef_rows[[paste0("pooled", suffix)]] <<- row
  # The rolling-mean actual-data plot keys its line style off the LOG model.
  if (log_distance) pooled_log_sig <<- significant

  grid$significant <- significant
  save_pooled_plot(grid, significant, MAX_DIST,
                   file.path(OUTPUT_DIR, sprintf("plot_phf1_p95_gradient_neurons%s.pdf", suffix)))

  fit_out <- grid %>% transmute(model = "pooled", distance_scale = scale_label,
                                dist_to_phf1_um, fitted, ci_lo, ci_hi, significant)
  write.table(fit_out, file.path(OUTPUT_DIR,
              sprintf("source_data_phf1_p95_gradient_neurons%s_fit.tsv", suffix)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  pooled %s: slope_per_sd=%s p=%.3g %s\n", scale_label,
              ifelse(is.na(row$slope_per_sd), "NA", sprintf("%.4f", row$slope_per_sd)),
              row$p, ifelse(significant, "*", "")))
  invisible(row)
}
for (mt in c("linear", "log", "gam")) run_pooled(mt)

##  ............................................................................
##  Rolling-mean actual-data figure (significance from the LOG model)       ####
## Model-free sliding-window mean of the observed p95_z vs distance (matches the
## module-score raw plots): points + rolling mean + 95% CI (= +/-1.96 SEM). The
## rolling-mean line is SOLID+thick when the LOG-distance LMM is significant,
## DASHED+thin otherwise -- i.e. significance is assigned by the log model.
roll_grid <- seq(0, MAX_DIST, length.out = 200)
roll_df   <- roll_mean(dat$dist_to_phf1_um, dat$p95_z, roll_grid, HALF_WINDOW)
roll_df$dist_to_phf1_um <- roll_grid

# PHF1+ neuron reference (mean +/- 95% CI). Reported, NOT drawn: PHF1+ intensity
# sits far above the non-PHF1 range and would squash the gradient off the panel.
phf1_ref <- mean_ci(md_all$p95_z[md_all$celltype %in% neuron_cts & md_all$PHF1])

roll_sig <- isTRUE(pooled_log_sig)
roll_lty <- if (roll_sig) "solid" else "dashed"
roll_lwd <- if (roll_sig) 0.8 else 0.4

pts <- dat[, c("dist_to_phf1_um", "p95_z")]
if (POINT_MAX > 0 && nrow(pts) > POINT_MAX) pts <- pts[sample.int(nrow(pts), POINT_MAX), ]
cat(sprintf("Rolling mean: window=%g um, %s points%s; log-model significant=%s\n",
            WINDOW_UM, if (POINT_MAX > 0) format(nrow(pts), big.mark = ",") else "no",
            if (POINT_MAX > 0 && nrow(dat) > POINT_MAX) sprintf(" (subsampled from %s)",
              format(nrow(dat), big.mark = ",")) else "", roll_sig))

p_roll <- ggplot()
if (POINT_MAX > 0)
  p_roll <- p_roll + geom_point(data = pts, aes(dist_to_phf1_um, p95_z),
                                colour = POOLED_COL, alpha = 0.12, size = 0.3, stroke = 0)
p_roll <- p_roll +
  geom_ribbon(data = roll_df, aes(dist_to_phf1_um, ymin = roll_mean - 1.96 * sem,
              ymax = roll_mean + 1.96 * sem), fill = POOLED_COL, alpha = 0.25) +
  geom_line(data = roll_df, aes(dist_to_phf1_um, roll_mean),
            colour = POOLED_COL, linetype = roll_lty, linewidth = roll_lwd) +
  labs(x = XLAB, y = YLAB) +
  coord_cartesian(xlim = c(0, MAX_DIST)) +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin  = margin(t = 4, r = 6, b = 4, l = 5))
g <- ggplot2::ggplotGrob(p_roll)
pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
g$widths[pcol] <- grid::unit(1.4, "in")
ggsave(file.path(OUTPUT_DIR, "plot_phf1_p95_gradient_neurons_rollmean.pdf"),
       g, width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
       height = 2.1, units = "in", device = "pdf")

write.table(roll_df %>% transmute(model = "pooled_rollmean", window_um = WINDOW_UM,
              dist_to_phf1_um, roll_mean, sem, ci_lo = roll_mean - 1.96 * sem,
              ci_hi = roll_mean + 1.96 * sem, n_window, log_model_significant = roll_sig),
            file.path(OUTPUT_DIR, "source_data_phf1_p95_gradient_neurons_rollmean.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

##  ............................................................................
##  Per-subtype linear LMM (sensitivity; BH across subtypes)                ####
sub_stats <- list(); sub_fit <- list()
for (ct in levels(dat$celltype)) {
  d <- dat[dat$celltype == ct, ]
  d$sample_id <- droplevels(d$sample_id)
  nd_don <- dplyr::n_distinct(d$sample_id)
  if (nrow(d) < MIN_CELLS_SUBTYPE || nd_don < MIN_DONORS) {
    cat(sprintf("  %-28s skipped (%d cells, %d donors)\n", ct, nrow(d), nd_don)); next
  }
  d$Age_s <- as.numeric(scale(d$Age)); d$PMI_s <- as.numeric(scale(d$PMI))
  fe_form <- ~ dist_scaled + Sex + Age_s + PMI_s
  # Fit BOTH scales (matching the module-score analysis): linear (headline, plotted) and
  # log (consistency; dist_scaled = log(dist)/sd(log(dist))). One coef row per scale.
  for (sc in c("linear", "log")) {
    log_distance <- sc == "log"
    scale_label  <- if (log_distance) "log_um" else "um"
    dt      <- if (log_distance) log(d$dist_to_phf1_um) else d$dist_to_phf1_um
    dist_sd <- sd(dt); d$dist_scaled <- dt / dist_sd
    obs <- fit_lmer(p95_z ~ dist_scaled + Sex + Age_s + PMI_s + (1 | sample_id), d)
    if (is.null(obs)) { cat(sprintf("  %-28s LMM (%s) failed\n", ct, scale_label)); next }
    if (!log_distance) {                     # fitted plot stays LINEAR (unchanged)
      grid <- predict_grid_lmer(obs$model, d, dist_sd, FALSE, MAX_DIST, fe_form,
                                include_celltype = FALSE)
      grid$celltype <- ct
      sub_fit[[ct]] <- grid
    }
    pc <- clamp_p(obs$p)
    sub_stats[[paste(ct, scale_label)]] <- tibble(
      model = ct, distance_scale = scale_label,
      slope_per_sd = obs$estimate, se = obs$se,
      ci_lo_per_sd = obs$estimate - 1.96 * obs$se,
      ci_hi_per_sd = obs$estimate + 1.96 * obs$se,
      slope_per_unit = obs$estimate / dist_sd,
      stat = obs$t, df = obs$df, p = pc$p, p_at_floor = pc$floored,
      singular = obs$singular, dist_sd = dist_sd,
      r2_marginal = r2m(obs$model), sigma_resid = stats::sigma(obs$model),
      outcome_sd = stats::sd(d$p95_z), n_cells = nrow(d), n_donors = nd_don)
  }
}

if (length(sub_stats)) {
  sub_stats <- bind_rows(sub_stats)
  # BH within each distance scale (each scale is its own family of subtype tests).
  sub_stats <- sub_stats %>% group_by(distance_scale) %>%
    mutate(padj = p.adjust(p, method = "BH"), significant = padj < SIG_ALPHA) %>%
    ungroup() %>% as.data.frame()
  # Export the full subtype rows (now incl. n_cells, n_donors, padj + R2/sigma/outcome_sd/p_at_floor).
  coef_rows[["subtypes"]] <- sub_stats %>%
    dplyr::select(model, distance_scale, slope_per_sd, se, ci_lo_per_sd, ci_hi_per_sd,
                  slope_per_unit, stat, df, p, p_at_floor, padj, singular, dist_sd,
                  r2_marginal, sigma_resid, outcome_sd, n_cells, n_donors, significant)

  sub_fit_df <- bind_rows(sub_fit) %>%
    left_join(sub_stats %>% filter(distance_scale == "um") %>%   # plot uses the LINEAR fit
                dplyr::select(celltype = model, significant), by = "celltype") %>%
    mutate(celltype = factor(celltype, levels = neuron_cts),
           sig = factor(ifelse(significant, "padj<0.05", "ns"), levels = c("padj<0.05", "ns")))

  p_sub <- ggplot(sub_fit_df, aes(dist_to_phf1_um, fitted, colour = celltype, group = celltype)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = celltype), alpha = 0.10, colour = NA) +
    geom_line(aes(linetype = sig, linewidth = sig)) +
    scale_colour_manual(values = celltype_palette, name = NULL, drop = FALSE) +
    scale_fill_manual(values = celltype_palette, guide = "none") +
    scale_linetype_manual(values = setNames(c("solid", "dashed"), c("padj<0.05", "ns")),
                          name = NULL, drop = FALSE, limits = c("padj<0.05", "ns")) +
    scale_linewidth_manual(values = setNames(c(0.8, 0.4), c("padj<0.05", "ns")),
                           name = NULL, drop = FALSE, limits = c("padj<0.05", "ns")) +
    guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2),
           linewidth = guide_legend(order = 2)) +
    labs(x = XLAB, y = YLAB) +
    coord_cartesian(xlim = c(0, MAX_DIST)) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.text = element_text(size = 6), legend.key.width = grid::unit(20, "pt"),
          legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
  g <- ggplot2::ggplotGrob(p_sub)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  g$widths[pcol] <- grid::unit(1.4, "in")
  ggsave(file.path(OUTPUT_DIR, "plot_phf1_p95_gradient_by_subtype.pdf"),
         g, width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
         height = 2.1, units = "in", device = "pdf")

  write.table(sub_fit_df %>% transmute(celltype = as.character(celltype), distance_scale = "um",
                                       dist_to_phf1_um, fitted, ci_lo, ci_hi, significant),
              file.path(OUTPUT_DIR, "source_data_phf1_p95_gradient_by_subtype_fit.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  sub_stats <- NULL
  cat("  No subtype passed the cell/donor guards; skipping subtype panel.\n")
}

##  ............................................................................
##  Per-subtype ROLLING-MEAN figure (75 um window; NO fit)                  ####
## Mirrors the module-score integrated_raw plot: one model-free sliding-window
## mean of observed p95_z per neuron subtype, coloured by celltype, with 95% CI
## ribbons (+/-1.96 SEM). No raw points, no fitted line, no significance styling.
MIN_WINDOW_N <- 10   # blank a grid point where its window holds < this many cells (avoids noisy tails)
sub_roll <- lapply(levels(dat$celltype), function(ct) {
  d <- dat[dat$celltype == ct, ]
  if (nrow(d) < MIN_CELLS_SUBTYPE) { cat(sprintf("  rollmean %-28s skipped (%d cells)\n", ct, nrow(d))); return(NULL) }
  r <- roll_mean(d$dist_to_phf1_um, d$p95_z, roll_grid, HALF_WINDOW)
  r$dist_to_phf1_um <- roll_grid
  r$roll_mean[r$n_window < MIN_WINDOW_N] <- NA
  r$sem[r$n_window < MIN_WINDOW_N]       <- NA
  r$celltype <- ct
  r
})
sub_roll_df <- bind_rows(sub_roll)

# Draw/legend order: put the first neuron_order subtype (Exc-IT-L2-3-CBLN2-HOPX)
# ON TOP and step backwards through neuron_order -> reverse the factor levels so it
# is drawn last (front); guide reverse keeps the legend reading in natural order
# (L2-3 first). FRONT_CT also gets its own single-cluster panel below.
FRONT_CT <- neuron_cts[1]

# Shared builder so the combined and single-cluster panels match exactly.
build_subtype_rollmean <- function(df, level_order, legend, panel_in = 1.6) {
  df$celltype <- factor(as.character(df$celltype), levels = level_order)
  p <- ggplot(df, aes(dist_to_phf1_um, roll_mean, colour = celltype, group = celltype)) +
    geom_ribbon(aes(ymin = roll_mean - 1.96 * sem, ymax = roll_mean + 1.96 * sem, fill = celltype),
                alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = celltype_palette, name = NULL, drop = FALSE) +
    scale_fill_manual(values = celltype_palette, guide = "none", drop = FALSE) +
    labs(x = XLAB, y = YLAB) +
    coord_cartesian(xlim = c(0, MAX_DIST)) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.text = element_text(size = 6), legend.key.width = grid::unit(14, "pt"),
          legend.key.height = grid::unit(9, "pt"), legend.position = if (legend) "right" else "none",
          legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
  if (legend) p <- p + guides(colour = guide_legend(reverse = TRUE))  # legend reads natural order
  g <- ggplot2::ggplotGrob(p)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  g$widths[pcol] <- grid::unit(panel_in, "in")
  g
}

if (nrow(sub_roll_df)) {
  # Combined: all subtypes, FRONT_CT drawn on top (reversed levels).
  g_all <- build_subtype_rollmean(sub_roll_df, level_order = rev(neuron_cts), legend = TRUE)
  ggsave(file.path(OUTPUT_DIR, "plot_phf1_p95_gradient_by_subtype_rollmean.pdf"),
         g_all, width = grid::convertWidth(sum(g_all$widths), "in", valueOnly = TRUE),
         height = 2.1, units = "in", device = "pdf")
  write.table(sub_roll_df %>% transmute(celltype = as.character(celltype), window_um = WINDOW_UM,
                dist_to_phf1_um, roll_mean, sem, ci_lo = roll_mean - 1.96 * sem,
                ci_hi = roll_mean + 1.96 * sem, n_window),
              file.path(OUTPUT_DIR, "source_data_phf1_p95_gradient_by_subtype_rollmean.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  wrote per-subtype rolling-mean plot (%d subtypes; %s on top).\n",
              dplyr::n_distinct(sub_roll_df$celltype), FRONT_CT))

  # Single-cluster version: FRONT_CT only.
  front_df <- sub_roll_df[as.character(sub_roll_df$celltype) == FRONT_CT, ]
  if (nrow(front_df)) {
    front_safe <- gsub("[^A-Za-z0-9_-]", "_", FRONT_CT)
    g_one <- build_subtype_rollmean(front_df, level_order = FRONT_CT, legend = FALSE)
    ggsave(file.path(OUTPUT_DIR, sprintf("plot_phf1_p95_gradient_%s_rollmean.pdf", front_safe)),
           g_one, width = grid::convertWidth(sum(g_one$widths), "in", valueOnly = TRUE),
           height = 2.1, units = "in", device = "pdf")
    write.table(front_df %>% transmute(celltype = as.character(celltype), window_um = WINDOW_UM,
                  dist_to_phf1_um, roll_mean, sem, ci_lo = roll_mean - 1.96 * sem,
                  ci_hi = roll_mean + 1.96 * sem, n_window),
                file.path(OUTPUT_DIR, sprintf("source_data_phf1_p95_gradient_%s_rollmean.tsv", front_safe)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  wrote single-cluster rolling-mean plot for %s.\n", FRONT_CT))
  }
}

##  ............................................................................
##  Supplementary source data: per-donor x distance-ring median z           ####
ring_long <- dat %>%
  mutate(dist_ring = cut(dist_to_phf1_um, breaks = RING_BREAKS, labels = RING_LAB,
                         include.lowest = TRUE, right = TRUE)) %>%
  filter(!is.na(dist_ring)) %>%
  group_by(sample_id, dist_ring) %>%
  summarise(median_p95_z = median(p95_z), mean_p95_z = mean(p95_z), n_cells = dplyr::n(),
            .groups = "drop") %>%
  mutate(ring_mid_um = RING_MID[as.character(dist_ring)]) %>%
  dplyr::select(sample_id, dist_ring, ring_mid_um, median_p95_z, mean_p95_z, n_cells)
write.table(ring_long, file.path(OUTPUT_DIR, "source_data_phf1_p95_gradient_neurons_persample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

##  ............................................................................
##  Structured coefficient table (supp)                                     ####
coeff_tab <- bind_rows(lapply(coef_rows, function(x) x))
write.table(coeff_tab, file.path(OUTPUT_DIR, "stats_phf1_p95_gradient_coeffs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

##  ............................................................................
##  Stats log (triple)                                                      ####
sink(file.path(OUTPUT_DIR, "stats_phf1_p95_gradient_neurons.txt"))
cat("PHF1 p95 intensity gradient in NON-PHF1 (manual) neurons vs distance to PHF1+ neuron\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Outcome: phf1_intensity_p95 (per-cell 95th-pct masked PHF1 DN),\n")
cat("  within-sample ROBUST z = (x - median)/MAD per sample_id (over all cells).\n")
cat("Population: neurons (celltype in neuron_order), manual PHF1 == FALSE,\n")
cat(sprintf("  0 < dist_to_phf1_um <= %g um. Cells=%d, donors=%d.\n\n", MAX_DIST, n_cells, n_donors))
cat("Models:\n")
cat("  pooled linear (headline): p95_z ~ dist_scaled + celltype + Sex + Age_s + PMI_s + (1|sample_id)\n")
cat("  pooled log:               as above, dist_scaled = log(dist)/sd(log(dist))\n")
cat(sprintf("  pooled GAM:               p95_z ~ s(dist,k=%d) + celltype + Sex + Age_s + PMI_s + s(sample_id,bs='re')\n", GAM_K))
cat("  per-subtype linear + log: p95_z ~ dist_scaled + Sex + Age_s + PMI_s + (1|sample_id)\n")
cat("    fit on both scales (um, log_um), one coef row each; BH WITHIN each scale. Plot = linear.\n")
cat("  Age/PMI z-scaled; dist_scaled = distance/SD (per model). slope_per_unit = slope_per_sd/dist_sd.\n")
cat("  Every LMM row also carries r2_marginal (Nakagawa, fixed effects), sigma_resid, and\n")
cat("    outcome_sd = sd(p95_z) in the modelled set (slope relative to cell-to-cell variation).\n")
cat(sprintf("  p clamped to P_FLOOR=%.3g on underflow; p_at_floor records it (e.g. pooled GAM p=0).\n", P_FLOOR))
cat("  A NEGATIVE, significant distance slope => graded PHF1 halo around PHF1+ neurons.\n")
cat(sprintf("  Rolling-mean actual-data plot: model-free sliding window (%g um) of observed p95_z;\n", WINDOW_UM))
cat("    line SOLID+thick if the LOG-distance LMM is significant, DASHED+thin otherwise\n")
cat(sprintf("    (log-model significant = %s). 95%% CI = +/-1.96 SEM.\n\n", isTRUE(pooled_log_sig)))
cat("Reference intensity by neuron PHF1 status (descriptive, NOT modelled):\n")
print(as.data.frame(ref_tab))
cat(sprintf("\n  PHF1+ neuron median p95_z = %.3f (the 'distance 0' ceiling; PHF1- neuron curve\n", phf1_ref_z))
cat("  should sit below this and be highest near distance 0 if a halo exists).\n")
if (!is.null(phf1_ref))
  cat(sprintf("  PHF1+ neuron mean p95_z = %.3f (95%% CI %.3f, %.3f; n=%d) -- reference band value,\n",
              phf1_ref$mean, phf1_ref$ci_lo, phf1_ref$ci_hi, phf1_ref$n))
cat("  reported not drawn (it sits far above the non-PHF1 range and would squash the panel).\n\n")
cat("Pooled model coefficients:\n")
print(as.data.frame(coeff_tab[coeff_tab$model == "pooled", ]))
if (!is.null(sub_stats)) {
  cat("\nPer-subtype coefficients (BH within distance scale; singular flag carried):\n")
  print(as.data.frame(sub_stats %>% dplyr::select(model, distance_scale, slope_per_sd,
        ci_lo_per_sd, ci_hi_per_sd, t = stat, df, p, p_at_floor, padj, significant, singular,
        r2_marginal, sigma_resid, outcome_sd, n_cells, n_donors)))
}
cat("\nDistance rings (um):", paste(RING_LAB, collapse = ", "), "\n")
cat("\nsessionInfo():\n"); print(sessionInfo())
sink()

cat("\nDone. Outputs in", OUTPUT_DIR, "\n")
