#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# imc_phf1_glia_distance.R
#
# Figure panels: 6C - upstream step
# Writes the IMC GFAP rolling mean and Set 3 fit
# (plots/imc_phf1_glia_distance/) that imc_cosmx_gfap_distance_panel.R reads.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# imc_phf1_glia_distance.R
#
# IMC (EC project) GLIAL response over continuous distance to the nearest PHF1+ neuron,
# AD cases only. Four panels:
#   1. CD68 intensity in POOLED microglia            (2 clusters combined)
#   2. GFAP intensity in POOLED astrocytes           (2 clusters combined)
#   3. Proportion of cells in each MICROGLIAL cluster over distance
#   4. Proportion of cells in each ASTROCYTE  cluster over distance
#
# Statistics and plotting deliberately mirror the CosMx module-score figures
# (R/plot_modulescore_vs_phf1_distance_modelp.R): 75 um sliding-window mean on a 200-point grid
# pooled across donors, ribbon = mean +/- 1.96 x within-window SEM, log-distance LMM slope with a
# donor random intercept, BH within the panel, and significance encoded as linetype+linewidth
# rather than as printed text.
#
# CLUSTERING INPUTS. CD68, GFAP, Iba1 and S100b are all `type` markers, and
# rowData(spe)$use_channel (the PCA input feeding Harmony -> Phenograph) is exactly those 22 type
# markers (08-phenotyping.Rmd:559-568, 07-batch_correction.Rmd:217-229), so the glial clusters
# were defined using the markers scored here. The mixing decomposition at the end of this script
# relates the pooled-intensity panels (1-2) to the cluster-mix panels (3-4).
#
# PHF1+ glia are kept in the intensity analysis.

suppressPackageStartupMessages({
  library(SpatialExperiment); library(dplyr); library(tibble); library(tidyr)
  library(ggplot2); library(lmerTest); library(RANN); library(grid)
})

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")
source("R/imc_utils.R")   # cat_re_variance(): every IMC log reports its donor RE variance

out_dir <- "plots/imc_phf1_glia_distance"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

SIG_ALPHA   <- 0.05
# 300 um cap, chosen from the measured anisotropy of the ROIs. ROIs are cortical strips (median
# 478 um wide x 2525 tall), so beyond some radius "far" is only reachable along the cortical
# axis. The measured decay is gradual, not a cliff: mean |dx|/dist is 0.59 at 150-200 um, 0.57
# at 200-250, 0.53 at 250-300, and only collapses to 0.31 beyond 400 um. Depth drift over
# 0-300 um is ~10 of 100 scaled_Y units with the SD steady (~25-27), so the far band is not
# collapsing onto one lamina. At 300 um all 44 donors and 129 of 131 ROIs still contribute.
MAX_DIST    <- 300
WINDOW_UM   <- 75      # rolling-mean window WIDTH, matches plot_modulescore_vs_phf1_distance_modelp.R
HALF_WINDOW <- WINDOW_UM / 2
N_GRID      <- 200
MIN_WINDOW_N <- 10     # do not draw a "mean" over fewer than 10 cells (as
                       # phf1_intensity_gradient_neurons.R / nn3_over_phf1_distance.R do)
# must reach MAX_DIST, else cut() silently NAs the outermost cells out of the ring tables
# steps no coarser than 50 um all the way to the cap, matching brks in the marker/morphology
# scripts -- a 100 um-wide final ring would dominate the supplementary per-sample table.
RING_BREAKS <- c(0, 25, 50, 100, 150, 200, 250, MAX_DIST)
SIG_BASIS   <- "padj<0.05"

# Distinct qualitative palette, positional -- NOT celltype_palette (matches MODULE_COLOURS in
# plot_modulescore_vs_phf1_distance_modelp.R).
MODULE_COLOURS <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                    "#A65628", "#F781BF", "#666666", "#66C2A5", "#FC8D62")

MICRO <- c("Microglia cluster 1 (IBA1)", "Microglia cluster 2 (IBA1, CD68)")
ASTRO <- c("Astrocyte cluster 1 (GFAP, S100b)", "Astrocyte cluster 2 (S100b)")

localwd <- "<IMC_ROOT>/"
spe <- readRDS(paste0(localwd, "spe.rds"))
stopifnot("PHF1_Otsu" %in% colnames(colData(spe)),
          "Matched_4G8_40" %in% colnames(colData(spe)))
stopifnot(all(c(MICRO, ASTRO) %in% unique(as.character(colData(spe)$celltype_clusters))))
stopifnot(all(c("CD68", "GFAP") %in% rownames(spe)))

# -------------------------------------------------------------------
# Cohort + distance backbone (shared with the other IMC distance scripts)
# -------------------------------------------------------------------
cd <- as.data.frame(colData(spe)); cd$cell_id <- colnames(spe)
sc <- spatialCoords(spe); cd$x <- sc[, "Pos_X"]; cd$y <- sc[, "Pos_Y"]

cv_removal <- cd %>% distinct(patient_id, sample_id) %>% count(patient_id) %>%
  filter(n < 3) %>% pull(patient_id)

base <- cd %>%
  filter(!celltype_clusters %in% c("Artefact cluster", "Unassigned cluster"),
         !patient_id %in% cv_removal, BraakGroup != "Braak_0_1", !is.na(PHF1_Otsu)) %>%
  mutate(patient_id = as.character(patient_id), sample_id = as.character(sample_id),
         celltype_clusters = as.character(celltype_clusters),
         Sex = factor(Sex)) %>%
  group_by(sample_id) %>%
  mutate(scaled_Y = 100 * (y - min(y)) / (max(y) - min(y))) %>%
  ungroup() %>%
  # donor covariates z-scored on the full pre-cap set, as the CosMx script does
  # z-score all continuous covariates: glmer warns "Rescale variables?" if scaled_Y is left on
  # its native 0-100 range alongside unit-scale predictors (large eigenvalue ratio).
  mutate(Age_s = as.numeric(scale(Age)), PMI_s = as.numeric(scale(PMI)),
         depth_s = as.numeric(scale(scaled_Y)),
         # amyloid proximity, defined as in the other IMC scripts: nucleus within a 40 um-dilated
         # 4G8+ plaque mask, so the glial panels carry the same with/without-amyloid pair as the
         # marker and morphology panels.
         plaque = factor(ifelse(Matched_4G8_40, "plaque-prox", "plaque-distal"),
                         levels = c("plaque-distal", "plaque-prox")))

anchors <- base %>% filter(PHF1_Otsu == "PHF1_pos",
                           grepl("neuron", celltype_clusters, ignore.case = TRUE))

# distance from EVERY cell to the nearest PHF1+ neuron in its own ROI (the proportion panels
# need all cells as the denominator; the intensity panels subset afterwards)
base$dist_um <- NA_real_
for (s in unique(base$sample_id)) {
  i <- which(base$sample_id == s)
  a <- anchors[anchors$sample_id == s, c("x", "y"), drop = FALSE]
  if (!nrow(a)) next
  base$dist_um[i] <- RANN::nn2(as.matrix(a), as.matrix(base[i, c("x", "y")]), k = 1)$nn.dists[, 1]
}

# dist > 0 BEFORE the log: a zero distance would give -Inf -> sd = NaN. PHF1+ neurons are their
# own anchors so they sit at exactly 0.
md <- base %>% filter(!is.na(dist_um), dist_um > 0, dist_um <= MAX_DIST)
DIST_SD <- sd(log(md$dist_um))
stopifnot(is.finite(DIST_SD), DIST_SD > 0)
md$dist_scaled <- log(md$dist_um) / DIST_SD      # plain log, /SD, NOT centred (canonical)

md$CD68 <- as.numeric(assay(spe, "exprs")["CD68", md$cell_id])
md$GFAP <- as.numeric(assay(spe, "exprs")["GFAP", md$cell_id])
md$is_micro <- md$celltype_clusters %in% MICRO
md$is_astro <- md$celltype_clusters %in% ASTRO

GRID <- seq(0, MAX_DIST, length.out = N_GRID)

# Lineage composition by distance bin. The proportion panels use ALL cells as the denominator, so
# they are ZERO-SUM: any lineage rising with distance is mirrored by another falling. This table is
# printed into the proportion stats logs so that property is visible rather than implicit.
COMPOSITION <- md %>%
  mutate(lineage = case_when(
           grepl("neuron", celltype_clusters, ignore.case = TRUE) ~ "neuron",
           celltype_clusters %in% MICRO                          ~ "microglia",
           celltype_clusters %in% ASTRO                          ~ "astrocyte",
           TRUE                                                  ~ "other"),
         bin = cut(dist_um, RING_BREAKS, include.lowest = TRUE, right = TRUE)) %>%
  filter(!is.na(bin)) %>%
  count(bin, lineage) %>% group_by(bin) %>%
  mutate(prop = n / sum(n), n_bin = sum(n)) %>% ungroup() %>%
  select(bin, lineage, prop, n_bin) %>%
  tidyr::pivot_wider(names_from = lineage, values_from = prop)

# -------------------------------------------------------------------
# roll_mean(): copied from plot_modulescore_vs_phf1_distance_modelp.R
# Sliding-window mean of y over x on a grid. Returns mean, +/-1 SEM, window n.
# -------------------------------------------------------------------
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

# -------------------------------------------------------------------
# Shared figure builder (CosMx module-score spec)
# -------------------------------------------------------------------
XLAB <- expression("Distance to PHF1+ neuron (" * mu * "m)")

# A guide key is built from LAYER DATA, so a level that never occurs (e.g. every series here is
# significant) yields a key with its label but NO GLYPH -- the "ns" dashed line is simply absent.
# drop = FALSE and override.aes both fail to fix that on their own. Carrying one undrawn row per
# missing level gives the key its glyph; nothing is plotted because roll_mean is NA. Ported from
# R/phf1_channel_intensity_exc.R; keep the two in step.
pad_levels <- function(d, col = "sig") {
  miss <- setdiff(c(SIG_BASIS, "ns"), as.character(d[[col]]))
  if (!length(miss)) return(d)
  filler <- d[rep(1L, length(miss)), , drop = FALSE]
  filler[[col]] <- factor(miss, levels = c(SIG_BASIS, "ns"))
  for (v in intersect(c("roll_mean", "sem"), names(filler))) filler[[v]] <- NA_real_
  dplyr::bind_rows(d, filler)
}

make_panel <- function(roll_df, ylab, legend_title, cols) {
  roll_df <- pad_levels(roll_df, "sig")
  # group on series AND sig: sig is constant within a series, so this never splits a real line,
  # but it does guarantee the padded row is its own group and cannot join one.
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
    labs(x = XLAB, y = ylab) +
    coord_cartesian(xlim = c(0, MAX_DIST)) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.text = element_text(size = 6),
          legend.key.width = grid::unit(18, "pt"), legend.key.height = grid::unit(8, "pt"),
          legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
  # fixed physical panel width so x-axes align across panels regardless of legend size
  g <- ggplot2::ggplotGrob(p)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  g$widths[pcol] <- grid::unit(1.6, "in")
  g
}

# total_w_in = NULL keeps the CosMx single-panel behaviour (panel forced to 1.6 in, total width
# floats with the legend). Given a value, the PANEL is shrunk so the TOTAL figure width equals it
# exactly -- used for the intensity panels so their x-axis matches one half of
# plot_modulescore_raw_COMBINED_Micro_Astro.pdf (360 pt / 2 panels = 180 pt = 2.5 in).
save_panel <- function(g, slug, total_w_in = NULL, total_h_in = 2.1) {
  if (!is.null(total_w_in)) {
    pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
    furniture <- sum(grid::convertWidth(g$widths[-pcol], "in", valueOnly = TRUE))
    g$widths[pcol] <- grid::unit(max(total_w_in - furniture, 0.6), "in")
  }
  ggsave(file.path(out_dir, sprintf("plot_%s.pdf", slug)), g,
         width = if (is.null(total_w_in))
                   grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE) else total_w_in,
         height = total_h_in, units = "in", device = "pdf")
}
# Canvas is matched to two DIFFERENT existing figures, deliberately -- kept as expressions rather
# than magic numbers so the provenance survives:
#   WIDTH  from plot_modulescore_raw_COMBINED_Micro_Astro.pdf (12.7 x 5 cm for TWO panels, so one
#          panel is 6.35 cm wide)
#   HEIGHT from plots/dist_to_phf1_glia_vs_null/*.pdf (6.6 x 4.5 cm)
PANEL_TOTAL_W_IN <- (12.7 / 2) / 2.54     # 2.500 in = 180 pt
PANEL_TOTAL_H_IN <- 4.5 / 2.54            # 1.772 in = 128 pt

wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE,
                                row.names = FALSE)

# per-donor x ring means, supplementary (not drawn), matching the CosMx _persample.tsv
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

CAVEATS <- function() {
  cat("\nNOTES:\n")
  cat(" - CD68, GFAP, Iba1 and S100b are `type` markers and were inputs to the glial clustering\n")
  cat("   (rowData(spe)$use_channel, 08-phenotyping.Rmd:559-568, 07-batch_correction.Rmd:217-229).\n")
  cat("   stats_mixing_decomposition.txt relates the intensity and cluster-mix panels.\n")
  cat(" - PHF1+ glia are retained in the intensity analysis (586 PHF1+ glia exist object-wide).\n")
  cat(sprintf(" - ROI geometry: cortical strips (median 478 x 2525 um), hence the %d um cap.\n",
              MAX_DIST))
  cat("   Distance is computed WITHIN ROI.\n")
  cat(" - Microglia cluster 2 has <10 cells in 31 of 132 ROIs. The pooled rolling window tolerates\n")
  cat("   this; a per-ROI analysis would not.\n")
  cat(" - Cells are nested in ROIs within donors; the donor level is modelled (n = 44 donors).\n")
}

# ===================================================================
# PANELS 1-2: pooled intensity
# ===================================================================
intensity_panel <- function(slug, cells, value_col, ylab, legend_title, series_label, title) {
  d <- md[cells, ]
  d$value_raw <- d[[value_col]]
  d <- d[is.finite(d$value_raw), ]
  # WITHIN-DONOR Z-SCORE, computed on the analysed population. Donors differ in their distance
  # distributions and in global staining brightness, so without it a between-donor association
  # can dominate the pooled rolling curve. z-scoring per donor removes both the level and the
  # spread differences, so the drawn curve and the tested coefficient describe the same
  # within-donor gradient. Precedent: the per-sample robust z in
  # plot_phf1_intensity_validation.R.
  d <- d %>% group_by(patient_id) %>%
    mutate(value = as.numeric(scale(value_raw))) %>% ungroup()
  d <- d[is.finite(d$value), ]

  # Sex / Age / PMI ARE INCLUDED (Set 3).
  #   The z-score standardises the OUTCOME within donor. It does nothing to the PREDICTOR -- donors
  #   still differ in their distance distributions. With every donor mean ybar_j = 0, the normal
  #   equation for a donor-level covariate A in y ~ d + A is
  #       beta_A = -beta_d * sum_j n_j A_j dbar_j / sum_j n_j A_j^2
  #   so beta_A vanishes only if A is uncorrelated with donor-mean distance. Donors with heavier
  #   tangle burden have shorter distances throughout, so that sum is not zero and the covariates
  #   are retained.
  # Reported by cat_re_variance() below: with the outcome z-scored within donor, the donor
  # random-intercept variance is estimated at zero (~1e-30), so the fit is equivalent to pooled
  # least squares on the standardised outcome. The no-covariate fit is retained as m_nc and both
  # estimates are logged.
  # Cortical depth (depth_s) is not in the headline, matching the CosMx distance models; it is
  # tested as a sensitivity term in R/imc_covariate_sensitivity.R.
  # The headline carries no amyloid term; the with-plaque fit is retained as a reported
  # sensitivity.
  m <- tryCatch(lmerTest::lmer(
    value ~ dist_scaled + Sex + Age_s + PMI_s + (1 | patient_id),
    data = d, REML = TRUE), error = function(e) NULL)
  m_pl <- tryCatch(lmerTest::lmer(
    value ~ dist_scaled + plaque + Sex + Age_s + PMI_s + (1 | patient_id),
    data = d, REML = TRUE), error = function(e) NULL)
  # no-covariate fit, kept as a reported sensitivity so the covariate shift is visible
  m_nc <- tryCatch(lmerTest::lmer(
    value ~ dist_scaled + (1 | patient_id),
    data = d, REML = TRUE), error = function(e) NULL)
  co <- if (is.null(m)) NULL else summary(m)$coefficients
  est <- se <- dfv <- tv <- pv <- NA_real_
  if (!is.null(co) && "dist_scaled" %in% rownames(co)) {
    est <- co["dist_scaled", "Estimate"]; se <- co["dist_scaled", "Std. Error"]
    dfv <- co["dist_scaled", "df"];       tv <- co["dist_scaled", "t value"]
    pv  <- co["dist_scaled", "Pr(>|t|)"]
  }
  co_pl <- if (is.null(m_pl)) NULL else summary(m_pl)$coefficients
  est_pl <- pv_pl <- beta_plq <- pv_plq <- NA_real_
  if (!is.null(co_pl) && "dist_scaled" %in% rownames(co_pl)) {
    est_pl <- co_pl["dist_scaled", "Estimate"]; pv_pl <- co_pl["dist_scaled", "Pr(>|t|)"]
  }
  if (!is.null(co_pl) && "plaqueplaque-prox" %in% rownames(co_pl)) {
    beta_plq <- co_pl["plaqueplaque-prox", "Estimate"]
    pv_plq   <- co_pl["plaqueplaque-prox", "Pr(>|t|)"]
  }
  # pre-covariate estimate, for the shift audit
  co_nc <- if (is.null(m_nc)) NULL else summary(m_nc)$coefficients
  est_nc <- pv_nc <- NA_real_
  if (!is.null(co_nc) && "dist_scaled" %in% rownames(co_nc)) {
    est_nc <- co_nc["dist_scaled", "Estimate"]; pv_nc <- co_nc["dist_scaled", "Pr(>|t|)"]
  }
  # donor-level covariate coefficients: the empirical test of the algebra in the comment above.
  # If these are all ~0 the covariates really were inert here; if not, they were absorbing
  # between-donor structure that the singular random intercept could not.
  cov_rows <- if (is.null(co)) character(0)
              else intersect(rownames(co), grep("^(Sex|Age_s|PMI_s)", rownames(co), value = TRUE))
  # single test per panel -> BH is a no-op; kept explicit so the basis is unambiguous
  padj <- p.adjust(pv, method = "BH")
  sig  <- isTRUE(padj < SIG_ALPHA)

  roll <- blank_sparse(roll_mean(d$dist_um, d$value, GRID, HALF_WINDOW)) %>%
    mutate(dist_um = GRID, series = series_label, sig = sig_factor(sig))
  cols <- setNames(MODULE_COLOURS[1], series_label)

  save_panel(make_panel(roll, ylab, legend_title, cols), slug,
             PANEL_TOTAL_W_IN, PANEL_TOTAL_H_IN)
  wr(roll %>% mutate(window_um = WINDOW_UM, significant = sig) %>%
       select(series, window_um, significant, dist_um, roll_mean, sem, n_window),
     sprintf("source_data_%s_rollmean.tsv", slug))
  wr(ring_table(d, "value_raw", series_label), sprintf("source_data_%s_persample.tsv", slug))

  sink(file.path(out_dir, sprintf("stats_%s.txt", slug)))
  cat(title, "\n"); cat(strrep("=", nchar(title)), "\n", sep = "")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("LMM (log distance):\n")
  cat("  HEADLINE   : value_z ~ dist_to_phf1_um_scaled + Sex + Age_s + PMI_s + (1|patient_id)\n")
  cat("               (NO amyloid term)\n")
  cat("  sensitivity: value_z ~ dist_to_phf1_um_scaled + plaque + Sex + Age_s + PMI_s + (1|patient_id)\n")
  cat("  sensitivity: value_z ~ dist_to_phf1_um_scaled + (1|patient_id)   (no covariates)\n")
  cat(sprintf("    dist_scaled beta  HEADLINE (covariates, no plaque) = %+.5f (p = %.3g)\n", est, pv))
  cat(sprintf("    dist_scaled beta  covariates + plaque              = %+.5f (p = %.3g)\n", est_pl, pv_pl))
  cat(sprintf("    dist_scaled beta  NO covariates                    = %+.5f (p = %.3g)\n", est_nc, pv_nc))
  cat(sprintf("    -> covariate shift in beta: %+.5f (%.1f%% of the no-covariate estimate)\n",
              est - est_nc, if (isTRUE(est_nc != 0)) 100 * (est - est_nc) / abs(est_nc) else NA_real_))
  cat(sprintf("    plaque-prox main effect       = %+.5f (p = %.3g)\n", beta_plq, pv_plq))
  cat("    plaque = Matched_4G8_40, nucleus within a 40 um-dilated 4G8+ plaque mask. If the\n")
  cat("    distance beta survives the amyloid term, the gradient is tau-associated rather than\n")
  cat("    a generic proximity-to-pathology effect.\n")
  cat("  Outcome =", value_col, "from assay 'exprs' = asinh(counts/1), then Z-SCORED WITHIN DONOR\n")
  cat("    over the analysed population, so the slope is in within-donor SD units.\n")
  cat("  WHY: donors differ in their distance distributions and in global staining brightness,\n")
  cat("    so without within-donor z-scoring a between-donor association can dominate the pooled\n")
  cat("    rolling curve. z-scoring per donor removes the between-donor level and spread, so the\n")
  cat("    drawn curve and the tested coefficient describe the same within-donor gradient.\n")
  cat("    Precedent: the per-sample robust z in plot_phf1_intensity_validation.R.\n")
  cat("  Sex/Age/PMI ARE INCLUDED (Set 3). The z-score standardises the OUTCOME, not the\n")
  cat("    PREDICTOR, and donors still differ in their distance distributions. With every donor\n")
  cat("    mean = 0, the normal equation for a donor-level covariate A in y ~ d + A gives\n")
  cat("        beta_A = -beta_d * sum_j n_j A_j dbar_j / sum_j n_j A_j^2\n")
  cat("    so beta_A is zero only if A is uncorrelated with donor-mean distance. Donors with\n")
  cat("    heavier tangle burden have shorter distances throughout, so it is not. The fitted\n")
  cat("    covariate coefficients below are the empirical check on that:\n")
  if (length(cov_rows)) {
    for (r in cov_rows)
      cat(sprintf("      %-12s %+.5f  (p = %.3g)\n", r, co[r, "Estimate"], co[r, "Pr(>|t|)"]))
  } else cat("      (none estimable)\n")
  cat("    A coefficient far from 0 means the covariate is absorbing between-donor structure.\n")
  cat("    With the outcome z-scored within donor, the donor random-intercept variance is\n")
  cat("    estimated at zero (see the random-effects block below), so the fit is equivalent to\n")
  cat("    pooled least squares on the standardised outcome.\n")
  cat("  dist_to_phf1_um_scaled = log(dist)/sd(log(dist)), plain log, /SD, NOT centred\n")
  cat(sprintf("    (canonical transform; dist_sd = %.4f, computed after the cap).\n", DIST_SD))
  cat("  Cortical depth: not in the headline, so that every IMC distance model matches the CosMx\n")
  cat("    distance models; it is tested as a sensitivity term in R/imc_covariate_sensitivity.R.\n")
  cat(sprintf("  Distance cap: cells with 0 < dist_to_phf1_um <= %d um.\n", MAX_DIST))
  cat(sprintf("  Cells modelled: %d | donors: %d | ROIs: %d\n",
              nrow(d), dplyr::n_distinct(d$patient_id), dplyr::n_distinct(d$sample_id)))
  cat("\n")
  cat("Rolling curve: sliding-window mean, window WIDTH", WINDOW_UM, "um (half-window",
      HALF_WINDOW, "um),\n")
  cat("  grid:", N_GRID, "points over 0 -", MAX_DIST, "um, POOLED across donors.\n")
  cat("  Ribbon = mean +/- 1.96 * SEM, SEM = sd(window)/sqrt(n_window) -- a CELL-level SEM.\n")
  cat(sprintf("  Windows with n_window < %d are blanked (innermost windows are truncated at 0).\n\n",
              MIN_WINDOW_N))
  cat("Significance basis: BH-adjusted model p; padj <", SIG_ALPHA, "-> solid+thick line.\n")
  cat("  Single test in this panel, so padj == p.\n\n")
  cat("=== RESULT (per s.d. of log distance; POSITIVE = higher FURTHER from a tangle) ===\n")
  cat(sprintf("  estimate %+.5f  SE %.5f  df %.1f  t %+.3f  p %.4g  padj %.4g  -> %s\n",
              est, se, dfv, tv, pv, padj, if (sig) "significant" else "n.s."))
  cat(sprintf("  95%% CI [%+.5f, %+.5f]\n", est - 1.96 * se, est + 1.96 * se))
  cat(sprintf("  Re-expressed per s.d. CLOSER to a tangle: %+.5f\n", -est))
  cat_re_variance(m, "Random effects (donor), HEADLINE model")
  if (!is.null(m)) { cat("\n=== Full model ===\n"); print(summary(m)) }
  CAVEATS()
  cat("\n=== sessionInfo() ===\n"); print(sessionInfo())
  sink()

  tibble(slug = slug, series = series_label, estimate = est, SE = se, df = dfv,
         t = tv, pval = pv, padj = padj, significant = sig, n_cells = nrow(d))
}

res_cd68 <- intensity_panel(
  "cd68_microglia_distance", which(md$is_micro), "CD68",
  "CD68 (within-donor z)", "Celltype", "Microglia",
  "CD68 in pooled microglia vs distance to nearest PHF1+ neuron (IMC, AD cases)")

res_gfap <- intensity_panel(
  "gfap_astrocytes_distance", which(md$is_astro), "GFAP",
  "GFAP (within-donor z)", "Celltype", "Astrocytes",
  "GFAP in pooled astrocytes vs distance to nearest PHF1+ neuron (IMC, AD cases)")

# ===================================================================
# PANELS 3-4: per-cluster proportion of ALL cells
# ===================================================================
# Denominator is the LINEAGE, not all cells: "what % of microglia are cluster 1 vs 2 over
# distance". With exactly two clusters the curves are complements summing to 1, so this is ONE
# test (the mixing ratio), not two -- the reciprocal model would give OR_1 = 1/OR_2 and the same
# p. The single model is fitted on P(second cluster | lineage cell) and its verdict is applied to
# both drawn lines.
proportion_panel <- function(slug, clusters, legend_title, denom_label, title) {
  d <- md %>% filter(celltype_clusters %in% clusters) %>%
    mutate(is_second = as.integer(celltype_clusters == clusters[2]))

  # no amyloid term in the headline; the with-plaque fit is retained as a reported
  # sensitivity, matching the intensity panels above.
  m <- tryCatch(lme4::glmer(
    is_second ~ dist_scaled + Sex + Age_s + PMI_s + (1 | patient_id),
    data = d, family = binomial,
    control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
    error = function(e) NULL)
  m_pl <- tryCatch(lme4::glmer(
    is_second ~ dist_scaled + plaque + Sex + Age_s + PMI_s + (1 | patient_id),
    data = d, family = binomial,
    control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
    error = function(e) NULL)
  co_pl <- if (is.null(m_pl)) NULL else summary(m_pl)$coefficients
  est_pl <- pv_pl <- NA_real_
  if (!is.null(co_pl) && "dist_scaled" %in% rownames(co_pl)) {
    est_pl <- co_pl["dist_scaled", "Estimate"]; pv_pl <- co_pl["dist_scaled", "Pr(>|z|)"]
  }
  co <- if (is.null(m)) NULL else summary(m)$coefficients
  est <- se <- zv <- pv <- NA_real_
  if (!is.null(co) && "dist_scaled" %in% rownames(co)) {
    est <- co["dist_scaled", "Estimate"]; se <- co["dist_scaled", "Std. Error"]
    zv  <- co["dist_scaled", "z value"];  pv <- co["dist_scaled", "Pr(>|z|)"]
  }
  padj <- p.adjust(pv, method = "BH")     # single test in the panel -> no-op, kept explicit
  sig  <- isTRUE(padj < SIG_ALPHA)

  st <- tibble(
    lineage = denom_label, reference = clusters[1], modelled = clusters[2],
    estimate = est, SE = se, z = zv, pval = pv, padj = padj, significant = sig,
    OR = exp(est), OR_lo = exp(est - 1.96 * se), OR_hi = exp(est + 1.96 * se),
    n_cells = nrow(d),
    n_ref = sum(d$is_second == 0), n_modelled = sum(d$is_second == 1),
    singular = if (is.null(m)) NA else lme4::isSingular(m))

  roll <- bind_rows(lapply(clusters, function(cl) {
    ind <- as.integer(d$celltype_clusters == cl)
    blank_sparse(roll_mean(d$dist_um, ind, GRID, HALF_WINDOW)) %>%
      mutate(dist_um = GRID, series = cl)
  })) %>%
    mutate(series = factor(series, levels = clusters), sig = sig_factor(sig))
  cols <- setNames(MODULE_COLOURS[seq_along(clusters)], clusters)

  save_panel(make_panel(roll, paste0("Proportion of ", denom_label), legend_title, cols), slug)
  wr(roll %>% mutate(window_um = WINDOW_UM, significant = sig == SIG_BASIS) %>%
       select(series, window_um, significant, dist_um, roll_mean, sem, n_window),
     sprintf("source_data_%s_rollmean.tsv", slug))
  wr(bind_rows(lapply(clusters, function(cl)
       ring_table(d %>% mutate(.v = as.integer(celltype_clusters == cl)), ".v", cl))),
     sprintf("source_data_%s_persample.tsv", slug))

  sink(file.path(out_dir, sprintf("stats_%s.txt", slug)))
  cat(title, "\n"); cat(strrep("=", nchar(title)), "\n", sep = "")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("Mixed-effects logistic (the default for compositional questions):\n")
  cat("  HEADLINE   : is_", clusters[2],
      " ~ dist_to_phf1_um_scaled + Sex + Age_s + PMI_s + (1|patient_id)\n", sep = "")
  cat("  sensitivity: same WITH a plaque term.\n")
  cat(sprintf("    OR per s.d. of log distance,  HEADLINE (no plaque) = %.4f (p = %.3g)\n",
              exp(est), pv))
  cat(sprintf("    OR per s.d. of log distance,  with plaque          = %.4f (p = %.3g)\n",
              exp(est_pl), pv_pl))
  cat("  Cortical depth is not in the headline, matching the CosMx distance models.\n")
  cat("  DENOMINATOR = ", denom_label, " only, i.e. 'what % of ", denom_label,
      " are each cluster'.\n", sep = "")
  cat("  This is the MIXING RATIO within the lineage, so it is NOT affected by shifts in other\n")
  cat("  lineages -- unlike a proportion-of-all-cells denominator, which is zero-sum against\n")
  cat("  neurons (and the anchor is itself a neuron).\n")
  cat("  ONE TEST, not two: with two clusters the drawn curves are complements summing to 1, so\n")
  cat("  the reciprocal model would return OR = 1/OR and an identical p. Reference cluster =\n")
  cat("  ", clusters[1], "; modelled = ", clusters[2], ".\n", sep = "")
  cat(sprintf("  dist_sd = %.4f. Cap: 0 < dist <= %d um.\n", DIST_SD, MAX_DIST))
  cat(sprintf("  Cells: %d (%s %d / %s %d) | donors: %d | ROIs: %d\n", nrow(d),
              clusters[1], st$n_ref, clusters[2], st$n_modelled,
              dplyr::n_distinct(d$patient_id), dplyr::n_distinct(d$sample_id)))
  cat("\n")
  cat("Rolling curve: window WIDTH", WINDOW_UM, "um, grid", N_GRID, "points, pooled across\n")
  cat("  donors; ribbon = mean +/- 1.96 * SEM of the 0/1 membership indicator; windows with\n")
  cat(sprintf("  n_window < %d blanked.\n\n", MIN_WINDOW_N))
  cat("=== RESULT (per s.d. of log distance) ===\n")
  cat(sprintf("  OR %.4f [%.4f, %.4f]  z %+.3f  p %.4g  padj %.4g  -> %s\n",
              st$OR, st$OR_lo, st$OR_hi, zv, pv, padj, if (sig) "significant" else "n.s."))
  cat("  OR > 1 = ", clusters[2], " makes up a LARGER share of ", denom_label,
      " FURTHER from a tangle.\n", sep = "")
  cat_re_variance(m, "Random effects (donor), HEADLINE model")
  if (!is.null(m)) { cat("\n=== Full model ===\n"); print(summary(m)) }
  cat("\n=== Lineage composition of ALL cells by distance bin (context) ===\n")
  cat("Not the denominator used above; included so a lineage-level shift is visible. The anchor\n")
  cat("is itself a neuron.\n")
  print(as.data.frame(COMPOSITION), row.names = FALSE, digits = 3)
  CAVEATS()
  cat("\n=== sessionInfo() ===\n"); print(sessionInfo())
  sink()
  st %>% mutate(slug = slug)
}

res_micro_prop <- proportion_panel("microglia_cluster_proportion_distance", MICRO,
  "Microglia", "microglia",
  "Microglial cluster mix vs distance to nearest PHF1+ neuron (IMC, AD cases)")
res_astro_prop <- proportion_panel("astrocyte_cluster_proportion_distance", ASTRO,
  "Astrocytes", "astrocytes",
  "Astrocyte cluster mix vs distance to nearest PHF1+ neuron (IMC, AD cases)")

# ===================================================================
# MIXING DECOMPOSITION: how much of the pooled marker gradient reflects the cluster mix?
# The clusters were defined using these markers (see CLUSTERING INPUTS), so a pooled-intensity
# gradient can include a contribution from a shift in cluster mix. This estimates it:
#   d(pooled marker)/d(s.d. log dist) predicted by mixing alone = dp1 * (mean1 - mean2)
# where dp1 is the shift in cluster-1 share implied by the fitted mixing OR.
# ===================================================================
decompose <- function(marker, clusters, lineage, int_res, prop_res) {
  v  <- md[[marker]]
  i1 <- md$celltype_clusters == clusters[1]; i2 <- md$celltype_clusters == clusters[2]
  m1 <- mean(v[i1], na.rm = TRUE); m2 <- mean(v[i2], na.rm = TRUE); gap <- m1 - m2
  p2 <- mean(md$celltype_clusters[i1 | i2] == clusters[2])
  lor <- log(prop_res$OR); dp2 <- lor * p2 * (1 - p2); dp1 <- -dp2
  # int_res$estimate is now on the WITHIN-DONOR Z scale while `gap` is on the arcsinh scale, so
  # convert the gap to donor-z units using the mean within-donor SD of the marker.
  sd_within <- md %>% filter(celltype_clusters %in% clusters) %>% group_by(patient_id) %>%
    summarise(s = sd(.data[[marker]]), .groups = "drop") %>% pull(s) %>% mean(na.rm = TRUE)
  pred <- dp1 * gap / sd_within; obs <- int_res$estimate
  tibble(lineage = lineage, marker = marker,
         mean_cluster1 = m1, mean_cluster2 = m2, gap = gap,
         mean_within_donor_sd = sd_within, gap_in_donor_z = gap / sd_within, share_cluster2 = p2,
         mix_OR = prop_res$OR, d_share_cluster1_per_sd = dp1,
         predicted_pooled_slope = pred, observed_pooled_slope = obs,
         pct_explained_by_mixing = 100 * pred / obs,
         # the ratio is only interpretable when there is a gradient and a mix shift to decompose
         interpretable = int_res$significant & prop_res$significant,
         pooled_sig = int_res$significant, mix_sig = prop_res$significant)
}
dec <- bind_rows(
  decompose("CD68", MICRO, "microglia",  res_cd68, res_micro_prop),
  decompose("GFAP", ASTRO, "astrocytes", res_gfap, res_astro_prop)
)
write.table(dec, file.path(out_dir, "stats_mixing_decomposition.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(out_dir, "stats_mixing_decomposition.txt"))
cat("How much of the pooled marker gradient reflects a change in cluster mix?\n")
cat("========================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("The glial clusters were defined using these markers (rowData(spe)$use_channel = the 22 type\n")
cat("markers, which include CD68 and GFAP), so a pooled-intensity gradient over distance can\n")
cat("include a contribution from a shift in the lineage's cluster mix. This decomposition\n")
cat("estimates that contribution.\n\n")
cat("Method: the fitted mixing OR implies a shift dp1 in the cluster-1 share per s.d. of log\n")
cat("distance (dp1 = -log(OR) * p2 * (1 - p2) at the observed share p2). Multiplying dp1 by the\n")
cat("between-cluster marker gap (mean1 - mean2) gives the pooled-marker slope expected from\n")
cat("mixing ALONE. Comparing that with the observed pooled slope gives the share explained.\n\n")
print(as.data.frame(dec), row.names = FALSE, digits = 4)
cat("\nINTERPRETATION\n")
for (i in seq_len(nrow(dec))) {
  if (dec$interpretable[i]) {
    cat(sprintf("  %s / %s: mixing explains %.0f%% of the observed pooled gradient.\n",
                dec$lineage[i], dec$marker[i], dec$pct_explained_by_mixing[i]))
  } else {
    cat(sprintf("  %s / %s: no decomposition -- pooled gradient %s, mix shift %s.\n",
                dec$lineage[i], dec$marker[i],
                if (dec$pooled_sig[i]) "significant" else "n.s.",
                if (dec$mix_sig[i]) "significant" else "n.s."))
    cat("    The percentage requires both a pooled gradient and a mix shift.\n")
  }
}
cat("\n=== sessionInfo() ===\n"); print(sessionInfo())
sink()

cat("\n-- mixing decomposition --\n")
print(as.data.frame(dec %>% select(lineage, marker, mix_OR, predicted_pooled_slope,
                                   observed_pooled_slope, pct_explained_by_mixing)),
      row.names = FALSE, digits = 3)

message(sprintf("Done. %d cells, %d donors, %d ROIs, dist_sd = %.4f",
                nrow(md), dplyr::n_distinct(md$patient_id),
                dplyr::n_distinct(md$sample_id), DIST_SD))
cat("\n-- intensity panels (per s.d. of LOG distance; negative = higher NEAR tangles) --\n")
print(as.data.frame(bind_rows(res_cd68, res_gfap) %>%
                      select(series, estimate, pval, padj, significant, n_cells)),
      row.names = FALSE, digits = 4)
cat("\n-- cluster-mix panels: OR per s.d. of log distance for the MODELLED cluster's share of\n   its own lineage; OR < 1 = that cluster is a larger share NEAR tangles --\n")
print(as.data.frame(bind_rows(res_micro_prop, res_astro_prop) %>%
                      select(lineage, modelled, OR, OR_lo, OR_hi, padj, significant, n_cells)),
      row.names = FALSE, digits = 4)
