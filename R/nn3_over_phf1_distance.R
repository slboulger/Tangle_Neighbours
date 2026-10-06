#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# nn3_over_phf1_distance.R
#
# Figure panels: S3 (all) - second step of run_nn3.sh
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# nn3_over_phf1_distance.R
#
# ANALYSIS B (NON-CELL-AUTONOMOUS): among PHF1- neurons, does 3-nearest-neighbour
# spacing vary with DISTANCE TO THE NEAREST PHF1+ NEURON, and does that differ
# between neuronal subtypes?
#
# This is the spatial-gradient counterpart to R/nn3_phf1_vs_neg.R (which asks the
# cell-autonomous question, PHF1+ vs PHF1-). Same metric, same cells, same
# neighbour pool -- it reads the cache that script writes, so the two analyses
# describe exactly the same neurons and cannot drift.
#
# The Zwang et al. 2024 prediction, if tangles drive local neuron loss, is that
# spacing should be WIDER near a tangle. See the geometry section below.
#
# ---------------------------------------------------------------------------
# GEOMETRY AND THE PERMUTATION NULL
# ---------------------------------------------------------------------------
# dist_to_phf1_um is itself a nearest-neighbour distance to a SUBSET OF NEURONS,
# so a cell in a sparse neighbourhood is far from the nearest PHF1+ neuron by
# construction: spacing and distance-to-tangle share geometry
# (rho(nn3_um, dist_to_phf1_um) is printed in the stats log). The coefficient is
# therefore calibrated against the PHF1 label-permutation null in
# R/null_replica_nn3.R (docs/MODELS.md, Permutation nulls): nn3_um is invariant to
# PHF1 relabelling, so the shuffled distance columns give the geometric baseline.
#
# ---------------------------------------------------------------------------
# MODEL. Mirrors the canonical linear-distance spec in
# deg_dream_linear_distance_phf1.r where it applies:
#   dist_scaled = log(dist_to_phf1_um) / sd(log(dist_to_phf1_um))
# natural log (NEVER log1p, per the canonical distance transform), divided
# by SD, NOT centred, and dist_sd computed within each subtype exactly as the DEG
# script computes it per celltype object. Distance is strictly > 0 for every
# modelled PHF1-negative cell, so no offset is needed and a non-positive value is
# a data error -- stop() on it rather than absorbing it.
#
#   log(nn3_um) ~ dist_scaled + edge_dist_um + Sex + Age_s + PMI_s + (1|sample_id)   HEADLINE
#               + depth_rel                                                          reported
#
# COVARIATES. The outcome is image-derived, so the specification is Set 2 plus
# edge_dist_um (docs/MODELS.md, Documented exceptions), both stated in the stats log:
#   * nUMI_log and percent_neg are not included. They are library-size and
#     negative-probe QC terms for COUNT outcomes. nn3_um is a purely geometric
#     measurement with no library-size dependence, so including them would add
#     noise without controlling anything.
#   * edge_dist_um is ADDED as the measure-specific technical covariate (FOV edge
#     censoring inflates nn3_um; it is adjusted for rather than excluded).
#
# SIGN CONVENTION: beta_closer = -beta_dist, so POSITIVE means spacing is WIDER
# NEARER a tangle (the Zwang et al. direction).
#
# Outputs, under plots/nn3_neuron_spacing/:
#   plot_nn3_over_distance_by_subtype.pdf     rolling-mean curves, one per subtype
#   plot_nn3_over_distance_faceted.pdf        rolling mean + fitted model + Wald CI
#   plot_nn3_distance_coef_by_subtype.pdf     coefficient forest
#   source_data_nn3_over_distance_roll.tsv    exact drawn rolling-mean rows
#   source_data_nn3_over_distance_fit.tsv     exact drawn fitted-curve rows
#   source_data_nn3_over_distance_rings.tsv   per-sample x subtype x ring means
#   source_data_nn3_distance_coef.tsv
#   stats_nn3_over_distance.txt
#   stats_nn3_distance_coef_effectsize.tsv
#   stats_nn3_over_distance_ring_rolling_check.tsv
#
# Requires R/nn3_phf1_vs_neg.R to have been run first (writes the cache).
# Run: Rscript R/nn3_over_phf1_distance.R    (~1 min)

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
  library(lmerTest)
})

hpc <- "<PROJECT_ROOT>/phf1_v2"
loc <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(hpc)) hpc else loc)
source("R/palettes.R")
source("R/nn3_utils.r")

out_dir <- "plots/nn3_neuron_spacing"
cache   <- file.path(out_dir, "source_data_nn3_cells.tsv")
if (!file.exists(cache))
  stop("Cache not found: ", cache, "\nRun R/nn3_phf1_vs_neg.R first.")

MAX_DIST_UM <- 1000    # canonical cap, matches run_deg_linear_distance.sh
WINDOW_UM   <- 75      # rolling-window WIDTH, matches plot_modulescore_vs_phf1_distance_modelp.R
N_GRID      <- 200
MIN_CELLS   <- 50L     # canonical guards from deg_dream_linear_distance_phf1.r
MIN_DONORS  <- 3L
FDR         <- 0.05
RING_BREAKS <- c(0, 50, 100, 200, 300, 500, 700, 1000)

wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE)

# -------------------------------------------------------------------
# Load the cache and build the modelling frame
# -------------------------------------------------------------------
md <- read.delim(cache, sep = "\t", header = TRUE, check.names = FALSE,
                 colClasses = c(Braak = "character", Sex = "character",
                                sample_id = "character", celltype = "character"))

d0 <- md %>%
  filter(!phf1_pos, !is.na(dist_to_phf1_um), dist_to_phf1_um <= MAX_DIST_UM,
         !is.na(nn3_um)) %>%
  mutate(celltype = factor(celltype, levels = c(neuron_order, "Unassigned Neuron")),
         Braak    = factor(Braak, levels = braak_levels),
         Sex      = factor(Sex),
         log_nn3  = log(nn3_um),
         Age_s    = as.numeric(scale(as.numeric(Age))),
         PMI_s    = as.numeric(scale(as.numeric(PMI))))
if (any(d0$dist_to_phf1_um <= 0))
  stop("Non-positive dist_to_phf1_um encountered; log() requires dist > 0.")
stopifnot(nrow(d0) > 0, !any(is.na(d0$celltype)))

subtypes <- levels(d0$celltype)

# -------------------------------------------------------------------
# Rolling-window mean, verbatim from plot_modulescore_vs_phf1_distance_modelp.R
# so the descriptive overlay is drawn the same way as the existing distance
# figures. Model-free; the SEM is over cells and is descriptive only.
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

grid_um <- seq(0, MAX_DIST_UM, length.out = N_GRID)
half    <- WINDOW_UM / 2

# -------------------------------------------------------------------
# Per-subtype fits
# -------------------------------------------------------------------
fit_one <- function(ct) {
  d <- d0 %>% filter(celltype == ct)
  base <- tibble(celltype = ct, n_cells = nrow(d),
                 n_donors = dplyr::n_distinct(d$sample_id),
                 fitted = FALSE, dist_sd = NA_real_,
                 beta_closer = NA_real_, SE = NA_real_, df = NA_real_,
                 CI.L = NA_real_, CI.R = NA_real_, pval = NA_real_,
                 cohens_d = NA_real_, beta_closer_depth_adj = NA_real_,
                 pval_depth_adj = NA_real_, beta_closer_unadj = NA_real_,
                 pval_unadj = NA_real_, beta_closer_nn3rd = NA_real_,
                 pval_nn3rd = NA_real_, singular = NA,
                 geo_mean_nn3_um = if (nrow(d)) exp(mean(d$log_nn3)) else NA_real_,
                 rho_nn3_dist = NA_real_, rho_nn3_local = NA_real_)
  if (nrow(d) < MIN_CELLS || base$n_donors < MIN_DONORS) return(base)

  # dist_sd computed WITHIN subtype, as deg_dream_linear_distance_phf1.r does per
  # celltype object. Divide by SD only -- not centred.
  dt <- log(d$dist_to_phf1_um)
  ds <- stats::sd(dt)
  d$dist_scaled <- dt / ds
  base$dist_sd <- ds
  base$rho_nn3_dist  <- suppressWarnings(cor(d$nn3_um, d$dist_to_phf1_um,
                                             method = "spearman", use = "complete.obs"))
  base$rho_nn3_local <- suppressWarnings(cor(d$nn3_um, d$n_local,
                                             method = "spearman", use = "complete.obs"))

  fitm <- function(f) tryCatch(lmerTest::lmer(f, data = d, REML = TRUE),
                               error = function(e) NULL)
  grab <- function(m, term = "dist_scaled") {
    if (is.null(m)) return(rep(NA_real_, 5))
    co <- summary(m)$coefficients
    if (!term %in% rownames(co)) return(rep(NA_real_, 5))
    c(co[term, "Estimate"], co[term, "Std. Error"], co[term, "df"],
      co[term, "t value"], co[term, "Pr(>|t|)"])
  }

  m_head  <- fitm(log_nn3 ~ dist_scaled + edge_dist_um + Sex + Age_s + PMI_s + (1 | sample_id))
  m_depth <- fitm(log_nn3 ~ dist_scaled + edge_dist_um + Sex + Age_s + PMI_s + depth_rel + (1 | sample_id))
  m_unadj <- fitm(log_nn3 ~ dist_scaled + (1 | sample_id))
  m_3rd   <- fitm(log(nn3rd_um) ~ dist_scaled + edge_dist_um + Sex + Age_s + PMI_s + (1 | sample_id))
  if (is.null(m_head)) return(base)

  a  <- grab(m_head); ad <- grab(m_depth); au <- grab(m_unadj); a3 <- grab(m_3rd)
  crit <- stats::qt(0.975, ifelse(is.na(a[3]), Inf, a[3]))
  sd_tot <- { vc <- as.data.frame(lme4::VarCorr(m_head)); sqrt(sum(vc$vcov[is.na(vc$var2)])) }

  base$fitted      <- TRUE
  base$beta_closer <- -a[1]; base$SE <- a[2]; base$df <- a[3]; base$pval <- a[5]
  base$CI.L <- -a[1] - crit * a[2]; base$CI.R <- -a[1] + crit * a[2]
  base$cohens_d <- -a[1] / sd_tot
  base$beta_closer_depth_adj <- -ad[1]; base$pval_depth_adj <- ad[5]
  base$beta_closer_unadj     <- -au[1]; base$pval_unadj     <- au[5]
  base$beta_closer_nn3rd     <- -a3[1]; base$pval_nn3rd     <- a3[5]
  # A singular fit means the sample_id variance is estimated at 0. Flagged rather
  # than dropped, and listed in the stats log.
  base$singular <- lme4::isSingular(m_head)
  attr(base, "model") <- m_head
  attr(base, "data")  <- d
  base
}

res_list <- lapply(subtypes, fit_one)
res <- bind_rows(res_list) %>% mutate(padj = NA_real_)
# BH over fitted subtypes only; p.adjust()'s default n counts NAs and would over-correct.
res$padj[res$fitted] <- p.adjust(res$pval[res$fitted], method = "BH")
res <- res %>%
  mutate(sig = factor(ifelse(!is.na(padj) & padj < FDR, "padj<0.05", "ns"),
                      levels = c("padj<0.05", "ns")),
         # A per-SD log change back-transformed to a % change in spacing.
         pct_closer   = 100 * (exp(beta_closer) - 1),
         pct_closer_L = 100 * (exp(CI.L) - 1),
         pct_closer_R = 100 * (exp(CI.R) - 1),
         # Same coefficient expressed per 100 um nearer at the median distance,
         # which is what a reader can actually picture.
         um_per_sd_closer = geo_mean_nn3_um * (exp(beta_closer) - 1))
wr(res %>% select(-sig), "source_data_nn3_distance_coef.tsv")
wr(res %>% select(celltype, n_cells, n_donors, dist_sd, beta_closer, SE, CI.L, CI.R,
                  cohens_d, pval, padj, pct_closer, pct_closer_L, pct_closer_R,
                  um_per_sd_closer, geo_mean_nn3_um, beta_closer_unadj, pval_unadj,
                  beta_closer_depth_adj, pval_depth_adj, beta_closer_nn3rd, pval_nn3rd,
                  rho_nn3_dist, rho_nn3_local, singular),
   "stats_nn3_distance_coef_effectsize.tsv")

# -------------------------------------------------------------------
# Descriptive rolling means (drawn) and per-sample ring means (source data)
# -------------------------------------------------------------------
roll_df <- bind_rows(lapply(subtypes, function(ct) {
  d <- d0 %>% filter(celltype == ct)
  if (nrow(d) < MIN_CELLS) return(NULL)
  rm_ <- roll_mean(d$dist_to_phf1_um, d$nn3_um, grid_um, half)
  tibble(celltype = ct, dist_to_phf1_um = grid_um,
         roll_mean = rm_$roll_mean, sem = rm_$sem, n_window = rm_$n_window)
})) %>%
  filter(n_window >= 10) %>%      # do not draw a "mean" over fewer than 10 cells
  left_join(res %>% select(celltype, sig, padj), by = "celltype") %>%
  mutate(celltype = factor(celltype, levels = subtypes))
wr(roll_df, "source_data_nn3_over_distance_roll.tsv")

.rl <- head(RING_BREAKS, -1); .ru <- tail(RING_BREAKS, -1)
RING_LAB <- sprintf("%g-%g", .rl, .ru)
RING_MID <- setNames((.rl + .ru) / 2, RING_LAB)
rings <- d0 %>%
  mutate(dist_ring = cut(dist_to_phf1_um, breaks = RING_BREAKS, labels = RING_LAB,
                         include.lowest = TRUE, right = TRUE)) %>%
  group_by(celltype, sample_id, Braak, dist_ring) %>%
  summarise(mean_nn3_um = mean(nn3_um), sd_nn3_um = sd(nn3_um),
            median_nn3_um = median(nn3_um), n_cells = dplyr::n(), .groups = "drop") %>%
  mutate(ring_mid_um = RING_MID[as.character(dist_ring)])
wr(rings, "source_data_nn3_over_distance_rings.tsv")

# Agreement check between the two descriptive summaries: the cell-weighted ring
# mean must land within 1.96 SEM of the rolling mean evaluated at the ring midpoint.
# Only the innermost ring is allowed to disagree -- there the 75 um window is
# truncated at distance 0 and so covers a different span than the 0-50 um ring.
ring_check <- rings %>%
  group_by(celltype, ring_mid_um) %>%
  summarise(ring_cell_mean = sum(mean_nn3_um * n_cells) / sum(n_cells),
            n_cells = sum(n_cells), .groups = "drop") %>%
  rowwise() %>%
  mutate(rl = { r <- roll_df[roll_df$celltype == celltype, ]
                if (!nrow(r)) NA_integer_ else which.min(abs(r$dist_to_phf1_um - ring_mid_um)) },
         roll_at_mid = { r <- roll_df[roll_df$celltype == celltype, ]
                         if (is.na(rl)) NA_real_ else r$roll_mean[rl] },
         roll_sem    = { r <- roll_df[roll_df$celltype == celltype, ]
                         if (is.na(rl)) NA_real_ else r$sem[rl] }) %>%
  ungroup() %>%
  mutate(abs_diff = abs(ring_cell_mean - roll_at_mid),
         within_ci = abs_diff <= 1.96 * roll_sem)
ring_bad <- ring_check %>% filter(!is.na(within_ci), !within_ci, ring_mid_um > min(RING_MID))
if (nrow(ring_bad))
  warning(sprintf("Ring and rolling means disagree beyond 1.96 SEM in %d non-innermost bin(s); see stats log.",
                  nrow(ring_bad)))
wr(ring_check %>% select(-rl), "stats_nn3_over_distance_ring_rolling_check.tsv")

# -------------------------------------------------------------------
# Fitted model curves. 200-point Wald-CI grid, covariates at their means and the
# reference Sex level, random effect excluded -- the predict_grid() pattern from
# plot_modulescore_vs_phf1_distance_modelp.R. x stays RAW microns for plotting;
# the grid starts at min(dist) because log() is undefined at 0.
# -------------------------------------------------------------------
fit_df <- bind_rows(lapply(seq_along(subtypes), function(i) {
  r <- res_list[[i]]
  if (!isTRUE(r$fitted)) return(NULL)
  m <- attr(r, "model"); d <- attr(r, "data")
  g  <- seq(min(d$dist_to_phf1_um), MAX_DIST_UM, length.out = N_GRID)
  nd <- data.frame(dist_scaled = log(g) / r$dist_sd,
                   edge_dist_um = mean(d$edge_dist_um),
                   Sex = factor(levels(d$Sex)[1], levels = levels(d$Sex)),
                   Age_s = mean(d$Age_s), PMI_s = mean(d$PMI_s))
  X <- model.matrix(~ dist_scaled + edge_dist_um + Sex + Age_s + PMI_s, nd)
  beta <- lme4::fixef(m)
  X <- X[, names(beta), drop = FALSE]
  ft <- as.numeric(X %*% beta)
  V  <- as.matrix(vcov(m)); se <- sqrt(rowSums((X %*% V) * X))
  tibble(celltype = subtypes[i], dist_to_phf1_um = g,
         fitted_nn3_um = exp(ft), ci_lo = exp(ft - 1.96 * se), ci_hi = exp(ft + 1.96 * se))
})) %>%
  left_join(res %>% select(celltype, sig, padj), by = "celltype") %>%
  mutate(celltype = factor(celltype, levels = subtypes))
wr(fit_df, "source_data_nn3_over_distance_fit.tsv")

# -------------------------------------------------------------------
# Figures
# -------------------------------------------------------------------
XLAB    <- expression("Distance to PHF1+ neuron (" * mu * "m)")
# Two-line y label: the single-line version is wider than these panels are tall and
# gets clipped. mu via plotmath so it survives the default pdf device.
YLAB    <- expression(atop("Mean distance to 3", "nearest neurons (" * mu * "m)"))
sig_lv  <- c("padj<0.05", "ns")
lt_vals <- setNames(c("solid", "dashed"), sig_lv)
lw_vals <- setNames(c(0.8, 0.4), sig_lv)

# 1. Overlay of the model-free rolling means, one line per subtype.
# NO ribbons here: 10 overlapping +/-1.96 SEM bands render as an unreadable wash.
# The per-subtype ribbons are in the faceted figure below, where they are legible.
p1 <- ggplot(roll_df, aes(dist_to_phf1_um, roll_mean, colour = celltype, group = celltype)) +
  geom_line(aes(linetype = sig, linewidth = sig)) +
  scale_colour_manual(values = celltype_palette, name = NULL, drop = FALSE) +
  scale_linetype_manual(values = lt_vals, name = NULL, drop = FALSE, limits = sig_lv) +
  scale_linewidth_manual(values = lw_vals, name = NULL, drop = FALSE, limits = sig_lv) +
  guides(colour = guide_legend(order = 1),
         linetype = guide_legend(order = 2), linewidth = guide_legend(order = 2)) +
  labs(x = XLAB, y = YLAB) +
  coord_cartesian(xlim = c(0, MAX_DIST_UM)) +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.key.width = grid::unit(20, "pt"),
        legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
# Panel pinned to a constant physical width so the x-axis aligns with the other
# distance figures in this project regardless of legend size.
g1 <- ggplot2::ggplotGrob(p1)
pc <- unique(g1$layout$l[grepl("^panel", g1$layout$name)])
g1$widths[pc] <- grid::unit(1.6, "in")
ggsave(file.path(out_dir, "plot_nn3_over_distance_by_subtype.pdf"), g1,
       width = grid::convertWidth(sum(g1$widths), "in", valueOnly = TRUE),
       height = 2.6, units = "in", device = "pdf")

# 2. Faceted: rolling mean (grey) with the fitted model and its Wald CI on top.
p2 <- ggplot(fit_df, aes(dist_to_phf1_um)) +
  geom_line(data = roll_df, aes(y = roll_mean), colour = "grey60", linewidth = 0.3) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = celltype), alpha = 0.20, colour = NA) +
  geom_line(aes(y = fitted_nn3_um, colour = celltype, linetype = sig, linewidth = sig)) +
  facet_wrap(~ celltype, ncol = 3) +
  scale_colour_manual(values = celltype_palette, guide = "none") +
  scale_fill_manual(values = celltype_palette, guide = "none") +
  scale_linetype_manual(values = lt_vals, name = NULL, drop = FALSE, limits = sig_lv) +
  scale_linewidth_manual(values = lw_vals, name = NULL, drop = FALSE, limits = sig_lv) +
  labs(x = XLAB, y = YLAB) +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        strip.background = element_blank(), strip.text = element_text(size = 5.5),
        legend.position = "top", legend.key.width = grid::unit(20, "pt"),
        plot.margin = margin(4, 4, 4, 4, unit = "pt"))
ggsave(file.path(out_dir, "plot_nn3_over_distance_faceted.pdf"), p2,
       width = 12.0, height = 10.0, units = "cm", device = "pdf")

# 3. Coefficient forest. Positive = spacing WIDER nearer a tangle (Zwang direction).
fp <- res %>% filter(fitted) %>%
  mutate(celltype = factor(celltype, levels = rev(subtypes)),
         sig_lab = ifelse(!is.na(padj) & padj < FDR, "padj < 0.05", "n.s."))
p3 <- ggplot(fp, aes(pct_closer, celltype, colour = sig_lab)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
  geom_errorbar(aes(xmin = pct_closer_L, xmax = pct_closer_R), orientation = "y",
                width = 0, linewidth = 0.4) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = c("padj < 0.05" = "#BD0026", "n.s." = "grey65"), name = NULL) +
  labs(x = "Change in spacing per s.d. closer\nto a PHF1+ neuron (%)", y = NULL) +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "top",
        plot.margin = margin(4, 10, 4, 4, unit = "pt"))
ggsave(file.path(out_dir, "plot_nn3_distance_coef_by_subtype.pdf"), p3,
       width = 9.6, height = 5.8, units = "cm", device = "pdf")

# -------------------------------------------------------------------
# Stats log
# -------------------------------------------------------------------
rho_all <- cor(d0$nn3_um, d0$dist_to_phf1_um, method = "spearman", use = "complete.obs")

sink(file.path(out_dir, "stats_nn3_over_distance.txt"))
cat("3-NN neuronal spacing vs distance to the nearest PHF1+ neuron, by neuronal subtype\n")
cat("=================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Cells:  PHF1-negative neurons only, dist_to_phf1_um <=", MAX_DIST_UM, "um.\n")
cat("Source: ", cache, " (written by R/nn3_phf1_vs_neg.R), so this analysis and the\n", sep = "")
cat("        cell-autonomous one describe exactly the same neurons and metric.\n")
cat("Metric: nn3_um = mean Euclidean distance to the 3 nearest other neurons, um,\n")
cat("        over the 10-label neuron pool (see R/nn3_utils.r).\n\n")
cat("MODEL (per subtype; outcome log(nn3_um))\n")
cat("  HEADLINE : log(nn3_um) ~ dist_scaled + edge_dist_um + Sex + Age_s + PMI_s + (1|sample_id)\n")
cat("  depth-adj: + depth_rel                                     (reported, not headline)\n")
cat("  unadj    : ~ dist_scaled + (1|sample_id)                   (reported)\n")
cat("  nn3rd    : outcome log(nn3rd_um), the 3rd-nearest reading of the metric\n")
cat("  dist_scaled = log(dist_to_phf1_um) / sd(log(dist_to_phf1_um)), natural log (never log1p),\n")
cat("  divided by SD, NOT centred; dist_sd computed WITHIN each subtype. This mirrors\n")
cat("  deg_dream_linear_distance_phf1.r exactly, per the canonical distance transform.\n")
cat("  COVARIATES (Set 2 + edge_dist_um): nUMI_log and percent_neg are not included\n")
cat("  (library-size / negative-probe QC terms for COUNT outcomes; nn3_um is geometric and has\n")
cat("  no library-size dependence), and edge_dist_um is added (FOV edge censoring inflates\n")
cat("  nn3_um; it is adjusted for rather than excluded).\n")
cat(sprintf("  BH across the %d fitted subtypes; significance padj < %.2f.\n",
            sum(res$fitted), FDR))
cat("  SIGN: beta_closer = -beta_dist, so POSITIVE = spacing WIDER nearer a tangle, which is\n")
cat("  the direction Zwang et al. 2024 would predict if tangles drive local neuron loss.\n\n")
cat("CELLS\n")
cat(sprintf("  PHF1- neurons in window: %d | donors: %d\n",
            nrow(d0), dplyr::n_distinct(d0$sample_id)))
cat("  dist_to_phf1_um:\n"); print(summary(d0$dist_to_phf1_um))
cat("  nn3_um:\n"); print(summary(d0$nn3_um))
cat("  per subtype:\n")
print(as.data.frame(res %>% select(celltype, n_cells, n_donors, fitted, dist_sd,
                                   geo_mean_nn3_um)),
      row.names = FALSE, digits = 4)
cat(sprintf("\n  Subtypes not fitted (< %d cells or < %d donors): %s\n", MIN_CELLS, MIN_DONORS,
            paste(res$celltype[!res$fitted], collapse = ", ")))

cat("\n*** GEOMETRY ***\n")
cat("dist_to_phf1_um is ITSELF a nearest-neighbour distance to a subset of neurons, so a neuron\n")
cat("in a sparse neighbourhood is far from the nearest PHF1+ neuron by construction; spacing and\n")
cat("distance-to-tangle share geometry.\n")
cat(sprintf("  rho(nn3_um, dist_to_phf1_um) over all PHF1- neurons in the window = %+.3f\n", rho_all))
cat("  per-subtype rho values are in the results table (rho_nn3_dist).\n")
cat("The coefficient is calibrated against PHF1 label permutation in R/null_replica_nn3.R:\n")
cat("nn3_um is invariant to PHF1 relabelling, so the 1000 shuffled distance columns in\n")
cat("deg/null_label_deg/engine/ give the geometric baseline directly.\n")

cat("\n=== RESULTS (headline model, sorted by p) ===\n")
print(as.data.frame(res %>% filter(fitted) %>% arrange(pval) %>%
        select(celltype, n_cells, beta_closer, CI.L, CI.R, pct_closer, pct_closer_L,
               pct_closer_R, um_per_sd_closer, cohens_d, pval, padj)),
      row.names = FALSE, digits = 3)
cat("\n--- Sensitivity / robustness (same units as beta_closer) ---\n")
print(as.data.frame(res %>% filter(fitted) %>% arrange(pval) %>%
        select(celltype, beta_closer, beta_closer_unadj, pval_unadj,
               beta_closer_depth_adj, pval_depth_adj, beta_closer_nn3rd, pval_nn3rd,
               rho_nn3_dist, rho_nn3_local, singular)),
      row.names = FALSE, digits = 3)
sing <- res$celltype[which(res$singular)]
if (length(sing)) {
  cat("\nSINGULAR FITS (sample_id variance estimated at 0):\n  ")
  cat(paste(sing, collapse = "\n  "), "\n")
}

cat("\n--- Ring vs rolling-mean agreement check ---\n")
cat("The cell-weighted ring mean should land within 1.96 SEM of the rolling mean at the ring\n")
cat("midpoint. The innermost ring is exempt: there the 75 um window is truncated at distance 0\n")
cat("and covers a different span than the 0-50 um ring.\n")
cat(sprintf("  bins checked: %d | outside 1.96 SEM: %d | outside, excluding the innermost ring: %d\n",
            sum(!is.na(ring_check$within_ci)), sum(!ring_check$within_ci, na.rm = TRUE),
            nrow(ring_bad)))
if (nrow(ring_bad)) print(as.data.frame(ring_bad %>% select(-rl)), row.names = FALSE, digits = 4)

cat("\n=== Effect sizes ===\n")
cat("Unstandardised effect with 95% CI = beta_closer (log um per s.d. of log distance), also\n")
cat("given as a % change in spacing and in microns at each subtype's geometric mean. Standardised\n")
cat("measure = cohens_d (beta / sqrt(sum of VarCorr variances)). Per subtype:\n")
for (i in which(res$fitted)) {
  cat(sprintf("  %-26s %+6.2f%% [%+6.2f, %+6.2f]  %+5.2f um  d = %+.3f  padj = %.3g\n",
              res$celltype[i], res$pct_closer[i], res$pct_closer_L[i], res$pct_closer_R[i],
              res$um_per_sd_closer[i], res$cohens_d[i], res$padj[i]))
}
cat("\nDistance is a WITHIN-donor predictor, so the donor random intercept is the correct\n")
cat("specification; the permutation null (R/null_replica_nn3.R) calibrates the coefficient.\n")

cat("\nNOTES:\n")
cat(" - Rolling-mean ribbons are +/-1.96 SEM over CELLS and are descriptive only.\n")
cat(sprintf(" - Rolling means are suppressed where a window holds < 10 cells, so curves can stop\n"))
cat("   short of the 1000 um cap for sparse subtypes. Ring means with their n are in\n")
cat("   source_data_nn3_over_distance_rings.tsv.\n")
cat(" - Sex, Age and PMI are donor-level covariates included alongside a donor random\n")
cat("   intercept, for consistency with the canonical model.\n")
cat("\n=== sessionInfo() ===\n"); print(sessionInfo())
sink()

message(sprintf("Done. %d PHF1- neurons, %d subtypes fitted, %d with padj < %.2f. rho(nn3, dist) = %+.3f",
                nrow(d0), sum(res$fitted),
                sum(res$padj < FDR, na.rm = TRUE), FDR, rho_all))
print(as.data.frame(res %>% filter(fitted) %>% arrange(pval) %>%
        select(celltype, pct_closer, pct_closer_L, pct_closer_R, padj, rho_nn3_dist)),
      row.names = FALSE, digits = 3)
