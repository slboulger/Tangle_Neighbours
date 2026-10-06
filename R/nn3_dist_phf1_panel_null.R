#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# nn3_dist_phf1_panel_null.R
#
# Figure panels: S3 (all), S3A, S3B
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# nn3_dist_phf1_panel_null.R
#
# THE DISTANCE PANEL DRAWN RELATIVE TO THE GEOMETRIC NULL. Exc-IT-L2-3-CBLN2-HOPX.
#
# Same grammar as R/nn3_dist_phf1_panel.R (75 um rolling mean, PHF1+ marginal strip, line
# style = distance verdict, point shape = contrast verdict) with the y axis changed from raw
# spacing to DEVIATION FROM THE PERMUTATION NULL. Reads the cache written by
# R/nn3_null_curves.R and fits NOTHING.
#
# WHY. dist_to_phf1_um is itself a nearest-neighbour distance over the same point pattern nn3
# measures, so the raw curve largely reflects that geometry. The assembled figure is
# therefore (a) the raw-scale OVERLAY -- observed and null curves together in microns --
# (b) the same comparison as a deviation, and (c) the coefficient against its null
# distribution. The bare raw curve is still written as plot_nn3_dist_phf1_raw.pdf for the
# companion panel in R/nn3_dist_phf1_panel.R, but it is not a panel of the composite.
#
# ---------------------------------------------------------------------------
# THREE DESIGN POINTS
# ---------------------------------------------------------------------------
# 1. UNITS ARE MICRONS, NOT PER CENT. A per-cent y axis would be confused with the per-cent
#    excess of the coefficient, which is a difference of SLOPES over a tenfold distance span.
#    The curve's height is a difference of LEVELS. Two different quantities must not wear the
#    same units in the same figure.
#
# 2. A CONSTANT OFFSET IS NOT A GRADIENT. Following from (1): a curve sitting at a flat
#    +1 um everywhere would have an excess coefficient of zero. The interpretable feature is
#    the RISE toward short distances, so the far-field level (>= 700 um) is drawn as its own
#    dashed rule and the reader is pointed at the rise above it, not at the height above zero.
#
# 3. THE POINTWISE BAND IS NOT A TEST. Over 200 grid points some excursion is expected. The
#    pointwise 95% band is drawn pale and is descriptive; the GLOBAL studentised-sup envelope
#    (k * null_sd, k chosen so 95% of null curves lie entirely inside) is drawn as dashed
#    bounds and carries the p-value. Only the global one is quoted.
#
# ---------------------------------------------------------------------------
# LEVEL OF THE STATISTICS
# ---------------------------------------------------------------------------
# The global envelope test, like the coefficient excess, is a pooled permutation statistic:
# the curve pools 9 donors and each labelling yields one pooled curve. The per-donor excess,
# each donor against its own null, is written as its own file (plot_nn3_perdonor_excess.pdf).
#
# Run: Rscript R/nn3_dist_phf1_panel_null.R [cap_um]     (seconds; nothing is fitted)

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(qs); library(patchwork)
})

hpc <- "<PROJECT_ROOT>/phf1_v2"
loc <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(hpc)) hpc else loc)
source("R/palettes.R")
source("R/nn3_null_hist_panel.R")

.cli <- commandArgs(trailingOnly = TRUE)
MAX_DIST <- if (length(.cli) >= 1) as.numeric(.cli[1]) else 1000
FDR <- 0.05
OBS_COL   <- "#BD0026"
BAND_PT   <- "grey86"    # pointwise null band  -- descriptive
BAND_GL   <- "grey55"    # global envelope      -- the test
SIG_BASIS <- "padj<0.05"
SIG_LEGEND <- "Log-distance"
CON_LEGEND <- "PHF1+ vs PHF1-"
CON_SHAPES <- c(16, 1)
sig_lv <- c(SIG_BASIS, "ns")

in_dir  <- "plots/nn3_null_deviation"
cur_f   <- file.path(in_dir, sprintf("null_curves_%dum.qs", MAX_DIST))
if (!file.exists(cur_f)) stop("Missing: ", cur_f, "\nRun R/nn3_null_curves.R first.")
cu <- qs::qread(cur_f)
CELLTYPE <- cu$celltype

out_dir <- sprintf("plots/nn3_dist_phf1_panel_null_%dum", MAX_DIST)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE)

XLAB <- expression("Distance to PHF1+ neuron (" * mu * "m)")
YLAB_A <- expression(atop("Mean distance to 3", "nearest neurons (" * mu * "m)"))
YLAB_B <- expression(atop("Spacing minus geometric", "null (" * mu * "m)"))

## ---------------------------------------------------------------------------
## Panel (a) -- the raw curve, for the SEM ribbon only. Same cells, same window.
## ---------------------------------------------------------------------------
md <- read.delim("plots/nn3_neuron_spacing/source_data_nn3_cells.tsv", sep = "\t",
                 check.names = FALSE, colClasses = c(sample_id = "character",
                 celltype = "character", Sex = "character", Braak = "character")) %>%
  filter(celltype == CELLTYPE, !is.na(nn3_um)) %>% mutate(phf1_pos = as.logical(phf1_pos))
mdn <- md %>% filter(!phf1_pos, !is.na(dist_to_phf1_um), dist_to_phf1_um <= MAX_DIST)
hw <- cu$window_um / 2
sem_a <- vapply(cu$grid, function(g) {
  i <- which(mdn$dist_to_phf1_um >= g - hw & mdn$dist_to_phf1_um <= g + hw)
  if (length(i) < 2) NA_real_ else stats::sd(mdn$nn3_um[i]) / sqrt(length(i))
}, numeric(1))

# The distance verdict on panel (a) is the PUBLISHED cell-level one, unchanged, so (a) is the
# same figure as R/nn3_dist_phf1_panel.R and the pair is directly comparable.
pubc <- read.delim("plots/nn3_neuron_spacing/source_data_nn3_distance_coef.tsv", sep = "\t") %>%
  filter(celltype == CELLTYPE)
pubk <- read.delim("plots/nn3_neuron_spacing/source_data_nn3_phf1_vs_neg_by_celltype.tsv",
                   sep = "\t") %>% filter(celltype == CELLTYPE)
sig_a <- isTRUE(pubc$padj < FDR); con_a <- isTRUE(pubk$padj < FDR)

sem_strip <- stats::sd(md$nn3_um[md$phf1_pos]) / sqrt(sum(md$phf1_pos))
A_roll <- tibble(dist = cu$grid, mean = cu$obs_mean, sem = sem_a, n = cu$obs_n) %>%
  filter(!is.na(mean)) %>%
  mutate(sig = factor(ifelse(sig_a, SIG_BASIS, "ns"), levels = sig_lv))
A_ref <- tibble(x_pos = -0.08 * MAX_DIST, mean = cu$strip_obs,
                lo = cu$strip_obs - 1.96 * sem_strip, hi = cu$strip_obs + 1.96 * sem_strip,
                con_sig = factor(ifelse(con_a, SIG_BASIS, "ns"), levels = sig_lv))

## ---------------------------------------------------------------------------
## Panel (b) -- deviation from the null
## ---------------------------------------------------------------------------
sig_b <- cu$p_global < FDR
strip_dev_lo <- unname(quantile(cu$strip_null - mean(cu$strip_null), 0.025))
strip_dev_hi <- unname(quantile(cu$strip_null - mean(cu$strip_null), 0.975))
con_b <- cu$p_strip < FDR

B_roll <- tibble(dist = cu$grid, dev = cu$dev_obs,
                 pt_lo = cu$null_q025, pt_hi = cu$null_q975,
                 gl_lo = -cu$env_k * cu$null_sd, gl_hi = cu$env_k * cu$null_sd) %>%
  filter(!is.na(dev)) %>%
  mutate(sig = factor(ifelse(sig_b, SIG_BASIS, "ns"), levels = sig_lv))
B_ref <- tibble(x_pos = -0.08 * MAX_DIST, dev = cu$strip_dev,
                lo = strip_dev_lo, hi = strip_dev_hi,
                con_sig = factor(ifelse(con_b, SIG_BASIS, "ns"), levels = sig_lv))

## ---------------------------------------------------------------------------
## Shared panel furniture
## ---------------------------------------------------------------------------
brk <- pretty(c(0, MAX_DIST)); brk <- brk[brk >= 0 & brk <= MAX_DIST]
XLO <- -0.14 * MAX_DIST

base_layers <- function(p, ref, ylo_lab = TRUE) {
  p +
    annotate("segment", x = -0.025 * MAX_DIST, xend = -0.025 * MAX_DIST,
             y = -Inf, yend = Inf, linewidth = 0.25, colour = "grey70") +
    scale_shape_manual(values = setNames(CON_SHAPES, sig_lv), name = CON_LEGEND,
                       drop = FALSE, limits = sig_lv) +
    scale_linetype_manual(values = setNames(c("solid", "dashed"), sig_lv),
                          name = SIG_LEGEND, drop = FALSE, limits = sig_lv) +
    scale_linewidth_manual(values = setNames(c(0.8, 0.4), sig_lv),
                           name = SIG_LEGEND, drop = FALSE, limits = sig_lv) +
    scale_x_continuous(breaks = brk) +
    guides(linetype = guide_legend(order = 1), linewidth = guide_legend(order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(colour = "grey25", size = 1.3))) +
    coord_cartesian(xlim = c(XLO, MAX_DIST), clip = "off") +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "right", legend.direction = "vertical", legend.box = "vertical",
          legend.text = element_text(size = 6), legend.title = element_text(size = 6),
          legend.key.width = grid::unit(16, "pt"), legend.key.height = grid::unit(7, "pt"),
          legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
          legend.spacing.y = grid::unit(1, "pt"),
          plot.margin = margin(t = 14, r = 6, b = 4, l = 5))
}

pa <- ggplot(pad_levels(A_roll, "sig", sig_lv), aes(dist, mean, group = sig)) +
  geom_ribbon(aes(ymin = mean - 1.96 * sem, ymax = mean + 1.96 * sem),
              fill = OBS_COL, alpha = 0.15, colour = NA) +
  geom_line(aes(linetype = sig, linewidth = sig), colour = OBS_COL) +
  geom_errorbar(data = pad_levels(A_ref, "con_sig", sig_lv, keep = "x_pos"),
                aes(x = x_pos, y = mean, ymin = lo, ymax = hi), inherit.aes = FALSE,
                width = 0, linewidth = 0.4, colour = OBS_COL) +
  geom_point(data = pad_levels(A_ref, "con_sig", sig_lv, keep = "x_pos"),
             aes(x = x_pos, y = mean, shape = con_sig), inherit.aes = FALSE,
             size = 1.3, stroke = 0.5, colour = OBS_COL) +
  annotate("text", x = A_ref$x_pos[1], y = Inf, vjust = -0.45, size = 2.6,
           fontface = "bold", colour = "grey15", label = "PHF1+") +
  labs(x = XLAB, y = YLAB_A)
pa <- base_layers(pa, A_ref)

pb <- ggplot(pad_levels(B_roll, "sig", sig_lv), aes(dist, dev, group = sig)) +
  # pointwise band: pale, descriptive only
  geom_ribbon(aes(ymin = pt_lo, ymax = pt_hi), fill = BAND_PT, colour = NA) +
  # global envelope: the actual test
  geom_line(aes(y = gl_lo), linetype = 2, linewidth = 0.3, colour = BAND_GL) +
  geom_line(aes(y = gl_hi), linetype = 2, linewidth = 0.3, colour = BAND_GL) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3, colour = "grey35") +
  # the far-field level: the line the RISE is read against (see header note 2)
  geom_hline(yintercept = cu$far_dev, linetype = 2, linewidth = 0.25, colour = OBS_COL,
             alpha = 0.6) +
  geom_line(aes(linetype = sig, linewidth = sig), colour = OBS_COL) +
  geom_errorbar(data = pad_levels(B_ref, "con_sig", sig_lv, keep = "x_pos"),
                aes(x = x_pos, y = dev, ymin = lo, ymax = hi), inherit.aes = FALSE,
                width = 0, linewidth = 0.4, colour = OBS_COL) +
  geom_point(data = pad_levels(B_ref, "con_sig", sig_lv, keep = "x_pos"),
             aes(x = x_pos, y = dev, shape = con_sig), inherit.aes = FALSE,
             size = 1.3, stroke = 0.5, colour = OBS_COL) +
  annotate("text", x = B_ref$x_pos[1], y = Inf, vjust = -0.45, size = 2.6,
           fontface = "bold", colour = "grey15", label = "PHF1+") +
  labs(x = XLAB, y = YLAB_B)
pb <- base_layers(pb, B_ref)

## ---------------------------------------------------------------------------
## Panel (b2) -- the same comparison on the RAW axis: both curves in microns.
##
## (b) plots observed MINUS null, which is the quantity the test is on but throws away the
## values. (b2) keeps the actual spacing on the y axis and puts the null's own curve and
## envelope behind it, on the same axis as the observed curve. Same numbers, and the cache
## stores the null quantiles AS DEVIATIONS, so the raw band is null_mean + those.
## ---------------------------------------------------------------------------
B2_roll <- tibble(dist = cu$grid, obs = cu$obs_mean, null = cu$null_mean,
                  pt_lo = cu$null_mean + cu$null_q025, pt_hi = cu$null_mean + cu$null_q975,
                  gl_lo = cu$null_mean - cu$env_k * cu$null_sd,
                  gl_hi = cu$null_mean + cu$env_k * cu$null_sd) %>%
  filter(!is.na(obs), !is.na(null)) %>%
  mutate(sig = factor(ifelse(sig_b, SIG_BASIS, "ns"), levels = sig_lv))
B2_ref <- tibble(x_pos = -0.08 * MAX_DIST, obs = cu$strip_obs, null = mean(cu$strip_null),
                 lo = cu$strip_obs - 1.96 * sem_strip, hi = cu$strip_obs + 1.96 * sem_strip,
                 con_sig = factor(ifelse(con_b, SIG_BASIS, "ns"), levels = sig_lv))

pb2 <- ggplot(pad_levels(B2_roll, "sig", sig_lv), aes(dist, obs, group = sig)) +
  geom_ribbon(aes(ymin = pt_lo, ymax = pt_hi), fill = BAND_PT, colour = NA) +
  geom_line(aes(y = gl_lo), linetype = 2, linewidth = 0.3, colour = BAND_GL) +
  geom_line(aes(y = gl_hi), linetype = 2, linewidth = 0.3, colour = BAND_GL) +
  geom_line(aes(y = null), linewidth = 0.5, colour = "grey30", linetype = 2) +
  geom_line(aes(linetype = sig, linewidth = sig), colour = OBS_COL) +
  # both strips: observed in red, the null PHF1+ mean in grey beside it
  geom_errorbar(data = pad_levels(B2_ref, "con_sig", sig_lv, keep = "x_pos"),
                aes(x = x_pos, y = obs, ymin = lo, ymax = hi), inherit.aes = FALSE,
                width = 0, linewidth = 0.4, colour = OBS_COL) +
  geom_point(data = pad_levels(B2_ref, "con_sig", sig_lv, keep = "x_pos"),
             aes(x = x_pos, y = obs, shape = con_sig), inherit.aes = FALSE,
             size = 1.3, stroke = 0.5, colour = OBS_COL) +
  geom_point(data = B2_ref, aes(x = x_pos, y = null), inherit.aes = FALSE,
             shape = 95, size = 2.6, colour = "grey30") +
  annotate("text", x = B2_ref$x_pos[1], y = Inf, vjust = -0.45, size = 2.6,
           fontface = "bold", colour = "grey15", label = "PHF1+") +
  labs(x = XLAB, y = YLAB_A)
pb2 <- base_layers(pb2, B2_ref)

## ---------------------------------------------------------------------------
## Third panel of the composite -- the coefficient against its null distribution.
## Built by the SHARED builder in R/nn3_null_hist_panel.R so this and the standalone
## plot_nn3_null_deviation_CBLN2.pdf cannot drift apart.
## ---------------------------------------------------------------------------
.perm <- read.delim(file.path(in_dir, "source_data_nn3_null_deviation_CBLN2.tsv"), sep = "\t")
.sm   <- read.delim("plots/nn3_neuron_spacing/stats_nn3_null_replica_summary.tsv", sep = "\t")
.sm   <- .sm[.sm$analysis == "B_distance" & .sm$celltype == CELLTYPE, ]
stopifnot(nrow(.perm) > 0, nrow(.sm) == 1L)
phist <- nn3_null_hist(beta_null = .perm$beta_closer_null,
                       beta_obs = .sm$beta_obs_fast, null_median = .sm$null_median,
                       null_q025 = .sm$null_q025, null_q975 = .sm$null_q975,
                       dist_sd = pubc$dist_sd, obs_col = OBS_COL, null_col = "grey78")

## ---------------------------------------------------------------------------
## The per-donor panel -- written as its own file (plot_nn3_perdonor_excess.pdf). It is not
## in the composite; the composite is the null comparison end to end.
## ---------------------------------------------------------------------------
pd <- cu$per_donor %>%
  mutate(Braak = factor(as.character(Braak), levels = braak_levels),
         sample_id = factor(sample_id, levels = rev(sample_id[order(excess)])))
mt <- cu$meta
pc <- ggplot(pd, aes(excess, sample_id)) +
  # each donor against its OWN null spread, recentred on that donor's null median
  geom_errorbar(aes(xmin = null_q025 - null_median, xmax = null_q975 - null_median),
                orientation = "y", width = 0, linewidth = 2.0, colour = "grey86") +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.3, colour = "grey35") +
  geom_vline(xintercept = mt$mean_excess, linetype = 2, linewidth = 0.4, colour = OBS_COL) +
  annotate("rect", xmin = mt$CI.L, xmax = mt$CI.R, ymin = -Inf, ymax = Inf,
           fill = OBS_COL, alpha = 0.10) +
  geom_point(aes(fill = Braak), shape = 21, size = 1.9, stroke = 0.3, colour = "grey20") +
  scale_fill_manual(values = braak_palette, name = "Braak", drop = FALSE) +
  labs(x = expression("Excess over the geometric null (log " * mu * "m)"), y = NULL) +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "right", legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.height = grid::unit(7, "pt"),
        plot.margin = margin(4, 8, 4, 4, unit = "pt"))

## ---------------------------------------------------------------------------
## Save: each panel alone (project standard) and the assembled figure
## ---------------------------------------------------------------------------
sv <- function(p, f, w, h) ggsave(file.path(out_dir, f), p, width = w, height = h,
                                  units = "cm", device = "pdf")
sv(pa, "plot_nn3_dist_phf1_raw.pdf",      10.0, 6.0)
sv(pb, "plot_nn3_dist_phf1_nullrel.pdf",  10.0, 6.0)
sv(pb2, "plot_nn3_dist_phf1_nulloverlay.pdf", 10.0, 6.0)
sv(pc, "plot_nn3_perdonor_excess.pdf",     9.6, 6.0)
# Composite: raw-scale overlay, then the deviation, then the coefficient against its null.
# The bare raw curve is written separately but is not a panel here. Panel tags are NOT
# added: they are applied downstream.
sv(pb2 / pb / phist, "plot_nn3_null_relative_figure.pdf", 10.0, 17.0)

## ---------------------------------------------------------------------------
## Source data + stats log
## ---------------------------------------------------------------------------
wr(tibble(dist_to_phf1_um = cu$grid, obs_mean = cu$obs_mean, obs_sem = sem_a,
          n_window = cu$obs_n, null_mean = cu$null_mean, null_sd = cu$null_sd,
          dev_obs = cu$dev_obs, pointwise_lo = cu$null_q025, pointwise_hi = cu$null_q975,
          global_lo = -cu$env_k * cu$null_sd, global_hi = cu$env_k * cu$null_sd,
          # the same bands on the RAW axis, as drawn in the overlay panel
          raw_pointwise_lo = cu$null_mean + cu$null_q025,
          raw_pointwise_hi = cu$null_mean + cu$null_q975,
          raw_global_lo = cu$null_mean - cu$env_k * cu$null_sd,
          raw_global_hi = cu$null_mean + cu$env_k * cu$null_sd),
   "source_data_nn3_null_relative_curve.tsv")
wr(tibble(panel = c("a", "b"), what = c("raw", "null-relative"),
          strip_value = c(cu$strip_obs, cu$strip_dev),
          lo = c(A_ref$lo, strip_dev_lo), hi = c(A_ref$hi, strip_dev_hi),
          basis = c("cell-level SEM", "null 2.5-97.5% of PHF1+ means"),
          p = c(pubk$padj, cu$p_strip)),
   "source_data_nn3_null_relative_strip.tsv")
wr(cu$per_donor, "source_data_nn3_perdonor_excess.tsv")

sink(file.path(out_dir, "stats_nn3_null_relative.txt"))
ttl <- sprintf("nn3 over PHF1 distance, relative to the geometric null -- %s (%d um)",
               CELLTYPE, MAX_DIST)
cat(ttl, "\n"); cat(strrep("=", nchar(ttl)), "\n"); cat("Date:", format(Sys.time()), "\n\n")
cat("NOTHING IS FITTED HERE. All values come from ", cur_f, ",\n", sep = "")
cat("written by R/nn3_null_curves.R, plus the published coefficient tables for panel (a)'s\n")
cat("line style and point shape.\n\n")
cat("ASSEMBLED FIGURE plot_nn3_null_relative_figure.pdf\n")
cat("PANEL (a) RAW-SCALE OVERLAY (plot_nn3_dist_phf1_nulloverlay.pdf): observed (red) and the\n")
cat("   mean of the 1000 null rolling means (grey dashed), both in microns, with the two\n")
cat("   bands drawn around the NULL curve.\n")
cat("PANEL (b) the same comparison as a DEVIATION: observed minus null mean, in um. Same\n")
cat("   numbers as (a) on a magnified axis; this is the scale the test is computed on.\n")
cat("   pale band  pointwise 2.5-97.5% of the null curves. DESCRIPTIVE ONLY -- over 200 grid\n")
cat("              points some excursion is expected by chance.\n")
cat(sprintf("   dashed     global studentised-sup envelope, k = %.3f null s.d. THIS carries the test.\n",
            cu$env_k))
cat(sprintf("   red dashed far-field level (>= %d um) = %+.3f um. A CONSTANT OFFSET IS NOT A\n",
            cu$far_from, cu$far_dev))
cat("              GRADIENT: read the RISE above this line, not the height above zero.\n")
cat("THIRD PANEL the 1000 null coefficients with the observed value marked (shared builder\n")
cat("   R/nn3_null_hist_panel.R; identical to plot_nn3_null_deviation_CBLN2.pdf).\n")
cat("The per-donor excess is written on its own as plot_nn3_perdonor_excess.pdf.\n")
cat("The bare raw curve (plot_nn3_dist_phf1_raw.pdf) is written but is not a panel of the\n")
cat("assembled figure.\n\n")

cat("=== PANEL (b): GLOBAL ENVELOPE TEST (cell-level) ===\n")
cat(sprintf("  sup_obs %.3f | null 95%% %.3f | p = %.4g -> line drawn %s\n",
            cu$sup_obs, cu$env_k, cu$p_global, if (sig_b) "solid" else "dashed"))
# WHERE the curve leaves the envelope, and in WHICH DIRECTION, reported alongside the
# global p.
.ok  <- is.finite(B_roll$dev)
.out <- .ok & (B_roll$dev > B_roll$gl_hi | B_roll$dev < B_roll$gl_lo)
if (any(.out)) {
  .r <- range(B_roll$dist[.out]); .i <- which(.out)[which.max(abs(B_roll$dev[.out]))]
  cat(sprintf("  EXCURSION: %d of %d grid points outside, %.0f-%.0f um, %s the null;\n",
              sum(.out), sum(.ok), .r[1], .r[2],
              paste(unique(ifelse(B_roll$dev[.out] > 0, "ABOVE", "BELOW")), collapse = "/")))
  cat(sprintf("             most extreme %+.2f um at %.0f um.\n",
              B_roll$dev[.i], B_roll$dist[.i]))
} else cat("  No grid point leaves the global envelope.\n")
cat(sprintf("  PHF1+ strip: %+.3f um [%+.3f, %+.3f], p = %.4g -> point drawn %s\n",
            cu$strip_dev, strip_dev_lo, strip_dev_hi, cu$p_strip,
            if (con_b) "filled" else "hollow"))

cat("\n=== PER-DONOR EXCESS (plot_nn3_perdonor_excess.pdf) ===\n")
print(as.data.frame(cu$per_donor %>%
        select(sample_id, Braak, n_cells, n_phf1_pos, beta_obs, null_median, excess, p_emp)),
      row.names = FALSE, digits = 4)
cat(sprintf("\n  mean excess %+.5f [%+.5f, %+.5f], t(%d) = %.2f, p = %.4g\n",
            mt$mean_excess, mt$CI.L, mt$CI.R, mt$df, mt$t, mt$p_t))
cat(sprintf("  %d of %d donors positive (sign test p = %.3f); cell-weighted excess %+.5f\n",
            mt$n_positive, mt$n_donors, mt$p_sign, mt$cell_weighted_excess))

cat("\n=== Pooled vs per-donor excess ===\n")
cat(sprintf("  pooled            : excess %+.5f, BH-adjusted empirical p = %.4g  (published)\n",
            mt$pooled_excess_published, mt$pooled_padj_emp))
cat(sprintf("  donor-level       : excess %+.5f [%+.5f, %+.5f], p = %.3f\n",
            mt$mean_excess, mt$CI.L, mt$CI.R, mt$p_t))

cat("\nNOTES:\n")
cat(" - The permutation does not balance edge_dist_um between labelled and unlabelled sets;\n")
cat("   the regression adjusts for it, the drawn curves do not.\n")
cat(" - Empirical p uses (1+r)/(1+n), bounded below by 2/1001.\n")
cat("\n=== sessionInfo() ===\n"); print(sessionInfo())
sink()

cat("Wrote", length(list.files(out_dir)), "files to", out_dir, "\n")
cat(sprintf("  global envelope p = %.4g | strip %+.3f um p = %.4g | donor-level excess %+.5f p = %.3f\n",
            cu$p_global, cu$strip_dev, cu$p_strip, mt$mean_excess, mt$p_t))
