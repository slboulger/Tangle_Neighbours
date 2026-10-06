#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# imc_cosmx_gfap_distance_panel.R
#
# Figure panels: 6C
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# imc_cosmx_gfap_distance_panel.R
#
# Dual-axis panel: IMC GFAP protein (left axis) and CosMx GFAP transcripts (right axis) in
# astrocytes vs distance to the nearest PHF1+ neuron.
#
# NO NEW MODEL IS FITTED HERE. Each line takes its significance from its canonical fit:
#   IMC    : plots/imc_phf1_glia_distance/ (R/imc_phf1_glia_distance.R, Set 3, 300 um
#            cap).
#            The drawn curve is read verbatim from that script's rollmean source data, so it is
#            the published curve, not a recomputation.
#   CosMx  : deg/de_linear_distance/Astro/ (R/deg_dream_linear_distance_phf1.r,
#            Set 1, 1000 um cap) -- the GFAP row of the Astro distance-DEG table.
# The x-axis therefore runs to 1000 um and each line is drawn ONLY over the window its model was
# fitted on (IMC stops at 300 um), so each line carries the verdict of the model fitted on that
# window.
#
# CosMx curve = GFAP counts per cell, CENTRED WITHIN DONOR and re-offset to the cohort mean.
# Donors contribute unequally along the distance axis, so the curve is centred within donor to
# match the model's donor random intercept (the same reason the IMC side uses a within-donor z).
# Donor-centred library size is flat over the same range, so the curve is not a library-size
# artefact. Units are kept (transcripts per cell) so the right axis is readable as abundance.
#
# The two y-axes are joined by a linear map chosen so the two drawn curves (with ribbons) span
# the same vertical extent. The relative vertical scaling of the two lines is therefore
# ARBITRARY: compare direction and shape across modalities, never height.

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(qs); library(dplyr); library(tibble); library(ggplot2)
})

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")

# The IMC curve is read verbatim from imc_phf1_glia_distance.R's output, so run that first.
IMC_DIR   <- "plots/imc_phf1_glia_distance"
IMC_ROLL  <- file.path(IMC_DIR, "source_data_gfap_astrocytes_distance_rollmean.tsv")
IMC_STATS <- file.path(IMC_DIR, "stats_gfap_astrocytes_distance.txt")
COSMX_SCE <- "celltype_sce_neighbours/Astro_sce_neighbours.qs"
COSMX_DEG <- "deg/de_linear_distance/Astro/Astro_dist_to_phf1_um_scaled.tsv"

out_dir <- "plots/imc_cosmx_gfap_distance_panel"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
slug <- "gfap_imc_cosmx_distance"

SIG_ALPHA    <- 0.05
SIG_BASIS    <- "padj<0.05"
COSMX_CAP    <- 1000   # canonical dream cap (--max_dist_um default)
WINDOW_UM    <- 75     # rolling-mean window width, as the IMC panel and the module-score figures
HALF_WINDOW  <- WINDOW_UM / 2
N_GRID       <- 200
MIN_WINDOW_N <- 10
GENE         <- "GFAP"

# Colourblind-safe pair. IMC GFAP protein = process cyan, the hue of
# the GFAP channel in the representative IMC ROI (plots/imc_gfap_images/, drawn pure #00FFFF,
# which is illegible as a line on white). CosMx GFAP transcripts = Okabe-Ito vermillion.
# Cyan vs vermillion stays separable under protan/deutan/tritan vision. Axis TITLES carry the
# series colour; axis numbers and the linetype key stay black (the key applies to both lines).
COL_IMC   <- "#00AEEF"
COL_COSMX <- "#D55E00"
FIG_W_CM  <- 7
FIG_H_CM  <- 4.5

# -------------------------------------------------------------------
# IMC: published curve + headline result
# -------------------------------------------------------------------
imc <- read.delim(IMC_ROLL, check.names = FALSE)
stopifnot(all(c("dist_um", "roll_mean", "sem", "n_window", "significant") %in% names(imc)))
imc_sig <- unique(as.character(imc$significant))
stopifnot(length(imc_sig) == 1, imc_sig %in% c(SIG_BASIS, "ns"))
imc_stats_txt <- readLines(IMC_STATS)
imc_result <- imc_stats_txt[grep("^=== RESULT", imc_stats_txt) + 0:2]
imc_n_line <- grep("Cells modelled", imc_stats_txt, value = TRUE)

# -------------------------------------------------------------------
# CosMx: the dream model's cells, GFAP counts, donor-centred
# -------------------------------------------------------------------
deg <- read.delim(COSMX_DEG, check.names = FALSE)
g <- deg[deg$gene == GENE, ]
stopifnot(nrow(g) == 1)

sce <- qs::qread(COSMX_SCE)
# identical cell filter to deg_dream_linear_distance_phf1.r stage 1 (PHF1+ cells carry NA)
sce <- sce[, !is.na(sce$dist_to_phf1_um) & sce$dist_to_phf1_um <= COSMX_CAP]
stopifnot(all(!sce$PHF1), all(sce$dist_to_phf1_um > 0))
if ("n_cells" %in% names(g)) stopifnot(ncol(sce) == g$n_cells)

cx <- tibble(donor = as.character(sce$sample_id),
             dist_um = sce$dist_to_phf1_um,
             count = as.numeric(counts(sce)[GENE, ]),
             lib = as.numeric(colSums(counts(sce))))
grand_mean <- mean(cx$count)
cx <- cx %>% group_by(donor) %>%
  mutate(count_c = count - mean(count) + grand_mean,
         lib_c = lib - mean(lib)) %>% ungroup()
dist_sd <- sd(log(cx$dist_um))
if ("dist_sd" %in% names(g)) stopifnot(abs(dist_sd - g$dist_sd) < 1e-8)
cx$dist_scaled <- log(cx$dist_um) / dist_sd

cx_padj <- g$padj
cx_sig  <- if (isTRUE(cx_padj < SIG_ALPHA)) SIG_BASIS else "ns"

roll_mean <- function(x, y, grid, half_window) {
  out <- lapply(grid, function(gp) {
    idx <- which(x >= gp - half_window & x <= gp + half_window); n <- length(idx)
    if (n < 1) return(c(roll_mean = NA_real_, sem = NA_real_, n_window = 0))
    m <- mean(y[idx]); s <- if (n > 1) stats::sd(y[idx]) / sqrt(n) else NA_real_
    c(roll_mean = m, sem = s, n_window = n)
  })
  r <- as.data.frame(do.call(rbind, out))
  r$roll_mean[r$n_window < MIN_WINDOW_N] <- NA
  r$sem[r$n_window < MIN_WINDOW_N] <- NA
  r
}
GRID <- seq(0, COSMX_CAP, length.out = N_GRID)
cx_roll <- roll_mean(cx$dist_um, cx$count_c, GRID, HALF_WINDOW) %>% mutate(dist_um = GRID)

# donor-level view (descriptive): per-donor OLS slope of donor-centred counts on scaled log
# distance, sign count and dz
per_donor <- cx %>% group_by(donor) %>%
  summarise(n_cells = n(),
            slope_counts_per_sd = unname(coef(lm(count_c ~ dist_scaled))[2]),
            .groups = "drop")
dz <- mean(per_donor$slope_counts_per_sd) / sd(per_donor$slope_counts_per_sd)

# -------------------------------------------------------------------
# Axis map: right (CosMx counts) -> left (IMC z), matching the drawn extents incl. ribbons
# -------------------------------------------------------------------
ext <- function(r) range(c(r$roll_mean - 1.96 * r$sem, r$roll_mean + 1.96 * r$sem), na.rm = TRUE)
imc_ext <- ext(imc); cx_ext <- ext(cx_roll)
map_b <- diff(imc_ext) / diff(cx_ext)
map_a <- imc_ext[1] - map_b * cx_ext[1]
to_left  <- function(v) map_a + map_b * v
to_right <- function(v) (v - map_a) / map_b

plot_df <- bind_rows(
  imc %>% transmute(modality = "IMC", dist_um, value = roll_mean, sem, n_window,
                    sig = imc_sig),
  cx_roll %>% transmute(modality = "CosMx", dist_um, value = roll_mean, sem, n_window,
                        sig = cx_sig)
) %>%
  mutate(y  = ifelse(modality == "IMC", value, to_left(value)),
         lo = ifelse(modality == "IMC", value - 1.96 * sem, to_left(value - 1.96 * sem)),
         hi = ifelse(modality == "IMC", value + 1.96 * sem, to_left(value + 1.96 * sem)),
         modality = factor(modality, levels = c("IMC", "CosMx")),
         sig = factor(sig, levels = c(SIG_BASIS, "ns")))

# Pad any absent linetype level with one undrawn row so the key still shows its glyph
# (see pad_levels() in imc_phf1_glia_distance.R).
miss <- setdiff(c(SIG_BASIS, "ns"), as.character(plot_df$sig))
if (length(miss)) {
  filler <- plot_df[rep(1L, length(miss)), ]
  filler$sig <- factor(miss, levels = c(SIG_BASIS, "ns"))
  filler[, c("y", "lo", "hi")] <- NA_real_
  plot_df <- bind_rows(plot_df, filler)
}

p <- ggplot(plot_df, aes(dist_um, y, colour = modality, fill = modality,
                         group = interaction(modality, sig))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(aes(linetype = sig, linewidth = sig)) +
  scale_colour_manual(values = c(IMC = COL_IMC, CosMx = COL_COSMX), guide = "none") +
  scale_fill_manual(values = c(IMC = COL_IMC, CosMx = COL_COSMX), guide = "none") +
  scale_linetype_manual(values = setNames(c("solid", "dashed"), c(SIG_BASIS, "ns")),
                        name = "Log-model FDR", drop = FALSE, limits = c(SIG_BASIS, "ns")) +
  scale_linewidth_manual(values = setNames(c(0.8, 0.4), c(SIG_BASIS, "ns")),
                         name = "Log-model FDR", drop = FALSE, limits = c(SIG_BASIS, "ns")) +
  guides(linetype = guide_legend(override.aes = list(colour = "black")),
         linewidth = guide_legend()) +
  scale_y_continuous(name = "IMC GFAP (within-donor z)",
                     sec.axis = sec_axis(~ to_right(.), name = "CosMx GFAP (counts/cell)")) +
  scale_x_continuous(breaks = seq(0, COSMX_CAP, 250)) +
  labs(x = expression("Distance to PHF1+ neuron (" * mu * "m)")) +
  coord_cartesian(xlim = c(0, COSMX_CAP)) +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title.y.left  = element_text(colour = COL_IMC),
        axis.title.y.right = element_text(colour = COL_COSMX),
        legend.position = "inside",
        legend.position.inside = c(0.42, 0.03), legend.justification = c(0, 0),
        legend.background = element_blank(),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6),
        legend.key.width = grid::unit(14, "pt"), legend.key.height = grid::unit(7, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(t = 4, r = 4, b = 4, l = 4))

ggsave(file.path(out_dir, sprintf("plot_%s.pdf", slug)), p,
       width = FIG_W_CM, height = FIG_H_CM, units = "cm", device = "pdf")

# -------------------------------------------------------------------
# Source data (long: one row per drawn grid point per modality)
# -------------------------------------------------------------------
write.table(plot_df %>% filter(!is.na(dist_um)) %>%
              transmute(modality, unit = ifelse(modality == "IMC", "within_donor_z",
                                                "GFAP_counts_per_cell_donor_centred"),
                        window_um = WINDOW_UM, significant = sig, dist_um,
                        roll_mean = value, sem, n_window, y_left_axis = y),
            file.path(out_dir, sprintf("source_data_%s.tsv", slug)),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(per_donor, file.path(out_dir, sprintf("source_data_%s_cosmx_perdonor.tsv", slug)),
            sep = "\t", quote = FALSE, row.names = FALSE)

eff <- tibble(
  modality = c("IMC", "CosMx"),
  outcome = c("GFAP protein, within-donor z", "GFAP transcripts (dream logFC, log2)"),
  estimate_per_sd_log_dist = c(as.numeric(sub(".*estimate ([+-][0-9.]+).*", "\\1", imc_result[2])),
                               g$logFC),
  ci_lo = c(as.numeric(sub(".*\\[([+-][0-9.]+),.*", "\\1", imc_result[3])), g$CI.L),
  ci_hi = c(as.numeric(sub(".*, ([+-][0-9.]+)\\].*", "\\1", imc_result[3])), g$CI.R),
  padj = c(as.numeric(sub(".*padj ([0-9.eE+-]+).*", "\\1", imc_result[2])), g$padj),
  cap_um = c(300, COSMX_CAP),
  source = c(IMC_STATS, COSMX_DEG))
write.table(eff, file.path(out_dir, sprintf("stats_%s_effectsize.tsv", slug)),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------------------------------------------------
# Stats log
# -------------------------------------------------------------------
sink(file.path(out_dir, sprintf("stats_%s.txt", slug)))
cat("GFAP in astrocytes vs distance to nearest PHF1+ neuron: IMC protein + CosMx transcripts\n")
cat("=====================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("NO MODEL IS FITTED IN THIS SCRIPT. Each line's linetype comes from its canonical fit.\n")
cat("Direction convention for both: POSITIVE = higher FURTHER from a tangle.\n\n")

cat("--- IMC (left axis, cyan) ---\n")
cat("Source:", IMC_STATS, "\n")
cat("Model (Set 3): GFAP_z ~ dist_scaled + Sex + Age_s + PMI_s + (1|patient_id), 300 um cap\n")
cat(imc_n_line, "\n")
cat(paste(imc_result, collapse = "\n"), "\n")
cat("Curve read verbatim from", IMC_ROLL, "\n")
cat("Donor RE is singular by construction (within-donor z); see the source log's RE block.\n\n")

cat("--- CosMx (right axis, vermillion) ---\n")
cat("Source:", COSMX_DEG, "(GFAP row)\n")
cat("Model (Set 1):", g$model, "\n")
cat(sprintf("  cap %d um | cells %d | donors %d | dist_sd %.6f (matches the DEG table)\n",
            COSMX_CAP, nrow(cx), dplyr::n_distinct(cx$donor), dist_sd))
cat(sprintf("  logFC (log2) per s.d. log distance %+.5f, 95%% CI [%+.5f, %+.5f]\n",
            g$logFC, g$CI.L, g$CI.R))
cat(sprintf("  t %+.3f  p %.3g  padj %.3g  -> %s\n", g$t, g$pval, g$padj,
            if (cx_sig == SIG_BASIS) "solid" else "dashed"))
cat(sprintf("  Fold change per s.d. log distance: %.4f [%.4f, %.4f]\n",
            2^g$logFC, 2^g$CI.L, 2^g$CI.R))
cat("Curve: GFAP counts per cell, centred within donor and re-offset to the cohort mean\n")
cat(sprintf("  (%.3f counts/cell), matching the model's donor random intercept.\n",
            grand_mean))
cat("  Rolling window", WINDOW_UM, "um,",
    N_GRID, "grid points over 0 -", COSMX_CAP, "um,\n")
cat(sprintf("  ribbon = mean +/- 1.96 x cell-level SEM, windows with < %d cells blanked.\n",
            MIN_WINDOW_N))
cat(sprintf("  Donor-centred library size, OLS slope per s.d. log distance: %+.2f counts\n",
            unname(coef(lm(cx$lib_c ~ cx$dist_scaled))[2])))
cat(sprintf("  (mean library %.0f), i.e. the curve is not a library-size artefact.\n\n",
            mean(cx$lib)))

cat("Donor-level view (9 donors):\n")
print(as.data.frame(per_donor), row.names = FALSE, digits = 4)
cat(sprintf("  donors with positive slope: %d / %d; donor dz = %.3f\n",
            sum(per_donor$slope_counts_per_sd > 0), nrow(per_donor), dz))
cat("  These are slopes of the DRAWN quantity (donor-centred raw counts, 0-1000 um), with no\n")
cat("  nUMI_log / percent_neg adjustment, so they describe the curve, not the dream model's\n")
cat("  per-donor analogue. The dream logFC + CI is the effect size.\n\n")

cat("WINDOW: each line is drawn over the window its model was fitted on (IMC 0-300 um,\n")
cat("  CosMx 0-1000 um).\n\n")

cat("Dual-axis map: left = a + b * right, a =", signif(map_a, 6), ", b =", signif(map_b, 6), "\n")
cat("  chosen so both drawn curves (incl. ribbons) span the same vertical extent. The relative\n")
cat("  vertical scaling is ARBITRARY: compare direction and shape, not height.\n")
cat(sprintf("\nFigure: %.1f x %.1f cm, default pdf device.\n", FIG_W_CM, FIG_H_CM))
cat("\n=== sessionInfo() ===\n"); print(sessionInfo())
sink()

message("Done: ", file.path(out_dir, sprintf("plot_%s.pdf", slug)))
