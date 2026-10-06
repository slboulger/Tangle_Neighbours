#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# phf1_morphology_dist_phf1_panel.R
#
# Figure panels: S7B
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# phf1_morphology_dist_phf1_panel.R
#
# CosMx replicate of R/imc_phf1_morphology_dist_phf1_panel.R
# (plots/imc_phf1_morphology_dist_phf1_panel/plot_nucleus_area_dist_phf1.pdf), measure for
# measure and element for element, on Exc-IT-L2-3-CBLN2-HOPX instead of IMC excitatory
# neurons: a rolling mean over distance to the nearest PHF1+ neuron, with the PHF1+ cells
# themselves as a MARGINAL STRIP at the left edge of the same panel, so the distance gradient
# and the cell-autonomous contrast are read off one axis. One figure per measure, plus a
# _nostrip variant.
#
# ENCODINGS, copied from the reference:
#   line style   log-distance model verdict     solid 0.8 = padj<0.05, dashed 0.4 = ns
#   point shape  PHF1+ vs PHF1- contrast        filled 16 = padj<0.05, hollow 1 = ns
#   The two legend NAMES must differ or ggplot merges the guides and the encodings collide.
#
# COMPLETE LEGEND EVEN WHEN NOTHING IS SIGNIFICANT. A guide key is built from LAYER DATA, so
# a factor level that never occurs renders its label with NO glyph -- the reference figure has
# a blank "ns" key for exactly this reason, and a panel where every measure is n.s. would
# instead lose its "padj<0.05" key. drop = FALSE and override.aes do not fix it. This script
# carries one UNDRAWN row per missing level (value NA), so both keys always draw. The line
# `group` includes the sig factor, or the filler joins a real line and ggplot refuses to draw
# one line with two linetypes.
#
# SCALE: RAW native units (um^2, um, unitless), as the reference plots raw values
# (CENTRE_ON_DONOR = FALSE). Both curve and strip are absolute, so the PHF1+ level and the
# PHF1- field read off one axis. Unlike the channel intensities, morphology is a physical
# measurement and IS comparable across donors, so no within-donor normalisation is needed.
#
# READ THE CURVE AND THE LINE STYLE AS TWO SEPARATE THINGS. The drawn curve POOLS donors; the
# model behind the line style is WITHIN donor. Both slopes are computed and printed, and the
# script WARNS on a sign clash, so a figure whose curve slopes against its own significance
# encoding cannot be produced silently. Set CENTRE_ON_DONOR <- TRUE for the donor-centred
# version of exactly these panels.
#
# The strip is a REFERENCE, not a distance-zero data point. PHF1+ neurons ARE the anchors, so
# they have no "distance to the nearest OTHER PHF1+ neuron" (dist_to_phf1_um is NA for them by
# design); they are drawn off-axis to the left of a separator rule, as in the reference.
#
# STATISTICS -- both refitted here so the figure cannot drift from its own annotations:
#   distance   value ~ dist_scaled + Sex + Age_s + PMI_s + (1 | sample_id), PHF1- cells,
#              matching R/phf1_distance_morphology_exc.R's headline exactly.
#   contrast   value ~ grp + Sex + Age_s + PMI_s + (1 | sample_id), all cells, matching
#              R/phf1_morphology_exc.R's headline. NOTE this DIFFERS from the IMC panel,
#              whose contrast headline is a donor-mean paired Wilcoxon -- that script has 44
#              donors, this cohort has 9 with 1-60 PHF1+ cells each, so the cell-level LMM is
#              the CosMx headline. The donor-paired Wilcoxon and dz are
#              computed and reported alongside.
#   BH across the 5 measures within each family.
#
# Nuclear measures require a plausible nucleus -- see R/phf1_morphology_filters.R.
#
# Distance cap is the single positional argument:
#   Rscript R/phf1_morphology_dist_phf1_panel.R 1000
#   Rscript R/phf1_morphology_dist_phf1_panel.R 300     # matches the IMC family
#
# Outputs the standard TRIPLE under plots/phf1_morphology_dist_phf1_panel_<cap>um/.

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(qs)
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
  library(lme4); library(lmerTest); library(emmeans)
})

hpc <- "<PROJECT_ROOT>/phf1_v2"
loc <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(hpc)) hpc else loc)
source("R/palettes.R")
source("R/phf1_morphology_filters.R")

CELLTYPE <- "Exc-IT-L2-3-CBLN2-HOPX"
.cli <- commandArgs(trailingOnly = TRUE)
MAX_DIST <- if (length(.cli) >= 1) as.numeric(.cli[1]) else 1000
if (!is.finite(MAX_DIST) || MAX_DIST <= 0) stop("max_dist must be a positive number")

FDR          <- 0.05
WINDOW_UM    <- 75
HALF_WINDOW  <- WINDOW_UM / 2
GRID_N       <- 200
MIN_WINDOW_N <- 10

SIG_BASIS  <- "padj<0.05"
SIG_LEGEND <- "Log-distance"        # linetype/linewidth: the distance-model verdict
CON_LEGEND <- "PHF1+ vs PHF1-"      # point shape: the contrast verdict. MUST differ.
CON_SHAPES <- c(16, 1)              # filled = significant, hollow = n.s.
SERIES_COL <- "#BD0026"             # the PHF1 red used across the morphology family

XLAB <- expression("Distance to PHF1+ neuron (" * mu * "m)")

# TOTAL figure width fixed at 6.14 cm, the family standard; the panel takes whatever the axes
# and the two guide boxes leave. Same caveat as the reference: the strip variant carries ~14%
# more x data range in the same physical width, so a micron is ~12% shorter there. Compare
# curves within a variant, not across the two.
TOTAL_W_IN  <- 6.14 / 2.54
PANEL_H_MIN <- 1.66

# TRUE plots each cell relative to its own donor's PHF1-negative mean, which makes the curve
# depict the same within-donor quantity the LMM tests and the strip height exactly the
# donor-paired difference. FALSE (default) plots raw values, as the reference figure does.
CENTRE_ON_DONOR <- FALSE

out_dir <- sprintf("plots/phf1_morphology_dist_phf1_panel_%dum", MAX_DIST)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE)

MEASURES <- tibble::tribble(
  ~col,           ~label,                 ~unit,      ~nuclear,
  "nucarea",      "Nucleus area",         "um^2",     TRUE,
  "nucaspect",    "Nucleus aspect ratio", "unitless", TRUE,
  "circularity",  "Cell circularity",     "unitless", FALSE,
  "eccentricity", "Cell eccentricity",    "unitless", FALSE,
  "perimeter",    "Cell perimeter",       "um",       FALSE
)
YLAB_EXPR <- if (CENTRE_ON_DONOR) list(
  nucarea      = expression(Delta * " nucleus area (" * mu * "m"^2 * ")"),
  nucaspect    = expression(Delta * " nucleus aspect ratio"),
  circularity  = expression(Delta * " cell circularity"),
  eccentricity = expression(Delta * " cell eccentricity"),
  perimeter    = expression(Delta * " cell perimeter (" * mu * "m)")
) else list(
  nucarea      = expression("Nucleus area (" * mu * "m"^2 * ")"),
  nucaspect    = "Nucleus aspect ratio",
  circularity  = "Cell circularity",
  eccentricity = "Cell eccentricity",
  perimeter    = expression("Cell perimeter (" * mu * "m)")
)

emmeans::emm_options(lmer.df = "satterthwaite", lmerTest.limit = 1e5)

## ---------------------------------------------------------------------------
## Cells -- identical construction to the two source scripts
## ---------------------------------------------------------------------------
sce_path <- file.path("celltype_sce_neighbours", paste0(CELLTYPE, "_sce_neighbours.qs"))
if (!file.exists(sce_path)) stop("Not found: ", sce_path)
sce <- qread(sce_path)
cd  <- as.data.frame(colData(sce)); cd$cell_id <- colnames(sce)
scale_err <- assert_px_scale(cd)
nuc_audit <- nucleus_filter_audit(cd$NucArea, group = "all CBLN2")

exc <- cd %>%
  transmute(
    cell_id, sample_id = as.character(sample_id),
    Braak = factor(as.character(Braak), levels = braak_levels),
    Sex = factor(as.character(Sex)), Age = as.numeric(Age), PMI = as.numeric(PMI),
    phf1_pos = as.logical(PHF1 %in% c(TRUE, "TRUE", "True", 1, "1")),
    dist_um  = as.numeric(dist_to_phf1_um),
    nucarea      = ifelse(nucleus_ok(NucArea), NucArea * PX2_UM2, NA_real_),
    nucaspect    = ifelse(nucleus_ok(NucArea), NucAspectRatio, NA_real_),
    circularity  = Circularity,
    eccentricity = Eccentricity,
    perimeter    = Perimeter * PX_UM
  ) %>%
  mutate(grp   = factor(ifelse(phf1_pos, "PHF1+", "PHF1-"), levels = c("PHF1-", "PHF1+")),
         Age_s = as.numeric(scale(Age)), PMI_s = as.numeric(scale(PMI)))
stopifnot(!anyNA(exc$Sex), !anyNA(exc$Age_s), !anyNA(exc$PMI_s))

n_cells <- nrow(exc); n_pos <- sum(exc$phf1_pos); n_don <- n_distinct(exc$sample_id)
if (n_don != 9) stop("Expected 9 donors, got ", n_don)

# Donor reference means, over each donor's PHF1-NEGATIVE cells, for the centred variant.
donor_ref <- exc %>% filter(grp == "PHF1-") %>% group_by(sample_id) %>%
  summarise(across(all_of(MEASURES$col), ~ mean(.x, na.rm = TRUE), .names = "ref_{.col}"),
            .groups = "drop")
exc <- exc %>% inner_join(donor_ref, by = "sample_id")
for (mc in MEASURES$col) {
  exc[[paste0("c_", mc)]] <-
    if (CENTRE_ON_DONOR) exc[[mc]] - exc[[paste0("ref_", mc)]] else exc[[mc]]
}

# Distance frame: PHF1-negative cells only (PHF1+ cells ARE the anchors).
md <- exc %>% filter(!phf1_pos, !is.na(dist_um), dist_um <= MAX_DIST)
if (any(md$dist_um <= 0)) stop("non-positive distance; log() would fail")
DIST_SD <- stats::sd(log(md$dist_um))
md$dist_scaled <- log(md$dist_um) / DIST_SD
md$Sex <- droplevels(md$Sex)

cat(sprintf("%s | %d cells (%d PHF1+) | %d donors | distance frame %d cells | dist_sd %.5f | cap %d um\n",
            CELLTYPE, n_cells, n_pos, n_don, nrow(md), DIST_SD, MAX_DIST))

## ---------------------------------------------------------------------------
## Statistics -- both refitted here, matching the two source scripts' headlines
## ---------------------------------------------------------------------------
dist_one <- function(mc) {
  d <- md %>% transmute(value = .data[[mc]], dist_scaled, sample_id, Sex, Age_s, PMI_s)
  d <- d[is.finite(d$value), ]
  m <- lmerTest::lmer(value ~ dist_scaled + Sex + Age_s + PMI_s + (1 | sample_id),
                      data = d, REML = TRUE)
  co <- summary(m)$coefficients["dist_scaled", ]
  crit <- stats::qt(0.975, co["df"])
  vc <- as.data.frame(lme4::VarCorr(m)); sd_tot <- sqrt(sum(vc$vcov[is.na(vc$var2)]))
  re <- vc$vcov[vc$grp == "sample_id" & is.na(vc$var2)]
  pl <- summary(stats::lm(value ~ dist_scaled, data = d))$coefficients["dist_scaled", ]
  tibble(measure = mc,
         beta_closer = -co["Estimate"], SE = co["Std. Error"], df = co["df"],
         CI.L = -co["Estimate"] - crit * co["Std. Error"],
         CI.R = -co["Estimate"] + crit * co["Std. Error"],
         cohens_d = -co["Estimate"] / sd_tot, pval = co["Pr(>|t|)"],
         n_cells = nrow(d), n_donors = n_distinct(d$sample_id),
         re_var_donor = if (length(re)) re[1] else NA_real_,
         re_singular = lme4::isSingular(m),
         pooled_beta_closer = -pl[1], pooled_p = pl[4])
}
res_dist <- bind_rows(lapply(MEASURES$col, dist_one)) %>%
  mutate(padj = p.adjust(pval, method = "BH")) %>%
  left_join(MEASURES, by = c("measure" = "col"))

contrast_one <- function(mc) {
  d <- exc %>% transmute(value = .data[[mc]], grp, sample_id, Sex, Age_s, PMI_s)
  d <- d[is.finite(d$value), ]
  m <- lmerTest::lmer(value ~ grp + Sex + Age_s + PMI_s + (1 | sample_id),
                      data = d, REML = TRUE)
  pr <- as.data.frame(summary(pairs(emmeans(m, ~ grp)), infer = TRUE))
  vc <- as.data.frame(lme4::VarCorr(m)); sd_tot <- sqrt(sum(vc$vcov[is.na(vc$var2)]))
  # pairs() returns (PHF1-) - (PHF1+); re-orient as PHF1+ minus PHF1-
  est <- -pr$estimate[1]; lo <- -pr$upper.CL[1]; hi <- -pr$lower.CL[1]
  # donor-paired supporting statistics (dz, paired Wilcoxon, direction count)
  dm <- d %>% group_by(sample_id, grp) %>% summarise(m = mean(value), .groups = "drop") %>%
    pivot_wider(names_from = grp, values_from = m)
  dm <- dm[stats::complete.cases(dm), ]
  diff <- dm[["PHF1+"]] - dm[["PHF1-"]]
  pw <- suppressWarnings(wilcox.test(dm[["PHF1+"]], dm[["PHF1-"]], paired = TRUE))
  tibble(measure = mc, lmm_diff = est, CI.L = lo, CI.R = hi,
         cohens_d = est / sd_tot, pval = pr$p.value[1],
         mean_PHF1_neg = mean(d$value[d$grp == "PHF1-"]),
         pct_of_PHF1_neg = 100 * est / mean(d$value[d$grp == "PHF1-"]),
         dz = mean(diff) / stats::sd(diff), paired_diff = mean(diff),
         n_donors = nrow(dm), n_same_direction = sum(sign(diff) == sign(mean(diff))),
         wilcox_p = pw$p.value,
         n_cells = nrow(d), n_PHF1_pos = sum(d$grp == "PHF1+"))
}
res_con <- bind_rows(lapply(MEASURES$col, contrast_one)) %>%
  mutate(padj = p.adjust(pval, method = "BH")) %>%
  left_join(MEASURES, by = c("measure" = "col"))

sign_clash <- with(res_dist, sign(beta_closer) != sign(pooled_beta_closer) & padj < FDR)
if (any(sign_clash))
  warning(sprintf(paste("SIGN CLASH: pooled curve and within-donor model disagree for %s.",
                        "The drawn line slopes against its own significance encoding.",
                        "Set CENTRE_ON_DONOR <- TRUE, or read the curve as descriptive only."),
                  paste(res_dist$label[sign_clash], collapse = ", ")))

## ---------------------------------------------------------------------------
## Rolling mean + PHF1+ strip
## ---------------------------------------------------------------------------
GRID <- seq(0, MAX_DIST, length.out = GRID_N)
roll_mean <- function(x, y, grid, hw) {
  as.data.frame(do.call(rbind, lapply(grid, function(g) {
    idx <- which(x >= g - hw & x <= g + hw); n <- length(idx)
    if (n < 1) return(c(roll_mean = NA_real_, sem = NA_real_, n_window = 0))
    c(roll_mean = mean(y[idx]),
      sem = if (n > 1) stats::sd(y[idx]) / sqrt(n) else NA_real_, n_window = n)
  })))
}
sig_factor <- function(x) factor(ifelse(x, SIG_BASIS, "ns"), levels = c(SIG_BASIS, "ns"))

roll_df <- bind_rows(lapply(MEASURES$col, function(mc) {
  cc <- paste0("c_", mc)
  d <- md %>% filter(is.finite(.data[[cc]]))
  r <- roll_mean(d$dist_um, d[[cc]], GRID, HALF_WINDOW)
  r$roll_mean[r$n_window < MIN_WINDOW_N] <- NA
  r$sem[r$n_window < MIN_WINDOW_N] <- NA
  sg <- isTRUE(res_dist$padj[res_dist$measure == mc] < FDR)
  r %>% mutate(measure = mc, dist_um = GRID, sig = sig_factor(sg), significant = sg)
})) %>% left_join(MEASURES, by = c("measure" = "col"))

ref_df <- bind_rows(lapply(MEASURES$col, function(mc) {
  cc <- paste0("c_", mc)
  s <- exc[[cc]][exc$grp == "PHF1+" & is.finite(exc[[cc]])]
  n <- length(s); sem <- stats::sd(s) / sqrt(n)
  tibble(measure = mc, ref_mean = mean(s), ci_lo = mean(s) - 1.96 * sem,
         ci_hi = mean(s) + 1.96 * sem, n = n,
         con_sig = sig_factor(isTRUE(res_con$padj[res_con$measure == mc] < FDR)))
})) %>% left_join(MEASURES, by = c("measure" = "col"))

## ---------------------------------------------------------------------------
## Panels
## ---------------------------------------------------------------------------
# One undrawn row per absent level, so BOTH keys of BOTH guides always render -- including
# on a panel where nothing is significant, or where everything is.
pad_levels <- function(d, col) {
  miss <- setdiff(c(SIG_BASIS, "ns"), as.character(d[[col]]))
  if (!length(miss)) return(d)
  filler <- d[rep(1L, length(miss)), , drop = FALSE]
  filler[[col]] <- factor(miss, levels = c(SIG_BASIS, "ns"))
  for (v in intersect(c("roll_mean", "sem", "ref_mean", "ci_lo", "ci_hi"), names(filler)))
    filler[[v]] <- NA_real_
  dplyr::bind_rows(d, filler)
}

save_fixed_panel <- function(p, path, panel_h = PANEL_H_MIN) {
  g <- ggplot2::ggplotGrob(p)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  prow <- unique(g$layout$t[grepl("^panel", g$layout$name)])
  g$heights[prow] <- grid::unit(panel_h, "in")
  # Measure everything that is NOT the panel column (axis, labels, both guide boxes,
  # margins) by zeroing the panel, then give the panel exactly the remainder of the total.
  g$widths[pcol] <- grid::unit(0, "in")
  nonpanel_w <- grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE)
  panel_w <- TOTAL_W_IN - nonpanel_w
  if (!is.finite(panel_w) || panel_w < 0.35) {
    warning(sprintf("%s: only %.2f in left for the panel at a %.2f in total width",
                    basename(path), panel_w, TOTAL_W_IN))
    panel_w <- max(panel_w, 0.35)
  }
  g$widths[pcol] <- grid::unit(panel_w, "in")
  total_h <- tryCatch(grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE),
                      error = function(e) NA_real_)
  if (!is.finite(total_h)) total_h <- panel_h + 0.44
  ggsave(path, g, width = TOTAL_W_IN, height = total_h, units = "in", device = "pdf")
}

build_panel <- function(mc, with_strip = TRUE) {
  rd <- pad_levels(roll_df %>% filter(measure == mc, !is.na(roll_mean)), "sig")
  rf <- pad_levels(ref_df  %>% filter(measure == mc), "con_sig")
  ylab <- YLAB_EXPR[[mc]]

  p <- ggplot(rd, aes(dist_um, roll_mean, group = sig))
  if (CENTRE_ON_DONOR)
    p <- p + geom_hline(yintercept = 0, linetype = 3, linewidth = 0.25, colour = "grey60")
  p <- p +
    geom_ribbon(aes(ymin = roll_mean - 1.96 * sem, ymax = roll_mean + 1.96 * sem),
                fill = SERIES_COL, alpha = 0.15, colour = NA) +
    geom_line(aes(linetype = sig, linewidth = sig), colour = SERIES_COL)

  if (with_strip) {
    x_pos <- -0.07 * MAX_DIST
    p <- p +
      annotate("segment", x = -0.025 * MAX_DIST, xend = -0.025 * MAX_DIST,
               y = -Inf, yend = Inf, linewidth = 0.25, colour = "grey70") +
      geom_errorbar(data = rf, aes(x = x_pos, y = ref_mean, ymin = ci_lo, ymax = ci_hi),
                    inherit.aes = FALSE, width = 0, linewidth = 0.4, colour = SERIES_COL) +
      geom_point(data = rf, aes(x = x_pos, y = ref_mean, shape = con_sig),
                 inherit.aes = FALSE, size = 1.3, stroke = 0.5, colour = SERIES_COL) +
      scale_shape_manual(values = setNames(CON_SHAPES, c(SIG_BASIS, "ns")),
                         name = CON_LEGEND, drop = FALSE, limits = c(SIG_BASIS, "ns")) +
      annotate("text", x = x_pos, y = Inf, vjust = -0.45,
               size = 2.6, fontface = "bold", colour = "grey15", label = "PHF1+")
  }

  # n = 3, not pretty()'s default: the right-hand guide box leaves only ~3 cm of panel, so
  # the default collides (7 ticks at a 300 um cap, and "800"/"1000" overlap at 1000).
  brk <- pretty(c(0, MAX_DIST), n = 3); brk <- brk[brk >= 0 & brk <= MAX_DIST]
  xlo <- if (with_strip) -0.14 * MAX_DIST else 0

  p +
    scale_linetype_manual(values = setNames(c("solid", "dashed"), c(SIG_BASIS, "ns")),
                          name = SIG_LEGEND, drop = FALSE, limits = c(SIG_BASIS, "ns")) +
    scale_linewidth_manual(values = setNames(c(0.8, 0.4), c(SIG_BASIS, "ns")),
                           name = SIG_LEGEND, drop = FALSE, limits = c(SIG_BASIS, "ns")) +
    scale_x_continuous(breaks = brk) +
    guides(linetype = guide_legend(order = 1), linewidth = guide_legend(order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(colour = "grey25", size = 1.3))) +
    labs(x = XLAB, y = ylab) +
    coord_cartesian(xlim = c(xlo, MAX_DIST), clip = "off") +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "right", legend.direction = "vertical", legend.box = "vertical",
          legend.text = element_text(size = 6), legend.title = element_text(size = 6),
          legend.key.width = grid::unit(16, "pt"), legend.key.height = grid::unit(7, "pt"),
          legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
          legend.spacing.y = grid::unit(1, "pt"),
          # headroom for the "PHF1+" label; applied to BOTH variants so the strip and
          # _nostrip figures are not vertically offset from each other on the page
          plot.margin = margin(t = 14, r = 6, b = 4, l = 5))
}

for (mc in MEASURES$col) {
  save_fixed_panel(build_panel(mc, TRUE),
                   file.path(out_dir, sprintf("plot_%s_dist_phf1.pdf", mc)))
  save_fixed_panel(build_panel(mc, FALSE),
                   file.path(out_dir, sprintf("plot_%s_dist_phf1_nostrip.pdf", mc)))
}

## ---------------------------------------------------------------------------
## Source data + stats log, mirroring the IMC panel directory
## ---------------------------------------------------------------------------
wr(roll_df %>% select(measure, label, unit, dist_um, roll_mean, sem, n_window,
                      sig, significant),
   "source_data_morphology_dist_phf1_rollmean.tsv")
wr(ref_df %>% select(measure, label, ref_mean, ci_lo, ci_hi, n, con_sig),
   "source_data_morphology_dist_phf1_ref.tsv")
wr(res_con %>% select(measure, label, unit, lmm_diff, CI.L, CI.R, cohens_d,
                      mean_PHF1_neg, pct_of_PHF1_neg, dz, paired_diff, n_donors,
                      n_same_direction, wilcox_p, pval, padj, n_cells, n_PHF1_pos),
   "source_data_morphology_dist_phf1_contrast_stats.tsv")
wr(res_dist %>% select(measure, label, unit, beta_closer, SE, df, CI.L, CI.R, cohens_d,
                       pval, padj, n_cells, n_donors, pooled_beta_closer, pooled_p,
                       re_var_donor, re_singular),
   "source_data_morphology_dist_phf1_distance_stats.tsv")

sink(file.path(out_dir, "stats_morphology_dist_phf1.txt"))
ttl <- sprintf("CosMx morphology over PHF1 distance with the PHF1+ strip -- %s (%d um)",
               CELLTYPE, MAX_DIST)
cat(ttl, "\n"); cat(strrep("=", nchar(ttl)), "\n"); cat("Date:", format(Sys.time()), "\n\n")
cat("Source: ", sce_path, "\n", sep = "")
cat("CosMx replicate of R/imc_phf1_morphology_dist_phf1_panel.R, element for element.\n")
cat("Companion analyses: plots/phf1_morphology/ (contrast) and\n")
cat(sprintf("plots/phf1_distance_morphology_%dum/ (distance).\n", MAX_DIST))
cat("\n")

cat(sprintf("SCALE: %s.\n", if (CENTRE_ON_DONOR)
  "each cell relative to its own donor's PHF1-negative mean (CENTRE_ON_DONOR = TRUE)" else
  "RAW native units (um^2 / um / unitless), as the reference figure (CENTRE_ON_DONOR = FALSE)"))
cat("  Morphology is a physical measurement and IS comparable across donors, so unlike the\n")
cat("  channel intensities no within-donor normalisation is applied.\n")
cat("  The drawn curve POOLS donors; the model behind the line style is WITHIN donor. Those\n")
cat("  are different quantities. Both slopes are printed below and the script warns on a\n")
cat("  sign clash.\n\n")

cat("POOLED (donor-ignoring) slope of the drawn curve, per s.d. CLOSER:\n")
print(as.data.frame(res_dist %>% select(label, beta_closer, padj, pooled_beta_closer, pooled_p)),
      row.names = FALSE, digits = 4)
if (any(sign_clash)) cat("  *** SIGN CLASH -- see warning; read the curve as descriptive ***\n")

cat("\nSTRIP: PHF1+ neurons are the ANCHORS and have no distance to a nearest OTHER PHF1+\n")
cat("  neuron (dist_to_phf1_um is NA for them by design), so they are drawn off-axis left of\n")
cat("  a separator rule. It is a REFERENCE LEVEL, not a distance-zero data point.\n")
cat("ENCODINGS: line style = log-distance verdict (solid padj<0.05, dashed ns);\n")
cat("  point shape = PHF1+ vs PHF1- verdict (filled padj<0.05, hollow ns).\n")
cat("  Both keys of both guides always render, including when a level does not occur.\n")

cat_nucleus_filter_note(nuc_audit)

cat("\nCELLS\n")
cat(sprintf("  contrast : %d cells (%d PHF1+, %d PHF1-) from %d donors\n",
            n_cells, n_pos, n_cells - n_pos, n_don))
cat(sprintf("  distance : %d PHF1-negative cells from %d donors, window 0-%d um\n",
            nrow(md), n_distinct(md$sample_id), MAX_DIST))
cat(sprintf("  dist_sd  : %.5f (log distance, within this celltype, after the cap)\n", DIST_SD))
cat("  n per measure differs where the nucleus filter applies -- see n_cells in the tables.\n")
cat("  PHF1+ cells per donor:\n")
print(as.data.frame(exc %>% filter(phf1_pos) %>% count(sample_id, Braak, name = "n_PHF1_pos")),
      row.names = FALSE)

cat("\nDISTANCE MODEL (PHF1-negative cells; matches R/phf1_distance_morphology_exc.R)\n")
cat("  value ~ dist_scaled + Sex + Age_s + PMI_s + (1 | sample_id), REML.\n")
cat("  No FOV-edge term, matching the canonical CosMx distance models\n")
cat("  and the IMC panel. dist_scaled = log(dist)/sd(log(dist)); per s.d. CLOSER.\n\n")
print(as.data.frame(res_dist %>%
        select(label, beta_closer, CI.L, CI.R, cohens_d, pval, padj, n_cells, n_donors,
               re_var_donor, re_singular)),
      row.names = FALSE, digits = 4)

cat("\nCONTRAST (PHF1+ vs PHF1-; matches R/phf1_morphology_exc.R headline)\n")
cat("  value ~ grp + Sex + Age_s + PMI_s + (1 | sample_id), cell-level LMM.\n")
cat("  This DIFFERS from the IMC panel, whose headline is a donor-mean paired Wilcoxon: that\n")
cat("  cohort has 44 donors, this one has 9 with 1-60 PHF1+ cells each, so the cell-level\n")
cat("  LMM is the CosMx headline. The donor-paired statistics are reported alongside.\n\n")
print(as.data.frame(res_con %>%
        select(label, lmm_diff, CI.L, CI.R, pct_of_PHF1_neg, cohens_d, dz,
               n_same_direction, n_donors, wilcox_p, pval, padj, n_cells, n_PHF1_pos)),
      row.names = FALSE, digits = 4)

cat("\n=== PHF1+ reference strip (mean +/- 95% CI, drawn scale) ===\n")
print(as.data.frame(ref_df %>% select(label, ref_mean, ci_lo, ci_hi, n, con_sig)),
      row.names = FALSE, digits = 4)

cat("\nNOTES:\n")
cat(" - Effect sizes: Cohen's d and the donor-paired dz are reported alongside p.\n")
cat(" - Library-size-adjusted fits are in plots/phf1_morphology/ and\n")
cat(sprintf("   plots/phf1_distance_morphology_%dum/.\n", MAX_DIST))
cat("\nExact p-values and n are in source_data_morphology_dist_phf1_{contrast,distance}_stats.tsv\n")
cat("at full precision; the tables above are rounded for reading only.\n")
cat("\n=== sessionInfo() ===\n"); print(sessionInfo())
sink()

cat("\nWrote", length(list.files(out_dir)), "files to", out_dir, "\n")
print(as.data.frame(res_dist %>% select(measure, beta_closer, padj, n_cells)),
      row.names = FALSE, digits = 3)
print(as.data.frame(res_con %>% select(measure, lmm_diff, padj, n_cells, n_PHF1_pos)),
      row.names = FALSE, digits = 3)
