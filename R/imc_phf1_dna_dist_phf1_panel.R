#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# imc_phf1_dna_dist_phf1_panel.R
#
# Figure panels: 5B
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# imc_phf1_dna_dist_phf1_panel.R
#
# DNA density and total DNA (191Ir; IMC, excitatory neurons) drawn in the style of
# R/plot_reactome_stress_death_vs_phf1_distance.R: a rolling mean over distance to the nearest
# PHF1+ neuron, with the PHF1+ cells themselves shown as a MARGINAL STRIP at the left edge of the
# same panel, so the distance gradient and the cell-autonomous contrast are read off one axis.
# One figure per measure (plus a _nostrip variant). The panel design, encodings, 6.14 cm total
# width and vertical legend are IDENTICAL to R/imc_phf1_morphology_dist_phf1_panel.R and must not
# drift from it -- the two are meant to sit together on a page.
#
# OUTCOMES:
#   dna_density = asinh(mean 191Ir counts in the mask)  -- DNA per unit area (compaction).
#                 Identical to assay(spe,"exprs")["191Ir", ].
#   dna_total   = mean 191Ir counts * area              -- integrated DNA signal per cell, i.e.
#                 the sum of pixel intensities over the mask (area IS the pixel count).
# Both are needed because they can move in opposite directions: PHF1+ nuclei are +26.7% larger,
# so the same DNA spread over more area lowers density while leaving total DNA unchanged.
#
# ENCODINGS, copied from the reference figure:
#   line style   log-distance model verdict     solid 0.8 = padj<0.05, dashed 0.4 = ns
#   point shape  PHF1+ vs PHF1- contrast        filled 16 = padj<0.05, hollow 1 = ns
#   The two legend NAMES must differ or ggplot merges the guides and the encodings collide.
#
# SCALE: RAW arcsinh 191Ir intensity / raw 191Ir x nucleus area, matching the reference figure,
#        which plots raw module scores.
# Both curve and strip are absolute, so the PHF1+ level and the PHF1- field read off one axis.
#
# READ THE CURVE AND THE LINE STYLE AS TWO SEPARATE THINGS. The drawn curve POOLS donors; the
# model behind the line style is WITHIN donor. Those are different quantities, because donors
# differ in mean nucleus size AND in tangle burden (hence in typical distance). Both slopes are
# therefore computed and printed, and the script WARNS on a sign clash, so a figure whose curve
# slopes against its own significance encoding cannot be produced silently.
# Set CENTRE_ON_DONOR <- TRUE for the donor-centred version of exactly these panels.
#
# The strip is a REFERENCE, not a distance-zero data point. PHF1+ neurons ARE the anchors, so they
# have no "distance to the nearest OTHER PHF1+ neuron"; they are drawn off-axis to the left of a
# separator rule, exactly as in the reference figure.
#
# STATISTICS -- both refitted here so the figure cannot drift from its own annotations:
#   distance   value ~ dist_z + Sex + Age_s + PMI_s + (1 | patient_id), PHF1- cells
#              (Set 3, docs/MODELS.md).
#   contrast   paired Wilcoxon on donor means, PHF1+ vs PHF1-.
#   BH across the 2 measures within each family.
#
# Outputs the standard TRIPLE under plots/imc_phf1_dna_dist_phf1_panel/.

suppressPackageStartupMessages({
  library(SpatialExperiment); library(dplyr); library(tibble); library(tidyr)
  library(ggplot2); library(lmerTest); library(RANN)
})

hpc <- "<PROJECT_ROOT>/phf1_v2"
loc <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(hpc)) hpc else loc)
source("R/palettes.R")
source("R/imc_utils.R")   # re_variance_summary(), re_total_sd(): donor RE variance in every log

out_dir <- "plots/imc_phf1_dna_dist_phf1_panel"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

FDR          <- 0.05
MAX_DIST     <- 300
EDGE_BUFFER  <- 0
MIN_ANCHORS  <- 1
WINDOW_UM    <- 75
HALF_WINDOW  <- WINDOW_UM / 2
GRID_N       <- 200
MIN_WINDOW_N <- 10

SIG_BASIS  <- "padj<0.05"
SIG_LEGEND <- "Log-distance"        # linetype/linewidth: the distance-model verdict
CON_LEGEND <- "PHF1+ vs PHF1-"      # point shape: the contrast verdict. MUST differ from above.
CON_SHAPES <- c(16, 1)              # filled = significant, hollow = n.s.
SERIES_COL <- "#BD0026"             # the PHF1 red used across the IMC morphology family

XLAB <- expression("Distance to PHF1+ neuron (" * mu * "m)")

# TOTAL figure width is fixed at 6.14 cm (the width used across the IMC morphology family), so
# the panel takes whatever the axes and the two guide boxes leave. This DIFFERS from the reference
# script, which instead pins panel width and lets the total float in order to keep microns-per-inch
# identical between the strip and _nostrip variants. With the total fixed, that guarantee is gone:
# the strip variant carries ~14% more x data range in the same physical width, so a micron is ~12%
# shorter there. Compare curves within a variant, not across the two.
TOTAL_W_IN  <- 6.14 / 2.54
PANEL_H_MIN <- 1.66

# TRUE plots each cell relative to its own donor's PHF1-negative mean, which makes the curve
# depict the same within-donor quantity the LMM tests and the strip height exactly the
# donor-paired difference. FALSE (default) plots raw values, as the reference figure does.
CENTRE_ON_DONOR <- FALSE

MEASURES <- tibble::tribble(
  ~col,          ~label,                 ~ylab,
  "dna_density", "DNA density (191Ir)",  "delta_density",
  "dna_total",   "Total DNA (191Ir x Area)", "delta_total"
)
# Axis titles carry the delta only. "vs the donor\'s PHF1- mean" is what the delta IS, and it is
# stated in the stats log -- spelled out on the axis it overruns the 1.66 in panel height.
YLAB_EXPR <- if (CENTRE_ON_DONOR) list(
  delta_density = expression(Delta * " DNA intensity (arcsinh)"),
  delta_total   = expression(Delta * " Total DNA (intensity x area)")
) else list(
  delta_density = expression("DNA intensity (arcsinh)"),
  delta_total   = expression("Total DNA (intensity x area)")
)

localwd <- "<IMC_ROOT>/"
spe <- readRDS(paste0(localwd, "spe.rds"))
stopifnot("PHF1_Otsu" %in% colnames(colData(spe)),
          "Matched_4G8_40" %in% colnames(colData(spe)),
          "191Ir" %in% rownames(spe), "area" %in% colnames(colData(spe)))
# 191Ir is NOT a clustering input, unlike CD68/GFAP -- asserted, not assumed
stopifnot(!isTRUE(rowData(spe)["191Ir", "use_channel"]))

# -------------------------------------------------------------------
# Cells -- identical backbone to the other IMC panels
# -------------------------------------------------------------------
cd <- as.data.frame(colData(spe))
cd$cell_id <- colnames(spe)
sc <- spatialCoords(spe)
cd$x <- sc[, "Pos_X"]; cd$y <- sc[, "Pos_Y"]

min_rois_per_patient <- 3
cv_removal <- cd %>% distinct(patient_id, sample_id) %>% count(patient_id) %>%
  filter(n < min_rois_per_patient) %>% pull(patient_id)

base <- cd %>%
  filter(!celltype_clusters %in% c("Artefact cluster", "Unassigned cluster"),
         !patient_id %in% cv_removal,
         BraakGroup != "Braak_0_1",
         !is.na(PHF1_Otsu)) %>%
  mutate(patient_id = as.character(patient_id), sample_id = as.character(sample_id),
         celltype_clusters = as.character(celltype_clusters),
         Sex = factor(Sex),
         grp = factor(ifelse(PHF1_Otsu == "PHF1_pos", "PHF1+", "PHF1-"),
                      levels = c("PHF1-", "PHF1+")))
stopifnot(!anyNA(base$Sex), !anyNA(base$Age), !anyNA(base$PMI))
base <- base %>%
  mutate(Age_s = as.numeric(scale(Age)), PMI_s = as.numeric(scale(PMI)),
         plaque = factor(ifelse(Matched_4G8_40, "plaque-prox", "plaque-distal"),
                         levels = c("plaque-distal", "plaque-prox"))) %>%
  group_by(sample_id) %>%
  mutate(edge_dist = pmin(x - min(x), max(x) - x, y - min(y), max(y) - y)) %>%
  ungroup()

# --- the two DNA outcomes -------------------------------------------------------------------
# counts is steinbock's MEAN ion count per pixel inside the mask, and area is the mask's pixel
# count, so counts * area is the SUM over the mask. area px^2 == um^2 at 1 um/px.
ir_counts <- as.numeric(assay(spe, "counts")["191Ir", ])[match(base$cell_id, colnames(spe))]
stopifnot(all(ir_counts > 0), all(base$area > 0))
base$dna_density <- asinh(ir_counts)
# dna_total on the RAW multiplicative scale.
# UNITS: total 191Ir ION COUNTS per cell. assay(spe,"counts") holds steinbock's MEAN ion count per
# pixel within the mask (median 29.9, non-integer -- it is a mean, not a sum), and `area` is the
# pixel count at 1 um/pixel, so mean-counts-per-pixel x n-pixels = total counts and the pixel
# dimension cancels. NOT counts x um^2. Median ~1712 counts/cell.
# The coefficient is in raw counts per s.d. of log-distance, not a proportional change.
base$dna_total   <- ir_counts * base$area
stopifnot(all(MEASURES$col %in% colnames(base)))

anchors <- base %>% filter(PHF1_Otsu == "PHF1_pos",
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

exc <- base %>% filter(grepl("^Excitatory", celltype_clusters))

# --- centre every cell on its own donor's PHF1-NEGATIVE mean (see header) --------------------
donor_ref <- exc %>% filter(grp == "PHF1-") %>%
  group_by(patient_id) %>%
  summarise(across(all_of(MEASURES$col), ~ mean(.x, na.rm = TRUE), .names = "ref_{.col}"),
            .groups = "drop")
exc <- exc %>% inner_join(donor_ref, by = "patient_id")
for (mc in MEASURES$col) {
  exc[[paste0("c_", mc)]] <- if (CENTRE_ON_DONOR) exc[[mc]] - exc[[paste0("ref_", mc)]] else exc[[mc]]
}

# --- distance frame: PHF1- excitatory neurons only (PHF1+ cells are the anchors) -------------
tg <- exc %>% filter(grp == "PHF1-")
tg$dist_um <- nearest_dist(tg, anchors)
md <- tg %>%
  filter(!is.na(dist_um), dist_um <= MAX_DIST, edge_dist >= EDGE_BUFFER,
         sample_id %in% roi_ok) %>%
  mutate(dist_z = log(dist_um) / sd(log(dist_um)))
if (any(md$dist_um <= 0)) stop("non-positive distance; log() would fail")
DIST_SD <- sd(log(md$dist_um))

# -------------------------------------------------------------------
# Statistics
# -------------------------------------------------------------------
dist_one <- function(mc) {
  d <- md %>% transmute(value = .data[[mc]], dist_z, patient_id, Sex, Age_s, PMI_s)
  d <- d[is.finite(d$value), ]
  m <- lmerTest::lmer(value ~ dist_z + Sex + Age_s + PMI_s + (1 | patient_id),
                      data = d, REML = TRUE)
  co <- summary(m)$coefficients["dist_z", ]
  crit <- stats::qt(0.975, co["df"])
  sd_tot <- re_total_sd(m)                                   # R/imc_utils.R
  rv <- re_variance_summary(m); rv <- rv[rv$group == "patient_id", , drop = FALSE]
  tibble(measure = mc,
         beta_closer = -co["Estimate"], SE = co["Std. Error"], df = co["df"],
         CI.L = -co["Estimate"] - crit * co["Std. Error"],
         CI.R = -co["Estimate"] + crit * co["Std. Error"],
         cohens_d = -co["Estimate"] / sd_tot,
         pval = co["Pr(>|t|)"], n_cells = nrow(d),
         re_var_donor = if (nrow(rv)) rv$vcov[1] else NA_real_,
         re_pct_donor = if (nrow(rv)) rv$pct_of_total[1] else NA_real_,
         re_singular  = if (nrow(rv)) rv$is_singular[1] else NA)
}
res_dist <- bind_rows(lapply(MEASURES$col, dist_one)) %>%
  mutate(padj = p.adjust(pval, method = "BH")) %>%
  left_join(MEASURES, by = c("measure" = "col"))

contrast_one <- function(mc) {
  d <- exc %>% transmute(value = .data[[mc]], grp, patient_id) %>% filter(is.finite(value))
  dm <- d %>% group_by(patient_id, grp) %>% summarise(m = mean(value), .groups = "drop") %>%
    pivot_wider(names_from = grp, values_from = m) %>%
    filter(is.finite(`PHF1+`), is.finite(`PHF1-`))
  diff <- dm[["PHF1+"]] - dm[["PHF1-"]]
  pw <- suppressWarnings(wilcox.test(dm[["PHF1+"]], dm[["PHF1-"]], paired = TRUE))
  pt <- t.test(dm[["PHF1+"]], dm[["PHF1-"]], paired = TRUE)
  tibble(measure = mc, paired_diff = mean(diff),
         paired_CI.L = pt$conf.int[1], paired_CI.R = pt$conf.int[2],
         dz = mean(diff) / sd(diff), n_donors = nrow(dm),
         n_same_direction = sum(sign(diff) == sign(mean(diff))),
         pct_of_PHF1_neg = 100 * mean(diff) / mean(dm[["PHF1-"]]),
         wilcox_p = pw$p.value)
}
res_con <- bind_rows(lapply(MEASURES$col, contrast_one)) %>%
  mutate(padj = p.adjust(wilcox_p, method = "BH")) %>%
  left_join(MEASURES, by = c("measure" = "col"))

# Pooled (donor-ignoring) slope, purely as a DIAGNOSTIC for the drawn curve.
pooled_slope <- function(mc) {
  d <- md %>% transmute(value = .data[[mc]], dist_z) %>% filter(is.finite(value))
  co <- summary(stats::lm(value ~ dist_z, data = d))$coefficients["dist_z", ]
  tibble(measure = mc, pooled_beta_closer = -co[1], pooled_p = co[4])
}
res_dist <- res_dist %>% left_join(bind_rows(lapply(MEASURES$col, pooled_slope)), by = "measure")
sign_clash <- with(res_dist, sign(beta_closer) != sign(pooled_beta_closer) & padj < FDR)
if (any(sign_clash)) {
  warning(sprintf(paste("SIGN CLASH: pooled curve and within-donor model disagree for %s.",
                        "The drawn line slopes against its own significance encoding.",
                        "Set CENTRE_ON_DONOR <- TRUE, or read the curve as descriptive only."),
                  paste(res_dist$label[sign_clash], collapse = ", ")))
}

# -------------------------------------------------------------------
# Rolling mean + PHF1+ strip
# -------------------------------------------------------------------
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

# PHF1+ strip: mean + 95% CI of the DONOR-CENTRED value among PHF1+ cells. Because the centring
# is on each donor's PHF1- mean, this height IS the donor-paired PHF1+ minus PHF1- difference.
ref_df <- bind_rows(lapply(MEASURES$col, function(mc) {
  cc <- paste0("c_", mc)
  s <- exc[[cc]][exc$grp == "PHF1+" & is.finite(exc[[cc]])]
  n <- length(s); sem <- stats::sd(s) / sqrt(n)
  cs <- isTRUE(res_con$padj[res_con$measure == mc] < FDR)
  tibble(measure = mc, ref_mean = mean(s), ci_lo = mean(s) - 1.96 * sem,
         ci_hi = mean(s) + 1.96 * sem, n = n, con_sig = sig_factor(cs))
})) %>% left_join(MEASURES, by = c("measure" = "col"))

# -------------------------------------------------------------------
# Panel builder + fixed-geometry save, both copied from the reference script
# -------------------------------------------------------------------
save_fixed_panel <- function(p, path, panel_h = PANEL_H_MIN) {
  g <- ggplot2::ggplotGrob(p)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  prow <- unique(g$layout$t[grepl("^panel", g$layout$name)])
  # Setting the panel row to an absolute unit is what makes sum(g$heights) convertible: by
  # default it is unit(1, "null"), which has no inch value.
  g$heights[prow] <- grid::unit(panel_h, "in")
  # Measure everything that is NOT the panel column (axis, labels, both guide boxes, margins)
  # by zeroing the panel, then give the panel exactly the remainder of the 6.14 cm total.
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

# A guide key is built from LAYER DATA, so a level that never occurs (e.g. every series is
# significant) yields a key with its label but NO GLYPH. drop = FALSE and override.aes both fail
# to fix that on their own. Carrying one undrawn row per missing level gives the key its glyph;
# nothing is plotted because the value is NA. Same approach as R/phf1_channel_intensity_exc.R;
# the aes group = sig below is what stops the filler row joining a real line.
pad_levels <- function(d, col) {
  miss <- setdiff(c(SIG_BASIS, "ns"), as.character(d[[col]]))
  if (!length(miss)) return(d)
  filler <- d[rep(1L, length(miss)), , drop = FALSE]
  filler[[col]] <- factor(miss, levels = c(SIG_BASIS, "ns"))
  for (v in intersect(c("roll_mean", "sem", "ref_mean", "ci_lo", "ci_hi"), names(filler)))
    filler[[v]] <- NA_real_
  dplyr::bind_rows(d, filler)
}

build_panel <- function(mc, with_strip = TRUE) {
  rd <- pad_levels(roll_df %>% filter(measure == mc), "sig")
  rf <- pad_levels(ref_df  %>% filter(measure == mc), "con_sig")
  ylab <- YLAB_EXPR[[MEASURES$ylab[MEASURES$col == mc]]]

  p <- ggplot(rd, aes(dist_um, roll_mean, group = sig))
  # a zero rule is meaningful only on the centred scale; on raw um^2 it is far off-panel
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

  # pretty() gives 0/50/.../300 here, which collides at this panel width -- the axis is only
  # ~4.8 cm wide once the legend is placed. Fixed 100 um steps instead.
  brk <- seq(0, MAX_DIST, by = 100)
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
    # VERTICAL legend on the right, as in the reference figure. Note the cost at the required
    # 6.14 cm total width: the guide box takes 1.24 in, leaving 1.18 in (3.00 cm) of panel.
    # Shrinking keys and text recovers almost nothing (1.20 in) because the guide TITLES set the
    # box width. A bottom legend would leave 1.89 in (4.81 cm) of panel instead -- switch
    # legend.position back to "bottom" if the data area matters more than the arrangement.
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

slug_of <- function(mc) if (mc == "dna_density") "dna_density" else "dna_total"
for (mc in MEASURES$col) {
  sl <- slug_of(mc)
  save_fixed_panel(build_panel(mc, TRUE),
                   file.path(out_dir, sprintf("plot_%s_dist_phf1.pdf", sl)))
  save_fixed_panel(build_panel(mc, FALSE),
                   file.path(out_dir, sprintf("plot_%s_dist_phf1_nostrip.pdf", sl)))
}

# -------------------------------------------------------------------
# Source data + stats
# -------------------------------------------------------------------
wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE)
wr(roll_df %>% mutate(window_um = WINDOW_UM) %>%
     select(measure, label, window_um, significant, dist_um, roll_mean, sem, n_window),
   "source_data_dna_dist_phf1_rollmean.tsv")
wr(ref_df %>% select(measure, label, ref_mean, ci_lo, ci_hi, n, con_sig),
   "source_data_dna_dist_phf1_ref.tsv")
wr(res_dist, "source_data_dna_dist_phf1_distance_stats.tsv")
wr(res_con,  "source_data_dna_dist_phf1_contrast_stats.tsv")

sink(file.path(out_dir, "stats_dna_dist_phf1.txt"))
cat("191Ir DNA density / total DNA over PHF1 distance with the PHF1+ strip (IMC, AD cases)\n")
cat(strrep("=", 86), "\n", sep = "")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat(sprintf("SCALE: %s (CENTRE_ON_DONOR = %s).\n",
            if (CENTRE_ON_DONOR) "each cell relative to its own donor's PHF1-negative mean"
            else "RAW arcsinh 191Ir intensity / raw 191Ir x area, as the reference figure plots raw module scores",
            CENTRE_ON_DONOR))
cat("  The drawn curve POOLS donors; the model behind the line style is WITHIN donor. Those are\n")
cat("  different quantities. Both slopes are printed below and the script warns on a sign clash.\n\n")
cat("POOLED (donor-ignoring) slope of the drawn curve, per s.d. CLOSER:\n")
print(as.data.frame(res_dist %>% select(label, beta_closer, padj, pooled_beta_closer, pooled_p)),
      row.names = FALSE, digits = 3)
cat("\n")
cat("STRIP: PHF1+ neurons are the ANCHORS and have no distance to a nearest OTHER PHF1+ neuron,\n")
cat("  so they are drawn off-axis left of a separator rule. It is a REFERENCE LEVEL, not a\n")
cat("  distance-zero data point.\n\n")
cat("ENCODINGS: line style = log-distance verdict (solid padj<0.05, dashed ns);\n")
cat("  point shape = PHF1+ vs PHF1- verdict (filled padj<0.05, hollow ns).\n\n")
cat("DISTANCE MODEL (PHF1- cells, Set 3):\n")
cat("  value ~ dist_z + Sex + Age_s + PMI_s + (1 | patient_id), REML.\n")
cat(sprintf("  dist_z = log(dist)/sd(log(dist)); dist_sd = %.5f; cap %d um; per s.d. CLOSER.\n\n",
            DIST_SD, MAX_DIST))
print(as.data.frame(res_dist %>% select(label, beta_closer, CI.L, CI.R, cohens_d, pval, padj,
                                        n_cells)), row.names = FALSE, digits = 3)
cat("\nCONTRAST (paired Wilcoxon on donor means):\n")
print(as.data.frame(res_con %>% select(label, paired_diff, paired_CI.L, paired_CI.R, dz,
                                       n_same_direction, n_donors, pct_of_PHF1_neg,
                                       wilcox_p, padj)), row.names = FALSE, digits = 3)
cat(sprintf("\nDRAWN CURVE: %d um sliding window (half-window %.1f um) over a %d-point grid;\n",
            WINDOW_UM, HALF_WINDOW, GRID_N))
cat(sprintf("  ribbon mean +/- 1.96 SEM; windows holding fewer than %d cells left blank.\n",
            MIN_WINDOW_N))
cat("\nNOTES\n-----\n")
cat("- The PHF1+ vs PHF1- contrast is reported at DONOR level (paired Wilcoxon, dz and CI).\n")
cat("- 191Ir IS THE SEGMENTATION CHANNEL: masks are grown from the DNA signal, so mask geometry\n")
cat("  and DNA intensity are not independent measurements. Measured coupling is weak, though --\n")
cat("  cor(log counts, log area) = 0.06 across all cells.\n")
cat("- total DNA is area-weighted BY CONSTRUCTION (cor with log area = 0.68), so never adjust it\n")
cat("  for area: that regresses a variable on its own component and collapses to the density fit.\n")
cat("- PHF1+ nuclei are +26.7% larger, so a density change can be geometric. Read the two\n")
cat("  outcomes together, not separately.\n")
cat("- No cortical-depth term (Set 3).\n")
cat("\n=== Random effects (donor), distance model per measure ===\n")
cat("A variance at ~0 with re_singular TRUE means the model does NOT adjust for donor and is\n")
cat("equivalent to pooling.\n")
print(as.data.frame(res_dist %>% select(label, re_var_donor, re_pct_donor, re_singular)),
      row.names = FALSE, digits = 4)
cat("\n"); print(sessionInfo())
sink()

cat(sprintf("\n== distance (per s.d. CLOSER, %s) ==\n",
            if (CENTRE_ON_DONOR) "donor-centred units" else "raw units"))
print(as.data.frame(res_dist %>% select(label, beta_closer, CI.L, CI.R, padj, n_cells)),
      row.names = FALSE, digits = 3)
cat("\n== Random effects (donor), distance model per measure ==\n")
cat("A variance at ~0 with re_singular TRUE means the model does NOT adjust for donor.\n")
print(as.data.frame(res_dist %>% select(label, re_var_donor, re_pct_donor, re_singular)),
      row.names = FALSE, digits = 4)
cat("\n== PHF1+ strip (donor-paired PHF1+ minus PHF1-) ==\n")
print(as.data.frame(res_con %>% select(label, paired_diff, paired_CI.L, paired_CI.R, dz,
                                       n_same_direction, n_donors, padj)),
      row.names = FALSE, digits = 3)
