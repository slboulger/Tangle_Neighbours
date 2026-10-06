#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# phf1_module_distance_intensity_adjust.R
#
# Figure panels: 4B, S4A, S4B, S4C, S4D, S4E
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# phf1_module_distance_intensity_adjust.R
#
# Is the ptau-signature-over-distance effect DEPENDENT ON or INDEPENDENT OF the
# sub-threshold PHF1 protein signal, within PHF1-negative neurons?
#
# The module-score-vs-distance LMM (plot_phf1_module_vs_phf1_distance_modelp.R) shows
# the PHF1 / Otero-Garcia signatures rise near PHF1+ neurons. Those same PHF1-negative
# neurons also carry a graded SUB-THRESHOLD PHF1 intensity that decays with distance
# (the "halo"; phf1_intensity_gradient_neurons.R). This adds per-cell PHF1 intensity
# (phf1_intensity_p95_z) to the distance model and asks whether the distance effect
# survives adjustment.
#
# STATS (LOG distance only), per celltype x signature, on PHF1-negative cells within
# MAX_DIST with complete covariates AND non-NA intensity (complete-case; the unadjusted
# model is refit on the SAME cells as the adjusted, so they are comparable):
#   unadjusted: score ~ dist_scaled + nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)
#   adjusted  : score ~ dist_scaled + phf1_intensity_p95_z + <covariates> + (1|sample_id)
#   a-path    : phf1_intensity_p95_z ~ dist_scaled + <covariates> + (1|sample_id)
#   (dist_scaled = log(distance)/SD; score higher NEAR PHF1 => NEGATIVE slope)
#   distance slope collapses toward 0 / n.s. when intensity is added -> DEPENDENT on the
#     sub-threshold protein halo (intensity mediates). Slope persists -> INDEPENDENT.
#
# VISUALISATION (raw + CI, matching the module raw plots): per celltype x signature,
# the OBSERVED module score vs distance (points + 75 um rolling mean + 95% CI) against
# the INTENSITY-ADJUSTED score (observed minus the fitted intensity-linear contribution,
# b_hat * (intensity - mean); same scale). If the adjusted curve flattens near tangles,
# the transcriptomic gradient is carried by the protein halo.
#
# Additive: does NOT touch the existing module-distance outputs.
#
# NOTES (also written to the stats log):
#   * COLLINEARITY: strong dist/intensity correlation makes the slopes hard to separate
#     (inflated SE) -- the correlation is reported.
#   * MEASUREMENT ERROR: phf1_intensity_p95_z is a noisy tau proxy, so regressing it out
#     under-adjusts.
#   * Associational adjustment on PHF1-negative cells.
#
# OUTPUT (standard triple) under OUTPUT_DIR:
#   plot_distance_intensity_adjust_<ct>_<sig>.pdf      -- observed vs adjusted, points + rolling mean + 95% CI
#   source_data_distance_intensity_adjust_<ct>_<sig>_rollmean.tsv  -- rolling-mean grids (observed + adjusted)
#   source_data_distance_intensity_adjust_<ct>_<sig>_cells.tsv     -- per-cell rows behind the points
#   stats_distance_intensity_adjust.txt                -- log-LMM unadj/adj/a-path, collinearity, verdicts, sessionInfo
#   stats_distance_intensity_adjust_summary.tsv        -- structured per celltype x signature
#
# EXTENSIONS (estimates only; no mechanistic interpretation printed) -- sub-threshold
# intensity is itself negatively correlated with distance, so these quantify how the
# distance and intensity contributions overlap:
#   (1) conditional intensity coefficient (95% CI, sign), semi-partial R2 for distance and
#       intensity, and a commonality decomposition (unique-dist, unique-intensity, shared =
#       UNASSIGNABLE) -- added as columns to the summary.
#   (2) distance slope stratified by within-donor PHF1-intensity decile (slope vs decile;
#       bottom decile also in the summary): *_by_intensity_decile.{tsv,pdf}
#   (3) distance model on PHF1 intensity as OUTCOME in non-neuronal (tangle-incapable)
#       celltypes (NONNEURONAL_CELLTYPES) -- optical/segmentation spillover bound:
#       stats_distance_intensity_spillover_nonneuronal.tsv
#   (4) refits excluding cells within EXCLUSION_RADII um (slope vs radius): *_by_exclusion_radius.{tsv,pdf}
#   (5) adjusted refit with intensity as donor-scaled rank / natural spline (df=3) /
#       top-decile indicator -- adjusted distance slope columns in the summary.
#
# Run INTERACTIVELY on the HPC (seu_PHF1.rds): edit the Settings block, then source()/step through.

##  ............................................................................
##  Packages + setup                                                        ####
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(splines)   # ns() for the spline sensitivity refit
})

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")   # fig_theme, neuron_order, ...

##  ............................................................................
##  Settings -- edit these, then run interactively (source / step through)   ####
SEU_PATH      <- "seu_PHF1.rds"
CELLTYPES     <- c("Exc-IT-L2-3-CBLN2-HOPX", "Exc-IT-L3-5-CHGA-IL1RAPL2")  # eligible neuron subtypes
MARKERS_DIR   <- "phf1_markers"                             # <ct>/<ct>_phf1_geneset.txt (PHF1 signature)
OTERO_RDS     <- "otero_signatures/otero_at8_signatures.rds"
OTERO_UP_SET  <- "otero_L23_up"                             # Otero-Garcia Exc1/2 UP
INTENSITY_COL <- "phf1_intensity_p95_z"                     # per-cell PHF1 intensity (robust-z within sample)
MAX_DIST      <- 1000                                       # cap (match the module-distance analysis)
WINDOW_UM     <- 75                                         # rolling-mean window width (match the raw plots)
OUTPUT_DIR    <- "plots/phf1_distance_intensity_adjust_1000um"
SEED          <- 42
N_POINTS      <- 6000                                       # display subsample of observed points
MIN_GENES_PRESENT <- 10
CTRL          <- 100
NBIN          <- 24

# --- extension settings ---
NONNEURONAL_CELLTYPES <- c("Astro", "Oligo", "OPC", "Micro", "Endo", "VLMC")  # tangle-incapable (spillover bound)
EXCLUSION_RADII       <- c(0, 50, 100, 150, 200)            # refit after excluding cells within r um
N_DECILES             <- 10                                 # within-donor PHF1-intensity strata
SPLINE_DF             <- 3                                  # natural-spline df for the intensity sensitivity refit
MIN_CELLS_FIT         <- 50                                 # skip a stratified/excluded fit below this many cells
MIN_DONORS_FIT        <- 2                                  # ... or below this many donors

# shared covariate + random-effects string (reused by every fit)
COV_TERMS <- "nUMI_log + percent_neg + Sex + Age_s + PMI_s + (1 | sample_id)"

HALF_WINDOW <- WINDOW_UM / 2
set.seed(SEED)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
XLAB <- expression("Distance to PHF1+ neuron (" * mu * "m)")

SIG_DEF <- list(
  phf1     = list(label = "PHF1 signature",         colour = "#D55E00"),   # Okabe-Ito vermillion
  otero_up = list(label = "Otero-Garcia Exc1/2 UP", colour = "#0072B2")    # Okabe-Ito blue (CB-safe pair)
)
OBS_LAB <- "Observed"
ADJ_LAB <- "Adjusted (PHF1 intensity removed)"

##  ............................................................................
##  Helpers                                                                  ####
pick_col <- function(meta, candidates, what) {
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0) stop(sprintf("None of the %s columns (%s) found.", what, paste(candidates, collapse=", ")))
  hit[1]
}
read_geneset <- function(markers_dir, ct) {
  ct_safe <- gsub("[^A-Za-z0-9_-]", "_", ct)
  f <- file.path(markers_dir, ct_safe, paste0(ct_safe, "_phf1_geneset.txt"))
  if (!file.exists(f)) return(NULL)
  g <- trimws(readLines(f)); unique(g[nzchar(g)])
}
score_sets <- function(seu_ct, present_sets, seed) {
  set.seed(seed)
  seu_ct <- tryCatch(
    AddModuleScore(seu_ct, features = present_sets, name = "SIG_", assay = "SCT",
                   ctrl = CTRL, nbin = NBIN, seed = seed),
    error = function(e) {
      cat("  AddModuleScore failed at nbin=", NBIN, " (", conditionMessage(e), ") - retry nbin=15\n", sep="")
      set.seed(seed)
      AddModuleScore(seu_ct, features = present_sets, name = "SIG_", assay = "SCT",
                     ctrl = CTRL, nbin = 15, seed = seed)
    })
  cols <- paste0("SIG_", seq_along(present_sets))
  sc <- seu_ct@meta.data[, cols, drop = FALSE]; colnames(sc) <- names(present_sets); sc
}
term_row <- function(m, term) {
  if (is.null(m)) return(NULL)
  cf <- coef(summary(m)); if (!term %in% rownames(cf)) return(NULL)
  est <- cf[term, "Estimate"]; se <- cf[term, "Std. Error"]; p <- cf[term, "Pr(>|t|)"]
  list(est = est, se = se, ci_lo = est - 1.96*se, ci_hi = est + 1.96*se, p = p, singular = lme4::isSingular(m))
}
fit_lmm <- function(form, data)
  tryCatch(lmerTest::lmer(form, data = data, REML = TRUE),
           error = function(e) { message("    lmer failed: ", conditionMessage(e)); NULL })
# safe field extractor from a term_row() list (NA if the fit/term is missing)
ge <- function(x, f) if (is.null(x) || is.null(x[[f]])) NA_real_ else x[[f]]
# Nakagawa marginal R2 (fixed-effects variance explained) for a fitted lmer
r2m <- function(m) {
  if (is.null(m)) return(NA_real_)
  vf <- stats::var(as.numeric(stats::model.matrix(m) %*% lme4::fixef(m)))
  vc <- as.data.frame(lme4::VarCorr(m))
  vr <- sum(vc$vcov[vc$grp != "Residual"])
  ve <- vc$vcov[vc$grp == "Residual"]
  as.numeric(vf / (vf + vr + ve))
}
# Fitted-model VIF (performance::check_collinearity) for named terms; NA if unavailable.
vif_terms <- function(m, terms) {
  out <- setNames(rep(NA_real_, length(terms)), terms)
  if (is.null(m) || !requireNamespace("performance", quietly = TRUE)) return(out)
  cc <- tryCatch(as.data.frame(performance::check_collinearity(m)), error = function(e) NULL)
  if (is.null(cc) || !"Term" %in% names(cc)) return(out)
  vcol <- if ("VIF" %in% names(cc)) "VIF" else grep("VIF", names(cc), value = TRUE)[1]
  for (t in terms) { i <- which(cc$Term == t); if (length(i)) out[t] <- cc[[vcol]][i[1]] }
  out
}
# slope (+95% CI) vs an ordered x variable; df needs columns slope, ci_lo, ci_hi + xcol
plot_slope_vs_x <- function(df, xcol, xlab, ylab, colour) {
  ggplot(df, aes(.data[[xcol]], slope)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0, linewidth = 0.4, colour = colour) +
    geom_line(linewidth = 0.3, colour = colour) +
    geom_point(size = 1.2, colour = colour) +
    labs(x = xlab, y = ylab) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
}
# model-free rolling mean of y over x on a grid (mean, +/- 1 SEM; 95% CI = +/- 1.96 SEM)
roll_mean <- function(x, y, grid, half_window) {
  as.data.frame(do.call(rbind, lapply(grid, function(g) {
    idx <- which(x >= g - half_window & x <= g + half_window); n <- length(idx)
    if (n < 1) return(c(roll_mean = NA_real_, sem = NA_real_, n_window = 0))
    m <- mean(y[idx]); s <- if (n > 1) stats::sd(y[idx]) / sqrt(n) else NA_real_
    c(roll_mean = m, sem = s, n_window = n)
  })))
}
# much-darker shade of a colour (multiply RGB toward black)
darken <- function(col, f = 0.45) {
  v <- as.numeric(grDevices::col2rgb(col)) * f
  grDevices::rgb(v[1], v[2], v[3], maxColorValue = 255)
}
# PHF1+ reference band (mean + 95% CI of the mean) for one score vector
phf1_reference <- function(s) {
  s <- s[!is.na(s)]; n <- length(s)
  if (n < 1) return(NULL)
  m <- mean(s); sem <- if (n > 1) stats::sd(s) / sqrt(n) else NA_real_
  list(mean = m, ci_lo = m - 1.96 * sem, ci_hi = m + 1.96 * sem, n = n)
}
# fix panel width so x-axes align across plots; same sizing as the module raw plots
save_fixed_panel <- function(p, path) {
  g <- ggplot2::ggplotGrob(p)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  g$widths[pcol] <- grid::unit(1.6, "in")
  ggsave(path, g, width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
         height = 2.1, units = "in", device = "pdf")
}

# EXACT layout of plot_integrated_raw_dist_*.pdf (build_integrated, mode="raw"): NO raw
# points; rolling-mean lines + 95% CI ribbons; dotted PHF1+ reference band; colour legend
# on the right; solid lines (no significance solid/dashed -- that lives in the stats file).
# Series = Observed (signature colour) and Adjusted (much darker shade), adjusted drawn
# BEHIND the observed. The PHF1+ reference is the OBSERVED reference (Observed colour).
make_adjust_raw_plot <- function(roll_obs, roll_adj, ref, x_cap, module_colour, dark_colour, ylab) {
  ser_levels <- c(OBS_LAB, ADJ_LAB)
  ser_cols   <- setNames(c(module_colour, dark_colour), ser_levels)
  roll_obs$series <- factor(OBS_LAB, levels = ser_levels)
  roll_adj$series <- factor(ADJ_LAB, levels = ser_levels)
  ref_band <- data.frame(dist_to_phf1_um = c(0, x_cap), ref_mean = ref$mean,
                         ci_lo = ref$ci_lo, ci_hi = ref$ci_hi,
                         series = factor(OBS_LAB, levels = ser_levels))
  ggplot() +
    # PHF1+ reference (observed), dotted -- matches the integrated ref styling
    geom_ribbon(data = ref_band, aes(dist_to_phf1_um, ymin = ci_lo, ymax = ci_hi, fill = series),
                alpha = 0.10, colour = NA) +
    geom_line(data = ref_band, aes(dist_to_phf1_um, ref_mean, colour = series),
              linetype = "dotted", linewidth = 0.4) +
    # ADJUSTED behind observed
    geom_ribbon(data = roll_adj, aes(dist_to_phf1_um, ymin = roll_mean - 1.96 * sem,
                                     ymax = roll_mean + 1.96 * sem, fill = series), alpha = 0.15, colour = NA) +
    geom_line(data = roll_adj, aes(dist_to_phf1_um, roll_mean, colour = series), linewidth = 0.7) +
    # OBSERVED on top
    geom_ribbon(data = roll_obs, aes(dist_to_phf1_um, ymin = roll_mean - 1.96 * sem,
                                     ymax = roll_mean + 1.96 * sem, fill = series), alpha = 0.15, colour = NA) +
    geom_line(data = roll_obs, aes(dist_to_phf1_um, roll_mean, colour = series), linewidth = 0.7) +
    scale_colour_manual(values = ser_cols, name = NULL, drop = FALSE, limits = ser_levels) +
    scale_fill_manual(values = ser_cols, guide = "none", drop = FALSE) +
    labs(x = XLAB, y = ylab) +
    coord_cartesian(xlim = c(0, x_cap)) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.text = element_text(size = 6), legend.key.width = grid::unit(14, "pt"),
          legend.key.height = grid::unit(9, "pt"), legend.position = "right",
          legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
}

# COMBINED per-celltype panel: exactly the plot_integrated_raw layout (both signatures
# overlaid, colour = key, right legend, rolling means + 95% CI, dotted PHF1+ refs, no
# points), PLUS each signature's intensity-ADJUSTED line in a much darker shade (drawn
# behind observed). key_cols = named colours ordered signature-by-signature (obs, adj).
make_combined_adjust_plot <- function(roll_obs_all, roll_adj_all, ref_all, key_cols, x_cap, ylab) {
  lv <- names(key_cols)
  roll_obs_all$key <- factor(roll_obs_all$key, levels = lv)
  roll_adj_all$key <- factor(roll_adj_all$key, levels = lv)
  ref_all$key      <- factor(ref_all$key,      levels = lv)
  ggplot() +
    geom_ribbon(data = ref_all, aes(dist_to_phf1_um, ymin = ci_lo, ymax = ci_hi, fill = key),
                alpha = 0.10, colour = NA) +
    geom_line(data = ref_all, aes(dist_to_phf1_um, ref_mean, colour = key),
              linetype = "dotted", linewidth = 0.4) +
    geom_ribbon(data = roll_adj_all, aes(dist_to_phf1_um, ymin = roll_mean - 1.96 * sem,
                                         ymax = roll_mean + 1.96 * sem, fill = key), alpha = 0.15, colour = NA) +
    geom_line(data = roll_adj_all, aes(dist_to_phf1_um, roll_mean, colour = key), linewidth = 0.7) +
    geom_ribbon(data = roll_obs_all, aes(dist_to_phf1_um, ymin = roll_mean - 1.96 * sem,
                                         ymax = roll_mean + 1.96 * sem, fill = key), alpha = 0.15, colour = NA) +
    geom_line(data = roll_obs_all, aes(dist_to_phf1_um, roll_mean, colour = key), linewidth = 0.7) +
    scale_colour_manual(values = key_cols, name = NULL, drop = FALSE, limits = lv) +
    scale_fill_manual(values = key_cols, guide = "none", drop = FALSE) +
    labs(x = XLAB, y = ylab) +
    coord_cartesian(xlim = c(0, x_cap)) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.text = element_text(size = 6), legend.key.width = grid::unit(14, "pt"),
          legend.key.height = grid::unit(9, "pt"), legend.position = "right",
          legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
}

# ADJUSTED-ONLY variant: PHF1+ dotted reference + ADJUSTED rolling mean + 95% CI only
# (observed rolling omitted). Same integrated_raw layout; works for one signature (1 key)
# or the combined panel (>=2 keys). key_cols maps the adjusted-series label(s) to colour(s).
make_adjusted_only_plot <- function(roll_adj_all, ref_all, key_cols, x_cap, ylab) {
  lv <- names(key_cols)
  roll_adj_all$key <- factor(roll_adj_all$key, levels = lv)
  ref_all$key      <- factor(ref_all$key,      levels = lv)
  ggplot() +
    geom_ribbon(data = ref_all, aes(dist_to_phf1_um, ymin = ci_lo, ymax = ci_hi, fill = key),
                alpha = 0.10, colour = NA) +
    geom_line(data = ref_all, aes(dist_to_phf1_um, ref_mean, colour = key),
              linetype = "dotted", linewidth = 0.4) +
    geom_ribbon(data = roll_adj_all, aes(dist_to_phf1_um, ymin = roll_mean - 1.96 * sem,
                                         ymax = roll_mean + 1.96 * sem, fill = key), alpha = 0.15, colour = NA) +
    geom_line(data = roll_adj_all, aes(dist_to_phf1_um, roll_mean, colour = key), linewidth = 0.7) +
    scale_colour_manual(values = key_cols, name = NULL, drop = FALSE, limits = lv) +
    scale_fill_manual(values = key_cols, guide = "none", drop = FALSE) +
    labs(x = XLAB, y = ylab) +
    coord_cartesian(xlim = c(0, x_cap)) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.text = element_text(size = 6), legend.key.width = grid::unit(14, "pt"),
          legend.key.height = grid::unit(9, "pt"), legend.position = "right",
          legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
}

# Forest plot of the distance-slope beta (per SD log-dist) with 95% CI: x = beta,
# y = series (PHF1 / PHF1 adj. / Otero / Otero adj.), circle + horizontal CI. Same
# theme + colours as the score-over-distance plots; y labels identify the series so
# no colour legend. key_cols orders the series (top-to-bottom via reversed levels).
make_forest_plot <- function(fdf, key_cols, xlab) {
  lv <- names(key_cols)
  fdf$series <- factor(fdf$series, levels = rev(lv))       # first entry at the TOP
  ggplot(fdf, aes(beta, series, colour = series)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_segment(aes(x = ci_lo, xend = ci_hi, y = series, yend = series), linewidth = 0.4) +
    geom_point(shape = 16, size = 2.2) +
    scale_colour_manual(values = key_cols, guide = "none", drop = FALSE) +
    labs(x = xlab, y = NULL) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.text.y = element_text(size = 6),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
}

##  ............................................................................
##  Load object, checks                                                      ####
cat("Loading Seurat object (this is large)...\n")
seu <- readRDS(SEU_PATH)
DefaultAssay(seu) <- "SCT"
stopifnot("celltype missing" = "celltype" %in% colnames(seu@meta.data),
          "PHF1 missing"      = "PHF1"     %in% colnames(seu@meta.data),
          "dist_to_phf1_um missing" = "dist_to_phf1_um" %in% colnames(seu@meta.data),
          "sample_id missing" = "sample_id" %in% colnames(seu@meta.data))
if (!INTENSITY_COL %in% colnames(seu@meta.data))
  stop(sprintf("Intensity column '%s' not on the object. Run add_phf1_intensity.R against %s first.",
               INTENSITY_COL, SEU_PATH))
sex_col <- pick_col(seu@meta.data, c("Sex"), "Sex")
age_col <- pick_col(seu@meta.data, c("Age"), "Age")
pmi_col <- pick_col(seu@meta.data, c("PMI", "PostMortemInterval", "PostMortem_Interval"), "PMI")
panel   <- rownames(seu[["SCT"]])
otero   <- readRDS(OTERO_RDS)
stopifnot("Otero UP set not found" = OTERO_UP_SET %in% names(otero))

##  ............................................................................
##  Score each celltype (both signatures)                                    ####
CT_DAT <- list()
for (ct in CELLTYPES) {
  if (!ct %in% seu$celltype) { cat("[skip]", ct, "not in celltype\n"); next }
  phf1_genes  <- intersect(read_geneset(MARKERS_DIR, ct) %||% character(0), panel)
  otero_genes <- intersect(unique(otero[[OTERO_UP_SET]]), panel)
  present <- list()
  if (length(phf1_genes)  >= MIN_GENES_PRESENT) present$phf1     <- phf1_genes
  if (length(otero_genes) >= MIN_GENES_PRESENT) present$otero_up <- otero_genes
  if (length(present) == 0) { cat("[skip]", ct, "no usable signature\n"); next }
  cat(sprintf("Scoring %s (%s)...\n", ct,
              paste(sprintf("%s=%d", names(present), lengths(present)), collapse=", ")))
  seu_ct <- subset(seu, celltype == ct); DefaultAssay(seu_ct) <- "SCT"
  sc <- score_sets(seu_ct, present, SEED)
  md <- seu_ct@meta.data
  id_cands <- c("cell_ID","cell_id","CellID"); id_col <- id_cands[id_cands %in% colnames(md)]
  cell_id <- if (length(id_col) > 0) as.character(md[[id_col[1]]]) else rownames(md)
  base <- data.frame(
    cell_id = cell_id, phf1_pos = as.logical(as.character(md$PHF1)),
    dist = as.numeric(md$dist_to_phf1_um), intensity = as.numeric(md[[INTENSITY_COL]]),
    percent_neg = as.numeric(md$percent.neg), nUMI_log = log2(as.numeric(md$nCount_RNA) + 1),
    sample_id = as.character(md$sample_id), Sex = as.character(md[[sex_col]]),
    Age = as.numeric(md[[age_col]]), PMI = as.numeric(md[[pmi_col]]), stringsAsFactors = FALSE)
  CT_DAT[[ct]] <- list(df = cbind(base, sc), sigs = names(present))
}

# Non-neuronal (tangle-incapable) celltypes: pull the same per-cell frame (NO module
# score needed -- intensity is the outcome). Extracted here, while seu is still loaded.
id_cands <- c("cell_ID", "cell_id", "CellID")
NN_DAT <- list()
for (ct in NONNEURONAL_CELLTYPES) {
  if (!ct %in% seu$celltype) { cat("[skip NN]", ct, "not present\n"); next }
  seu_ct <- subset(seu, celltype == ct); md <- seu_ct@meta.data
  id_col  <- id_cands[id_cands %in% colnames(md)]
  cell_id <- if (length(id_col) > 0) as.character(md[[id_col[1]]]) else rownames(md)
  NN_DAT[[ct]] <- data.frame(
    cell_id = cell_id, phf1_pos = as.logical(as.character(md$PHF1)),
    dist = as.numeric(md$dist_to_phf1_um), intensity = as.numeric(md[[INTENSITY_COL]]),
    percent_neg = as.numeric(md$percent.neg), nUMI_log = log2(as.numeric(md$nCount_RNA) + 1),
    sample_id = as.character(md$sample_id), Sex = as.character(md[[sex_col]]),
    Age = as.numeric(md[[age_col]]), PMI = as.numeric(md[[pmi_col]]), stringsAsFactors = FALSE)
}
rm(seu); invisible(gc())
if (length(CT_DAT) == 0) stop("No celltype scored.")

##  ............................................................................
##  Fit (log distance) + build observed/adjusted rolling curves             ####
summ_rows <- list()

for (ct in names(CT_DAT)) {
  df0 <- CT_DAT[[ct]]$df
  keep <- !df0$phf1_pos & !is.na(df0$dist) & df0$dist > 0 & df0$dist <= MAX_DIST &
          !is.na(df0$intensity) & !is.na(df0$percent_neg) & !is.na(df0$nUMI_log) &
          !is.na(df0$Age) & !is.na(df0$PMI) & !is.na(df0$Sex)
  d <- df0[keep, ]
  d$Sex <- droplevels(factor(d$Sex)); d$sample_id <- factor(d$sample_id)
  d$Age_s <- as.numeric(scale(d$Age)); d$PMI_s <- as.numeric(scale(d$PMI))
  dist_sd <- stats::sd(log(d$dist)); d$dist_scaled <- log(d$dist) / dist_sd
  n_na_int <- sum(!df0$phf1_pos & !is.na(df0$dist) & df0$dist > 0 & df0$dist <= MAX_DIST & is.na(df0$intensity))
  cat(sprintf("\n%s: PHF1-negative modelled cells (intensity non-NA) = %d (dropped %d NA-intensity), donors = %d\n",
              ct, nrow(d), n_na_int, dplyr::n_distinct(d$sample_id)))

  x_cap <- MAX_DIST; grid <- seq(0, x_cap, length.out = 200)
  cov_terms <- COV_TERMS

  # accumulators for the COMBINED per-celltype panel (both signatures, obs + adj)
  comb_obs <- list(); comb_adj <- list(); comb_ref <- list(); comb_cols <- c()
  comb_adjonly <- list(); comb_refadjonly <- list(); comb_adjonly_cols <- c()  # adjusted-only (plain label, dark colour)
  comb_forest <- list(); comb_forest_cols <- c()   # beta +/- CI per series (obs + adj) for the forest

  for (sig in CT_DAT[[ct]]$sigs) {
    sub <- d; sub$score <- sub[[sig]]
    m_unadj <- fit_lmm(as.formula(paste("score ~ dist_scaled +", cov_terms)), sub)
    m_adj   <- fit_lmm(as.formula(paste("score ~ dist_scaled + intensity +", cov_terms)), sub)
    m_apath <- fit_lmm(as.formula(paste("intensity ~ dist_scaled +", cov_terms)), sub)
    u <- term_row(m_unadj, "dist_scaled"); a <- term_row(m_adj, "dist_scaled")
    b <- term_row(m_adj, "intensity"); ap <- term_row(m_apath, "dist_scaled")
    if (is.null(u) || is.null(a) || is.null(b)) { cat("  ", sig, ": model failed\n"); next }
    pct_atten <- if (abs(u$est) < 1e-9) NA_real_ else 100 * (1 - a$est / u$est)
    collin <- suppressWarnings(cor(sub$dist_scaled, sub$intensity, use = "complete.obs"))
    vif <- vif_terms(m_adj, c("dist_scaled", "intensity"))   # (6) fitted-model VIF, adjusted fit
    # (7) verdict on the SIGNED attenuation: pct_atten<0 => adjusted slope STEEPER (suppression),
    # not attenuation. Positive => attenuated (weaker after adjustment).
    verdict <- if (is.na(a$p) || is.na(pct_atten)) "not testable"
      else if (pct_atten < -10) "distance effect STEEPENED by adjustment (suppression)"
      else if (pct_atten >= 50 && a$p >= 0.05) "distance effect largely DEPENDENT on intensity"
      else if (pct_atten <= 25 && a$p < 0.05) "distance effect largely INDEPENDENT of intensity"
      else "partial attenuation (see slopes/CI)"
    sig_col  <- SIG_DEF[[sig]]$colour

    # within-donor intensity transforms (used by the decile and sensitivity additions)
    sub <- sub %>% group_by(sample_id) %>%
      mutate(int_rank = dplyr::percent_rank(intensity),
             int_decile = dplyr::ntile(intensity, N_DECILES)) %>%
      ungroup() %>% as.data.frame()
    sub$int_topdecile <- as.integer(sub$int_decile == N_DECILES)

    ## --- (1) conditional intensity coefficient, semi-partial R2, commonality ---
    # reuse the same formula/fits; intensity-only and covariate-only models complete the
    # 2x2 needed for the variance partition (dist-only = m_unadj, full = m_adj).
    m_int <- fit_lmm(as.formula(paste("score ~ intensity +", cov_terms)), sub)
    m_cov <- fit_lmm(as.formula(paste("score ~", cov_terms)), sub)
    r2_full <- r2m(m_adj); r2_dist <- r2m(m_unadj); r2_int <- r2m(m_int); r2_cov <- r2m(m_cov)
    dR2_full <- r2_full - r2_cov; dR2_dist <- r2_dist - r2_cov; dR2_int <- r2_int - r2_cov
    sp_r2_dist <- dR2_full - dR2_int        # unique distance  (semi-partial R2)
    sp_r2_int  <- dR2_full - dR2_dist       # unique intensity (semi-partial R2)
    shared_unassign <- dR2_full - sp_r2_dist - sp_r2_int   # shared -> UNASSIGNABLE

    ## --- (5) adjusted refit with intensity as rank / natural spline / top-decile indicator ---
    a_rank <- term_row(fit_lmm(as.formula(paste("score ~ dist_scaled + int_rank +", cov_terms)), sub), "dist_scaled")
    a_spl  <- term_row(fit_lmm(as.formula(paste0("score ~ dist_scaled + ns(intensity, df = ", SPLINE_DF, ") + ", cov_terms)), sub), "dist_scaled")
    a_top  <- term_row(fit_lmm(as.formula(paste("score ~ dist_scaled + int_topdecile +", cov_terms)), sub), "dist_scaled")

    ## --- (2) distance slope stratified by within-donor PHF1-intensity decile ---
    dec_df <- bind_rows(lapply(sort(unique(sub$int_decile)), function(dec) {
      sdd <- sub[sub$int_decile == dec, ]
      if (nrow(sdd) < MIN_CELLS_FIT || dplyr::n_distinct(sdd$sample_id) < MIN_DONORS_FIT) return(NULL)
      tr <- term_row(fit_lmm(as.formula(paste("score ~ dist_scaled +", cov_terms)), sdd), "dist_scaled")
      if (is.null(tr)) return(NULL)
      data.frame(celltype = ct, signature = sig, decile = dec, slope = tr$est,
                 ci_lo = tr$ci_lo, ci_hi = tr$ci_hi, p = tr$p, n = nrow(sdd), singular = tr$singular)
    }))
    if (nrow(dec_df) > 0) {
      write.table(dec_df, file.path(OUTPUT_DIR, sprintf("stats_distance_intensity_adjust_%s_%s_by_intensity_decile.tsv", ct, sig)),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      save_fixed_panel(plot_slope_vs_x(dec_df, "decile", "Within-donor PHF1-intensity decile",
                                       "Distance slope (per SD log-dist)", sig_col),
                       file.path(OUTPUT_DIR, sprintf("plot_distance_intensity_adjust_%s_%s_by_intensity_decile.pdf", ct, sig)))
    }
    dec1 <- if (nrow(dec_df) > 0) dec_df[dec_df$decile == 1, , drop = FALSE] else dec_df  # bottom decile

    ## --- (4) refits after excluding cells within each exclusion radius ---
    excl_df <- bind_rows(lapply(EXCLUSION_RADII, function(r) {
      sr <- sub[sub$dist > r, ]
      if (nrow(sr) < MIN_CELLS_FIT || dplyr::n_distinct(sr$sample_id) < MIN_DONORS_FIT) return(NULL)
      sdl <- stats::sd(log(sr$dist)); sr$dist_scaled <- log(sr$dist) / sdl
      tr <- term_row(fit_lmm(as.formula(paste("score ~ dist_scaled + intensity +", cov_terms)), sr), "dist_scaled")
      if (is.null(tr)) return(NULL)
      data.frame(celltype = ct, signature = sig, exclusion_radius_um = r, n = nrow(sr), dist_sd = sdl,
                 slope_per_sd = tr$est, ci_lo = tr$ci_lo, ci_hi = tr$ci_hi, p = tr$p,
                 slope_per_logunit = tr$est / sdl, ci_lo_logunit = tr$ci_lo / sdl, ci_hi_logunit = tr$ci_hi / sdl)
    }))
    if (nrow(excl_df) > 0) {
      write.table(excl_df, file.path(OUTPUT_DIR, sprintf("stats_distance_intensity_adjust_%s_%s_by_exclusion_radius.tsv", ct, sig)),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      pe <- excl_df %>% transmute(exclusion_radius_um, slope = slope_per_logunit, ci_lo = ci_lo_logunit, ci_hi = ci_hi_logunit)
      save_fixed_panel(plot_slope_vs_x(pe, "exclusion_radius_um", expression("Exclusion radius (" * mu * "m)"),
                                       "Distance slope (per log-um)", sig_col),
                       file.path(OUTPUT_DIR, sprintf("plot_distance_intensity_adjust_%s_%s_by_exclusion_radius.pdf", ct, sig)))
    }

    # intensity-adjusted score (same scale): remove the fitted intensity-linear contribution
    sub$score_adj <- sub$score - b$est * (sub$intensity - mean(sub$intensity))

    # observed + adjusted rolling means over raw distance
    ro <- roll_mean(sub$dist, sub$score,     grid, HALF_WINDOW)
    ra <- roll_mean(sub$dist, sub$score_adj, grid, HALF_WINDOW)
    roll_obs <- data.frame(dist_to_phf1_um = grid, ro)
    roll_adj <- data.frame(dist_to_phf1_um = grid, ra)

    # PHF1+ reference band (observed module score of manual PHF1+ neurons)
    ref <- phf1_reference(df0[[sig]][df0$phf1_pos])

    sig_col  <- SIG_DEF[[sig]]$colour
    dark_col <- darken(sig_col, 0.45)
    # per-signature plot: exact plot_integrated_raw layout (no points, right legend, rolling
    # means + 95% CI, dotted PHF1+ ref). Stats (significance/attenuation) stay in the stats file.
    p <- make_adjust_raw_plot(roll_obs, roll_adj, ref, x_cap,
                              module_colour = sig_col, dark_colour = dark_col, ylab = "Module score")
    save_fixed_panel(p, file.path(OUTPUT_DIR, sprintf("plot_distance_intensity_adjust_%s_%s.pdf", ct, sig)))

    # accumulate for the combined per-celltype panel
    lab_o <- SIG_DEF[[sig]]$label; lab_a <- paste0(lab_o, " (adjusted)")
    comb_obs[[sig]] <- roll_obs %>% mutate(key = lab_o)
    comb_adj[[sig]] <- roll_adj %>% mutate(key = lab_a)
    comb_ref[[sig]] <- data.frame(dist_to_phf1_um = c(0, x_cap), ref_mean = ref$mean,
                                  ci_lo = ref$ci_lo, ci_hi = ref$ci_hi, key = lab_o)
    comb_cols[lab_o] <- sig_col; comb_cols[lab_a] <- dark_col
    # forest rows: observed (labelled "(unadjusted)" here only) + adjusted distance-slope
    # beta +/- 95% CI (per SD). Forest-specific labels/colours so the combined/adjusted-only
    # plots (which reuse lab_o) are unaffected.
    lab_u <- paste0(lab_o, " (unadjusted)")
    comb_forest[[paste(sig, "obs")]] <- data.frame(series = lab_u, beta = u$est, ci_lo = u$ci_lo, ci_hi = u$ci_hi)
    comb_forest[[paste(sig, "adj")]] <- data.frame(series = lab_a, beta = a$est, ci_lo = a$ci_lo, ci_hi = a$ci_hi)
    comb_forest_cols[lab_u] <- sig_col; comb_forest_cols[lab_a] <- dark_col
    # adjusted-only: PLAIN signature label (no "(adjusted)" suffix), dark colour. Keyed to
    # lab_o in its own accumulators so it doesn't clash with the full obs+adj combined.
    roll_adjonly <- roll_adj %>% mutate(key = lab_o)
    ref_adjonly  <- data.frame(dist_to_phf1_um = c(0, x_cap), ref_mean = ref$mean,
                               ci_lo = ref$ci_lo, ci_hi = ref$ci_hi, key = lab_o)
    comb_adjonly[[sig]] <- roll_adjonly; comb_refadjonly[[sig]] <- ref_adjonly
    comb_adjonly_cols[lab_o] <- dark_col
    # per-signature ADJUSTED-ONLY plot (adjusted line + PHF1+ dotted ref; observed omitted)
    save_fixed_panel(
      make_adjusted_only_plot(roll_adjonly, ref_adjonly, setNames(dark_col, lab_o), x_cap, "Module score"),
      file.path(OUTPUT_DIR, sprintf("plot_distance_intensity_adjust_%s_%s_adjustedonly.pdf", ct, sig)))

    roll_long <- rbind(
      data.frame(series = OBS_LAB, dist_to_phf1_um = grid, ro),
      data.frame(series = ADJ_LAB, dist_to_phf1_um = grid, ra))
    write.table(roll_long %>% mutate(celltype = ct, signature = sig, window_um = WINDOW_UM) %>%
                  dplyr::select(celltype, signature, series, window_um, dist_to_phf1_um, roll_mean, sem, n_window),
                file.path(OUTPUT_DIR, sprintf("source_data_distance_intensity_adjust_%s_%s_rollmean.tsv", ct, sig)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(sub %>% transmute(celltype = ct, signature = sig, cell_id,
                                  sample_id = as.character(sample_id), dist_to_phf1_um = dist,
                                  score, score_adj, phf1_intensity_p95_z = intensity),
                file.path(OUTPUT_DIR, sprintf("source_data_distance_intensity_adjust_%s_%s_cells.tsv", ct, sig)),
                sep = "\t", quote = FALSE, row.names = FALSE)

    summ_rows[[paste(ct, sig)]] <- data.frame(
      celltype = ct, signature = sig, sig_label = SIG_DEF[[sig]]$label, n_cells = nrow(sub),
      n_intensity_dropped = n_na_int, dist_slope_unadj = u$est, dist_p_unadj = u$p,
      dist_slope_adj = a$est, dist_p_adj = a$p, pct_attenuation = pct_atten,
      intensity_slope_adj = b$est, intensity_p_adj = b$p,           # b-path
      apath_dist_slope = if (is.null(ap)) NA_real_ else ap$est,     # a-path: dist->intensity
      apath_dist_p     = if (is.null(ap)) NA_real_ else ap$p,
      dist_intensity_cor = collin, singular_adj = a$singular, verdict = verdict,
      # (1) conditional intensity coefficient (95% CI, sign) + variance partition
      cond_intensity_coef = b$est, cond_intensity_ci_lo = b$ci_lo, cond_intensity_ci_hi = b$ci_hi,
      cond_intensity_sign = sign(b$est),
      sp_r2_dist = sp_r2_dist, sp_r2_intensity = sp_r2_int,
      commonality_unique_dist = sp_r2_dist, commonality_unique_intensity = sp_r2_int,
      commonality_shared_unassignable = shared_unassign,
      r2m_full = r2_full, r2m_covariates = r2_cov,
      # (5) adjusted distance slope under alternative intensity representations
      dist_slope_adj_rank = ge(a_rank, "est"), dist_p_adj_rank = ge(a_rank, "p"),
      dist_ci_lo_adj_rank = ge(a_rank, "ci_lo"), dist_ci_hi_adj_rank = ge(a_rank, "ci_hi"),
      dist_slope_adj_spline = ge(a_spl, "est"), dist_p_adj_spline = ge(a_spl, "p"),
      dist_ci_lo_adj_spline = ge(a_spl, "ci_lo"), dist_ci_hi_adj_spline = ge(a_spl, "ci_hi"),
      dist_slope_adj_topdecile = ge(a_top, "est"), dist_p_adj_topdecile = ge(a_top, "p"),
      dist_ci_lo_adj_topdecile = ge(a_top, "ci_lo"), dist_ci_hi_adj_topdecile = ge(a_top, "ci_hi"),
      # (2) bottom (lowest-intensity) within-donor decile, reported separately
      dist_slope_decile1 = if (nrow(dec1)) dec1$slope[1] else NA_real_,
      dist_ci_lo_decile1 = if (nrow(dec1)) dec1$ci_lo[1] else NA_real_,
      dist_ci_hi_decile1 = if (nrow(dec1)) dec1$ci_hi[1] else NA_real_,
      dist_p_decile1     = if (nrow(dec1)) dec1$p[1] else NA_real_,
      # (6) fitted-model VIF for the adjusted fit (vs the bivariate dist_intensity_cor)
      vif_dist_adj = unname(vif["dist_scaled"]), vif_intensity_adj = unname(vif["intensity"]),
      # (5) dist_sd + CI for the headline slopes, and every per-SD slope also per-log-unit
      # (slope / dist_sd), so slopes are comparable across celltypes and across scripts.
      dist_sd = dist_sd,
      dist_slope_unadj_ci_lo = u$ci_lo, dist_slope_unadj_ci_hi = u$ci_hi,
      dist_slope_adj_ci_lo   = a$ci_lo, dist_slope_adj_ci_hi   = a$ci_hi,
      dist_slope_unadj_perlog = u$est / dist_sd,
      dist_ci_lo_unadj_perlog = u$ci_lo / dist_sd, dist_ci_hi_unadj_perlog = u$ci_hi / dist_sd,
      dist_slope_adj_perlog   = a$est / dist_sd,
      dist_ci_lo_adj_perlog   = a$ci_lo / dist_sd, dist_ci_hi_adj_perlog   = a$ci_hi / dist_sd,
      apath_dist_slope_perlog = if (is.null(ap)) NA_real_ else ap$est / dist_sd,
      dist_slope_adj_rank_perlog      = ge(a_rank, "est") / dist_sd,
      dist_slope_adj_spline_perlog    = ge(a_spl,  "est") / dist_sd,
      dist_slope_adj_topdecile_perlog = ge(a_top,  "est") / dist_sd)
    cat(sprintf("  %s: dist slope %.4f -> %.4f (%.0f%% atten), adj p=%.2g | %s\n",
                sig, u$est, a$est, pct_atten, a$p, verdict))
  }

  # COMBINED per-celltype panel: both signatures, observed + adjusted (exact integrated_raw look)
  if (length(comb_obs) > 0) {
    roll_obs_all <- bind_rows(comb_obs); roll_adj_all <- bind_rows(comb_adj)
    ref_all      <- bind_rows(comb_ref)
    pc <- make_combined_adjust_plot(roll_obs_all, roll_adj_all, ref_all, comb_cols, x_cap, "Module score")
    save_fixed_panel(pc, file.path(OUTPUT_DIR, sprintf("plot_distance_intensity_adjust_combined_%s.pdf", ct)))
    write.table(rbind(roll_obs_all, roll_adj_all) %>% mutate(celltype = ct, window_um = WINDOW_UM) %>%
                  dplyr::select(celltype, key, window_um, dist_to_phf1_um, roll_mean, sem, n_window),
                file.path(OUTPUT_DIR, sprintf("source_data_distance_intensity_adjust_combined_%s_rollmean.tsv", ct)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    # ADJUSTED-ONLY combined (both signatures' adjusted lines + PHF1+ dotted refs; observed
    # omitted). Plain signature labels in the legend (no "(adjusted)" suffix), dark colours.
    pc_adj <- make_adjusted_only_plot(bind_rows(comb_adjonly), bind_rows(comb_refadjonly),
                                      comb_adjonly_cols, x_cap, "Module score")
    save_fixed_panel(pc_adj, file.path(OUTPUT_DIR, sprintf("plot_distance_intensity_adjust_combined_%s_adjustedonly.pdf", ct)))

    # FOREST: distance-slope beta +/- 95% CI (per SD log-dist), y = series (obs + adj),
    # same colours as the score-over-distance plots (obs = signature colour, adj = dark).
    fdf <- bind_rows(comb_forest)
    save_fixed_panel(make_forest_plot(fdf, comb_forest_cols,
                                      expression(paste(beta, ", per S.D. log-distance"))),
                     file.path(OUTPUT_DIR, sprintf("plot_distance_intensity_adjust_forest_%s.pdf", ct)))
    write.table(fdf %>% mutate(celltype = ct) %>% dplyr::select(celltype, series, beta, ci_lo, ci_hi),
                file.path(OUTPUT_DIR, sprintf("source_data_distance_intensity_adjust_forest_%s.tsv", ct)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  combined panel %s: %d signatures\n", ct, length(comb_obs)))
  }
}
summ_df <- bind_rows(summ_rows)

##  ............................................................................
##  (3) PHF1 intensity ~ distance in non-neuronal (tangle-incapable) celltypes ####
##  Bound on optical / segmentation spillover from an adjacent tangle. Same distance
##  model + covariates + random effects as above, intensity as the outcome.
nn_df <- bind_rows(lapply(names(NN_DAT), function(ct) {
  df0 <- NN_DAT[[ct]]
  keep <- !df0$phf1_pos & !is.na(df0$dist) & df0$dist > 0 & df0$dist <= MAX_DIST &
          !is.na(df0$intensity) & !is.na(df0$percent_neg) & !is.na(df0$nUMI_log) &
          !is.na(df0$Age) & !is.na(df0$PMI) & !is.na(df0$Sex)
  d <- df0[keep, ]
  if (nrow(d) < MIN_CELLS_FIT || dplyr::n_distinct(d$sample_id) < MIN_DONORS_FIT) return(NULL)
  d$Sex <- droplevels(factor(d$Sex)); d$sample_id <- factor(d$sample_id)
  d$Age_s <- as.numeric(scale(d$Age)); d$PMI_s <- as.numeric(scale(d$PMI))
  sdl <- stats::sd(log(d$dist)); d$dist_scaled <- log(d$dist) / sdl
  tr <- term_row(fit_lmm(as.formula(paste("intensity ~ dist_scaled +", COV_TERMS)), d), "dist_scaled")
  if (is.null(tr)) return(NULL)
  data.frame(celltype = ct, n_cells = nrow(d), donors = dplyr::n_distinct(d$sample_id), dist_sd = sdl,
             intensity_dist_slope_per_sd = tr$est, ci_lo = tr$ci_lo, ci_hi = tr$ci_hi, p = tr$p,
             intensity_dist_slope_per_logunit = tr$est / sdl, singular = tr$singular)
}))
if (!is.null(nn_df) && nrow(nn_df) > 0)
  write.table(nn_df, file.path(OUTPUT_DIR, "stats_distance_intensity_spillover_nonneuronal.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

##  ............................................................................
##  Stats tables                                                             ####
write.table(summ_df, file.path(OUTPUT_DIR, "stats_distance_intensity_adjust_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(OUTPUT_DIR, "stats_distance_intensity_adjust.txt"))
cat("Ptau-signature-over-distance: dependence on sub-threshold PHF1 intensity (LOG distance)\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Population: PHF1-NEGATIVE neurons within", MAX_DIST, "um, complete covariates AND non-NA",
    INTENSITY_COL, "(complete-case; unadjusted refit on the SAME cells as adjusted).\n")
cat("Models (per celltype x signature):\n")
cat("  unadjusted: score ~ dist_scaled + nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)\n")
cat("  adjusted  : score ~ dist_scaled +", INTENSITY_COL, "+ <covariates> + (1|sample_id)\n")
cat("  a-path    :", INTENSITY_COL, "~ dist_scaled + <covariates> + (1|sample_id)\n")
cat("  dist_scaled = log(distance)/SD; score higher NEAR PHF1 => NEGATIVE slope.\n")
cat("  pct_attenuation = 100*(1 - adjusted/unadjusted): large + adjusted n.s. => dependent on intensity.\n")
cat("Visualisation removes the fitted intensity-linear contribution (b_hat*(intensity-mean)) and\n")
cat("  compares OBSERVED vs ADJUSTED rolling means over raw distance (window", WINDOW_UM, "um).\n\n")
cat("Per celltype x signature:\n"); print(as.data.frame(summ_df))

cat("\n--- Extension outputs (estimates only) ---\n")
cat("Summary columns added: conditional intensity coefficient (cond_intensity_coef, 95% CI, sign);\n")
cat("  semi-partial R2 (sp_r2_dist, sp_r2_intensity); commonality (commonality_unique_dist,\n")
cat("  commonality_unique_intensity, commonality_shared_unassignable = shared/UNASSIGNABLE);\n")
cat("  adjusted distance slope with intensity as rank / ns(df=", SPLINE_DF,
    ") / top-decile indicator; bottom within-donor intensity decile (dist_slope_decile1).\n", sep = "")
cat("Per-signature TSV/PDF: *_by_intensity_decile.* (slope vs within-donor decile);\n")
cat("  *_by_exclusion_radius.* (slope vs exclusion radius", paste(EXCLUSION_RADII, collapse = ","), "um).\n")
cat("\n(3) PHF1 intensity ~ distance in non-neuronal (tangle-incapable) celltypes",
    "[spillover bound]:\n")
if (!is.null(nn_df) && nrow(nn_df) > 0) print(as.data.frame(nn_df)) else cat("  (no non-neuronal fits)\n")

cat("\nNOTES:\n")
cat(" - COLLINEARITY: see dist_intensity_cor; strong |cor| inflates SE and destabilises the split.\n")
cat(" - MEASUREMENT ERROR: phf1_intensity_p95_z is a noisy tau proxy -> regressing it out under-adjusts.\n")
cat(" - Associational adjustment on PHF1-negative cells.\n")
cat("\nsessionInfo():\n"); print(sessionInfo())
sink()

cat("\nAll outputs written to:", OUTPUT_DIR, "\n"); cat("Done.\n")
