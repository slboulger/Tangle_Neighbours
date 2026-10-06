#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# phf1_distance_morphology_exc.R
#
# Figure panels: S7B - companion analysis run by run_morphology_observed.sh
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# phf1_distance_morphology_exc.R
#
# NON-CELL-AUTONOMOUS question: does the cell/nucleus morphology of PHF1-NEGATIVE
# Exc-IT-L2-3-CBLN2-HOPX neurons vary with distance to the nearest PHF1+ neuron?
#
# CosMx port of R/imc_phf1_distance_morphology_exc.R. The cell-autonomous counterpart
# (PHF1+ vs PHF1- themselves) is R/phf1_morphology_exc.R -> plots/phf1_morphology/.
# The one-panel figure combining both is R/phf1_morphology_dist_phf1_panel.R.
# Local neuron density (n_local) is carried as a sensitivity covariate (density-adj).
#
# MODELS (per measure), sign = per s.d. CLOSER to a tangle (positive = larger near tangles)
#   HEADLINE   : value ~ dist_scaled + Sex + Age_s + PMI_s + (1|sample_id)
#   depth-adj  : + depth_rel                                (reported, not headline)
#   density-adj: + n_local                                  (reported, not headline)
#   umi-adj    : + nUMI_log                                 (reported, not headline)
#   unadj      : value ~ dist_scaled + (1|sample_id)        (reported)
# No FOV-edge term, for consistency: the canonical CosMx distance models
# (deg_dream_linear_distance_phf1.r, plot_modulescore_vs_phf1_distance_modelp.R), the
# channel-intensity family and the IMC panel this mirrors all omit it; only nn3 carries it, as
# a spacing-specific technical term.
# The headline is not density-adjusted, matching the IMC analysis.
#
# DISTANCE TRANSFORM -- canonical: plain natural log, divided by SD only,
# NOT centred, with dist_sd computed WITHIN this celltype after the cap. log1p is
# forbidden: dist_to_phf1_um is strictly > 0 for every modelled (PHF1-negative) cell,
# and the script stop()s rather than absorbing a non-positive value with an offset.
#
# MEASURES: NucArea, NucAspectRatio, Circularity, Eccentricity, Perimeter.
#   'Solidity' is EXCLUDED -- the AtoMx column of that name equals Area/Perimeter, not a
#   convex-hull solidity. Do not re-add it. See R/phf1_morphology_exc.R for the details.
#   Nuclear measures additionally require a PLAUSIBLE nucleus; see R/phf1_morphology_filters.R.
#
# Distance cap is the single positional argument (1000 canonical, 300 to match the IMC family):
#   Rscript R/phf1_distance_morphology_exc.R 1000
#   Rscript R/phf1_distance_morphology_exc.R 300
#
# Triple-output convention, output under plots/phf1_distance_morphology_<cap>um/.
# Read-only on all inputs.

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(qs)
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
  library(lme4); library(lmerTest)
})

hpc <- "<PROJECT_ROOT>/phf1_v2"
loc <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(hpc)) hpc else loc)
source("R/palettes.R")                 # fig_theme; braak_palette["6"] == "#BD0026"
source("R/phf1_morphology_filters.R")  # PX_UM, PX2_UM2, NUC_MIN_UM2, nucleus_ok(), audit

CELLTYPE <- "Exc-IT-L2-3-CBLN2-HOPX"
.cli <- commandArgs(trailingOnly = TRUE)
MAX_DIST_UM <- if (length(.cli) >= 1) as.numeric(.cli[1]) else 1000
if (!is.finite(MAX_DIST_UM) || MAX_DIST_UM <= 0) stop("max_dist_um must be a positive number")

WINDOW_UM    <- 75            # rolling-window WIDTH, matches the module-score figures
N_GRID       <- 200
MIN_WINDOW_N <- 10            # do not draw a "mean" over fewer than 10 cells
RING_BREAKS  <- c(0, 50, 100, 200, 300, 500, 700, 1000)
FDR          <- 0.05
SIG_BASIS    <- "padj<0.05"

out_dir <- sprintf("plots/phf1_distance_morphology_%dum", MAX_DIST_UM)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE)

## ---------------------------------------------------------------------------
## 1. Load morphology + PHF1 distance
## ---------------------------------------------------------------------------
sce_path <- file.path("celltype_sce_neighbours", paste0(CELLTYPE, "_sce_neighbours.qs"))
if (!file.exists(sce_path)) stop("Not found: ", sce_path)
sce <- qread(sce_path)
cd  <- as.data.frame(colData(sce)); cd$cell_id <- colnames(sce)
scale_err <- assert_px_scale(cd)

nuc_audit <- nucleus_filter_audit(cd$NucArea, group = "all CBLN2")
nuc_ok    <- nucleus_ok(cd$NucArea)

base <- cd %>%
  transmute(
    cell_id, sample_id = as.character(sample_id),
    Braak = factor(as.character(Braak), levels = braak_levels),
    Sex = factor(as.character(Sex)), Age = as.numeric(Age), PMI = as.numeric(PMI),
    phf1_pos = as.logical(PHF1 %in% c(TRUE, "TRUE", "True", 1, "1")),
    dist_to_phf1_um = as.numeric(dist_to_phf1_um),
    nCount_RNA = as.numeric(nCount_RNA),
    nUMI_log   = log2(as.numeric(nCount_RNA) + 1),
    # NUCLEAR measures: NA unless the nucleus is present AND plausibly sized
    nucarea      = ifelse(nuc_ok, NucArea * PX2_UM2, NA_real_),
    nucaspect    = ifelse(nuc_ok, NucAspectRatio,    NA_real_),
    # CELL measures: defined for every cell, including those with no segmented nucleus
    circularity  = Circularity,
    eccentricity = Eccentricity,
    perimeter    = Perimeter * PX_UM
  )

## ---------------------------------------------------------------------------
## 2. Spatial covariates -- REUSED, not recomputed.
##    plots/nn3_neuron_spacing/source_data_nn3_cells.tsv (written by R/nn3_phf1_vs_neg.R)
##    already carries n_local / depth_rel for the whole neuron pool, computed by
##    local_neuron_count() and add_cortical_depth() in R/nn3_utils.r. Reusing it keeps this
##    analysis and the nn3 spacing analysis on byte-identical covariates.
## ---------------------------------------------------------------------------
cache <- "plots/nn3_neuron_spacing/source_data_nn3_cells.tsv"
if (!file.exists(cache))
  stop("Covariate cache not found: ", cache, "\nRun R/nn3_phf1_vs_neg.R first.")
cov_tab <- read.delim(cache, sep = "\t", header = TRUE, check.names = FALSE,
                      colClasses = c(cell_id = "character")) %>%
  select(cell_id, n_local, depth_rel)

miss <- sum(!base$cell_id %in% cov_tab$cell_id)
if (miss > 0)
  stop(sprintf("%d / %d %s cells are absent from %s -- the cache is stale.",
               miss, nrow(base), CELLTYPE, cache))
base <- left_join(base, cov_tab, by = "cell_id")
cat(sprintf("Covariate cache join: %d / %d cells matched, NAs: n_local %d | depth_rel %d\n",
            sum(base$cell_id %in% cov_tab$cell_id), nrow(base),
            sum(is.na(base$n_local)), sum(is.na(base$depth_rel))))

## ---------------------------------------------------------------------------
## 3. Modelling set: PHF1-NEGATIVE cells within the distance window.
##    dist_to_phf1_um is NA for PHF1+ cells by design (they are the source set, not
##    query cells -- see R/phf1_distance_utils.r), so !is.na() already drops them.
## ---------------------------------------------------------------------------
md <- base %>%
  filter(!phf1_pos, !is.na(dist_to_phf1_um), dist_to_phf1_um <= MAX_DIST_UM,
         !is.na(n_local), !is.na(depth_rel),
         !is.na(Sex), !is.na(Age), !is.na(PMI))

if (any(md$dist_to_phf1_um <= 0))
  stop(sprintf(paste0("%d modelled cell(s) have dist_to_phf1_um <= 0. Distance to the ",
                      "nearest PHF1+ neuron must be strictly > 0 for every PHF1-negative ",
                      "cell; log() is the canonical transform and log1p is not permitted ",
                      "here. Fix upstream in R/label_phf1_neighbours.r rather ",
                      "than adding an offset."), sum(md$dist_to_phf1_um <= 0)), call. = FALSE)

dt      <- log(md$dist_to_phf1_um)
dist_sd <- stats::sd(dt)
md$dist_scaled <- dt / dist_sd
md$Age_s <- as.numeric(scale(md$Age))
md$PMI_s <- as.numeric(scale(md$PMI))
md$Sex   <- droplevels(md$Sex)

n_don <- n_distinct(md$sample_id)
cat(sprintf("Modelled: %d PHF1- cells | %d donors | dist_sd = %.4f | cap %d um\n",
            nrow(md), n_don, dist_sd, MAX_DIST_UM))
if (n_don != 9) stop("Expected 9 donors, got ", n_don)

rho_dist_density <- suppressWarnings(
  cor(log(md$dist_to_phf1_um), md$n_local, method = "spearman"))

## ---------------------------------------------------------------------------
## 4. Features
## ---------------------------------------------------------------------------
features <- tibble::tribble(
  ~slug,          ~label,                 ~unit,      ~nuclear,
  "nucarea",      "Nucleus area",         "um^2",     TRUE,
  "nucaspect",    "Nucleus aspect ratio", "unitless", TRUE,
  "circularity",  "Cell circularity",     "unitless", FALSE,
  "eccentricity", "Cell eccentricity",    "unitless", FALSE,
  "perimeter",    "Cell perimeter",       "um",       FALSE
)
stopifnot(all(features$slug %in% colnames(md)))

delta_label <- function(slug) {
  switch(slug,
    nucarea      = expression(Delta * " nucleus area vs donor mean (" * mu * "m"^2 * ")"),
    nucaspect    = expression(Delta * " nucleus aspect ratio vs donor mean"),
    circularity  = expression(Delta * " cell circularity vs donor mean"),
    eccentricity = expression(Delta * " cell eccentricity vs donor mean"),
    perimeter    = expression(Delta * " cell perimeter vs donor mean (" * mu * "m)"),
    slug)
}

## ---------------------------------------------------------------------------
## 5. Helpers -- rolling window, verbatim from plot_modulescore_vs_phf1_distance_modelp.R
##    so the descriptive overlay is drawn the same way as every other distance figure.
##    Model-free; the SEM is over cells and the band is descriptive.
## ---------------------------------------------------------------------------
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
  bad <- r$n_window < MIN_WINDOW_N
  r$roll_mean[bad] <- NA_real_; r$sem[bad] <- NA_real_
  r
}
# Within-donor CENTRING (not z-scoring), so the axis keeps native units. Pooling across
# donors can invert the within-donor trend (Simpson's paradox) because donors differ in
# both baseline morphology and tangle burden; the model handles that with (1|sample_id),
# and the drawn curve must be on the matching scale or it will disagree with the fit.
centre_within_donor <- function(d, col) {
  d %>% group_by(sample_id) %>%
    mutate(value_c = .data[[col]] - mean(.data[[col]], na.rm = TRUE)) %>% ungroup()
}
sig_factor <- function(x) factor(ifelse(x, SIG_BASIS, "ns"), levels = c(SIG_BASIS, "ns"))

## ---------------------------------------------------------------------------
## 6. Fit
## ---------------------------------------------------------------------------
fit_one <- function(mc) {
  d <- md %>%
    transmute(value = .data[[mc]], dist_scaled, dist_to_phf1_um, depth_rel, n_local,
              nUMI_log, nCount_RNA, Sex, Age_s, PMI_s, sample_id)
  d <- d[is.finite(d$value), ]

  grab <- function(m, term = "dist_scaled") {
    if (is.null(m)) return(rep(NA_real_, 5))
    co <- summary(m)$coefficients
    if (!term %in% rownames(co)) return(rep(NA_real_, 5))
    c(co[term, "Estimate"], co[term, "Std. Error"], co[term, "df"],
      co[term, "t value"], co[term, "Pr(>|t|)"])
  }
  fitm <- function(f) tryCatch(lmerTest::lmer(f, data = d, REML = TRUE),
                               error = function(e) NULL)

  m_head  <- fitm(value ~ dist_scaled + Sex + Age_s + PMI_s + (1 | sample_id))   # HEADLINE
  m_depth <- fitm(value ~ dist_scaled + Sex + Age_s + PMI_s + depth_rel + (1 | sample_id))
  m_dens  <- fitm(value ~ dist_scaled + Sex + Age_s + PMI_s + n_local + (1 | sample_id))
  m_umi   <- fitm(value ~ dist_scaled + Sex + Age_s + PMI_s + nUMI_log + (1 | sample_id))
  m_unadj <- fitm(value ~ dist_scaled + (1 | sample_id))

  a <- grab(m_head); ad <- grab(m_depth); an <- grab(m_dens)
  am <- grab(m_umi); au <- grab(m_unadj)
  crit <- stats::qt(0.975, ifelse(is.na(a[3]), Inf, a[3]))
  sd_tot <- if (is.null(m_head)) NA_real_ else {
    vc <- as.data.frame(lme4::VarCorr(m_head)); sqrt(sum(vc$vcov[is.na(vc$var2)]))
  }
  # Pooled (donor-ignoring) slope of the DRAWN curve, as a diagnostic only.
  pl <- summary(stats::lm(value ~ dist_scaled, data = d))$coefficients["dist_scaled", ]

  tibble(
    measure = mc,
    # SIGN: every beta is negated so positive = LARGER CLOSER to a tangle
    beta_closer = -a[1], SE = a[2], df = a[3], t = -a[4], pval = a[5],
    CI.L = -a[1] - crit * a[2], CI.R = -a[1] + crit * a[2],
    cohens_d = -a[1] / sd_tot,
    beta_closer_depth_adj = -ad[1], pval_depth_adj = ad[5],
    beta_closer_dens_adj  = -an[1], pval_dens_adj  = an[5],
    beta_closer_umi_adj   = -am[1], pval_umi_adj   = am[5],
    beta_closer_unadj     = -au[1], pval_unadj     = au[5],
    pooled_beta_closer = -pl[1], pooled_p = pl[4],
    mean_value = mean(d$value), n_cells = nrow(d),
    n_donors = n_distinct(d$sample_id),
    rho_nCount    = suppressWarnings(cor(d$value, d$nCount_RNA, method = "spearman")),
    cor_density   = suppressWarnings(cor(d$value, d$n_local, method = "spearman")),
    cor_proximity = suppressWarnings(cor(d$value, -log(d$dist_to_phf1_um), method = "spearman")),
    singular = if (is.null(m_head)) NA else lme4::isSingular(m_head)
  )
}

fits <- lapply(features$slug, fit_one)
names(fits) <- features$slug
res <- bind_rows(lapply(fits, function(x) x)) %>%
  left_join(features, by = c("measure" = "slug")) %>%
  mutate(padj = p.adjust(pval, method = "BH"),
         significant = !is.na(padj) & padj < FDR) %>%
  arrange(pval)

# refit objects for the drawn line
models <- lapply(features$slug, function(mc) {
  d <- md %>% transmute(value = .data[[mc]], dist_scaled, Sex, Age_s, PMI_s, sample_id)
  d <- d[is.finite(d$value), ]
  list(m = tryCatch(lmerTest::lmer(value ~ dist_scaled + Sex + Age_s + PMI_s + (1 | sample_id),
                                   data = d, REML = TRUE), error = function(e) NULL), d = d)
})
names(models) <- features$slug

## ---------------------------------------------------------------------------
## 7. Drawn data
## ---------------------------------------------------------------------------
GRID <- seq(min(md$dist_to_phf1_um), MAX_DIST_UM, length.out = N_GRID)
half <- WINDOW_UM / 2

roll_df <- bind_rows(lapply(features$slug, function(mc) {
  d <- md %>% filter(is.finite(.data[[mc]]))
  d <- centre_within_donor(d, mc)
  r <- blank_sparse(roll_mean(d$dist_to_phf1_um, d$value_c, GRID, half))
  tibble(measure = mc,
         label = res$label[res$measure == mc], unit = res$unit[res$measure == mc],
         window_um = WINDOW_UM, significant = sig_factor(res$significant[res$measure == mc]),
         dist_um = GRID, roll_mean = r$roll_mean, sem = r$sem, n_window = r$n_window)
}))

fit_df <- bind_rows(lapply(features$slug, function(mc) {
  m <- models[[mc]]$m; d <- models[[mc]]$d
  if (is.null(m)) return(NULL)
  nd <- data.frame(dist_scaled = log(GRID) / dist_sd,
                   Sex = factor(levels(d$Sex)[1], levels = levels(d$Sex)),
                   Age_s = mean(d$Age_s), PMI_s = mean(d$PMI_s))
  X    <- model.matrix(~ dist_scaled + Sex + Age_s + PMI_s, nd)
  beta <- lme4::fixef(m); X <- X[, names(beta), drop = FALSE]
  # The drawn line is a DEVIATION from the grid mean, to match the within-donor-centred
  # rolling mean. Its uncertainty must therefore be that of the contrast: centre the design
  # matrix before propagating vcov, or the intercept and covariate-mean uncertainty (which
  # cancels out of the plotted quantity) inflates the band by an order of magnitude.
  Xc <- sweep(X, 2, colMeans(X), "-")
  ft <- as.numeric(Xc %*% beta)
  V  <- as.matrix(vcov(m)); se <- sqrt(rowSums((Xc %*% V) * Xc))
  tibble(measure = mc, label = res$label[res$measure == mc],
         dist_um = GRID, fitted = ft, ci_lo = ft - 1.96 * se, ci_hi = ft + 1.96 * se,
         significant = sig_factor(res$significant[res$measure == mc]))
}))

rings <- bind_rows(lapply(features$slug, function(mc) {
  d <- md %>% filter(is.finite(.data[[mc]]))
  d <- centre_within_donor(d, mc)
  d$ring <- cut(d$dist_to_phf1_um, breaks = RING_BREAKS, include.lowest = TRUE)
  d %>% group_by(sample_id, ring) %>%
    summarise(mean_value = mean(.data[[mc]]), mean_value_centred = mean(value_c),
              sem = sd(.data[[mc]]) / sqrt(n()), n_cells = n(), .groups = "drop") %>%
    mutate(measure = mc, label = res$label[res$measure == mc],
           unit = res$unit[res$measure == mc])
}))

wr(res,     "source_data_morphology_distance_coef.tsv")
wr(roll_df, "source_data_morphology_distance_rollmean.tsv")
wr(fit_df,  "source_data_morphology_distance_fit.tsv")
wr(rings,   "source_data_morphology_distance_rings.tsv")
wr(res %>% select(measure, label, unit, beta_closer, SE, df, CI.L, CI.R, cohens_d,
                  pval, padj, n_cells, n_donors,
                  beta_closer_unadj, beta_closer_depth_adj, beta_closer_dens_adj,
                  beta_closer_umi_adj, pval_umi_adj, pooled_beta_closer, pooled_p,
                  rho_nCount, cor_density, cor_proximity),
   "stats_morphology_distance_coef_effectsize.tsv")

## ---------------------------------------------------------------------------
## 8. Figures
## ---------------------------------------------------------------------------
# (a) coefficient forest. x = Cohen's d, NOT raw beta: the measures are in um^2, um and
#     unitless, so raw betas are not comparable on one axis. Raw beta + CI is in the tables.
fp <- res %>%
  mutate(label = factor(label, levels = rev(label[order(cohens_d)])),
         sig = ifelse(significant, "padj < 0.05", "n.s."),
         d_lo = cohens_d * CI.L / beta_closer, d_hi = cohens_d * CI.R / beta_closer)
p1 <- ggplot(fp, aes(x = cohens_d, y = label, colour = sig)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
  geom_errorbar(aes(xmin = d_lo, xmax = d_hi), orientation = "y", width = 0, linewidth = 0.4) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = c("padj < 0.05" = "#BD0026", "n.s." = "grey65"), name = NULL,
                      drop = FALSE, limits = c("padj < 0.05", "n.s.")) +
  labs(x = "Cohen's d per s.d. closer to a tangle\n(adjusted for sex, age, PMI)", y = NULL) +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "top",
        plot.margin = margin(4, 8, 4, 4, unit = "pt"))
ggsave(file.path(out_dir, "plot_morphology_distance_coef.pdf"), p1,
       width = 8.8, height = 4.6, units = "cm", device = "pdf")

# (b) rolling-mean panels, one per measure, with the model fit overlaid
XLAB <- expression("Distance to nearest PHF1+ neuron (" * mu * "m)")
for (mc in features$slug) {
  rd <- roll_df %>% filter(measure == mc, !is.na(roll_mean))
  fd <- fit_df  %>% filter(measure == mc)
  pr <- ggplot() +
    geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3, colour = "grey50") +
    geom_ribbon(data = rd, aes(dist_um, ymin = roll_mean - 1.96 * sem,
                               ymax = roll_mean + 1.96 * sem),
                fill = "grey60", alpha = 0.25, colour = NA) +
    geom_line(data = rd, aes(dist_um, roll_mean, colour = "Rolling mean (75 um)"),
              linewidth = 0.4) +
    geom_ribbon(data = fd, aes(dist_um, ymin = ci_lo, ymax = ci_hi),
                fill = "#BD0026", alpha = 0.15, colour = NA) +
    geom_line(data = fd, aes(dist_um, fitted, colour = "Model fit",
                             linetype = significant, linewidth = significant)) +
    scale_colour_manual(values = c("Rolling mean (75 um)" = "grey40", "Model fit" = "#BD0026"),
                        name = NULL, breaks = c("Rolling mean (75 um)", "Model fit")) +
    scale_linetype_manual(values = setNames(c("solid", "dashed"), c(SIG_BASIS, "ns")),
                          name = "Model FDR", drop = FALSE, limits = c(SIG_BASIS, "ns")) +
    scale_linewidth_manual(values = setNames(c(0.8, 0.4), c(SIG_BASIS, "ns")),
                           name = "Model FDR", drop = FALSE, limits = c(SIG_BASIS, "ns")) +
    guides(colour = guide_legend(order = 1),
           linetype = guide_legend(order = 2), linewidth = guide_legend(order = 2)) +
    labs(x = XLAB, y = delta_label(mc)) +
    coord_cartesian(xlim = c(0, MAX_DIST_UM)) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.text = element_text(size = 6),
          legend.key.width = grid::unit(18, "pt"), legend.key.height = grid::unit(8, "pt"),
          legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
  g <- ggplot2::ggplotGrob(pr)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  g$widths[pcol] <- grid::unit(1.6, "in")   # pin the panel so x-axes align across figures
  ggsave(file.path(out_dir, sprintf("plot_morphology_distance_rollmean_%s.pdf", mc)), g,
         width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
         height = 2.4, units = "in", device = "pdf")
}

## ---------------------------------------------------------------------------
## 9. Stats logs
## ---------------------------------------------------------------------------
sink(file.path(out_dir, "stats_morphology_distance_coef.txt"))
ttl <- sprintf("Morphology vs distance to nearest PHF1+ neuron -- PHF1- %s (CosMx, %d um)",
               CELLTYPE, MAX_DIST_UM)
cat(ttl, "\n"); cat(strrep("=", nchar(ttl)), "\n")
cat("Date:", format(Sys.time()), "\n\n")
cat("Source: ", sce_path, "\n", sep = "")
cat("Covariates n_local / depth_rel reused from ", cache, "\n", sep = "")
cat("  (written by R/nn3_phf1_vs_neg.R via R/nn3_utils.r), so this analysis and the nn3\n")
cat("  spacing analysis describe the same cells on identical covariates.\n")
cat("This is the NON-cell-autonomous question. The cell-autonomous one is\n")
cat("plots/phf1_morphology/; the combined one-panel figure is\n")
cat(sprintf("plots/phf1_morphology_dist_phf1_panel_%dum/.\n\n", MAX_DIST_UM))

cat("MODELS (per measure), sign = per s.d. CLOSER to a tangle (positive = larger near tangles)\n")
cat("  HEADLINE   : value ~ dist_scaled + Sex + Age_s + PMI_s + (1|sample_id)\n")
cat("  depth-adj  : + depth_rel                              (reported, NOT the headline)\n")
cat("  density-adj: + n_local                                (reported, NOT the headline)\n")
cat("  umi-adj    : + nUMI_log                               (reported, NOT the headline)\n")
cat("  unadj      : value ~ dist_scaled + (1|sample_id)      (reported)\n")
cat("  No FOV-edge term: the canonical CosMx distance models, the channel-intensity\n")
cat("  family and the IMC panel all omit it; only nn3 carries it.\n")
cat("  BH across the", nrow(features), "measures; significance padj <", FDR, "\n\n")
cat(sprintf("  dist_scaled = log(dist_to_phf1_um) / sd(log(dist_to_phf1_um)); dist_sd = %.5f,\n",
            dist_sd))
cat("  computed WITHIN this celltype after the cap. Plain log, divided by SD only, not\n")
cat("  centred. log1p is forbidden and a non-positive distance stop()s.\n")
cat(sprintf("  Pixel scale asserted on the object: max |Area.um2/Area - %.9f| = %.3g\n",
            PX2_UM2, scale_err))

cat_nucleus_filter_note(nuc_audit)

cat("\nCELLS\n")
cat(sprintf("  modelled PHF1- cells: %d | donors: %d | window: 0-%d um\n",
            nrow(md), n_don, MAX_DIST_UM))
cat("  n per measure differs where the nucleus filter applies -- see n_cells below.\n")
cat("  dist_to_phf1_um:\n"); print(summary(md$dist_to_phf1_um))
cat("  per donor:\n")
print(as.data.frame(md %>% count(sample_id, Braak, name = "n_cells")), row.names = FALSE)

cat("\nLOCAL NEURON DENSITY\n")
cat(sprintf("  Spearman rho(log(distance), neurons within 50 um) = %.3f\n", rho_dist_density))
cat("  Per-measure cor_density / cor_proximity are tabulated below, and\n")
cat("  beta_closer_dens_adj is given for direct comparison. The HEADLINE model is not\n")
cat("  density-adjusted.\n")

cat("\n=== RESULTS (headline model, sorted by p) ===\n")
print(as.data.frame(res %>% select(label, unit, beta_closer, CI.L, CI.R, cohens_d,
                                   pval, padj, n_cells, n_donors, mean_value)),
      row.names = FALSE, digits = 4)
cat("\n--- Sensitivity (same units as beta_closer) ---\n")
print(as.data.frame(res %>% select(measure, beta_closer, beta_closer_unadj,
                                   beta_closer_depth_adj, beta_closer_dens_adj,
                                   beta_closer_umi_adj, pval_umi_adj,
                                   cor_density, cor_proximity, singular)),
      row.names = FALSE, digits = 4)
cat("\n--- Pooled (donor-ignoring) slope of the drawn curve, for comparison ---\n")
print(as.data.frame(res %>% select(measure, beta_closer, pooled_beta_closer, pooled_p)),
      row.names = FALSE, digits = 4)
clash <- res %>% filter(is.finite(pooled_beta_closer), significant,
                        sign(beta_closer) != sign(pooled_beta_closer))
if (nrow(clash)) {
  cat("  *** SIGN CLASH: the drawn curve slopes against its own significance encoding for: ",
      paste(clash$measure, collapse = ", "), "\n", sep = "")
  cat("  The curve pools donors; the model is within donor. Read the curve as descriptive.\n")
}

cat("\nLIBRARY SIZE\n")
cat("  CosMx segmentation is transcript-informed, so library size (rho_nCount) and the\n")
cat("  nUMI-adjusted slope are tabulated:\n")
print(as.data.frame(res %>% select(measure, rho_nCount, beta_closer,
                                   beta_closer_umi_adj, pval_umi_adj)),
      row.names = FALSE, digits = 4)
cat("  REPORTED, not the headline: cell size plausibly causes both size and transcript\n")
cat("  count, so adjusting can remove real signal.\n")

cat("\n=== Effect sizes ===\n")
for (i in seq_len(nrow(res))) {
  cat(sprintf("  %-14s %+10.4g [%+10.4g, %+10.4g] %-9s per s.d. closer  (%+6.2f%% of mean)  d = %+.4f  padj = %.3g  n = %d\n",
              res$measure[i], res$beta_closer[i], res$CI.L[i], res$CI.R[i], res$unit[i],
              100 * res$beta_closer[i] / res$mean_value[i], res$cohens_d[i], res$padj[i],
              res$n_cells[i]))
}
cat("\n  Effect sizes are d and the 95% CI over ~", nrow(md), " cells from 9 donors.\n", sep = "")
cat("  Where the CI is tight around zero, read it as an exclusion bound, and compare\n")
cat("  against the cell-autonomous effect in plots/phf1_morphology/.\n")

cat("\nNOTES:\n")
cat(" - AtoMx 'Solidity' is Area/Perimeter and is excluded; Circularity is 4*pi*A/P^2, so\n")
cat("   Circularity and Perimeter between them re-express cell Area.\n")
cat(" - Local neuron density is a sensitivity covariate (rho above); headline unadjusted.\n")
cat(" - No FOV-edge covariate: it is absent from the canonical CosMx distance models, so it\n")
cat("   is absent here too.\n")
cat(" - Nucleus measures are conditioned on a plausible nucleus (see the filter note above).\n")
cat(" - No Braak stratification is attempted.\n")
cat(" - The drawn curve is within-donor CENTRED; the model fit is recentred to match.\n")
cat("\n=== sessionInfo() ===\n"); print(sessionInfo())
sink()

writeLines(c(
  sprintf("Rolling-mean overlay for %s, %d um window on a %d-point grid over 0-%d um.",
          CELLTYPE, WINDOW_UM, N_GRID, MAX_DIST_UM),
  "Values are within-donor centred (value - donor mean), native units.",
  sprintf("Windows with fewer than %d cells are blanked, not drawn.", MIN_WINDOW_N),
  "SEM is over cells: descriptive only, not an inferential band.",
  "The model fit drawn alongside is the headline LMM at covariate means, recentred to the",
  "same delta scale; its ribbon is a Wald CI on that CENTRED contrast (the design matrix is",
  "centred before propagating vcov), so it excludes the intercept and covariate-mean",
  "uncertainty that cancels out of the plotted quantity, and excludes the random intercept.",
  sprintf("Source data: source_data_morphology_distance_rollmean.tsv (%d rows), _fit.tsv (%d rows).",
          nrow(roll_df), nrow(fit_df))),
  file.path(out_dir, "stats_morphology_distance_rollmean.txt"))

writeLines(c(
  sprintf("Per-sample x ring means for %s.", CELLTYPE),
  sprintf("Ring breaks (um): %s", paste(RING_BREAKS, collapse = ", ")),
  "Both raw and within-donor-centred means are given; n_cells is the ring occupancy.",
  "Supplementary tabulation only -- no test is run on the rings.",
  sprintf("Source data: source_data_morphology_distance_rings.tsv (%d rows).", nrow(rings))),
  file.path(out_dir, "stats_morphology_distance_binned.txt"))

cat("\nDone. Wrote", length(list.files(out_dir)), "files to", out_dir, "\n")
print(as.data.frame(res %>% select(measure, beta_closer, CI.L, CI.R, cohens_d, padj,
                                   n_cells, cor_density)),
      row.names = FALSE, digits = 3)
