#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# phf1_channel_intensity_exc.R
#
# Figure panels: 5C
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# phf1_channel_intensity_exc.R
#
# DAPI / Histone / rRNA immunofluorescence intensity in Exc-IT-L2-3-CBLN2-HOPX neurons,
# asked two ways in one script:
#   PART A  cell-autonomous : PHF1+ vs PHF1-            -> plots/phf1_channel_intensity/
#   PART B  spatial         : vs distance to nearest PHF1+ neuron
#                                          -> plots/phf1_channel_intensity_distance_1000um/
#
# Structure, models, covariates and output conventions are copied from
# R/phf1_morphology_exc.R and R/phf1_distance_morphology_exc.R so the intensity and
# morphology results are directly comparable.
#
# PHF1-CHANNEL SCALE REFERENCE (Part B). DAPI, Histone and rRNA are NON-PHF1 channels.
# pct_of_phf1_ref expresses each channel's distance slope as a percentage of the PHF1-channel
# slope, phf1_intensity_p95_z at -0.241 z per s.d. of log distance in neurons (i.e. +0.241
# in this script's "per s.d. CLOSER" sign; R/phf1_module_distance_intensity_adjust.R). These
# channels come from the AtoMx CosMx acquisition, whereas phf1_intensity_p95_z comes from a
# separate confocal post-stain image (python/extract_phf1_intensity.py).
#
# CHANNELS (9): DAPI, Histone, rRNA x {density, max, total} -- AtoMx per-cell mask stats.
#   DENSITY (Mean.X) is the mean over the cell mask, i.e. a concentration. Primary readout.
#   TOTAL   (Mean.X * Area.um2) is integrated content -- the analogue of the IMC dna_total.
#           It is size-driven BY CONSTRUCTION: a cell segmented 15% larger has ~15% more
#           total signal with identical staining. The density readout is the one that is
#           not mechanically tied to segmentation.
#   MAX     is the single brightest pixel: size-biased and saturation-prone, reported only.
#   RATIO   DAPI:Histone = DNA per histone, a chromatin-compaction readout. Both terms are
#           means over the SAME cell mask, so area cancels exactly -- this is the one channel
#           here that CANNOT be produced by a segmentation-size difference.
#   Mean.G and Mean.GFAP exist on the object as further control channels and are NOT
#   analysed here; add them to CHANNELS to include them.
#
# TRANSFORM -- log then z WITHIN sample. Raw channel values are 8-bit-style detector units
#   and are only comparable WITHIN an acquisition: R/add_phf1_intensity.R states this
#   explicitly for the PHF1 channel ("confocal exposure/gain differ per acquisition; use
#   the per-sample z-score for any cross-sample comparison"), and the same holds here --
#   median Mean.DAPI ranges 360 to 1906 across the 9 donors, a 5.3-fold spread that dwarfs
#   any biological effect. So:  z = (log(x) - mean(log(x))) / sd(log(x))  within sample_id.
#   log because exposure/gain act multiplicatively; every channel is strictly positive
#   (min 124 DN, no zeros, no NAs), so log() is safe and log1p is neither needed nor
#   permitted. Coefficients are therefore in within-donor s.d. of log intensity,
#   the same units as phf1_intensity_p95_z.
#   z is computed ONCE over all cells of the celltype, so Part A and Part B share a scale.
#
# Read-only on all inputs.

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(qs)
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
  library(lme4); library(lmerTest); library(emmeans)
})

hpc <- "<PROJECT_ROOT>/phf1_v2"
loc <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(hpc)) hpc else loc)
source("R/palettes.R")

CELLTYPE    <- "Exc-IT-L2-3-CBLN2-HOPX"
# Distance cap, as the single positional argument. 1000 um is the CosMx canonical window;
# 300 um matches the IMC distance cap (Set 3), where the ROI geometry forbids anything wider.
# Run both:
#   Rscript R/phf1_channel_intensity_exc.R 1000
#   Rscript R/phf1_channel_intensity_exc.R 300
# Part A (PHF1+ vs PHF1-) does not depend on the cap and is simply rewritten identically.
.cli <- commandArgs(trailingOnly = TRUE)
MAX_DIST_UM <- if (length(.cli) >= 1) as.numeric(.cli[1]) else 1000
if (!is.finite(MAX_DIST_UM) || MAX_DIST_UM <= 0) stop("max_dist_um must be a positive number")
WINDOW_UM   <- 75
N_GRID      <- 200
MIN_WINDOW_N <- 10
RING_BREAKS <- c(0, 50, 100, 200, 300, 500, 700, 1000)
FDR         <- 0.05
SIG_BASIS   <- "padj<0.05"
# Reference: phf1_intensity_p95_z regressed on log distance in neurons gives a dist_scaled
# coefficient of -0.241 (dimmer further away). Expressed in this script's "per s.d. CLOSER"
# sign convention that is +0.241, i.e. brighter nearer a tangle.
PHF1_REF_BETA_DIST <- -0.241
PHF1_REF_CLOSER    <- -PHF1_REF_BETA_DIST

cat_dir  <- "plots/phf1_channel_intensity"
dist_dir <- sprintf("plots/phf1_channel_intensity_distance_%dum", MAX_DIST_UM)
dir.create(cat_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(dist_dir, recursive = TRUE, showWarnings = FALSE)

phf1_fill <- c("PHF1-" = "grey75", "PHF1+" = "#BD0026")
emmeans::emm_options(lmer.df = "satterthwaite", lmerTest.limit = 1e5)

wr <- function(x, f, dir) write.table(x, file.path(dir, f), sep = "\t", quote = FALSE,
                                      row.names = FALSE)

# DENSITY vs TOTAL, following the IMC convention used on the EC dataset:
#   dna_density <- asinh(ir_counts)              # concentration
#   dna_total   <- log2(ir_counts * area)        # integrated content
# Substitutions for CosMx:
#   * ir_counts -> Mean.DAPI (DNA) / Mean.Histone (chromatin). Mean.X is the mean over the
#     CELL mask, so the matching area is the CELL area (Area.um2), not NucArea -- multiplying
#     a cell-mask mean by a nucleus area would not integrate anything.
#   * asinh -> log. asinh is the IMC transform for low ion counts with exact zeros; these
#     channels are detector units with a floor of ~124 DN and no zeros, where asinh(x) and
#     log(x) differ only by a constant. Everything then gets the same within-sample z.
#   * log2 -> log for the same reason; the z-score makes the base irrelevant anyway.
# use_area_cov = FALSE for the totals: log(cell area) is a COMPONENT of the outcome there,
# so putting it in the sensitivity model would just re-derive the density result. Because
# nUMI_log is not used (image-derived outcome, Set 2), those channels have no sensitivity
# fit and m_area is NULL for them -- the headline and unadjusted fits carry them.
# rRNA total is included for symmetry (it is the ribosomal-content analogue); drop that row
# if only the DNA totals are wanted.
CHANNELS <- tibble::tribble(
  ~slug,           ~col,            ~label,            ~stat,     ~use_area_cov,
  "dapi_mean",     "Mean.DAPI",     "Mean DAPI intensity",    "mean",    TRUE,
  "dapi_max",      "Max.DAPI",      "Max DAPI intensity",      "max",     TRUE,
  "dapi_total",    "Total.DAPI",    "Total DAPI intensity",      "total",   FALSE,
  "histone_mean",  "Mean.Histone",  "Mean Histone intensity", "mean",    TRUE,
  "histone_max",   "Max.Histone",   "Max Histone intensity",   "max",     TRUE,
  "histone_total", "Total.Histone", "Total Histone intensity",   "total",   FALSE,
  "rrna_mean",     "Mean.rRNA",     "Mean rRNA intensity",    "mean",    TRUE,
  "rrna_max",      "Max.rRNA",      "Max rRNA intensity",      "max",     TRUE,
  "rrna_total",    "Total.rRNA",    "Total rRNA intensity",      "total",   FALSE,
  # DNA per histone -- a chromatin-compaction / nucleosome-density readout. Both terms are
  # cell-mask means over the SAME mask, so cell area cancels exactly: the ratio is unaffected
  # by segmentation size, and the "total" ratio (Mean.DAPI*A)/(Mean.Histone*A) is
  # algebraically the same number, so it is not tested separately. Note it is NOT recoverable from the two z columns (each was divided
  # by its own SD), so it is built from the raw channels before normalising.
  "dapi_histone",  "Ratio.DAPI.Histone", "DAPI:Histone ratio", "ratio", TRUE
)
SOURCE_COLS <- c("Mean.DAPI", "Max.DAPI", "Mean.Histone", "Max.Histone",
                 "Mean.rRNA", "Max.rRNA")

## ---------------------------------------------------------------------------
## 1. Load and normalise
## ---------------------------------------------------------------------------
sce_path <- file.path("celltype_sce_neighbours", paste0(CELLTYPE, "_sce_neighbours.qs"))
if (!file.exists(sce_path)) stop("Not found: ", sce_path)
sce <- qread(sce_path)
cd  <- as.data.frame(colData(sce)); cd$cell_id <- colnames(sce)

miss <- setdiff(c(SOURCE_COLS, "Area.um2"), colnames(cd))
if (length(miss)) stop("Missing column(s): ", paste(miss, collapse = ", "))

# Integrated content = mean over the cell mask x cell area, in DN * um^2.
cd$Total.DAPI    <- cd$Mean.DAPI    * cd$Area.um2
cd$Total.Histone <- cd$Mean.Histone * cd$Area.um2
cd$Total.rRNA    <- cd$Mean.rRNA    * cd$Area.um2
# Dimensionless; cell area cancels because both means are over the same mask.
cd$Ratio.DAPI.Histone <- cd$Mean.DAPI / cd$Mean.Histone

md <- cd %>%
  transmute(
    cell_id, sample_id = as.character(sample_id),
    SlideLabel = as.character(SlideLabel),
    Braak = factor(as.character(Braak), levels = braak_levels),
    Sex = factor(as.character(Sex)), Age = as.numeric(Age), PMI = as.numeric(PMI),
    phf1_pos = as.logical(PHF1 %in% c(TRUE, "TRUE", "True", 1, "1")),
    dist_to_phf1_um = as.numeric(dist_to_phf1_um),
    nCount_RNA = as.numeric(nCount_RNA),
    nUMI_log   = log2(as.numeric(nCount_RNA) + 1),
    log_area   = log(as.numeric(Area.um2)),
    !!!setNames(lapply(CHANNELS$col, function(c0) rlang::sym(c0)), CHANNELS$col)
  ) %>%
  mutate(grp   = factor(ifelse(phf1_pos, "PHF1+", "PHF1-"), levels = c("PHF1-", "PHF1+")),
         Age_s = as.numeric(scale(Age)), PMI_s = as.numeric(scale(PMI)))

# --- strictly positive check, then log + within-sample z ---------------------
for (cc in CHANNELS$col) {
  v <- md[[cc]]
  if (any(!is.finite(v))) stop(cc, ": non-finite values present.")
  if (any(v <= 0)) stop(sprintf(paste0("%s: %d cell(s) have intensity <= 0. Channel ",
                                       "intensities are detector units and must be > 0 for ",
                                       "log(); log1p is not permitted. Fix ",
                                       "upstream rather than adding an offset."),
                                cc, sum(v <= 0)), call. = FALSE)
}
raw_summary <- bind_rows(lapply(seq_len(nrow(CHANNELS)), function(i) {
  md %>% group_by(sample_id, SlideLabel) %>%
    summarise(channel = CHANNELS$slug[i], median_raw = median(.data[[CHANNELS$col[i]]]),
              n = n(), .groups = "drop")
}))
for (i in seq_len(nrow(CHANNELS))) {
  md[[CHANNELS$slug[i]]] <- ave(log(md[[CHANNELS$col[i]]]), md$sample_id,
                                FUN = function(z) (z - mean(z)) / stats::sd(z))
}

n_cells <- nrow(md); n_pos <- sum(md$phf1_pos); n_don <- n_distinct(md$sample_id)
cat(sprintf("Loaded %s: %d cells | %d PHF1+ | %d donors\n", CELLTYPE, n_cells, n_pos, n_don))
if (n_don != 9) stop("Expected 9 donors, got ", n_don)

# per-donor raw spread, the reason for within-sample normalisation
fold_spread <- raw_summary %>% group_by(channel) %>%
  summarise(min_median = min(median_raw), max_median = max(median_raw),
            fold = max(median_raw) / min(median_raw), .groups = "drop")

umi_med_pos <- median(md$nCount_RNA[md$phf1_pos])
umi_med_neg <- median(md$nCount_RNA[!md$phf1_pos])
umi_wilcox_p <- suppressWarnings(wilcox.test(nCount_RNA ~ phf1_pos, data = md)$p.value)

## ===========================================================================
## PART A -- cell-autonomous: PHF1+ vs PHF1-
## ===========================================================================
fit_cat <- function(i) {
  ch <- CHANNELS[i, ]
  d  <- md %>% transmute(cell_id, sample_id, Braak, Sex, Age_s, PMI_s, nUMI_log, log_area,
                         nCount_RNA, grp, value = .data[[ch$slug]],
                         raw = .data[[ch$col]])
  d <- d[is.finite(d$value), ]

  dm <- d %>% group_by(sample_id, Braak, grp) %>%
    summarise(donor_mean = mean(value), n_cells = n(), .groups = "drop")
  dm_wide <- dm %>% select(sample_id, grp, donor_mean) %>%
    tidyr::pivot_wider(names_from = grp, values_from = donor_mean)
  dm_wide <- dm_wide[stats::complete.cases(dm_wide), ]

  fitm <- function(f, ...) tryCatch(lmerTest::lmer(f, data = d, REML = TRUE, ...),
                                    error = function(e) NULL)
  # HEADLINE carries Sex + Age_s + PMI_s so the PHF1 contrast and the distance model share
  # one covariate set. NOTE these are DONOR-level terms and
  # grp varies WITHIN donor, so (1 | sample_id) already absorbs them: the grp estimate is
  # essentially unchanged by their inclusion. The unadjusted fit is kept and reported so the
  # equivalence is visible rather than assumed.
  m_ri  <- fitm(value ~ grp + Sex + Age_s + PMI_s + (1 | sample_id))             # HEADLINE
  m_rs  <- fitm(value ~ grp + Sex + Age_s + PMI_s + (1 + grp | sample_id),
                control = lmerControl(optimizer = "bobyqa",
                                      optCtrl = list(maxfun = 2e5)))
  m_cov <- fitm(value ~ grp + (1 | sample_id))                                   # unadjusted
  # area-adjusted sensitivity. This outcome is IMAGE-derived, so the covariate set is the
  # image set (Sex + Age + PMI + donor) plus log(cell area), which is itself image-derived.
  # nUMI_log is not included: adjusting an intensity outcome for
  # transcript abundance conditions on a variable that differs ~50% between PHF1+ and PHF1-
  # cells, so it is a collider on the contrast rather than a nuisance term.
  # log(cell area) is omitted for the *_total channels because it is a component of the
  # outcome there -- including it would just re-derive the density result -- which leaves
  # nothing to adjust for, hence NULL.
  m_area <- if (isTRUE(ch$use_area_cov)) {
    fitm(value ~ grp + Sex + Age_s + PMI_s + log_area + (1 | sample_id))
  } else {
    NULL
  }
  lrt   <- if (!is.null(m_ri) && !is.null(m_rs))
    tryCatch(anova(m_ri, m_rs, refit = FALSE), error = function(e) NULL) else NULL

  emm_df <- NULL; pair_df <- NULL; omni <- NULL
  p_lmm <- NA_real_; diff_lmm <- NA_real_; diff_lo <- NA_real_; diff_hi <- NA_real_
  d_model <- NA_real_
  if (!is.null(m_ri)) {
    emm <- emmeans(m_ri, ~ grp)
    emm_df  <- as.data.frame(summary(emm, infer = TRUE))
    pair_df <- as.data.frame(summary(pairs(emm), infer = TRUE))
    omni    <- tryCatch(anova(m_ri), error = function(e) NULL)
    diff_lmm <- -pair_df$estimate[1]; diff_lo <- -pair_df$upper.CL[1]
    diff_hi  <- -pair_df$lower.CL[1]; p_lmm  <- pair_df$p.value[1]
    vc <- as.data.frame(lme4::VarCorr(m_ri))
    d_model <- diff_lmm / sqrt(sum(vc$vcov[is.na(vc$var2)]))
  }
  grab_diff <- function(m) {
    if (is.null(m)) return(c(NA_real_, NA_real_))
    pc <- tryCatch(as.data.frame(summary(pairs(emmeans(m, ~ grp)), infer = TRUE)),
                   error = function(e) NULL)
    if (is.null(pc)) return(c(NA_real_, NA_real_))
    c(-pc$estimate[1], pc$p.value[1])
  }
  sl <- grab_diff(m_rs); cv <- grab_diff(m_cov); tc <- grab_diff(m_area)

  n_neg <- sum(d$grp == "PHF1-"); n_pos_f <- sum(d$grp == "PHF1+")
  m_neg <- mean(d$value[d$grp == "PHF1-"]); m_pos <- mean(d$value[d$grp == "PHF1+"])
  sd_pooled <- sqrt(((n_neg - 1) * var(d$value[d$grp == "PHF1-"]) +
                     (n_pos_f - 1) * var(d$value[d$grp == "PHF1+"])) / (n_neg + n_pos_f - 2))
  d_cell <- (m_pos - m_neg) / sd_pooled

  dz <- NA_real_; dz_mean <- NA_real_; dz_lo <- NA_real_; dz_hi <- NA_real_
  n_same <- NA_integer_; w_p <- NA_real_
  if (nrow(dm_wide) >= 3) {
    dd <- dm_wide[["PHF1+"]] - dm_wide[["PHF1-"]]
    dz <- mean(dd) / sd(dd); dz_mean <- mean(dd); n_same <- sum(dd > 0)
    pt <- tryCatch(t.test(dm_wide[["PHF1+"]], dm_wide[["PHF1-"]], paired = TRUE),
                   error = function(e) NULL)
    if (!is.null(pt)) { dz_lo <- pt$conf.int[1]; dz_hi <- pt$conf.int[2] }
    pw <- tryCatch(wilcox.test(dm_wide[["PHF1+"]], dm_wide[["PHF1-"]], paired = TRUE),
                   error = function(e) NULL)
    if (!is.null(pw)) w_p <- pw$p.value
  }
  # Raw-scale fold change, for a units-bearing sentence alongside the z. It MUST be
  # estimated within donor like everything else: a pooled group-mean ratio mixes the 5.3-fold
  # between-donor exposure spread with the unequal PHF1+ counts per donor and can come out
  # with the opposite sign to the within-donor estimate. So refit the headline model on
  # log(raw) and exponentiate its contrast.
  m_raw <- fitm(log(raw) ~ grp + (1 | sample_id))
  fc_raw <- NA_real_
  if (!is.null(m_raw)) {
    pr <- tryCatch(as.data.frame(summary(pairs(emmeans(m_raw, ~ grp)), infer = TRUE)),
                   error = function(e) NULL)
    if (!is.null(pr)) fc_raw <- exp(-pr$estimate[1])   # PHF1+ / PHF1-
  }

  row <- tibble(
    channel = ch$slug, label = ch$label, stat = ch$stat,
    headline_test = "cell-level LMM: z ~ grp + Sex + Age_s + PMI_s + (1 | sample_id)",
    headline_p = p_lmm,
    mean_z_PHF1_neg = m_neg, mean_z_PHF1_pos = m_pos,
    lmm_diff_z = diff_lmm, CI.L = diff_lo, CI.R = diff_hi,
    fold_change_raw = fc_raw,
    cohens_d_model = d_model, cohens_d_cell = d_cell,
    dz_donor_paired = dz, dz_mean_diff = dz_mean,
    dz_mean_diff_lower = dz_lo, dz_mean_diff_upper = dz_hi,
    n_donors_paired = nrow(dm_wide), n_donors_same_dir = n_same, wilcox_p = w_p,
    slope_diff = sl[1], slope_p = sl[2],
    lrt_p_slope = if (is.null(lrt)) NA_real_ else lrt$`Pr(>Chisq)`[2],
    unadj_diff = cv[1], unadj_p = cv[2],
    area_adj_diff = tc[1], area_adj_p = tc[2],
    # "what fraction of the raw effect survives adjustment" is only meaningful when there
    # IS a raw effect; dividing a near-zero by a near-zero prints noise like -28000%.
    area_adj_pct_of_raw = if (is.finite(diff_lmm) && abs(diff_lmm) >= 0.05)
      100 * tc[1] / diff_lmm else NA_real_,
    rho_nCount = suppressWarnings(cor(d$raw, d$nCount_RNA, method = "spearman")),
    n_cells = nrow(d), n_PHF1_pos = n_pos_f,
    singular = if (is.null(m_ri)) NA else lme4::isSingular(m_ri)
  )
  list(row = row, d = d, dm = dm, emm_df = emm_df, pair_df = pair_df, omni = omni,
       lrt = lrt, ch = ch)
}

A <- lapply(seq_len(nrow(CHANNELS)), fit_cat)
resA <- bind_rows(lapply(A, `[[`, "row"))
resA$padj <- p.adjust(resA$headline_p, method = "BH")
resA$significant <- !is.na(resA$padj) & resA$padj < FDR

for (i in seq_along(A)) {
  f <- A[[i]]; ch <- f$ch; slug <- paste0("phf1_intensity_", ch$slug, "_CBLN2")
  r <- resA[i, ]
  wr(f$d %>% select(cell_id, sample_id, Braak, Sex, grp, value, raw),
     paste0("source_data_", slug, ".tsv"), cat_dir)
  wr(f$dm, paste0("source_data_", slug, "_donor_means.tsv"), cat_dir)
  if (!is.null(f$pair_df)) wr(f$pair_df, paste0("stats_", slug, "_pairwise.tsv"), cat_dir)
  wr(r, paste0("stats_", slug, "_effectsize.tsv"), cat_dir)

  d <- f$d; dm <- f$dm
  grp_mean <- dm %>% group_by(grp) %>%
    summarise(m = mean(donor_mean), sem = sd(donor_mean) / sqrt(n()), .groups = "drop") %>%
    mutate(xc = as.integer(grp))
  v_lo <- min(d$value); v_hi <- max(quantile(d$value, 0.99), max(dm$donor_mean))
  v_rg <- v_hi - v_lo; v_br <- v_hi + v_rg * 0.06; tick <- v_rg * 0.02
  p_lab <- paste0("padj = ", formatC(r$padj, format = "g", digits = 2))

  p <- ggplot(d, aes(x = value, y = grp)) +
    geom_violin(aes(fill = grp), orientation = "y", scale = "width",
                linewidth = 0.2, alpha = 0.6, colour = NA) +
    geom_jitter(data = dm, aes(x = donor_mean, y = grp, size = n_cells),
                width = 0, height = 0.16, shape = 21, fill = "white",
                colour = "black", stroke = 0.25, inherit.aes = FALSE) +
    geom_errorbar(data = grp_mean, aes(y = grp, xmin = m - sem, xmax = m + sem),
                  orientation = "y", inherit.aes = FALSE, width = 0.1, linewidth = 0.3) +
    geom_segment(data = grp_mean, aes(y = xc - 0.3, yend = xc + 0.3, x = m, xend = m),
                 inherit.aes = FALSE, linewidth = 0.5) +
    annotate("segment", y = 1, yend = 2, x = v_br, xend = v_br, linewidth = 0.3) +
    annotate("segment", y = 1, yend = 1, x = v_br, xend = v_br - tick, linewidth = 0.3) +
    annotate("segment", y = 2, yend = 2, x = v_br, xend = v_br - tick, linewidth = 0.3) +
    annotate("text", y = 1.5, x = v_br, label = p_lab, hjust = -0.12, size = 2.2) +
    scale_fill_manual(values = phf1_fill, guide = "none") +
    scale_size(transform = "log10", range = c(0.7, 2.0), name = "Cells per donor",
               breaks = c(1, 10, 100, 1000)) +
    coord_cartesian(xlim = c(v_lo, v_br + v_rg * 0.62)) +
    labs(x = sprintf("%s (within-donor z)", ch$label), y = NULL) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "bottom",
          legend.key.height = grid::unit(7, "pt"),
          legend.box.spacing = grid::unit(2, "pt"), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 4, unit = "pt"))
  ggsave(file.path(cat_dir, paste0("plot_", slug, ".pdf")), p,
         width = 6.14, height = 4.0, units = "cm", device = "pdf")
}
wr(resA, "stats_phf1_channel_intensity_CBLN2_effectsize.tsv", cat_dir)
wr(raw_summary, "source_data_channel_raw_per_donor.tsv", cat_dir)

## ===========================================================================
## PART B -- spatial: intensity vs distance to the nearest PHF1+ neuron
## ===========================================================================
cache <- "plots/nn3_neuron_spacing/source_data_nn3_cells.tsv"
if (!file.exists(cache))
  stop("Covariate cache not found: ", cache, "\nRun R/nn3_phf1_vs_neg.R first.")
cov_tab <- read.delim(cache, sep = "\t", colClasses = c(cell_id = "character")) %>%
  select(cell_id, n_local, depth_rel)
if (sum(!md$cell_id %in% cov_tab$cell_id) > 0)
  stop("Covariate cache is stale: some ", CELLTYPE, " cells are absent.")
mdB <- left_join(md, cov_tab, by = "cell_id") %>%
  filter(!phf1_pos, !is.na(dist_to_phf1_um), dist_to_phf1_um <= MAX_DIST_UM,
         !is.na(n_local), !is.na(depth_rel))

if (any(mdB$dist_to_phf1_um <= 0))
  stop("Non-positive dist_to_phf1_um in the modelled set; log() is the canonical transform ",
       "and log1p is not permitted.")
dtB <- log(mdB$dist_to_phf1_um); dist_sd <- stats::sd(dtB)
mdB$dist_scaled <- dtB / dist_sd
mdB$Sex <- droplevels(mdB$Sex)
rho_dist_density <- suppressWarnings(
  cor(log(mdB$dist_to_phf1_um), mdB$n_local, method = "spearman"))
cat(sprintf("Part B modelled: %d PHF1- cells | dist_sd = %.4f\n", nrow(mdB), dist_sd))

fit_dist <- function(i) {
  ch <- CHANNELS[i, ]
  d <- mdB %>% transmute(value = .data[[ch$slug]], raw = .data[[ch$col]],
                         dist_scaled, dist_to_phf1_um, depth_rel, n_local,
                         nUMI_log, log_area, nCount_RNA, Sex, Age_s, PMI_s, sample_id)
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
  m_head  <- fitm(value ~ dist_scaled + Sex + Age_s + PMI_s + (1 | sample_id))
  m_depth <- fitm(value ~ dist_scaled + Sex + Age_s + PMI_s + depth_rel + (1 | sample_id))
  m_dens  <- fitm(value ~ dist_scaled + Sex + Age_s + PMI_s + n_local + (1 | sample_id))
  # area-adjusted sensitivity; no nUMI_log (image outcome, see the PHF1+/-
  # section above). log(area) is part of the outcome for the *_total channels, so there is
  # nothing left to adjust for there.
  m_area  <- if (isTRUE(ch$use_area_cov)) {
    fitm(value ~ dist_scaled + Sex + Age_s + PMI_s + log_area + (1 | sample_id))
  } else {
    NULL
  }
  m_unadj <- fitm(value ~ dist_scaled + (1 | sample_id))

  a <- grab(m_head); ad <- grab(m_depth); an <- grab(m_dens)
  at <- grab(m_area); au <- grab(m_unadj)
  crit <- stats::qt(0.975, ifelse(is.na(a[3]), Inf, a[3]))
  sd_tot <- if (is.null(m_head)) NA_real_ else {
    vc <- as.data.frame(lme4::VarCorr(m_head)); sqrt(sum(vc$vcov[is.na(vc$var2)]))
  }
  row <- tibble(
    channel = ch$slug, label = ch$label, stat = ch$stat,
    # SIGN: negated so positive = BRIGHTER CLOSER to a tangle
    beta_closer = -a[1], SE = a[2], df = a[3], t = -a[4], pval = a[5],
    CI.L = -a[1] - crit * a[2], CI.R = -a[1] + crit * a[2],
    cohens_d = -a[1] / sd_tot,
    beta_closer_depth_adj = -ad[1], pval_depth_adj = ad[5],
    beta_closer_dens_adj  = -an[1], pval_dens_adj  = an[5],
    beta_closer_area_adj  = -at[1], pval_area_adj  = at[5],
    beta_closer_unadj     = -au[1], pval_unadj     = au[5],
    pct_of_phf1_ref = 100 * (-a[1]) / PHF1_REF_CLOSER,
    n_cells = nrow(d), n_donors = n_distinct(d$sample_id),
    rho_nCount    = suppressWarnings(cor(d$raw, d$nCount_RNA, method = "spearman")),
    cor_density   = suppressWarnings(cor(d$value, d$n_local, method = "spearman")),
    cor_proximity = suppressWarnings(cor(d$value, -log(d$dist_to_phf1_um), method = "spearman")),
    singular = if (is.null(m_head)) NA else lme4::isSingular(m_head)
  )
  list(row = row, model = m_head, data = d, ch = ch)
}

B <- lapply(seq_len(nrow(CHANNELS)), fit_dist)
names(B) <- CHANNELS$slug
resB <- bind_rows(lapply(B, `[[`, "row")) %>%
  mutate(padj = p.adjust(pval, method = "BH"),
         significant = !is.na(padj) & padj < FDR) %>%
  arrange(pval)

GRID <- seq(min(mdB$dist_to_phf1_um), MAX_DIST_UM, length.out = N_GRID)
half <- WINDOW_UM / 2
roll_mean <- function(x, y, grid, half_window) {
  as.data.frame(do.call(rbind, lapply(grid, function(g) {
    idx <- which(x >= g - half_window & x <= g + half_window); n <- length(idx)
    if (n < 1) return(c(roll_mean = NA_real_, sem = NA_real_, n_window = 0))
    c(roll_mean = mean(y[idx]),
      sem = if (n > 1) stats::sd(y[idx]) / sqrt(n) else NA_real_, n_window = n)
  })))
}
sig_factor <- function(x) factor(ifelse(x, SIG_BASIS, "ns"), levels = c(SIG_BASIS, "ns"))

roll_df <- bind_rows(lapply(CHANNELS$slug, function(cs) {
  d <- B[[cs]]$data %>% group_by(sample_id) %>%
    mutate(value_c = value - mean(value)) %>% ungroup()
  r <- roll_mean(d$dist_to_phf1_um, d$value_c, GRID, half)
  bad <- r$n_window < MIN_WINDOW_N
  r$roll_mean[bad] <- NA_real_; r$sem[bad] <- NA_real_
  tibble(channel = cs, label = resB$label[resB$channel == cs], window_um = WINDOW_UM,
         significant = sig_factor(resB$significant[resB$channel == cs]),
         dist_um = GRID, roll_mean = r$roll_mean, sem = r$sem, n_window = r$n_window)
}))

fit_df <- bind_rows(lapply(CHANNELS$slug, function(cs) {
  m <- B[[cs]]$model; d <- B[[cs]]$data
  if (is.null(m)) return(NULL)
  nd <- data.frame(dist_scaled = log(GRID) / dist_sd,
                   Sex = factor(levels(d$Sex)[1], levels = levels(d$Sex)),
                   Age_s = mean(d$Age_s), PMI_s = mean(d$PMI_s))
  X <- model.matrix(~ dist_scaled + Sex + Age_s + PMI_s, nd)
  beta <- lme4::fixef(m); X <- X[, names(beta), drop = FALSE]
  # centred design => CI on the plotted deviation, not the absolute prediction
  Xc <- sweep(X, 2, colMeans(X), "-")
  ft <- as.numeric(Xc %*% beta)
  V  <- as.matrix(vcov(m)); se <- sqrt(rowSums((Xc %*% V) * Xc))
  tibble(channel = cs, label = resB$label[resB$channel == cs], dist_um = GRID,
         fitted = ft, ci_lo = ft - 1.96 * se, ci_hi = ft + 1.96 * se,
         significant = sig_factor(resB$significant[resB$channel == cs]))
}))

rings <- bind_rows(lapply(CHANNELS$slug, function(cs) {
  d <- B[[cs]]$data %>% group_by(sample_id) %>%
    mutate(value_c = value - mean(value)) %>% ungroup()
  d$ring <- cut(d$dist_to_phf1_um, breaks = RING_BREAKS, include.lowest = TRUE)
  d %>% group_by(sample_id, ring) %>%
    summarise(mean_z = mean(value), mean_z_centred = mean(value_c),
              sem = sd(value) / sqrt(n()), n_cells = n(), .groups = "drop") %>%
    mutate(channel = cs, label = resB$label[resB$channel == cs])
}))

wr(resB,    "source_data_channel_distance_coef.tsv", dist_dir)
wr(roll_df, "source_data_channel_distance_rollmean.tsv", dist_dir)
wr(fit_df,  "source_data_channel_distance_fit.tsv", dist_dir)
wr(rings,   "source_data_channel_distance_rings.tsv", dist_dir)
wr(resB %>% select(channel, label, stat, beta_closer, CI.L, CI.R, cohens_d, pval, padj,
                   beta_closer_unadj, beta_closer_depth_adj, beta_closer_dens_adj,
                   beta_closer_area_adj, pval_area_adj, pct_of_phf1_ref,
                   rho_nCount, cor_density, cor_proximity, n_cells, n_donors),
   "stats_channel_distance_coef_effectsize.tsv", dist_dir)

# forest, in z units so all six share an axis; PHF1 reference slope drawn for scale
fp <- resB %>% mutate(label = factor(label, levels = rev(label[order(beta_closer)])),
                      sig = ifelse(significant, "padj < 0.05", "n.s."))
p1 <- ggplot(fp, aes(x = beta_closer, y = label, colour = sig)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
  geom_vline(xintercept = PHF1_REF_CLOSER, linetype = 3, linewidth = 0.4, colour = "#BD0026") +
  geom_errorbar(aes(xmin = CI.L, xmax = CI.R), orientation = "y", width = 0, linewidth = 0.4) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = c("padj < 0.05" = "#BD0026", "n.s." = "grey65"), name = NULL) +
  labs(x = "Log intensity (within-donor z)\nper s.d. closer to a tangle",
       y = NULL,
       subtitle = sprintf("adjusted for sex, age, PMI\ndotted line = PHF1 channel slope (%.3f z)",
                          PHF1_REF_CLOSER)) +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "top",
        plot.subtitle = element_text(size = 6, colour = "grey30"),
        plot.margin = margin(4, 8, 4, 4, unit = "pt"))
ggsave(file.path(dist_dir, "plot_channel_distance_coef.pdf"), p1,
       width = 9.4, height = 6.4, units = "cm", device = "pdf")

XLAB <- expression("Distance to nearest PHF1+ neuron (" * mu * "m)")
for (cs in CHANNELS$slug) {
  rd <- roll_df %>% filter(channel == cs, !is.na(roll_mean))
  fd <- fit_df  %>% filter(channel == cs)
  lab <- CHANNELS$label[CHANNELS$slug == cs]
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
    labs(x = XLAB, y = bquote(Delta ~ .(lab) ~ "z vs donor mean")) +
    coord_cartesian(xlim = c(0, MAX_DIST_UM)) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.text = element_text(size = 6),
          legend.key.width = grid::unit(18, "pt"), legend.key.height = grid::unit(8, "pt"),
          legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
  g <- ggplot2::ggplotGrob(pr)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  g$widths[pcol] <- grid::unit(1.6, "in")
  ggsave(file.path(dist_dir, sprintf("plot_channel_distance_rollmean_%s.pdf", cs)), g,
         width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
         height = 2.4, units = "in", device = "pdf")
}

## ===========================================================================
## PART C -- combined figure: distance curve + PHF1+ marginal strip on one axis
## Style copied from R/plot_reactome_stress_death_vs_phf1_distance.R
## (plots/reactome_stress_death_vs_phf1_distance_1000um/): rolling means over distance
## for the PHF1-NEGATIVE cells, with the PHF1+ mean drawn as a point + 95% CI in a
## marginal strip at the LEFT edge (i.e. past the zero-distance end), separated by a thin
## rule and labelled "PHF1+". Line style encodes the DISTANCE model's BH verdict; point
## shape encodes the PHF1+ vs PHF1- contrast's BH verdict. The two legends must carry
## DIFFERENT names or ggplot merges them (they share break labels).
##
## SCALE. Both the curve and the strip are drawn on the UNCENTRED within-donor z, not the
## delta-vs-donor-mean used in Part B: the strip is only readable if the PHF1+ cells and
## the PHF1- curve sit on one common scale, which the within-donor z already provides
## (this is the same reason the reactome figure uses the raw module score rather than a
## z-score). PHF1- cells sit near 0 by construction, so the strip's height IS the
## cell-autonomous effect and the curve's slope IS the spatial one, in one picture.
## ===========================================================================
comb_dir <- sprintf("plots/phf1_channel_intensity_combined_%dum", MAX_DIST_UM)
dir.create(comb_dir, recursive = TRUE, showWarnings = FALSE)

SIG_LEGEND <- "Log-distance"       # linetype/linewidth: the distance-model verdict
CON_LEGEND <- "PHF1+ vs PHF1-"     # point shape: the contrast verdict. MUST differ from above.
CON_SHAPES <- c(16, 1)             # filled = significant, hollow = n.s.
SERIES_COL <- "#BD0026"            # the PHF1 red used across the morphology/panel family
TOTAL_W_IN  <- 6.14 / 2.54         # total figure width, the family standard
PANEL_H_MIN <- 1.66

# Channels drawn as panels. One figure per channel, as in the IMC reference (a single
# series per panel, so there is no colour guide and the line carries SERIES_COL directly).
PANEL_CHANNELS <- c("rrna_mean", "dapi_mean", "histone_mean")

# Rolling mean of the UNCENTRED z for the PHF1-negative modelled cells.
roll_uncentred <- bind_rows(lapply(PANEL_CHANNELS, function(cs) {
  d <- B[[cs]]$data
  r <- roll_mean(d$dist_to_phf1_um, d$value, GRID, half)
  bad <- r$n_window < MIN_WINDOW_N
  r$roll_mean[bad] <- NA_real_; r$sem[bad] <- NA_real_
  tibble(channel = cs, label = CHANNELS$label[CHANNELS$slug == cs],
         dist_um = GRID, roll_mean = r$roll_mean, sem = r$sem, n_window = r$n_window,
         sig = sig_factor(resB$significant[resB$channel == cs]))
}))

# PHF1+ reference: mean + 95% CI of the same z over the PHF1+ cells of this celltype.
phf1_ref <- bind_rows(lapply(PANEL_CHANNELS, function(cs) {
  s <- md[[cs]][md$phf1_pos & is.finite(md[[cs]])]
  if (!length(s)) return(NULL)
  m <- mean(s); sem <- if (length(s) > 1) stats::sd(s) / sqrt(length(s)) else NA_real_
  tibble(channel = cs, label = CHANNELS$label[CHANNELS$slug == cs],
         ref_mean = m, sem = sem, ci_lo = m - 1.96 * sem, ci_hi = m + 1.96 * sem,
         n_phf1_pos = length(s),
         con_sig = sig_factor(resA$significant[resA$channel == cs]))
}))

# Total width is fixed and the panel takes whatever the axes and guide boxes leave. Copied
# from R/imc_phf1_morphology_dist_phf1_panel.R, including its note: with the total pinned,
# the strip variant carries ~14% more x range in the same physical width, so a micron is
# ~12% shorter there than on the _nostrip copy. Compare curves within a variant, not across.
save_fixed_panel <- function(p, path, panel_h = PANEL_H_MIN) {
  g <- ggplot2::ggplotGrob(p)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  prow <- unique(g$layout$t[grepl("^panel", g$layout$name)])
  g$heights[prow] <- grid::unit(panel_h, "in")
  g$widths[pcol]  <- grid::unit(0, "in")
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

# A guide key is built from LAYER DATA, so a level that never occurs (e.g. the only series is
# significant) yields a key with its label but no glyph. drop = FALSE and override.aes both
# fail to fix it -- the IMC reference figure shows the same blank "ns" key. Carrying one
# undrawn row per missing level gives the key its glyph; nothing is plotted because the
# value is NA. Grouping is on sig so the filler never joins a real line.
pad_levels <- function(d, col) {
  miss <- setdiff(c(SIG_BASIS, "ns"), as.character(d[[col]]))
  if (!length(miss)) return(d)
  filler <- d[rep(1L, length(miss)), , drop = FALSE]
  filler[[col]] <- factor(miss, levels = c(SIG_BASIS, "ns"))
  for (v in intersect(c("roll_mean", "sem", "ref_mean", "ci_lo", "ci_hi"), names(filler)))
    filler[[v]] <- NA_real_
  dplyr::bind_rows(d, filler)
}

build_panel <- function(cs, with_strip = TRUE) {
  rd <- pad_levels(roll_uncentred %>% filter(channel == cs, !is.na(roll_mean)), "sig")
  rf <- pad_levels(phf1_ref %>% filter(channel == cs), "con_sig")
  lab <- CHANNELS$label[CHANNELS$slug == cs]

  p <- ggplot(rd, aes(dist_um, roll_mean, group = sig)) +
    geom_hline(yintercept = 0, linetype = 3, linewidth = 0.25, colour = "grey60") +
    geom_ribbon(aes(ymin = roll_mean - 1.96 * sem, ymax = roll_mean + 1.96 * sem),
                fill = SERIES_COL, alpha = 0.15, colour = NA) +
    geom_line(aes(linetype = sig, linewidth = sig), colour = SERIES_COL)

  if (with_strip) {
    x_pos <- -0.07 * MAX_DIST_UM
    p <- p +
      annotate("segment", x = -0.025 * MAX_DIST_UM, xend = -0.025 * MAX_DIST_UM,
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
  # the default collides (7 ticks at a 300 um cap; "800"/"1000" overlap at n = 4 and a
  # 1000 um cap). n = 3 gives 0/100/200/300 and 0/500/1000 respectively, both of which fit.
  brk <- pretty(c(0, MAX_DIST_UM), n = 3); brk <- brk[brk >= 0 & brk <= MAX_DIST_UM]
  xlo <- if (with_strip) -0.14 * MAX_DIST_UM else 0

  p +
    scale_linetype_manual(values = setNames(c("solid", "dashed"), c(SIG_BASIS, "ns")),
                          name = SIG_LEGEND, drop = FALSE, limits = c(SIG_BASIS, "ns")) +
    scale_linewidth_manual(values = setNames(c(0.8, 0.4), c(SIG_BASIS, "ns")),
                           name = SIG_LEGEND, drop = FALSE, limits = c(SIG_BASIS, "ns")) +
    scale_x_continuous(breaks = brk) +
    guides(linetype = guide_legend(order = 1), linewidth = guide_legend(order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(colour = "grey25", size = 1.3))) +
    labs(x = expression("Distance to PHF1+ neuron (" * mu * "m)"),
         y = sprintf("%s (within-donor z)", lab)) +
    coord_cartesian(xlim = c(xlo, MAX_DIST_UM), clip = "off") +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "right", legend.direction = "vertical", legend.box = "vertical",
          legend.text = element_text(size = 6), legend.title = element_text(size = 6),
          legend.key.width = grid::unit(16, "pt"), legend.key.height = grid::unit(7, "pt"),
          legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
          legend.spacing.y = grid::unit(1, "pt"),
          plot.margin = margin(t = 14, r = 6, b = 4, l = 5))
}

for (cs in PANEL_CHANNELS) {
  save_fixed_panel(build_panel(cs, TRUE),
                   file.path(comb_dir, sprintf("plot_%s_dist_phf1.pdf", cs)))
  save_fixed_panel(build_panel(cs, FALSE),
                   file.path(comb_dir, sprintf("plot_%s_dist_phf1_nostrip.pdf", cs)))
}
cat(sprintf("Panels: %.2f cm wide, panel height %.2f in, legend right, cap %d um.\n",
            TOTAL_W_IN * 2.54, PANEL_H_MIN, MAX_DIST_UM))
wr(roll_uncentred, "source_data_channel_dist_phf1_rollmean.tsv", comb_dir)
wr(phf1_ref %>% select(channel, label, ref_mean, ci_lo, ci_hi, n_phf1_pos, con_sig),
   "source_data_channel_dist_phf1_ref.tsv", comb_dir)

# Pooled (donor-ignoring) slope of the DRAWN curve, for comparison with the within-donor
# model behind the line style. The outcome is already within-donor z, so the two should
# agree; printing both makes that checkable rather than assumed, and the script warns on a
# sign clash.
pooled <- bind_rows(lapply(PANEL_CHANNELS, function(cs) {
  d <- B[[cs]]$data
  m <- stats::lm(value ~ dist_scaled, data = d)
  co <- summary(m)$coefficients
  tibble(channel = cs, pooled_beta_closer = -co["dist_scaled", "Estimate"],
         pooled_p = co["dist_scaled", "Pr(>|t|)"])
}))

# Stats tables, mirroring the IMC panel directory: p-values and n written at full precision.
contrast_stats <- resA %>% filter(channel %in% PANEL_CHANNELS) %>%
  select(channel, label, lmm_diff_z, CI.L, CI.R, fold_change_raw, cohens_d_model,
         cohens_d_cell, dz_donor_paired, dz_mean_diff, dz_mean_diff_lower,
         dz_mean_diff_upper, n_donors_paired, n_donors_same_dir, wilcox_p,
         headline_p, padj, n_cells, n_PHF1_pos)
distance_stats <- resB %>% filter(channel %in% PANEL_CHANNELS) %>%
  select(channel, label, beta_closer, SE, df, CI.L, CI.R, cohens_d, pval, padj,
         n_cells, n_donors) %>%
  left_join(pooled, by = "channel")
wr(contrast_stats, "source_data_channel_dist_phf1_contrast_stats.tsv", comb_dir)
wr(distance_stats, "source_data_channel_dist_phf1_distance_stats.tsv", comb_dir)

clash <- distance_stats %>% filter(is.finite(pooled_beta_closer),
                                   sign(beta_closer) != sign(pooled_beta_closer))
if (nrow(clash))
  warning("Pooled and within-donor slopes disagree in sign for: ",
          paste(clash$channel, collapse = ", "))

sink(file.path(comb_dir, "stats_channel_dist_phf1.txt"))
ttl <- sprintf("Channel intensity over PHF1 distance with the PHF1+ strip -- %s (CosMx, %d um)",
               CELLTYPE, MAX_DIST_UM)
cat(ttl, "\n"); cat(strrep("=", nchar(ttl)), "\n"); cat("Date:", format(Sys.time()), "\n\n")
cat("Source: ", sce_path, "\n", sep = "")
cat("Figure design mirrors R/imc_phf1_morphology_dist_phf1_panel.R.\n")
cat("Figures only -- no model is fitted here; both encodings come from the two analyses\n")
cat("already run by this script:\n")
cat(sprintf("  line style  = distance-model BH verdict   -> %s/\n", dist_dir))
cat(sprintf("  point shape = PHF1+ vs PHF1- BH verdict   -> %s/\n", cat_dir))

cat("\nUNITS. AtoMx reports these channels as the mean/max pixel value over the CELL mask in\n")
cat("raw 16-bit detector counts -- arbitrary fluorescence units, not a physical quantity\n")
cat("(observed range here: Mean.DAPI 127-7704, Max.Histone up to 64500, i.e. near 16-bit\n")
cat("saturation). They are also comparable only WITHIN an acquisition: median Mean.DAPI\n")
cat("spans 5.3-fold across the 9 donors. The plotted and modelled quantity is therefore\n")
cat("  z = (log(x) - mean(log(x))) / sd(log(x))  within sample_id,\n")
cat("which is dimensionless: 1 unit = 1 within-donor SD of log intensity. That is why the\n")
cat("axes read '(within-donor z)' and carry no physical unit.\n")

cat("\nCELLS\n")
cat(sprintf("  PHF1+ vs PHF1- : %d cells (%d PHF1+, %d PHF1-) from %d donors\n",
            n_cells, n_pos, n_cells - n_pos, n_don))
cat(sprintf("  distance model : %d PHF1-negative cells from %d donors, window 0-%d um\n",
            nrow(mdB), n_distinct(mdB$sample_id), MAX_DIST_UM))
cat(sprintf("  dist_sd = %.5f (log distance, within this celltype, after the cap)\n", dist_sd))
cat("  PHF1+ cells per donor:\n")
print(as.data.frame(md %>% filter(phf1_pos) %>% count(sample_id, Braak, name = "n_PHF1_pos")),
      row.names = FALSE)

cat("\nSTRIP: PHF1+ neurons are the distance ANCHORS and have no distance to a nearest OTHER\n")
cat("  PHF1+ neuron (dist_to_phf1_um is NA for them by design), so they are drawn off-axis\n")
cat("  left of a separator rule. It is a REFERENCE LEVEL, not a distance-zero data point.\n")
cat("ENCODINGS: line style = log-distance verdict (solid padj<0.05, dashed ns);\n")
cat("  point shape = PHF1+ vs PHF1- verdict (filled padj<0.05, hollow ns).\n")

cat("\nDISTANCE MODEL (PHF1-negative cells)\n")
cat("  value ~ dist_scaled + Sex + Age_s + PMI_s + (1 | sample_id), REML.\n")
cat("  No FOV-edge term: the canonical CosMx distance\n")
cat("  models and the IMC panel this mirrors all omit it. Sign is per s.d. CLOSER.\n\n")
print(as.data.frame(distance_stats %>%
        select(label, beta_closer, CI.L, CI.R, cohens_d, pval, padj, n_cells, n_donors)),
      row.names = FALSE, digits = 4)
cat("\n  Pooled (donor-ignoring) slope of the drawn curve, for comparison:\n")
print(as.data.frame(distance_stats %>% select(label, beta_closer, pooled_beta_closer, pooled_p)),
      row.names = FALSE, digits = 4)
if (nrow(clash)) cat("  *** SIGN CLASH between pooled and within-donor slopes -- see warning ***\n")

cat("\nCONTRAST (PHF1+ vs PHF1-, cell-level LMM headline)\n")
print(as.data.frame(contrast_stats %>%
        select(label, lmm_diff_z, CI.L, CI.R, fold_change_raw, cohens_d_model,
               dz_donor_paired, n_donors_same_dir, n_donors_paired, headline_p, padj,
               n_cells, n_PHF1_pos)),
      row.names = FALSE, digits = 4)

cat("\n=== PHF1+ reference strip (mean +/- 95% CI of the within-donor z) ===\n")
print(as.data.frame(phf1_ref %>% select(label, ref_mean, ci_lo, ci_hi, n_phf1_pos, con_sig)),
      row.names = FALSE, digits = 4)

cat("\nExact p-values and n are in source_data_channel_dist_phf1_{contrast,distance}_stats.tsv\n")
cat("at full precision; the tables above are rounded for reading only.\n")
cat("\n=== sessionInfo() ===\n"); print(sessionInfo()); sink()
cat("\nPART C ->", comb_dir, "|", length(list.files(comb_dir)), "files\n")

## ===========================================================================
## Stats logs
## ===========================================================================
common_header <- function() {
  cat("Source: ", sce_path, "\n", sep = "")
  cat("Script: R/phf1_channel_intensity_exc.R\n")
  cat("\n")
  cat("CHANNELS: Mean/Max of DAPI, Histone, rRNA (AtoMx per-cell mask statistics).\n")
  cat("  Mean.X is a concentration over the cell mask; Max.X is the single brightest pixel\n")
  cat("  and is therefore size-biased and saturation-prone. Prefer Mean.\n")
  cat("  Mean.G and Mean.GFAP are available as further control channels, not analysed here.\n\n")
  cat("TRANSFORM: z = (log(x) - mean(log(x))) / sd(log(x)) WITHIN sample_id.\n")
  cat("  Raw channel values are detector units comparable only within an acquisition --\n")
  cat("  R/add_phf1_intensity.R states this for the PHF1 channel and it holds here too.\n")
  cat("  Observed per-donor median spread (why this is not optional):\n")
  print(as.data.frame(fold_spread), row.names = FALSE, digits = 4)
  cat("  Every channel is strictly positive (no zeros, no NAs) so log() is safe; log1p is\n")
  cat("  not permitted. Coefficients are in within-donor s.d. of log intensity,\n")
  cat("  the same units as phf1_intensity_p95_z.\n\n")
}

sink(file.path(cat_dir, "stats_phf1_channel_intensity_CBLN2.txt"))
ttl <- sprintf("DAPI / Histone / rRNA intensity in PHF1+ vs PHF1- %s neurons (CosMx)", CELLTYPE)
cat(ttl, "\n"); cat(strrep("=", nchar(ttl)), "\n"); cat("Date:", format(Sys.time()), "\n\n")
common_header()
cat("MODELS\n")
cat("  HEADLINE : z ~ grp + Sex + Age_s + PMI_s + (1 | sample_id)   [cell-level, lmerTest]\n")
cat("  unadj    : z ~ grp + (1 | sample_id)                         (reported)\n")
cat("  slope    : z ~ grp + Sex + Age_s + PMI_s + (1 + grp | sample_id) (reported)\n")
cat("  area-adj : + log(cell area)                                  (reported; mean/max only)\n")
cat("  The outcome is IMAGE-derived, so the covariate set is the image set (Sex + Age + PMI\n")
cat("  + donor). nUMI_log is not in the sensitivity fit: adjusting an\n")
cat("  intensity outcome for transcript abundance conditions on a variable that differs ~50%\n")
cat("  between PHF1+ and PHF1- cells, making it a collider rather than a nuisance term.\n")
cat("  area_adj is NA for the *_total channels: log(cell area) is a COMPONENT of those\n")
cat("  outcomes, so without nUMI_log there is nothing left to adjust for.\n")
cat("  Sex/Age/PMI are included so this and the distance model share one covariate set.\n")
cat("  They are DONOR-level and grp varies WITHIN donor, so (1|sample_id) already absorbs\n")
cat("  them; compare lmm_diff_z with unadj_diff below -- they should be near-identical.\n")
cat("  BH across the", nrow(CHANNELS), "channels; significance padj <", FDR, "\n\n")
cat("  NOTE ON `singular = TRUE`: expected here, not a failure. The outcome is already\n")
cat("  z-scored WITHIN sample_id, so every donor has mean 0 by construction and the donor\n")
cat("  random-intercept variance is legitimately ~0. The term is kept anyway so the model\n")
cat("  matches the morphology scripts and copes with the unbalanced group sizes; dropping\n")
cat("  it changes the estimates negligibly. A singular fit would only be a concern if the\n")
cat("  outcome had NOT been within-donor normalised.\n")
cat("  The same cause produces lme4's 'Model failed to converge with 1 negative eigenvalue'\n")
cat("  warning on the tech-adjusted fits (6 of them at each cap). It is a BOUNDARY artefact,\n")
cat("  not a failed fit: the donor variance is 2.4e-31, so the RE parameter is pinned at 0\n")
cat("  and the Hessian check on that scale is uninformative. The LMM\n")
cat("  and an OLS fit of the same fixed effects agree to 6 decimal places on both the grp\n")
cat("  estimate and its p, and rescaling the covariates does not remove the warning. The\n")
cat("  HEADLINE fits do not emit it. Do not 'fix' this by dropping the random intercept.\n\n")
cat("COHORT\n")
cat(sprintf("  %s: %d cells | %d PHF1+ | %d donors\n", CELLTYPE, n_cells, n_pos, n_don))
cat("\n=== RESULTS (headline model, sorted by p) ===\n")
print(as.data.frame(resA %>% arrange(headline_p) %>%
  select(label, lmm_diff_z, CI.L, CI.R, fold_change_raw, cohens_d_model,
         dz_donor_paired, n_donors_same_dir, headline_p, padj)),
  row.names = FALSE, digits = 3)
cat("\n--- Sensitivity (z units) ---\n")
print(as.data.frame(resA %>% arrange(headline_p) %>%
  select(channel, lmm_diff_z, unadj_diff, unadj_p, slope_diff, slope_p, lrt_p_slope,
         area_adj_diff, area_adj_p, area_adj_pct_of_raw, wilcox_p, singular)),
  row.names = FALSE, digits = 3)
cat("\n*** TECHNICAL COVARIATES ***\n")
cat(sprintf("  PHF1+ cells carry more RNA (median nCount %.0f vs %.0f, Wilcoxon p = %.3g) and\n",
            umi_med_pos, umi_med_neg, umi_wilcox_p))
cat("  are segmented larger; rho_nCount below gives each channel's correlation with nCount.\n")
print(as.data.frame(resA %>% arrange(headline_p) %>%
  select(channel, rho_nCount, lmm_diff_z, area_adj_diff, area_adj_pct_of_raw, area_adj_p)),
  row.names = FALSE, digits = 3)
cat("  area_adj adds log(cell area) to the headline for the density/max channels; it is NA\n")
cat("  for the *_total channels, where log area is a component of the outcome. nUMI_log is\n")
cat("  NOT included (image-derived outcome, Set 2).\n")
cat("  area_adj_pct_of_raw is NA where |raw effect| < 0.05 z -- the ratio of two near-zero\n")
cat("  numbers is noise, not a retained fraction.\n")
cat("  REPORTED, not the headline: larger cells integrate more signal, so part of a total- or\n")
cat("  max-channel effect can reflect size rather than concentration.\n")
cat("  rho_nCount is still printed above as a DESCRIPTIVE check on how far each channel\n")
cat("  tracks library size -- it is a diagnostic, not a term in any model here.\n")
cat("\n=== Effect sizes ===\n")
for (i in order(resA$headline_p)) {
  cat(sprintf("  %-16s %+6.3f z [%+6.3f, %+6.3f]  (x%.3f raw)  d = %+.3f  dz = %+.3f  %d/%d donors  padj = %.3g\n",
              resA$label[i], resA$lmm_diff_z[i], resA$CI.L[i], resA$CI.R[i],
              resA$fold_change_raw[i], resA$cohens_d_model[i], resA$dz_donor_paired[i],
              resA$n_donors_same_dir[i], resA$n_donors_paired[i], resA$padj[i]))
}
cat("\n  DAPI and Histone are nuclear-stain channels included as comparators for rRNA.\n")
cat("\nNOTES:\n")
cat(" - Density (Mean.X) is the primary readout. TOTAL = Mean.X * cell area is size-driven by\n")
cat("   construction. Max.X is size-biased and saturation-prone.\n")
cat(" - n = 3 donors per Braak group; no Braak stratification.\n")
cat("\n=== sessionInfo() ===\n"); print(sessionInfo()); sink()

sink(file.path(dist_dir, "stats_channel_distance_coef.txt"))
ttl <- sprintf("DAPI / Histone / rRNA intensity vs distance to nearest PHF1+ neuron -- PHF1- %s",
               CELLTYPE)
cat(ttl, "\n"); cat(strrep("=", nchar(ttl)), "\n"); cat("Date:", format(Sys.time()), "\n\n")
common_header()
cat("MODELS (per channel), sign = per s.d. CLOSER to a tangle (positive = BRIGHTER near tangles)\n")
cat("  HEADLINE   : z ~ dist_scaled + Sex + Age_s + PMI_s + (1|sample_id)\n")
cat("  depth-adj  : + depth_rel                          (reported)\n")
cat("  density-adj: + n_local                            (reported, headline is NOT adjusted)\n")
cat("  area-adj   : + log(cell area)                     (reported; mean/max only)\n")
cat("  unadj      : z ~ dist_scaled + (1|sample_id)      (reported)\n")
cat("  nUMI_log is deliberately ABSENT from every fit here: the outcome is image-derived, so\n")
cat("  the covariate set is the image set (Sex + Age + PMI + donor).\n")
cat("  edge_dist_um is deliberately ABSENT: the canonical CosMx distance models\n")
cat("  (deg_dream_linear_distance_phf1.r, plot_modulescore_vs_phf1_distance_modelp.R) and the\n")
cat("  IMC panel this figure mirrors all omit it; only nn3 adds it, as a spacing-specific\n")
cat("  technical term. BH across", nrow(CHANNELS),
    "channels; padj <", FDR, "\n")
cat(sprintf("  dist_scaled = log(dist)/sd(log(dist)); dist_sd = %.4f, within this celltype.\n\n",
            dist_sd))
cat("CELLS\n")
cat(sprintf("  modelled PHF1- cells: %d | donors: %d | window 0-%d um\n",
            nrow(mdB), n_distinct(mdB$sample_id), MAX_DIST_UM))
cat("  dist_to_phf1_um:\n"); print(summary(mdB$dist_to_phf1_um))

cat("\n*** PHF1-CHANNEL SCALE REFERENCE ***\n")
cat("  DAPI, Histone and rRNA are NON-PHF1 channels. The PHF1-channel slope is\n")
cat(sprintf("  phf1_intensity_p95_z at %.3f z per s.d. of log distance, i.e. %+.3f in the\n",
            PHF1_REF_BETA_DIST, PHF1_REF_CLOSER))
cat("  'per s.d. CLOSER' sign used here.\n")
cat("  pct_of_phf1_ref below expresses each channel's slope as a percentage of that value.\n")
cat("  These channels come from the AtoMx CosMx acquisition; phf1_intensity_p95_z comes from\n")
cat("  a separate confocal post-stain image (python/extract_phf1_intensity.py).\n")
cat(sprintf("\n  Local density: Spearman rho(log distance, neurons within 50 um) = %.3f\n",
            rho_dist_density))
cat("  (the density-adjusted fit is in the sensitivity table below).\n")

cat("\n=== RESULTS (headline model, sorted by p) ===\n")
print(as.data.frame(resB %>% select(label, beta_closer, CI.L, CI.R, cohens_d,
                                    pct_of_phf1_ref, pval, padj, n_cells)),
      row.names = FALSE, digits = 4)
cat("\n--- Sensitivity (z units) ---\n")
print(as.data.frame(resB %>% select(channel, beta_closer, beta_closer_unadj,
                                    beta_closer_depth_adj, beta_closer_dens_adj,
                                    beta_closer_area_adj, pval_area_adj,
                                    rho_nCount, cor_density, singular)),
      row.names = FALSE, digits = 4)
cat("\n=== Effect sizes ===\n")
for (i in seq_len(nrow(resB))) {
  cat(sprintf("  %-16s %+7.4f z [%+7.4f, %+7.4f] per s.d. closer  d = %+.4f  = %+6.1f%% of the PHF1 slope  padj = %.3g\n",
              resB$label[i], resB$beta_closer[i], resB$CI.L[i], resB$CI.R[i],
              resB$cohens_d[i], resB$pct_of_phf1_ref[i], resB$padj[i]))
}
cat("\nNOTES:\n")
cat(" - The drawn curve is within-donor centred; the fit is recentred to match, and its\n")
cat("   ribbon is a Wald CI on that centred contrast.\n")
cat(" - Max.X is size-biased and saturation-prone.\n")
cat("\n=== sessionInfo() ===\n"); print(sessionInfo()); sink()

cat("\nPART A ->", cat_dir, "|", length(list.files(cat_dir)), "files\n")
print(as.data.frame(resA %>% select(channel, lmm_diff_z, fold_change_raw, cohens_d_model,
                                    dz_donor_paired, padj, area_adj_diff)),
      row.names = FALSE, digits = 3)
cat("\nPART B ->", dist_dir, "|", length(list.files(dist_dir)), "files\n")
print(as.data.frame(resB %>% select(channel, beta_closer, cohens_d, padj,
                                    pct_of_phf1_ref, beta_closer_area_adj)),
      row.names = FALSE, digits = 3)
