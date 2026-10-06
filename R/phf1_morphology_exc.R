#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# phf1_morphology_exc.R
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
# phf1_morphology_exc.R
#
# CELL-AUTONOMOUS question: do PHF1+ (tangle-bearing) Exc-IT-L2-3-CBLN2-HOPX neurons
# differ in cell/nucleus morphology from PHF1- neurons of the same subtype?
#
# CosMx counterpart of the IMC nucleus-morphology contrast. The non-cell-autonomous
# counterpart -- morphology as a function of distance to the nearest tangle -- is
# R/phf1_distance_morphology_exc.R.
#
# HEADLINE MODEL: value ~ grp + (1 | sample_id), cell-level, lmerTest.
#   The random intercept is CORRECTLY CALIBRATED for a within-donor predictor -- this
#   was verified empirically by permuting PHF1 within donor.
#   It answers "do PHF1+ cells differ from PHF1- cells WITHIN a donor?", replication
#   unit = cell. The random-SLOPE fit (1 + grp | sample_id) answers the different
#   question "is the average DONOR's effect non-zero?", replication unit = donor; it is
#   reported alongside with its LRT, and is NOT a bug fix for the intercept model.
#   The IMC script uses a donor-mean paired Wilcoxon as its headline because it has 44
#   donors. Here there are 9, three of which contribute <= 3 PHF1+ cells, so the
#   donor-level test is reported as a supporting statistic rather than the headline.
#
# MEASURES: NucArea, NucAspectRatio, Circularity, Eccentricity, Perimeter.
#   'Solidity' is deliberately EXCLUDED. The AtoMx column of that name is not a
#   convex-hull solidity: it equals Area/Perimeter exactly (verified on all 26,585
#   CBLN2 cells; range 3.3-72.7, not 0-1). Do not re-add it.
#   Circularity is likewise exactly 4*pi*Area/Perimeter^2, so Circularity and Perimeter
#   between them re-express cell Area -- they are not independent shape measures.
#
# UNITS: Area.um2/Area = 0.014467506 for every cell => 1 px = 0.120281 um.
#   NucArea (px^2) * 0.014467506 -> um^2;  Perimeter (px) * 0.120281 -> um.
#
# Triple-output convention, output under plots/phf1_morphology/.
# Read-only on all inputs.

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(qs)
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
  library(lme4); library(lmerTest); library(emmeans)
})

hpc <- "<PROJECT_ROOT>/phf1_v2"
loc <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(hpc)) hpc else loc)
source("R/palettes.R")   # fig_theme; braak_palette["6"] == "#BD0026" is the project PHF1 red
source("R/phf1_morphology_filters.R")  # PX_UM, PX2_UM2, NUC_MIN_UM2, nucleus_ok(), audit

out_dir <- "plots/phf1_morphology"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

CELLTYPE  <- "Exc-IT-L2-3-CBLN2-HOPX"
PX_UM     <- 0.120281        # microns per pixel
PX2_UM2   <- 0.014467506     # square microns per square pixel (Area.um2 / Area)
FDR       <- 0.05
SIG_BASIS <- "padj<0.05"

phf1_fill <- c("PHF1-" = "grey75", "PHF1+" = "#BD0026")

# emmeans defaults to Kenward-Roger and silently falls back to asymptotic z above 3000
# observations. There are ~26k cells here, so ask for Satterthwaite explicitly and raise
# the limit, otherwise the contrast df is reported as Inf.
emmeans::emm_options(lmer.df = "satterthwaite", lmerTest.limit = 1e5)

wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE)

## ---------------------------------------------------------------------------
## 1. Load
##    The per-celltype SCE (51 MB) is used in preference to seu_PHF1.rds (7.8 GB):
##    it is the canonical DEG object for this celltype and carries every column
##    needed here (all morphology metrics + PHF1 + covariates). Asserted below.
## ---------------------------------------------------------------------------
sce_path <- file.path("celltype_sce_neighbours",
                      paste0(CELLTYPE, "_sce_neighbours.qs"))
if (!file.exists(sce_path)) stop("Not found: ", sce_path)
sce <- qread(sce_path)

cd <- as.data.frame(colData(sce))
cd$cell_id <- colnames(sce)

# --- pixel-scale assertion: this is what fixes PX_UM / PX2_UM2 ---------------
scale_err <- max(abs(cd$Area.um2 / cd$Area - PX2_UM2), na.rm = TRUE)
if (!is.finite(scale_err) || scale_err > 1e-6)
  stop(sprintf("Area.um2/Area is not the constant %.9f (max deviation %.3g). ",
               PX2_UM2, scale_err),
       "The px->um conversion used for NucArea and Perimeter is invalid.")

# --- the algebraic identities, re-checked at runtime so they cannot rot ------
id_solidity   <- max(abs(cd$Solidity - cd$Area / cd$Perimeter), na.rm = TRUE)
id_circularity <- max(abs(cd$Circularity - 4 * pi * cd$Area / cd$Perimeter^2), na.rm = TRUE)

md <- cd %>%
  transmute(
    cell_id   = cell_id,
    sample_id = as.character(sample_id),
    celltype  = as.character(celltype),
    Braak     = factor(as.character(Braak), levels = braak_levels),
    Sex       = factor(as.character(Sex)),
    Age       = as.numeric(Age),
    PMI       = as.numeric(PMI),
    phf1_pos  = as.logical(PHF1 %in% c(TRUE, "TRUE", "True", 1, "1")),
    nCount_RNA = as.numeric(nCount_RNA),
    nUMI_log   = log2(as.numeric(nCount_RNA) + 1),
    # morphology, converted to physical units where the source column is in pixels
    nucarea      = ifelse(nucleus_ok(NucArea), NucArea * PX2_UM2, NA_real_),
    nucaspect    = ifelse(nucleus_ok(NucArea), NucAspectRatio, NA_real_),
    circularity  = Circularity,
    eccentricity = Eccentricity,
    perimeter    = Perimeter * PX_UM,
    # nuc_ok is the filter ACTUALLY applied to the nuclear measures (segmented AND
    # plausibly sized), so the bias check below tests the real conditioning rather than
    # the weaker "was anything segmented" version. See R/phf1_morphology_filters.R.
    nuc_ok       = nucleus_ok(NucArea)
  ) %>%
  mutate(grp   = factor(ifelse(phf1_pos, "PHF1+", "PHF1-"), levels = c("PHF1-", "PHF1+")),
         Age_s = as.numeric(scale(Age)),
         PMI_s = as.numeric(scale(PMI)))

n_cells <- nrow(md); n_pos <- sum(md$phf1_pos); n_don <- n_distinct(md$sample_id)
cat(sprintf("Loaded %s: %d cells | %d PHF1+ | %d donors\n",
            CELLTYPE, n_cells, n_pos, n_don))
if (n_don != 9)
  stop("Expected 9 donors, got ", n_don, " -- the cohort has changed; re-check the object.")
if (n_pos < 1) stop("No PHF1+ cells in ", CELLTYPE)

## ---------------------------------------------------------------------------
## 2. Nuclear-segmentation bias check.
##    26% of cells have NucArea == 0 (= no nucleus segmented). The two nucleus
##    measures are conditioned on successful segmentation, so if that conditioning
##    is itself PHF1-dependent, both nucleus results would reflect that selection.
## ---------------------------------------------------------------------------
seg_or <- NA_real_; seg_lo <- NA_real_; seg_hi <- NA_real_; seg_p <- NA_real_
m_seg <- tryCatch(
  lme4::glmer(nuc_ok ~ grp + (1 | sample_id), data = md, family = binomial),
  error = function(e) NULL, warning = function(w) NULL)
if (!is.null(m_seg)) {
  co <- summary(m_seg)$coefficients
  if ("grpPHF1+" %in% rownames(co)) {
    est <- co["grpPHF1+", "Estimate"]; se <- co["grpPHF1+", "Std. Error"]
    seg_or <- exp(est); seg_lo <- exp(est - 1.96 * se); seg_hi <- exp(est + 1.96 * se)
    seg_p  <- co["grpPHF1+", "Pr(>|z|)"]
  }
}
seg_tab <- md %>% group_by(grp) %>%
  summarise(n = n(), n_nuc_ok = sum(nuc_ok), pct_nuc_ok = 100 * mean(nuc_ok),
            .groups = "drop")

## Library size by PHF1 status, quantified once here and reported against every measure.
umi_med_pos <- median(md$nCount_RNA[md$phf1_pos])
umi_med_neg <- median(md$nCount_RNA[!md$phf1_pos])
umi_wilcox_p <- suppressWarnings(
  wilcox.test(nCount_RNA ~ phf1_pos, data = md)$p.value)

## ---------------------------------------------------------------------------
## 3. Features
## ---------------------------------------------------------------------------
features <- tibble::tribble(
  ~slug,          ~col,           ~unit,       ~needs_nucleus,
  "nucarea",      "nucarea",      "um^2",      TRUE,
  "nucaspect",    "nucaspect",    "unitless",  TRUE,
  "circularity",  "circularity",  "unitless",  FALSE,
  "eccentricity", "eccentricity", "unitless",  FALSE,
  "perimeter",    "perimeter",    "um",        FALSE
)
features$label <- c("Nucleus area", "Nucleus aspect ratio", "Cell circularity",
                    "Cell eccentricity", "Cell perimeter")
stopifnot(all(features$col %in% colnames(md)))

axis_label <- function(slug) {
  switch(slug,
    nucarea      = expression("Nucleus area (" * mu * "m"^2 * ")"),
    nucaspect    = "Nucleus aspect ratio",
    circularity  = "Cell circularity",
    eccentricity = "Cell eccentricity",
    perimeter    = expression("Cell perimeter (" * mu * "m)"),
    slug)
}

## ---------------------------------------------------------------------------
## 4. Fit every feature (pass 1). padj spans features, so nothing is written or
##    plotted until all fits are in.
## ---------------------------------------------------------------------------
fit_one <- function(i) {
  ft <- features[i, ]
  d  <- md %>%
    transmute(cell_id, sample_id, Braak, Sex, Age_s, PMI_s, nUMI_log, nCount_RNA,
              celltype, grp, value = .data[[ft$col]])
  d <- d[is.finite(d$value), ]

  # per-donor means -- the donor-level supporting statistics and the figure's dots
  dm <- d %>% group_by(sample_id, Braak, grp) %>%
    summarise(donor_mean = mean(value), n_cells = n(), .groups = "drop")
  dm_wide <- dm %>% select(sample_id, grp, donor_mean) %>%
    tidyr::pivot_wider(names_from = grp, values_from = donor_mean)
  dm_wide <- dm_wide[stats::complete.cases(dm_wide), ]

  fitm <- function(f, ...) tryCatch(lmerTest::lmer(f, data = d, REML = TRUE, ...),
                                    error = function(e) NULL)
  m_ri  <- fitm(value ~ grp + (1 | sample_id))                                   # HEADLINE
  m_rs  <- fitm(value ~ grp + (1 + grp | sample_id),
                control = lmerControl(optimizer = "bobyqa",
                                      optCtrl = list(maxfun = 2e5)))             # reported
  m_cov <- fitm(value ~ grp + Sex + Age_s + PMI_s + (1 | sample_id))             # reported
  m_umi <- fitm(value ~ grp + nUMI_log + (1 | sample_id))                        # reported
  lrt   <- if (!is.null(m_ri) && !is.null(m_rs))
    tryCatch(anova(m_ri, m_rs, refit = FALSE), error = function(e) NULL) else NULL

  emm_df <- NULL; pair_df <- NULL; omni <- NULL
  p_lmm <- NA_real_; diff_lmm <- NA_real_; diff_lo <- NA_real_; diff_hi <- NA_real_
  d_model <- NA_real_
  if (!is.null(m_ri)) {
    emm     <- emmeans(m_ri, ~ grp)
    emm_df  <- as.data.frame(summary(emm, infer = TRUE))
    pair_df <- as.data.frame(summary(pairs(emm), infer = TRUE))
    omni    <- tryCatch(anova(m_ri), error = function(e) NULL)
    # pairs() returns (PHF1-) - (PHF1+); re-orient as PHF1+ minus PHF1-
    diff_lmm <- -pair_df$estimate[1]
    diff_lo  <- -pair_df$upper.CL[1]
    diff_hi  <- -pair_df$lower.CL[1]
    p_lmm    <- pair_df$p.value[1]
    vc       <- as.data.frame(lme4::VarCorr(m_ri))
    d_model  <- diff_lmm / sqrt(sum(vc$vcov[is.na(vc$var2)]))   # drop covariance rows
  }

  # random-slope contrast, the donor-generalisation answer
  sl_d <- NA_real_; sl_lo <- NA_real_; sl_hi <- NA_real_; sl_p <- NA_real_
  if (!is.null(m_rs)) {
    ps <- tryCatch(as.data.frame(summary(pairs(emmeans(m_rs, ~ grp)), infer = TRUE)),
                   error = function(e) NULL)
    if (!is.null(ps)) {
      sl_d <- -ps$estimate[1]; sl_lo <- -ps$upper.CL[1]; sl_hi <- -ps$lower.CL[1]
      sl_p <- ps$p.value[1]
    }
  }
  # covariate-adjusted contrasts
  grab_diff <- function(m) {
    if (is.null(m)) return(c(NA_real_, NA_real_))
    pc <- tryCatch(as.data.frame(summary(pairs(emmeans(m, ~ grp)), infer = TRUE)),
                   error = function(e) NULL)
    if (is.null(pc)) return(c(NA_real_, NA_real_))
    c(-pc$estimate[1], pc$p.value[1])
  }
  cv <- grab_diff(m_cov); cov_d <- cv[1]; cov_p <- cv[2]
  uv <- grab_diff(m_umi); umi_d <- uv[1]; umi_p <- uv[2]
  # association of this measure with library size
  rho_umi <- suppressWarnings(cor(d$value, d$nCount_RNA, method = "spearman"))

  # --- effect sizes (always report one alongside the p) ----------------------
  n_neg <- sum(d$grp == "PHF1-"); n_pos_f <- sum(d$grp == "PHF1+")
  m_neg <- mean(d$value[d$grp == "PHF1-"]); m_pos <- mean(d$value[d$grp == "PHF1+"])
  sd_pooled <- sqrt(((n_neg - 1) * var(d$value[d$grp == "PHF1-"]) +
                     (n_pos_f - 1) * var(d$value[d$grp == "PHF1+"])) / (n_neg + n_pos_f - 2))
  d_cell <- (m_pos - m_neg) / sd_pooled

  # donor-level paired statistics (supporting, not the headline: n = 9)
  dz <- NA_real_; dz_mean <- NA_real_; dz_lo <- NA_real_; dz_hi <- NA_real_
  n_same <- NA_integer_; w_p <- NA_real_; pct_neg <- NA_real_
  if (nrow(dm_wide) >= 3) {
    dz_diff <- dm_wide[["PHF1+"]] - dm_wide[["PHF1-"]]
    dz      <- mean(dz_diff) / sd(dz_diff)
    dz_mean <- mean(dz_diff)
    n_same  <- sum(dz_diff > 0)
    pt <- tryCatch(t.test(dm_wide[["PHF1+"]], dm_wide[["PHF1-"]], paired = TRUE),
                   error = function(e) NULL)
    if (!is.null(pt)) { dz_lo <- pt$conf.int[1]; dz_hi <- pt$conf.int[2] }
    pw <- tryCatch(wilcox.test(dm_wide[["PHF1+"]], dm_wide[["PHF1-"]], paired = TRUE),
                   error = function(e) NULL)
    if (!is.null(pw)) w_p <- pw$p.value
    pct_neg <- 100 * dz_mean / mean(dm_wide[["PHF1-"]])
  }

  row <- tibble(
    metric = ft$slug, label = ft$label, unit = ft$unit,
    headline_test = "cell-level LMM: value ~ grp + (1 | sample_id)",
    headline_p = p_lmm,
    mean_PHF1_neg = m_neg, mean_PHF1_pos = m_pos,
    lmm_diff = diff_lmm, lmm_diff_lower = diff_lo, lmm_diff_upper = diff_hi,
    pct_of_PHF1_neg_lmm = 100 * diff_lmm / m_neg,
    cohens_d_model = d_model, cohens_d_cell = d_cell,
    dz_donor_paired = dz, dz_mean_diff = dz_mean,
    dz_mean_diff_lower = dz_lo, dz_mean_diff_upper = dz_hi,
    pct_of_PHF1_neg_donor = pct_neg,
    n_donors_paired = nrow(dm_wide), n_donors_same_dir = n_same, wilcox_p = w_p,
    slope_diff = sl_d, slope_diff_lower = sl_lo, slope_diff_upper = sl_hi,
    slope_p = sl_p,
    lrt_p_slope = if (is.null(lrt)) NA_real_ else lrt$`Pr(>Chisq)`[2],
    cov_adj_diff = cov_d, cov_adj_p = cov_p,
    umi_adj_diff = umi_d, umi_adj_p = umi_p,
    umi_adj_pct_of_raw = 100 * umi_d / diff_lmm,
    rho_nCount = rho_umi,
    n_cells = nrow(d), n_PHF1_pos = n_pos_f,
    singular = if (is.null(m_ri)) NA else lme4::isSingular(m_ri)
  )
  list(row = row, d = d, dm = dm, dm_wide = dm_wide, emm_df = emm_df,
       pair_df = pair_df, omni = omni, m_ri = m_ri, m_rs = m_rs, m_cov = m_cov,
       lrt = lrt, ft = ft)
}

fits <- lapply(seq_len(nrow(features)), fit_one)
res  <- bind_rows(lapply(fits, `[[`, "row"))
res$padj <- p.adjust(res$headline_p, method = "BH")
res$significant <- !is.na(res$padj) & res$padj < FDR
for (i in seq_along(fits)) fits[[i]]$row <- res[i, ]

## ---------------------------------------------------------------------------
## 5. Per-feature outputs (pass 2)
## ---------------------------------------------------------------------------
for (i in seq_along(fits)) {
  f <- fits[[i]]; ft <- f$ft; slug <- paste0("phf1_morph_", ft$slug, "_CBLN2")
  r <- res[i, ]

  wr(f$d %>% select(cell_id, sample_id, Braak, Sex, celltype, grp, value),
     paste0("source_data_", slug, ".tsv"))
  wr(f$dm, paste0("source_data_", slug, "_donor_means.tsv"))
  if (!is.null(f$pair_df)) wr(f$pair_df, paste0("stats_", slug, "_pairwise.tsv"))
  wr(r, paste0("stats_", slug, "_effectsize.tsv"))

  # ---- stats log ----
  sink(file.path(out_dir, paste0("stats_", slug, ".txt")))
  ttl <- sprintf("%s in PHF1+ vs PHF1- %s neurons (CosMx)", ft$label, CELLTYPE)
  cat(ttl, "\n"); cat(strrep("=", nchar(ttl)), "\n")
  cat("Date:", format(Sys.time()), "\n\n")
  cat("Source: ", sce_path, "\n", sep = "")
  cat("CosMx counterpart of the IMC nucleus-morphology contrast (cell-autonomous half).\n")
  cat("Non-cell-autonomous counterpart: plots/phf1_distance_morphology_1000um/.\n\n")

  cat("MEASURE: ", ft$col, " (", ft$unit, ")\n", sep = "")
  if (ft$col == "nucarea")   cat("  = NucArea (px^2) * ", PX2_UM2, "\n", sep = "")
  if (ft$col == "perimeter") cat("  = Perimeter (px) * ", PX_UM, "\n", sep = "")
  cat("\nMODELS\n")
  cat("  HEADLINE  : value ~ grp + (1 | sample_id)                  [cell-level, lmerTest]\n")
  cat("  slope     : value ~ grp + (1 + grp | sample_id)            (reported)\n")
  cat("  cov-adj   : value ~ grp + Sex + Age_s + PMI_s + (1|sample_id)  (reported)\n")
  cat("  donor-level: paired Wilcoxon / paired t on the 9 per-donor means (reported)\n")
  cat("  BH across the", nrow(features), "measures; significance padj <", FDR, "\n\n")
  cat("  The random intercept is CALIBRATED for this within-donor predictor (verified by\n")
  cat("  within-donor permutation). Intercept vs slope is a choice of INFERENTIAL TARGET\n")
  cat("  (cell vs donor), not a correctness fix.\n\n")

  cat("CELLS\n")
  cat(sprintf("  modelled: %d | PHF1+: %d | PHF1-: %d | donors: %d\n",
              r$n_cells, r$n_PHF1_pos, r$n_cells - r$n_PHF1_pos, r$n_donors_paired))
  cat("  value summary:\n"); print(summary(f$d$value))
  cat("\n  PHF1+ cells per donor:\n")
  print(as.data.frame(f$dm %>% filter(grp == "PHF1+") %>%
                        select(sample_id, Braak, n_cells)), row.names = FALSE)

  if (isTRUE(ft$needs_nucleus)) {
    cat_nucleus_filter_note()
    print(as.data.frame(seg_tab), row.names = FALSE)
    cat(sprintf("  glmer(nuc_ok ~ grp + (1|sample_id), binomial): OR(PHF1+) = %.3f [%.3f, %.3f], p = %.3g\n",
                seg_or, seg_lo, seg_hi, seg_p))
    cat("  An OR far from 1 would indicate PHF1-dependent selection on this measure.\n")
  }

  cat("\n=== Effect sizes ===\n")
  cat(sprintf("  LMM difference (PHF1+ minus PHF1-): %+.4g [%+.4g, %+.4g] %s = %+.2f%% of the PHF1- mean\n",
              r$lmm_diff, r$lmm_diff_lower, r$lmm_diff_upper, ft$unit, r$pct_of_PHF1_neg_lmm))
  cat(sprintf("  Cohen's d (model-based, diff / sqrt(sum VarCorr)) = %+.3f\n", r$cohens_d_model))
  cat(sprintf("  Cohen's d (cell-level, pooled SD)                 = %+.3f\n", r$cohens_d_cell))
  cat(sprintf("  Donor paired dz = %+.3f | mean donor diff %+.4g [%+.4g, %+.4g] = %+.2f%%\n",
              r$dz_donor_paired, r$dz_mean_diff, r$dz_mean_diff_lower,
              r$dz_mean_diff_upper, r$pct_of_PHF1_neg_donor))
  cat(sprintf("  Donors moving in the same direction: %d / %d | paired Wilcoxon p = %.3g\n",
              r$n_donors_same_dir, r$n_donors_paired, r$wilcox_p))
  cat(sprintf("  HEADLINE p = %.3g | padj = %.3g | %s\n",
              r$headline_p, r$padj, ifelse(r$significant, SIG_BASIS, "n.s.")))

  cat("\n=== Omnibus (Satterthwaite F, headline model) ===\n")
  if (!is.null(f$omni)) print(f$omni) else cat("  model did not converge\n")
  cat("\n=== EMMs (+/- 95% CI) ===\n")
  if (!is.null(f$emm_df)) print(f$emm_df, row.names = FALSE)
  cat("\n=== Contrast ===\n")
  if (!is.null(f$pair_df)) print(f$pair_df, row.names = FALSE)
  cat("\n=== Random-slope fit (donor-generalisation target) ===\n")
  cat(sprintf("  diff = %+.4g [%+.4g, %+.4g], p = %.3g\n",
              r$slope_diff, r$slope_diff_lower, r$slope_diff_upper, r$slope_p))
  cat("  LRT intercept-only vs slope: p =", format(r$lrt_p_slope, digits = 3), "\n")
  if (!is.null(f$lrt)) print(f$lrt)
  cat("\n=== Covariate-adjusted fit (Sex + Age + PMI) ===\n")
  cat(sprintf("  diff = %+.4g, p = %.3g\n", r$cov_adj_diff, r$cov_adj_p))

  cat("\n=== Library size ===\n")
  cat(sprintf("  Spearman rho(%s, nCount_RNA) = %+.3f\n", ft$col, r$rho_nCount))
  cat(sprintf("  median nCount PHF1+ %.0f vs PHF1- %.0f (Wilcoxon p = %.3g)\n",
              umi_med_pos, umi_med_neg, umi_wilcox_p))
  cat(sprintf("  Library-size-adjusted fit (value ~ grp + nUMI_log + (1|sample_id)):\n"))
  cat(sprintf("    diff = %+.4g, p = %.3g  -- %.0f%% of the unadjusted %+.4g\n",
              r$umi_adj_diff, r$umi_adj_p, r$umi_adj_pct_of_raw, r$lmm_diff))
  cat("  This is REPORTED, not the headline: nUMI_log is a collider on the PHF1 contrast (cell\n")
  cat("  size plausibly causes both), so adjusting for it can remove part of the effect under\n")
  cat("  test (docs/MODELS.md, Set 2).\n")

  cat("\nNOTES:\n")
  cat(" - AtoMx 'Solidity' is Area/Perimeter, not a solidity, and is EXCLUDED.\n")
  cat("   Circularity is exactly 4*pi*Area/Perimeter^2, so Circularity and Perimeter\n")
  cat("   between them re-express cell Area -- not independent shape measures.\n")
  cat(sprintf("   (runtime check: max|Solidity - Area/P| = %.3g, max|Circ - 4piA/P^2| = %.3g)\n",
              id_solidity, id_circularity))
  cat(" - These are AtoMx segmentation metrics.\n")
  cat(" -", n_pos, "PHF1+ cells in this subtype, spread 1-60 per donor; three donors\n")
  cat("   contribute <= 3.\n")
  cat(" - n = 3 donors per Braak group; no Braak stratification is attempted.\n")
  cat(" - The violin is clipped at the 99th percentile for display only; every cell is in\n")
  cat("   the model and in source_data_*.tsv.\n")

  cat("\n=== sessionInfo() ===\n"); print(sessionInfo())
  sink()

  # ---- figure: violin + per-donor means -------------------------------------
  d <- f$d; dm <- f$dm
  grp_mean <- dm %>% group_by(grp) %>%
    summarise(m = mean(donor_mean), sem = sd(donor_mean) / sqrt(n()), .groups = "drop") %>%
    mutate(xc = as.integer(grp))

  v_lo <- min(d$value, na.rm = TRUE)
  v_hi <- max(quantile(d$value, 0.99, na.rm = TRUE), max(dm$donor_mean, na.rm = TRUE))
  v_rg <- v_hi - v_lo
  v_br <- v_hi + v_rg * 0.06     # bracket position on the value axis
  tick <- v_rg * 0.02            # end ticks, pointing back toward the data
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
    # Donor dot size ~ cells contributing, so the near-singleton PHF1+ donors read as weak.
    # log10, not scale_size_area: PHF1+ donors carry 1-60 cells against ~2900 for PHF1-,
    # and an area-proportional scale shrinks every PHF1+ dot to invisibility.
    scale_size(transform = "log10", range = c(0.7, 2.0), name = "Cells per donor",
               breaks = c(1, 10, 100, 1000)) +
    coord_cartesian(xlim = c(v_lo, v_br + v_rg * 0.62)) +
    labs(x = axis_label(ft$slug), y = NULL) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "bottom",
          legend.key.height = grid::unit(7, "pt"),
          legend.box.spacing = grid::unit(2, "pt"),
          legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 4, unit = "pt"))

  ggsave(file.path(out_dir, paste0("plot_", slug, ".pdf")), p,
         width = 6.14, height = 4.0, units = "cm", device = "pdf")
}

## ---------------------------------------------------------------------------
## 6. Combined outputs (BH spans features, so these are the authoritative tables)
## ---------------------------------------------------------------------------
wr(res, "stats_phf1_morphology_CBLN2_effectsize.tsv")

sink(file.path(out_dir, "stats_phf1_morphology_CBLN2.txt"))
ttl <- sprintf("Cell/nucleus morphology in PHF1+ vs PHF1- %s neurons (CosMx) -- all measures",
               CELLTYPE)
cat(ttl, "\n"); cat(strrep("=", nchar(ttl)), "\n")
cat("Date:", format(Sys.time()), "\n\n")
cat("Source: ", sce_path, "\n", sep = "")
cat("Script: R/phf1_morphology_exc.R\n")
cat("CosMx counterpart of the IMC nucleus-morphology contrast.\n\n")

cat("MEASURES AND THEIR ALGEBRA -- READ FIRST\n")
cat("  AtoMx ships six segmentation metrics. 'Solidity' is EXCLUDED from this analysis:\n")
cat("  it is not a convex-hull solidity, it is exactly Area/Perimeter (range 3.3-72.7).\n")
cat("  Circularity is exactly 4*pi*Area/Perimeter^2. So Circularity and Perimeter between\n")
cat("  them re-express cell Area; do NOT read them as two independent shape findings.\n")
cat(sprintf("  Runtime identity checks: max|Solidity - Area/P| = %.3g | max|Circ - 4piA/P^2| = %.3g\n",
            id_solidity, id_circularity))
cat(sprintf("  Pixel scale: Area.um2/Area = %.9f exactly (max dev %.3g) => 1 px = %.6f um.\n",
            PX2_UM2, scale_err, PX_UM))
cat("  NucArea and Perimeter are converted to um^2 / um; the rest are unitless.\n\n")

cat("COHORT\n")
cat(sprintf("  %s: %d cells | %d PHF1+ | %d donors (3 x Braak 1, 3 x 4, 3 x 6)\n",
            CELLTYPE, n_cells, n_pos, n_don))
cat("  PHF1+ per donor:\n")
print(as.data.frame(md %>% filter(phf1_pos) %>% count(sample_id, Braak, name = "n_PHF1_pos")),
      row.names = FALSE)

cat_nucleus_filter_note(nucleus_filter_audit(cd$NucArea, group = "all CBLN2"))
cat("\nNUCLEAR SEGMENTATION (pct_nuc_ok = passes the full validity filter)\n")
print(as.data.frame(seg_tab), row.names = FALSE)
cat(sprintf("  glmer(nuc_ok ~ grp + (1|sample_id), binomial): OR(PHF1+) = %.3f [%.3f, %.3f], p = %.3g\n",
            seg_or, seg_lo, seg_hi, seg_p))
cat("  Nucleus measures (nucarea, nucaspect) are conditioned on this; an OR far from 1\n")
cat("  would indicate PHF1-dependent selection on them.\n")

cat("\n=== RESULTS (headline model, sorted by p) ===\n")
print(as.data.frame(res %>% arrange(headline_p) %>%
  select(metric, unit, lmm_diff, lmm_diff_lower, lmm_diff_upper, pct_of_PHF1_neg_lmm,
         cohens_d_model, cohens_d_cell, dz_donor_paired, n_donors_same_dir,
         headline_p, padj, n_cells)),
  row.names = FALSE, digits = 3)

cat("\n--- Sensitivity (same units as lmm_diff) ---\n")
print(as.data.frame(res %>% arrange(headline_p) %>%
  select(metric, lmm_diff, slope_diff, slope_p, lrt_p_slope, cov_adj_diff, cov_adj_p,
         wilcox_p, singular)), row.names = FALSE, digits = 3)

cat("\n=== Library size ===\n")
cat(sprintf("  median nCount_RNA PHF1+ %.0f vs PHF1- %.0f (Wilcoxon p = %.3g).\n",
            umi_med_pos, umi_med_neg, umi_wilcox_p))
print(as.data.frame(res %>% arrange(headline_p) %>%
  select(metric, rho_nCount, lmm_diff, umi_adj_diff, umi_adj_pct_of_raw, umi_adj_p)),
  row.names = FALSE, digits = 3)
cat("  umi_adj_* is value ~ grp + nUMI_log + (1|sample_id). REPORTED, not the headline:\n")
cat("  nUMI_log is a collider on the PHF1 contrast (cell size plausibly causes both size and\n")
cat("  transcript count), so adjusting for it can remove part of the effect under test\n")
cat("  (docs/MODELS.md, Set 2).\n")

cat("\n=== Effect sizes ===\n")
for (i in order(res$headline_p)) {
  cat(sprintf("  %-14s %+9.4g [%+9.4g, %+9.4g] %-9s (%+6.2f%%)  d = %+.3f  dz = %+.3f  %d/%d donors  padj = %.3g\n",
              res$metric[i], res$lmm_diff[i], res$lmm_diff_lower[i], res$lmm_diff_upper[i],
              res$unit[i], res$pct_of_PHF1_neg_lmm[i], res$cohens_d_model[i],
              res$dz_donor_paired[i], res$n_donors_same_dir[i], res$n_donors_paired[i],
              res$padj[i]))
}

cat("\nNOTES:\n")
cat(" - AtoMx segmentation metrics.\n")
cat(" - Library-size-adjusted fits are reported in the block above.\n")
cat(" - 185 PHF1+ cells, 1-60 per donor; three donors contribute <= 3.\n")
cat(" - n = 3 donors per Braak group; no Braak stratification attempted.\n")

cat("\n=== sessionInfo() ===\n"); print(sessionInfo())
sink()

cat("\nDone. Wrote", length(list.files(out_dir)), "files to", out_dir, "\n")
print(as.data.frame(res %>% select(metric, lmm_diff, cohens_d_model, dz_donor_paired,
                                   headline_p, padj)), row.names = FALSE, digits = 3)
