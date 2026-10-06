#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# nn3_phf1_vs_neg.R
#
# Figure panels: S3 (all) - first step of run_nn3.sh
# Also writes the per-neuron covariate cache
# (plots/nn3_neuron_spacing/source_data_nn3_cells.tsv) read by the 5C and S7B scripts.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# nn3_phf1_vs_neg.R
#
# ANALYSIS A (CELL-AUTONOMOUS): do PHF1+ (ptau+) neurons sit in SPARSER
# neighbourhoods than PHF1- neurons?
#
# Outcome: nn3_um, the mean distance to the 3 nearest other neurons (the
# Zwang et al. 2024 3-nearest-neighbour spacing metric; see R/nn3_utils.r for the
# full provenance and the 2D-vs-3D caveats). Computed over the 10-label neuron
# pool, identical to the PHF1+ distance-source set, so nn3_um and the canonical
# dist_to_phf1_um describe the same neuron population.
#
# HEADLINE TESTS (two, reported together):
#   cell level  : lmer(log(nn3_um) ~ PHF1 + celltype + edge_dist_um + (1|sample_id))
#                 emmeans contrast, back-transformed to a % difference
#   donor level : per donor, mean log(nn3_um) PHF1+ vs PHF1- -> paired Wilcoxon
#                 (n = 9) + paired t for the CI + Cohen's dz
#
# Per-donor PHF1+ counts span 4 to 118, so an unweighted paired test over 9 donors gives
# donors with few PHF1+ cells the same weight as the rest, while the LMM weights by
# information. Two further donor-level statistics are therefore reported: an
# inverse-variance weighted donor estimate (the donor-level analogue of what the LMM
# estimates) and the unweighted test restricted to donors with >= 20 PHF1+ neurons.
#
# REPORTED BUT NOT HEADLINE
#   + depth_rel            cortical-depth adjusted. PHF1+ neurons are layer-biased
#                          (255 of 390 are L2-3 / L3-5) and spacing varies with
#                          layer, so this is the most important sensitivity fit.
#   PHF1 x Braak           does any effect depend on stage (n = 3 donors/stage).
#   nn3rd_um               distance to the 3rd nearest rather than the mean of 3,
#                          the alternative reading of the Zwang et al. metric.
#   per subtype            fitted ONLY for subtypes that pass the project-wide
#                          eligibility criterion (>= 6 PHF1+ cells in >= 5 samples,
#                          R/eligible_neuron_subtypes.R), BH across those. On this
#                          cohort that is Exc-IT-L2-3-CBLN2-HOPX and
#                          Exc-IT-L3-5-CHGA-IL1RAPL2 only. The pooled model is NOT
#                          restricted.
#
# WHY edge_dist_um IS A COVARIATE AND NOT AN EXCLUSION. CosMx FOVs do not tile
# the tissue, so a neuron beside an unimaged region has invisible neighbours and
# an inflated nn3_um. The exposure is adjusted for rather than cells dropped, so
# nothing is lost. Verified mild --
# rho(nn3_um, edge_dist_um_raw) = +0.03, and only ~5% of neurons sit within
# 100 um of an unimaged region.
#
# WHY NO PERMUTATION NULL HERE. The 1000-labelling null is a SUPPLEMENTARY
# analysis in R/null_replica_nn3.R; the primary models carry model-based
# inference only.
#
# Outputs, under plots/nn3_neuron_spacing/:
#   source_data_nn3_cells.tsv                  SHARED CACHE -- one row per pool
#                                              neuron, read by nn3_over_phf1_distance.R
#                                              and null_replica_nn3.R so neither
#                                              re-reads the 327 MB SCE.
#   plot_nn3_phf1_vs_neg.pdf                   pooled, donor means overlaid
#   source_data_nn3_phf1_vs_neg.tsv            drawn cells
#   source_data_nn3_phf1_vs_neg_donor_means.tsv
#   stats_nn3_phf1_vs_neg.txt
#   stats_nn3_phf1_vs_neg_pairwise.tsv
#   stats_nn3_phf1_vs_neg_effectsize.tsv
#   plot_nn3_phf1_vs_neg_by_celltype.pdf       per-subtype forest
#   plot_nn3_phf1_vs_neg_donor_forest.pdf      per-donor contrast + inv-var pooled estimate
#   source_data_nn3_phf1_vs_neg_donor_forest.tsv
#   source_data_nn3_phf1_vs_neg_by_celltype.tsv
#   stats_nn3_phf1_vs_neg_by_celltype.txt
#
# Run: Rscript R/nn3_phf1_vs_neg.R    (~1 min, dominated by the SCE read)

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(qs)
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
  library(lmerTest); library(emmeans); library(RANN)
})

hpc <- "<PROJECT_ROOT>/phf1_v2"
loc <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(hpc)) hpc else loc)
source("R/palettes.R")
source("R/nn3_utils.r")
source("R/eligible_neuron_subtypes.R")

out_dir <- "plots/nn3_neuron_spacing"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

K            <- 3L      # the Zwang et al. metric
EDGE_CAP_UM  <- 300     # edge_dist_um cap; justified against the nn3 upper tail below
RASTER_UM    <- 10      # FOV-union raster resolution for edge_dist_um
FDR          <- 0.05
# Per-subtype fits are gated SOLELY by eligible_neuron_subtypes() (>= 6 PHF1+ cells
# in >= 5 samples) -- no local cell/donor threshold, so this comparison agrees with
# the nDEG / nPathway PHF1-filtered figures. See the per-subtype section below.

# Satterthwaite df: without this, emmeans prints df = Inf above ~3000 obs and the
# contrast CI is wrong. n here is ~63k.
emmeans::emm_options(lmer.df = "satterthwaite", lmerTest.limit = 1e5,
                     pbkrtest.limit = 1e5)

wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE)

# -------------------------------------------------------------------
# Load. sce_phf1_dist.qs is the all-cell SCE carrying the canonical
# dist_to_phf1_um (written by add_dist_to_full_sce.R, which reuses the values
# from celltype_sce_neighbours/ rather than recomputing them).
# -------------------------------------------------------------------
message("Reading sce_phf1_dist.qs ...")
sce <- qs::qread("sce_phf1_dist.qs")
cd  <- as.data.frame(colData(sce))
cd$cell_id <- colnames(sce)
rm(sce); gc()

req <- c("cell_id", "sample_id", "celltype", "PHF1", "x_slide_mm", "y_slide_mm",
         "fov", "x_FOV_px", "y_FOV_px", "dist_to_phf1_um", "Braak", "Sex", "Age", "PMI")
miss <- setdiff(req, colnames(cd))
if (length(miss)) stop("Missing colData columns: ", paste(miss, collapse = ", "))

cd$sample_id  <- as.character(cd$sample_id)
cd$celltype   <- as.character(cd$celltype)
cd$PHF1       <- as.character(cd$PHF1)
cd$x_slide_mm <- as.numeric(cd$x_slide_mm)
cd$y_slide_mm <- as.numeric(cd$y_slide_mm)

# Cortical depth normalised over ALL cells in a sample (not the neuron subset),
# so depth_rel means the same thing in every script that uses it.
cd <- add_cortical_depth(cd)

# -------------------------------------------------------------------
# Provenance assertion: the dist_to_phf1_um travelling with this object must be
# byte-identical to the canonical column in celltype_sce_neighbours/. Same check
# generate_phf1_null_labels.r makes, for the same reason -- so a new metric
# can never be silently pinned to stale distances.
# -------------------------------------------------------------------
neigh <- list.files("celltype_sce_neighbours", pattern = "_sce_neighbours\\.qs$",
                    full.names = TRUE)
if (length(neigh)) {
  chk_f <- neigh[which.min(file.size(neigh))]   # smallest file: cheapest valid check
  s_chk <- qs::qread(chk_f)
  chk   <- s_chk$dist_to_phf1_um; names(chk) <- colnames(s_chk)
  shared <- intersect(names(chk), cd$cell_id)
  cmp <- all.equal(unname(chk[shared]),
                   unname(setNames(cd$dist_to_phf1_um, cd$cell_id)[shared]))
  if (!isTRUE(cmp)) stop("dist_to_phf1_um does not match canonical ", basename(chk_f), ": ", cmp)
  message(sprintf("Provenance OK: dist_to_phf1_um matches %s for %d cells.",
                  basename(chk_f), length(shared)))
  rm(s_chk); gc()
} else {
  warning("celltype_sce_neighbours/ not found -- provenance check skipped.")
}

# -------------------------------------------------------------------
# Metric + technical covariate + construct-validity diagnostic
# -------------------------------------------------------------------
message("Computing nn", K, " spacing ...")
nn <- compute_nn3_um(cd, sample_col = "sample_id", neuron_pool = NEURON_POOL_10, k = K)

message("Computing FOV edge distance ...")
ed <- fov_edge_dist_um(cd, cap_um = EDGE_CAP_UM, raster_um = RASTER_UM)

message("Computing local neuron count (diagnostic only) ...")
lc <- local_neuron_count(cd, sample_col = "sample_id", neuron_pool = NEURON_POOL_10)

md <- cd %>%
  select(cell_id, sample_id, celltype, PHF1, x_slide_mm, y_slide_mm, fov,
         dist_to_phf1_um, depth_mm, depth_rel, Braak, Sex, Age, PMI) %>%
  inner_join(nn, by = "cell_id") %>%
  left_join(ed, by = "cell_id") %>%
  left_join(lc, by = "cell_id") %>%
  filter(!is.na(nn3_um)) %>%
  mutate(
    phf1_pos = PHF1 == "TRUE",
    grp      = factor(ifelse(phf1_pos, "PHF1+", "PHF1-"), levels = c("PHF1-", "PHF1+")),
    celltype = factor(celltype, levels = c(neuron_order, "Unassigned Neuron")),
    Braak    = factor(as.character(Braak), levels = braak_levels),
    Sex      = factor(as.character(Sex)),
    log_nn3  = log(nn3_um),
    Age_s    = as.numeric(scale(as.numeric(Age))),
    PMI_s    = as.numeric(scale(as.numeric(PMI)))
  )
stopifnot(!any(is.na(md$celltype)), nrow(md) > 0)

# Shared cache for Analysis B and the null replica.
wr(md %>% select(cell_id, sample_id, celltype, PHF1, phf1_pos, nn3_um, nn3rd_um,
                 n_coincident, n_pool_sample, edge_dist_um, edge_dist_um_raw,
                 n_local, depth_mm, depth_rel, dist_to_phf1_um,
                 Braak, Sex, Age, PMI, x_slide_mm, y_slide_mm, fov),
   "source_data_nn3_cells.tsv")

# Diagnostics that decide whether the metric and the edge covariate are sane.
rho_local <- cor(md$nn3_um, md$n_local,          method = "spearman", use = "complete.obs")
rho_edge  <- cor(md$nn3_um, md$edge_dist_um_raw, method = "spearman", use = "complete.obs")
rho_dist  <- cor(md$nn3_um, md$dist_to_phf1_um,  method = "spearman", use = "complete.obs")
q999      <- quantile(md$nn3_um, 0.999, na.rm = TRUE)

# -------------------------------------------------------------------
# Cell-level models
# -------------------------------------------------------------------
fit <- function(f, d = md) tryCatch(lmerTest::lmer(f, data = d, REML = TRUE),
                                    error = function(e) { message("fit failed: ", conditionMessage(e)); NULL })

m_head  <- fit(log_nn3 ~ grp + celltype + edge_dist_um + (1 | sample_id))
m_depth <- fit(log_nn3 ~ grp + celltype + edge_dist_um + depth_rel + (1 | sample_id))
m_min   <- fit(log_nn3 ~ grp + (1 | sample_id))
m_braak <- fit(log_nn3 ~ grp * Braak + celltype + edge_dist_um + (1 | sample_id))
m_3rd   <- fit(log(nn3rd_um) ~ grp + celltype + edge_dist_um + (1 | sample_id))

emm_of <- function(m) if (is.null(m)) NULL else emmeans::emmeans(m, ~ grp)
emm_h  <- emm_of(m_head)
emm_df  <- as.data.frame(summary(emm_h, infer = TRUE))
pair_df <- as.data.frame(summary(pairs(emm_h), infer = TRUE))   # PHF1- minus PHF1+
omni    <- anova(m_head)

# Re-orient every contrast to PHF1+ minus PHF1- so a positive number always means
# "PHF1+ neurons are more widely spaced". pairs() gives (PHF1-) - (PHF1+).
grab_diff <- function(m) {
  if (is.null(m)) return(rep(NA_real_, 4))
  p <- as.data.frame(summary(pairs(emmeans::emmeans(m, ~ grp)), infer = TRUE))
  c(-p$estimate[1], -p$upper.CL[1], -p$lower.CL[1], p$p.value[1])
}
d_head  <- grab_diff(m_head);  d_depth <- grab_diff(m_depth)
d_min   <- grab_diff(m_min);   d_3rd   <- grab_diff(m_3rd)

# -------------------------------------------------------------------
# Donor-level paired test (n = 9 donors)
# -------------------------------------------------------------------
donor <- md %>%
  group_by(sample_id, Braak, grp) %>%
  summarise(mean_log_nn3 = mean(log_nn3), geo_mean_nn3_um = exp(mean(log_nn3)),
            median_nn3_um = median(nn3_um), n_cells = dplyr::n(), .groups = "drop")

dw <- donor %>%
  select(sample_id, Braak, grp, mean_log_nn3) %>%
  pivot_wider(names_from = grp, values_from = mean_log_nn3) %>%
  filter(!is.na(`PHF1-`), !is.na(`PHF1+`)) %>%
  mutate(diff_log = `PHF1+` - `PHF1-`)

paired_w <- tryCatch(wilcox.test(dw$`PHF1+`, dw$`PHF1-`, paired = TRUE), error = function(e) NULL)
paired_t <- tryCatch(t.test(dw$`PHF1+`, dw$`PHF1-`, paired = TRUE),      error = function(e) NULL)
dz       <- mean(dw$diff_log) / sd(dw$diff_log)
n_same   <- sum(sign(dw$diff_log) == sign(mean(dw$diff_log)))

# -------------------------------------------------------------------
# DONOR WEIGHTING, and two further donor-level statistics.
#
# The unweighted paired test above gives all 9 donors equal weight, but their
# PHF1+ counts span 4 to 118. A donor contributing 4 PHF1+ neurons carries as
# much weight as one contributing 118. The cell-level LMM instead weights donors
# by information. Both of the following are therefore reported:
#   * inverse-variance weighted donor-level estimate (fixed-effect meta-analysis
#     over donors), the donor-level analogue of what the LMM actually estimates;
#   * the unweighted paired test restricted to donors with >= MIN_POS_DONOR PHF1+
#     neurons, i.e. those whose per-donor mean is estimated at all precisely.
# -------------------------------------------------------------------
MIN_POS_DONOR <- 20L

dstat <- md %>%
  group_by(sample_id, Braak) %>%
  summarise(n_pos = sum(phf1_pos), n_neg = sum(!phf1_pos),
            m_pos = mean(log_nn3[phf1_pos]),  m_neg = mean(log_nn3[!phf1_pos]),
            v_pos = var(log_nn3[phf1_pos]),   v_neg = var(log_nn3[!phf1_pos]),
            .groups = "drop") %>%
  filter(n_pos >= 2, n_neg >= 2) %>%
  mutate(diff_log = m_pos - m_neg,
         se_diff  = sqrt(v_pos / n_pos + v_neg / n_neg),
         w        = 1 / se_diff^2,
         pct_diff = 100 * (exp(diff_log) - 1))

iv     <- dstat %>% filter(is.finite(w), w > 0)
iv_est <- sum(iv$w * iv$diff_log) / sum(iv$w)
iv_se  <- sqrt(1 / sum(iv$w))
iv_ci  <- iv_est + c(-1.96, 1.96) * iv_se
iv_p   <- 2 * stats::pnorm(-abs(iv_est / iv_se))

dw_r     <- dw %>% semi_join(dstat %>% filter(n_pos >= MIN_POS_DONOR), by = "sample_id")
paired_w_r <- tryCatch(wilcox.test(dw_r$`PHF1+`, dw_r$`PHF1-`, paired = TRUE),
                       error = function(e) NULL)
paired_t_r <- tryCatch(t.test(dw_r$`PHF1+`, dw_r$`PHF1-`, paired = TRUE),
                       error = function(e) NULL)
dz_r     <- if (nrow(dw_r) > 1) mean(dw_r$diff_log) / sd(dw_r$diff_log) else NA_real_
n_same_r <- if (nrow(dw_r) > 0) sum(sign(dw_r$diff_log) == sign(mean(dw_r$diff_log))) else NA_integer_
wr(dstat, "stats_nn3_phf1_vs_neg_donor_contrasts.tsv")

# -------------------------------------------------------------------
# Effect sizes. The standard is an unstandardised effect WITH its 95% CI plus
# a standardised measure. On a log outcome, exp(diff) - 1 is the % difference.
# -------------------------------------------------------------------
vc_tot  <- { vc <- as.data.frame(lme4::VarCorr(m_head)); sqrt(sum(vc$vcov[is.na(vc$var2)])) }
d_model <- d_head[1] / vc_tot
sd_pool <- sd(md$log_nn3)
d_cell  <- d_head[1] / sd_pool
emm_neg <- emm_df$emmean[emm_df$grp == "PHF1-"]

es <- tibble(
  measure              = "nn3_um",
  n_cells              = nrow(md),
  n_phf1_pos           = sum(md$phf1_pos),
  n_donors             = dplyr::n_distinct(md$sample_id),
  # cell-level LMM, headline
  diff_log             = d_head[1],
  diff_log_CI.L        = d_head[2],
  diff_log_CI.R        = d_head[3],
  pct_diff             = 100 * (exp(d_head[1]) - 1),
  pct_diff_CI.L        = 100 * (exp(d_head[2]) - 1),
  pct_diff_CI.R        = 100 * (exp(d_head[3]) - 1),
  p_lmm                = d_head[4],
  geo_mean_phf1_neg_um = exp(emm_neg),
  um_diff              = exp(emm_neg) * (exp(d_head[1]) - 1),
  um_diff_CI.L         = exp(emm_neg) * (exp(d_head[2]) - 1),
  um_diff_CI.R         = exp(emm_neg) * (exp(d_head[3]) - 1),
  cohens_d_model       = d_model,
  cohens_d_cell        = d_cell,
  # donor level, paired
  dz_donor             = dz,
  diff_log_donor       = mean(dw$diff_log),
  diff_log_donor_CI.L  = if (!is.null(paired_t)) paired_t$conf.int[1] else NA_real_,
  diff_log_donor_CI.R  = if (!is.null(paired_t)) paired_t$conf.int[2] else NA_real_,
  pct_diff_donor       = 100 * (exp(mean(dw$diff_log)) - 1),
  p_wilcox_donor       = if (!is.null(paired_w)) paired_w$p.value else NA_real_,
  n_donors_paired      = nrow(dw),
  n_donors_same_dir    = n_same,
  # donor level, inverse-variance weighted (the donor-level analogue of the LMM)
  pct_diff_donor_iv     = 100 * (exp(iv_est) - 1),
  pct_diff_donor_iv_L   = 100 * (exp(iv_ci[1]) - 1),
  pct_diff_donor_iv_R   = 100 * (exp(iv_ci[2]) - 1),
  p_donor_iv            = iv_p,
  # donor level, unweighted but restricted to adequately-sampled donors
  pct_diff_donor_restr  = if (nrow(dw_r)) 100 * (exp(mean(dw_r$diff_log)) - 1) else NA_real_,
  dz_donor_restr        = dz_r,
  p_wilcox_donor_restr  = if (!is.null(paired_w_r)) paired_w_r$p.value else NA_real_,
  n_donors_restr        = nrow(dw_r),
  n_donors_restr_same_dir = n_same_r,
  # sensitivity fits
  pct_diff_depth_adj   = 100 * (exp(d_depth[1]) - 1), p_depth_adj = d_depth[4],
  pct_diff_unadj       = 100 * (exp(d_min[1])   - 1), p_unadj     = d_min[4],
  pct_diff_nn3rd       = 100 * (exp(d_3rd[1])   - 1), p_nn3rd     = d_3rd[4],
  # diagnostics
  rho_nn3_local_count  = rho_local,
  rho_nn3_edge_dist    = rho_edge,
  rho_nn3_dist_to_phf1 = rho_dist
)
wr(es, "stats_nn3_phf1_vs_neg_effectsize.tsv")
wr(pair_df, "stats_nn3_phf1_vs_neg_pairwise.tsv")
wr(as.data.frame(emm_df), "stats_nn3_phf1_vs_neg_emmeans.tsv")

# -------------------------------------------------------------------
# Per-subtype fits.
#
# THE ELIGIBILITY CRITERION IS THE ONLY GATE. A subtype is fitted if and only if
# eligible_neuron_subtypes() returns it -- >= 6 PHF1+ cells in >= 5 samples, read
# from the same count table the nDEG / nPathway PHF1-filtered figures use. That
# keeps this categorical comparison byte-consistent with the rest of the project
# rather than inventing a local threshold.
#
# On this cohort that admits exactly two subtypes:
#   Exc-IT-L2-3-CBLN2-HOPX      185 PHF1+
#   Exc-IT-L3-5-CHGA-IL1RAPL2    70 PHF1+
#
# Not fitted: "Unassigned Neuron" (89 PHF1+). It is absent from the helper's
# input table, which covers neuron_order only, and the helper's output defines
# eligibility here; an ambiguous label is not used for a subtype-level contrast.
# Its descriptive numbers are still reported below.
#
# The POOLED model above is NOT restricted: it keeps all 10 neuron-pool labels.
# It adjusts for celltype, so sparse subtypes shift it very
# little -- restricting it gives +6.45% [+1.02, +12.16] against +6.84%
# [+2.40, +11.47] pooled, i.e. the same magnitude with a wider interval.
# -------------------------------------------------------------------
elig_canon <- eligible_neuron_subtypes(getwd(), neuron_levels = neuron_order)
if (!length(elig_canon))
  stop("eligible_neuron_subtypes() returned nothing; cannot run the per-subtype comparison.")

ct_tab <- md %>% group_by(celltype) %>%
  summarise(n_cells = dplyr::n(), n_phf1_pos = sum(phf1_pos),
            n_donors_both = dplyr::n_distinct(sample_id[phf1_pos]), .groups = "drop")

fit_ct <- function(ct) {
  d <- md %>% filter(celltype == ct)
  row <- tibble(celltype = as.character(ct), n_cells = nrow(d),
                n_phf1_pos = sum(d$phf1_pos),
                n_donors_with_phf1 = dplyr::n_distinct(d$sample_id[d$phf1_pos]),
                geo_mean_neg_um = exp(mean(d$log_nn3[!d$phf1_pos])),
                geo_mean_pos_um = if (any(d$phf1_pos)) exp(mean(d$log_nn3[d$phf1_pos])) else NA_real_,
                diff_log = NA_real_, CI.L = NA_real_, CI.R = NA_real_, pval = NA_real_,
                fitted = FALSE)
  # Eligibility is the gate. Everything else is descriptive-only.
  if (!as.character(ct) %in% elig_canon) return(row)
  m <- fit(log_nn3 ~ grp + edge_dist_um + (1 | sample_id), d)
  if (is.null(m)) return(row)
  a <- grab_diff(m)
  row$diff_log <- a[1]; row$CI.L <- a[2]; row$CI.R <- a[3]; row$pval <- a[4]
  row$fitted <- TRUE
  row
}

ct_res <- bind_rows(lapply(levels(md$celltype), fit_ct)) %>%
  mutate(canonical_eligible = celltype %in% elig_canon,
         pct_diff   = 100 * (exp(diff_log) - 1),
         pct_diff_L = 100 * (exp(CI.L) - 1),
         pct_diff_R = 100 * (exp(CI.R) - 1),
         padj       = NA_real_)
# BH over the FITTED subtypes only. p.adjust()'s default n = length(p) counts NAs,
# which would over-correct here, so subset rather than pass NAs through.
ct_res$padj[ct_res$fitted] <- p.adjust(ct_res$pval[ct_res$fitted], method = "BH")
wr(ct_res, "source_data_nn3_phf1_vs_neg_by_celltype.tsv")

# -------------------------------------------------------------------
# Figures
# -------------------------------------------------------------------
# Two-line y label: the single-line version is wider than a 5 cm panel is tall and
# gets clipped. mu via plotmath so it survives the default pdf device.
nn3_lab     <- expression(atop("Mean distance to 3", "nearest neurons (" * mu * "m)"))
dist_breaks <- c(5, 10, 20, 30, 50, 100, 200, 400)
# Display window, not a filter: the boxplot hides its outliers, so let the axis
# match. Covers >99.9% of cells (99.9th percentile of nn3_um is ~97 um).
Y_VIEW      <- c(8, 130)

sd_cells <- md %>% select(cell_id, sample_id, Braak, celltype, grp, nn3_um, nn3rd_um,
                          edge_dist_um, depth_rel, n_local)
wr(sd_cells,  "source_data_nn3_phf1_vs_neg.tsv")
wr(donor,     "source_data_nn3_phf1_vs_neg_donor_means.tsv")

dm <- donor %>% mutate(Braak = factor(as.character(Braak), levels = braak_levels))

p1 <- ggplot(md, aes(grp, nn3_um)) +
  geom_boxplot(aes(fill = grp), outlier.shape = NA, linewidth = 0.3, width = 0.55,
               alpha = 0.5, colour = "grey30") +
  geom_line(data = dm, aes(grp, geo_mean_nn3_um, group = sample_id),
            linewidth = 0.25, colour = "grey55") +
  geom_point(data = dm, aes(grp, geo_mean_nn3_um, fill = Braak),
             shape = 21, size = 1.5, stroke = 0.25, colour = "grey20") +
  scale_fill_manual(values = c("PHF1-" = "grey75", "PHF1+" = "#BD0026", braak_palette),
                    breaks = braak_levels, name = "Braak") +
  scale_y_log10(breaks = dist_breaks, labels = dist_breaks) +
  coord_cartesian(ylim = Y_VIEW) +
  labs(x = NULL, y = nn3_lab) +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "right", plot.margin = margin(4, 4, 4, 6, unit = "pt"))
ggsave(file.path(out_dir, "plot_nn3_phf1_vs_neg.pdf"), p1,
       width = 6.6, height = 5.2, units = "cm", device = "pdf")

fp <- ct_res %>% filter(fitted) %>%
  mutate(celltype = factor(celltype, levels = rev(levels(md$celltype))),
         sig = ifelse(!is.na(padj) & padj < FDR, "padj < 0.05", "n.s."),
         lab = sprintf("n = %d", n_phf1_pos),
         lab_x = max(pct_diff_R, na.rm = TRUE))
if (nrow(fp)) {
  # Only eligible subtypes are drawn, so there is no eligible-vs-sparse shape
  # distinction to make any more: every point on this panel meets the >= 6 PHF1+
  # in >= 5 samples criterion. Subtypes that miss it are descriptive-only and live
  # in stats_nn3_phf1_vs_neg_by_celltype.txt.
  p2 <- ggplot(fp, aes(pct_diff, celltype, colour = sig)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
    geom_errorbar(aes(xmin = pct_diff_L, xmax = pct_diff_R), orientation = "y",
                  width = 0, linewidth = 0.4) +
    geom_point(size = 1.7) +
    geom_text(aes(x = lab_x, label = lab), hjust = 0, size = 1.9, colour = "grey35",
              nudge_x = 2) +
    scale_colour_manual(values = c("padj < 0.05" = "#BD0026", "n.s." = "grey65"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    labs(x = "Spacing in PHF1+ vs PHF1- neurons (%)", y = NULL) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "top",
          legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(4, 8, 4, 4, unit = "pt"))
  ggsave(file.path(out_dir, "plot_nn3_phf1_vs_neg_by_celltype.pdf"), p2,
         width = 9.4, height = 4.0, units = "cm", device = "pdf")
}

# Donor-level forest. Each donor's own PHF1+/PHF1- contrast is drawn with its own
# 95% CI and sized by PHF1+ count, over the inverse-variance pooled estimate.
dfor <- dstat %>%
  mutate(pct_L = 100 * (exp(diff_log - 1.96 * se_diff) - 1),
         pct_R = 100 * (exp(diff_log + 1.96 * se_diff) - 1),
         sample_id = factor(sample_id, levels = sample_id[order(diff_log)]))
p3 <- ggplot(dfor, aes(pct_diff, sample_id)) +
  annotate("rect", xmin = 100 * (exp(iv_ci[1]) - 1), xmax = 100 * (exp(iv_ci[2]) - 1),
           ymin = -Inf, ymax = Inf, fill = "#BD0026", alpha = 0.10) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
  geom_vline(xintercept = 100 * (exp(iv_est) - 1), linewidth = 0.4, colour = "#BD0026") +
  geom_errorbar(aes(xmin = pct_L, xmax = pct_R), orientation = "y", width = 0,
                linewidth = 0.35, colour = "grey40") +
  geom_point(aes(size = n_pos, fill = Braak), shape = 21, stroke = 0.25, colour = "grey20") +
  scale_fill_manual(values = braak_palette, name = "Braak", drop = FALSE) +
  scale_size_continuous(range = c(0.9, 2.8), name = "PHF1+ n") +
  labs(x = "Spacing in PHF1+ vs PHF1- neurons (%)", y = NULL) +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(4, 6, 4, 4, unit = "pt"))
ggsave(file.path(out_dir, "plot_nn3_phf1_vs_neg_donor_forest.pdf"), p3,
       width = 8.4, height = 5.4, units = "cm", device = "pdf")
wr(dfor %>% select(sample_id, Braak, n_pos, n_neg, diff_log, se_diff, pct_diff, pct_L, pct_R, w),
   "source_data_nn3_phf1_vs_neg_donor_forest.tsv")

# -------------------------------------------------------------------
# Stats logs
# -------------------------------------------------------------------
sink(file.path(out_dir, "stats_nn3_phf1_vs_neg.txt"))
cat("3-nearest-neuronal-neighbour spacing: PHF1+ vs PHF1- neurons (CosMx, 9 donors)\n")
cat("=============================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Source object: sce_phf1_dist.qs (all cells; canonical dist_to_phf1_um)\n")
cat("Metric:        nn3_um = mean Euclidean distance to the", K, "nearest other neurons, um.\n")
cat("Neighbour pool:", length(NEURON_POOL_10), "neuron labels (neuron_order + 'Unassigned Neuron'),\n")
cat("               identical to the PHF1+ distance-source set in run_label_phf1_neighbours.sh.\n\n")
cat("REFERENCE. Zwang et al. 2024 Cell Rep. (Hyman lab, with Bennett; PMC11441076) imaged\n")
cat("rTg4510 cortex longitudinally: tangle-FREE neurons died at >3x the rate of tangle-bearing\n")
cat("neurons and became more distant from their neighbours before dying (3-NN distance\n")
cat("33.9 +/- 1.9 um in dying vs 24.0 +/- 0.5 um in persisting neurons, p = 0.02, mixed model).\n")
cat(sprintf("Our pooled median nn3 is %.1f um, i.e. the metric lands in the same range despite\n",
            median(md$nn3_um)))
cat("being 2D. That is a level check only, not a replication: see NOTES.\n\n")
cat("MODELS (outcome log(nn3_um); contrasts re-oriented to PHF1+ minus PHF1-, so a\n")
cat("POSITIVE value means PHF1+ neurons are MORE WIDELY SPACED)\n")
cat("  HEADLINE cell  : log(nn3_um) ~ grp + celltype + edge_dist_um + (1 | sample_id)\n")
cat("  HEADLINE donor : per-donor mean log(nn3_um), PHF1+ vs PHF1-, paired Wilcoxon (n = 9),\n")
cat("                   reported three ways -- unweighted, inverse-variance weighted, and\n")
cat(sprintf("                   unweighted restricted to donors with >= %d PHF1+ neurons. See the\n",
            MIN_POS_DONOR))
cat("                   effect-size block for why all three are needed.\n")
cat("  depth-adjusted : + depth_rel                       (reported, not headline)\n")
cat("  unadjusted     : ~ grp + (1 | sample_id)           (reported)\n")
cat("  Braak          : ~ grp * Braak + celltype + edge_dist_um + (1 | sample_id)\n")
cat("  nn3rd sensitivity: outcome log(nn3rd_um), distance to the 3rd nearest neuron\n\n")
cat("CELLS\n")
cat(sprintf("  pool neurons: %d | PHF1+: %d | PHF1-: %d | donors: %d\n",
            nrow(md), sum(md$phf1_pos), sum(!md$phf1_pos),
            dplyr::n_distinct(md$sample_id)))
cat("  nn3_um:\n"); print(summary(md$nn3_um))
cat("  per-donor pool size and median nn3:\n")
print(as.data.frame(md %>% group_by(sample_id, Braak) %>%
        summarise(n_pool = dplyr::n(), n_phf1 = sum(phf1_pos),
                  median_nn3_um = round(median(nn3_um), 1), .groups = "drop")),
      row.names = FALSE)
cat("  PHF1+ neurons per subtype:\n")
print(as.data.frame(ct_tab), row.names = FALSE)
cat("\n=== CONSTRUCT VALIDITY / CONFOUND DIAGNOSTICS ===\n")
cat(sprintf("  rho(nn3_um, neurons within 50 um)   = %+.3f   (must be strongly NEGATIVE)\n", rho_local))
cat(sprintf("  rho(nn3_um, edge_dist_um_raw)       = %+.3f   (FOV edge censoring: mild)\n", rho_edge))
cat(sprintf("  rho(nn3_um, dist_to_phf1_um)        = %+.3f   <- THE GEOMETRIC COUPLING\n", rho_dist))
cat("  The last one matters for Analysis B, not for this one: dist_to_phf1_um is itself a\n")
cat("  nearest-neighbour distance to a neuron subset, so it is coupled to neuron sparseness by\n")
cat("  construction. This comparison does not use distance at all.\n")
cat(sprintf("  edge_dist_um cap = %d um vs the 99.9th percentile of nn3_um = %.1f um: the cap sits\n",
            EDGE_CAP_UM, q999))
cat("  well above the metric's upper tail, so censoring beyond it cannot affect nn3_um.\n")
cat(sprintf("  neurons within 100 um of an unimaged region: %.1f%%\n",
            100 * mean(md$edge_dist_um_raw < 100)))
cat(sprintf("  neurons with a coincident (identical-centroid) pool neuron: %d\n",
            sum(md$n_coincident > 0, na.rm = TRUE)))
cat("\n=== HEADLINE 1: cell-level LMM ===\n")
print(summary(m_head))
cat("\n--- Omnibus (Satterthwaite F) ---\n"); print(omni)
cat("\n--- EMMs (log scale) ---\n"); print(emm_df, row.names = FALSE)
cat("\n--- Contrast, PHF1- minus PHF1+ as emmeans reports it ---\n")
print(pair_df, row.names = FALSE)
cat("\n=== HEADLINE 2: donor-level paired (n =", nrow(dw), "donors) ===\n")
if (!is.null(paired_w)) print(paired_w)
if (!is.null(paired_t)) print(paired_t)
cat("\nPer-donor mean log(nn3_um):\n")
print(as.data.frame(dw %>% mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
cat("\nPer-donor contrast with its own precision (the input to the weighted estimate):\n")
print(as.data.frame(dstat %>% select(sample_id, Braak, n_pos, n_neg, diff_log, se_diff,
                                     pct_diff, w) %>%
                    mutate(across(where(is.numeric), ~signif(.x, 4)))), row.names = FALSE)
cat("\n=== Effect sizes ===\n")
cat(sprintf("  CELL-LEVEL LMM:            %+.2f%%  95%% CI [%+.2f%%, %+.2f%%]  p = %.3g\n",
            es$pct_diff, es$pct_diff_CI.L, es$pct_diff_CI.R, es$p_lmm))
cat(sprintf("    in microns:              %+.2f um 95%% CI [%+.2f, %+.2f] on a PHF1- geometric mean of %.1f um\n",
            es$um_diff, es$um_diff_CI.L, es$um_diff_CI.R, es$geo_mean_phf1_neg_um))
cat(sprintf("    Cohen's d (model-based)  = %+.3f ; (cell pooled SD) = %+.3f\n",
            es$cohens_d_model, es$cohens_d_cell))
cat(sprintf("  DONOR, unweighted (n = %d): %+.2f%% , paired dz = %+.3f , %d of %d donors same direction, p = %.3g\n",
            es$n_donors_paired, es$pct_diff_donor, es$dz_donor, es$n_donors_same_dir,
            es$n_donors_paired, es$p_wilcox_donor))
cat(sprintf("  DONOR, inv-var weighted:   %+.2f%%  95%% CI [%+.2f%%, %+.2f%%]  p = %.3g\n",
            es$pct_diff_donor_iv, es$pct_diff_donor_iv_L, es$pct_diff_donor_iv_R,
            es$p_donor_iv))
cat(sprintf("  DONOR, >= %d PHF1+ (n = %d): %+.2f%% , paired dz = %+.3f , %d of %d same direction, p = %.3g\n",
            MIN_POS_DONOR, es$n_donors_restr, es$pct_diff_donor_restr, es$dz_donor_restr,
            es$n_donors_restr_same_dir, es$n_donors_restr, es$p_wilcox_donor_restr))
cat("\nWEIGHTING. Per-donor PHF1+ counts span 4 to 118. The unweighted paired test gives a donor\n")
cat("with 4 PHF1+ neurons the same weight as one with 118; the LMM and the inverse-variance\n")
cat("estimate weight by information.\n")
cat(sprintf("Celltype-adjusted vs unadjusted cell-level fit: %+.2f%% vs %+.2f%%.\n",
            es$pct_diff, es$pct_diff_unadj))
cat("PHF1 status is a WITHIN-donor predictor, so a random intercept on sample_id is the correct\n")
cat("specification; the supplementary label-permutation null in R/null_replica_nn3.R checks its\n")
cat("calibration.\n")
cat("\n--- Sensitivity fits (all as % difference, PHF1+ vs PHF1-) ---\n")
cat(sprintf("  depth-adjusted        : %+.2f%%  p = %.3g\n", es$pct_diff_depth_adj, es$p_depth_adj))
cat(sprintf("  unadjusted            : %+.2f%%  p = %.3g\n", es$pct_diff_unadj, es$p_unadj))
cat(sprintf("  3rd-nearest (not mean): %+.2f%%  p = %.3g\n", es$pct_diff_nn3rd, es$p_nn3rd))
cat("\n=== Braak interaction (grp x Braak) ===\n")
if (!is.null(m_braak)) print(anova(m_braak))
cat("\nNOTES:\n")
cat(" - 2D, single section. A 2D nn3 over-estimates true 3D spacing because out-of-plane\n")
cat("   neighbours are invisible. Zwang et al. measured it in 3D z-stacks, so absolute values\n")
cat("   are not directly comparable even though they happen to be close.\n")
cat(" - CROSS-SECTIONAL. Zwang et al.'s finding is that neurons BECAME more distant before\n")
cat("   dying, observed longitudinally in the same cells. Nothing here can address change over\n")
cat("   time; only the standing spatial association is testable.\n")
cat(" - Neurons are transcriptomically typed, not pan-neuronally labelled. 40,681 'Unassigned'\n")
cat("   cells are excluded from the pool and certainly contain neurons, so nn3_um is biased\n")
cat("   upward. The bias should be near-identical across PHF1 status, so it inflates the level\n")
cat("   rather than the contrast.\n")
cat(" - 'Unassigned Neuron' IS in the pool (15,444 cells, ~25% of it) because that matches how\n")
cat("   dist_to_phf1_um was defined.\n")
cat(" - FOV edge censoring is adjusted for, not excluded; see the\n")
cat("   diagnostics above for how mild it is.\n")
cat(" - PHF1+ neurons are strongly layer- and subtype-biased, so the celltype term in the\n")
cat("   headline model is load-bearing. Compare it against the unadjusted fit above.\n")
cat(" - n = 3 donors per Braak stage; the interaction is underpowered and descriptive.\n")
cat(" - The supplementary permutation null is R/null_replica_nn3.R and is deliberately NOT\n")
cat("   part of this primary analysis.\n")
cat("\n=== sessionInfo() ===\n"); print(sessionInfo())
sink()

sink(file.path(out_dir, "stats_nn3_phf1_vs_neg_by_celltype.txt"))
cat("3-NN neuronal spacing, PHF1+ vs PHF1-, PER NEURONAL SUBTYPE\n")
cat("===========================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Model per subtype: log(nn3_um) ~ grp + edge_dist_um + (1 | sample_id)\n")
cat("Contrast oriented PHF1+ minus PHF1- and reported as a % difference.\n\n")
cat("SUBTYPE SELECTION. A subtype is fitted if and only if it passes the project-wide\n")
cat("eligibility criterion -- >= 6 PHF1+ cells in >= 5 samples, via\n")
cat("R/eligible_neuron_subtypes.R, reading the same count table the nDEG / nPathway\n")
cat("PHF1-filtered figures use. There is NO additional local threshold, so this comparison\n")
cat("selects subtypes identically to the rest of the project. BH is applied across the\n")
cat("fitted subtypes only. On this cohort the criterion admits:\n  ")
cat(paste(elig_canon, collapse = "\n  "), "\n\n")
cat("Every other subtype is DESCRIPTIVE-ONLY below (fitted = FALSE): counts and geometric\n")
cat("means are given and no contrast is estimated.\n\n")
cat("Two cases worth naming explicitly:\n")
cat(" - 'Unassigned Neuron' (89 PHF1+ cells) is NOT fitted: it is absent from the helper's input\n")
cat("   table (which covers neuron_order only), and the helper's output defines eligibility\n")
cat("   here; an ambiguous label is not used for a subtype-level contrast.\n")
cat(" - Exc-ET-L5-SPON1-FGD4 has 0 PHF1+ cells and cannot enter this comparison at all.\n\n")
cat("The POOLED model in stats_nn3_phf1_vs_neg.txt is deliberately NOT restricted -- it keeps\n")
cat("all 10 neuron-pool labels and adjusts for celltype. For reference, restricting it to the\n")
cat("eligible subtypes gives +6.45% [+1.02, +12.16] against +6.84% [+2.40, +11.47] pooled:\n")
cat("the same magnitude with a wider interval.\n\n")
cat("=== RESULTS ===\n")
print(as.data.frame(ct_res %>%
        select(celltype, n_cells, n_phf1_pos, n_donors_with_phf1, canonical_eligible,
               geo_mean_neg_um, geo_mean_pos_um, pct_diff, pct_diff_L, pct_diff_R,
               pval, padj, fitted)),
      row.names = FALSE, digits = 3)
cat("\n=== sessionInfo() ===\n"); print(sessionInfo())
sink()

message(sprintf("Done. %d pool neurons (%d PHF1+), %d donors. PHF1+ vs PHF1- = %+.2f%% [%+.2f, %+.2f], p = %.3g; donor dz = %+.3f",
                nrow(md), sum(md$phf1_pos), dplyr::n_distinct(md$sample_id),
                es$pct_diff, es$pct_diff_CI.L, es$pct_diff_CI.R, es$p_lmm, es$dz_donor))
