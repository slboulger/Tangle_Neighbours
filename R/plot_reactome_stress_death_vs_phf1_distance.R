#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_reactome_stress_death_vs_phf1_distance.R
#
# Figure panels: 3C, 3D, S2A, S9A, S9B, S9C, S9D, S9E, S9F
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# plot_reactome_stress_death_vs_phf1_distance.R
#
# Curated Reactome DEATH-COMMITMENT and STRESS-ADAPTATION module scores over distance
# to the nearest PHF1+ (ptau) neuron, in one excitatory neuron subtype and the two
# reactive glial types. Direct sibling of plot_phf1_module_vs_phf1_distance_modelp.R
# and plot_modulescore_vs_phf1_distance_modelp.R: SAME scoring (AddModuleScore on SCT),
# SAME canonical log-distance LMM, SAME model-FDR solid/dashed encoding, SAME figure
# style and output triple.
#
# WHY. The existing module-score-over-distance figures test DATA-DERIVED signatures --
# the celltype's own PHF1 marker set and the Otero-Garcia AT8 sets. They can say a
# gradient exists but not what process it is. This script adds a HYPOTHESIS-DRIVEN
# layer: Reactome pathways in three mechanistic groups, asking whether cells near
# tangle-bearing neurons are committing to die, mounting a stress response, or
# clearing damaged cargo -- and whether neurons and glia give the same answer.
#
# THE MODULE LIST IS DERIVED BY RULE, NOT HAND-PICKED. It is the
# Reactome Programmed Cell Death / Autophagy / Cellular-responses-to-stimuli subtrees
# filtered by CATALOGUE_RULES in R/reactome_filtered_modules.R -- panel size 25-80,
# panel coverage >= 70%, no surviving pair above Jaccard 0.8, one set per redundancy
# cluster, programme roots excluded. That yields 20 modules in three groups:
#
#   pcd        (6) : regulated necrosis, intrinsic apoptosis, apoptotic execution
#                    phase, RIPK1-mediated necrosis, BH3-only activation, extrinsic
#                    caspase activation
#   autophagy  (4) : selective autophagy, aggrephagy, mitophagy, late endosomal
#                    microautophagy
#   stress    (10) : heat stress, NFE2L2 nuclear events, UPR, hypoxia, HSP90 chaperone
#                    cycle, cytoprotection by HMOX1, heme signalling, ROS
#                    detoxification, oncogene-induced senescence, HSF1 transactivation
#
# Selecting by rule makes the choice checkable. Consequences of the rules (see
# R/reactome_filtered_modules.R): the ubiquitin-proteasome pathway (R-HSA-983168) is
# filed by Reactome outside all three roots, so no rule over these subtrees reaches it;
# the three UPR arms collapse into the parent UPR term; and sets are chosen for
# size/coverage/non-redundancy rather than for mechanistic interest.
#
# Gene sets come from the CATALOGUE, built in-process by build_reactome_catalogue()
# from the cached, provenance-stamped Reactome hierarchy + Enrichr GMT. The script is
# therefore self-contained -- it does not read results/reactome_programme_catalogue/ --
# yet returns byte-identical gene lists to it, because it calls the same builder on the
# same panel. Sets are panel-resolved WITH HGNC alias resolution (PARK2 -> PRKN,
# IL8 -> CXCL8, H2AFX -> H2AX, UFD1L -> UFD1); they are deliberately not re-derived by
# a plain GMT intersect, which would drop PRKN from PINK1-PRKN Mediated Mitophagy.
# The overlap audit is written to stats_reactome_geneset_overlap.tsv.
#
# ONE FIGURE PER CELLTYPE x GROUP (6 in total). Every module of the group is overlaid
# on one axis as a model-free rolling mean + 95% CI, on the RAW AddModuleScore scale:
# z-scoring flattens the modules onto each other, and the raw scale is what makes the
# PHF1+ reference readable (the glia script's reasoning, kept here).
#
# PHF1+ REFERENCE (NEURON celltypes only -- gated on membership of neuron_order, not on
# whether any cell happens to be flagged PHF1+, so stray PHF1+ glia from segmentation
# spillover can never put a strip on a glial panel). PHF1+ neurons define distance 0 and
# are EXCLUDED from every regression and rolling mean. They appear as a mean + 95% CI
# point per module in a MARGINAL STRIP at the LEFT edge, i.e. at the zero-distance end
# where they belong, sharing the y-axis with the curves so the gradient can be read
# against the level in the tangle-bearing cells themselves.
#
# The strip's POINT SHAPE carries the PHF1+ vs PHF1- contrast for that module:
# FILLED = significant (BH across the group's modules), HOLLOW = n.s. Shape rather than
# colour, so it stays readable in greyscale and for colourblind readers, and so it does
# not collide with the module colour. The contrast is a separate LMM,
#     score ~ phf1_pos + nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)
# fitted on ALL cells of the celltype with NO distance cap: capping would select the
# PHF1-negative comparison group on proximity to PHF1+ cells -- the very gradient the
# rest of the figure is about -- and bias the contrast toward the null. It is reported
# with a donor-paired Cohen's dz and the count of donors moving the same way.
#
# NOTES (also written into each stats log):
#   * These Reactome sets are EXTERNAL -- unlike the PHF1 marker set they are not
#     selected from the PHF1+/PHF1- contrast, so the PHF1+ strip is not elevated by
#     construction.
#   * Significance is the model FDR (BH within a celltype x group).
#   * AddModuleScore values from different gene sets are not on the same absolute
#     scale. Compare SHAPES across modules within a panel, not heights.
#
# DISTANCE TRANSFORM is the canonical one: natural log, scaled by SD only,
# not centred, with dist_sd computed WITHIN each celltype:
#     dist_to_phf1_um_scaled <- log(dist_to_phf1_um) / sd(log(dist_to_phf1_um))
# dist_sd is computed once per celltype, so a module's dist_scaled is identical in both
# groups. The script stop()s on a non-positive distance rather than absorbing it with an
# offset -- distance is strictly > 0 for every modelled (PHF1-negative) cell.
#
# OUTPUT (per celltype <ct> x group <grp>), under --output_dir:
#   plot_reactome_<grp>_dist_<ct>.pdf                    -- rolling means + 95% CI (+ PHF1+ strip)
#   plot_reactome_<grp>_dist_<ct>_nostrip.pdf            -- same figure, PHF1+ strip dropped (NEURONS only)
#   source_data_reactome_<grp>_dist_<ct>_rollmean.tsv    -- exact drawn rows (200 pts x module)
#   source_data_reactome_<grp>_dist_<ct>_phf1ref.tsv     -- PHF1+ strip rows (NEURONS only)
#   stats_reactome_<grp>_dist_<ct>.txt                   -- formula, LMM table, effect sizes, sessionInfo
#   stats_reactome_<grp>_dist_<ct>_lmm.tsv               -- structured per-module LMM table
#   stats_reactome_<grp>_dist_<ct>_effectsize.tsv        -- estimate + 95% CI + Cohen's d
#   PHF1+ vs PHF1- comparison (NEURONS only):
#     stats_reactome_<grp>_dist_<ct>_phf1contrast.tsv           -- LMM estimate + CI + Cohen's d + dz
#                                                                  + cell-level and donor-level Wilcoxon
#     plot_reactome_<grp>_phf1contrast_forest_<ct>.pdf          -- estimate + 95% CI per module, donor points behind
#     plot_reactome_<grp>_phf1contrast_bydonor_<ct>.pdf         -- a line per donor, PHF1- to PHF1+, Braak-coloured
#     source_data_reactome_<grp>_phf1contrast_bydonor_<ct>.tsv  -- exact drawn rows
#   Written once:
#     stats_reactome_geneset_overlap.tsv   -- per module: Reactome id/name, panel overlap, scored
#     stats_reactome_geneset_members.tsv   -- every gene of every set with an on_panel flag
#     stats_reactome_lmm_coeffs.tsv        -- all celltype x group units pooled
#
# Use run_reactome_stress_death_distance.sh

##  ............................................................................
##  Packages + setup                                                        ####
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(argparse)
})

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")                    # fig_theme, celltype_palette, neuron_order, ...
source("R/reactome_module_genesets.R")    # REACTOME_GMT_DEFAULT
source("R/reactome_filtered_modules.R")   # reactome_filtered_modules(), CATALOGUE_RULES,
                                          # MODULE_PALETTE. Derives the module list by
                                          # rule from the Reactome subtrees; see that
                                          # file for what the rules are and cost.

set.seed(42)

##  ............................................................................
##  Arguments                                                               ####
parser <- ArgumentParser()
parser$add_argument("--seu", default = "seu_PHF1.rds",
  help = "Path to seu_PHF1.rds [default: seu_PHF1.rds]")
parser$add_argument("--gmt", default = REACTOME_GMT_DEFAULT,
  help = "Cached Enrichr Reactome GMT [default: genesets/enrichr_gmt/Reactome_Pathways_2024.gmt]")
parser$add_argument("--output_dir",
  default = "plots/reactome_stress_death_vs_phf1_distance_1000um",
  help = "Output directory")
parser$add_argument("--seed", type = "integer", default = 42,
  help = "Seed (AddModuleScore + provenance) [default: 42]")
parser$add_argument("--max_dist_um", type = "double", default = 1000,
  help = "Restrict modelled cells to dist_to_phf1_um <= this (um); 0 = no cap [default: 1000]")
parser$add_argument("--window_um", type = "double", default = 75,
  help = "Rolling-mean window WIDTH (um) [default: 75]")
parser$add_argument("--overwrite", default = "yes", choices = c("yes", "no"),
  help = paste("Clear --output_dir before writing, so a re-run after a script change",
               "cannot leave orphaned files from an older version behind.",
               "'no' keeps whatever",
               "is already there [default: yes]"))
args <- parser$parse_args()

MAX_DIST <- if (args$max_dist_um > 0) args$max_dist_um else Inf
if (is.finite(MAX_DIST))
  cat(sprintf("Distance cap: modelling cells with dist_to_phf1_um <= %g um\n", MAX_DIST))

# Clear the output dir first (default). Re-running after a script change otherwise leaves
# orphaned files from the older version sitting alongside the new ones, indistinguishable
# except by timestamp -- and an orphan that the current script CANNOT produce reads as
# real output. Only ever removes files this script writes, and only inside output_dir.
if (args$overwrite == "yes" && dir.exists(args$output_dir)) {
  old <- list.files(args$output_dir,
                    pattern = "^(plot|source_data|stats)_.*\\.(pdf|tsv|txt)$", full.names = TRUE)
  if (length(old)) {
    cat("Clearing", length(old), "existing output file(s) from", args$output_dir, "\n")
    unlink(old)
  }
}
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

# Constants (mirror the siblings) --------------------------------------------
MIN_GENES_PRESENT <- 10                   # drop genesets with fewer panel genes
CTRL              <- 100
NBIN              <- 24
SIG_ALPHA         <- 0.05                 # significance threshold (model FDR)
HALF_WINDOW       <- args$window_um / 2   # rolling-mean half-window (um)
GRID_N            <- 200
SIG_BASIS         <- "padj<0.05"
# Legend names. CON_LEGEND MUST differ from SIG_LEGEND, or ggplot merges the two guides
# (they share break labels) and the linetype/shape encodings collide.
SIG_LEGEND        <- "Log-distance"       # linetype/linewidth: the distance-model verdict
                                          # (the "padj<0.05" key below it says what the test is)
CON_LEGEND        <- "PHF1+ vs PHF1-"     # point shape: the PHF1+/- contrast verdict
CON_SHAPES        <- c(16, 1)   # filled = significant contrast, hollow = n.s.

XLAB <- expression("Distance to PHF1+ neuron (" * mu * "m)")  # mu via plotmath (survives default pdf device)
YLAB <- "Module score"

# Qualitative palette for the modules: 12 colours from R/reactome_filtered_modules.R
# (Okabe-Ito extended with four Tol muted hues), because the rule-derived stress group
# carries 10 modules. The pale yellow stays LAST -- it is illegible as a thin line on white.
# reactome_filtered_modules() stop()s if any group outgrows this, so a colour can never
# be silently reused for two modules.
MODULE_COLOURS <- MODULE_PALETTE
##  ............................................................................
##  Module configuration                                                    ####
# DERIVED BY RULE, not hand-picked. The module list is the Reactome PCD / autophagy /
# stress subtrees filtered by CATALOGUE_RULES (R/reactome_filtered_modules.R): panel
# size 25-80, panel coverage >= 70%, no surviving pair above Jaccard 0.8, one set per
# redundancy cluster, programme roots themselves excluded. 20 modules in 3 groups.
#
# Two things the rules cannot give, both recorded in every stats log:
#   - the ubiquitin-proteasome pathway (R-HSA-983168) is not included: Reactome files
#     it outside all three roots, so no rule over these subtrees can reach it.
#   - the three UPR arms are not separable; ATF6 (10 genes) and IRE1alpha fall below
#     the size floor, so UPR appears as the parent term.
#
# REACTOME_PRIMARY / MODULE_GROUPS / MODULE_LABELS / GROUP_LABELS are assigned AFTER the
# Seurat object is loaded, because the rules are defined relative to the CosMx panel and
# the panel is rownames(seu[["SCT"]]). See "Resolve gene sets" below.

CELLTYPES <- c("Exc-IT-L2-3-CBLN2-HOPX", "Astro", "Micro", "Exc-IT-L3-5-CHGA-IL1RAPL2", "Inh-VIP")


##  ............................................................................
##  Helpers (copied verbatim from the sibling scripts unless noted)          ####

pick_col <- function(meta, candidates, what) {
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0)
    stop(sprintf("None of the %s columns (%s) found in seu meta.data",
                 what, paste(candidates, collapse = ", ")))
  hit[1]
}

# Intersect each set with the panel; keep those with >= min_genes; build overlap table.
prep_sets <- function(sets, panel, min_genes, collection, group) {
  present <- lapply(sets, function(g) intersect(g, panel))
  n_total <- vapply(sets, length, integer(1))
  n_pres  <- vapply(present, length, integer(1))
  scored  <- n_pres >= min_genes
  overlap <- tibble(
    collection = collection, group = group, set = names(sets),
    n_total = as.integer(n_total), n_present = as.integer(n_pres),
    frac_present = round(n_pres / pmax(n_total, 1L), 4), scored = scored)
  list(present = present[scored], overlap = overlap)
}

# AddModuleScore for one prefix; returns data.frame of scores renamed to set names.
score_group <- function(seu_ct, present_sets, prefix, seed) {
  set.seed(seed)
  seu_ct <- tryCatch(
    AddModuleScore(seu_ct, features = present_sets, name = prefix,
                   assay = "SCT", ctrl = CTRL, nbin = NBIN, seed = seed),
    error = function(e) {
      cat("    AddModuleScore failed at nbin=", NBIN, " (", conditionMessage(e),
          ") - retrying at nbin=15\n", sep = "")
      set.seed(seed)
      AddModuleScore(seu_ct, features = present_sets, name = prefix,
                     assay = "SCT", ctrl = CTRL, nbin = 15, seed = seed)
    })
  cols <- paste0(prefix, seq_along(present_sets))     # AddModuleScore appends list index
  stopifnot(
    "AddModuleScore column count != number of scored sets" =
      length(cols) == length(present_sets) && all(cols %in% colnames(seu_ct@meta.data))
  )
  sc <- seu_ct@meta.data[, cols, drop = FALSE]
  colnames(sc) <- names(present_sets)                 # load-bearing index->name mapping
  sc
}

# Subset Seurat to a celltype and score every module in ONE AddModuleScore call, so all
# modules share one set of expression-matched control bins. Scoring is on the FULL
# celltype population (PHF1+ and PHF1-, all distances) -- the cap and the PHF1-negative
# filter are applied downstream and must never touch the control-gene background.
extract_celltype_data <- function(seu, ct, present_sets, seed, sex_col, age_col, pmi_col) {
  seu_ct <- subset(seu, celltype == ct)
  DefaultAssay(seu_ct) <- "SCT"
  md <- seu_ct@meta.data
  cat(sprintf("  %s cells: %d\n", ct, nrow(md)))

  id_cands <- c("cell_ID", "cell_id", "CellID")
  id_col   <- id_cands[id_cands %in% colnames(md)]
  cell_id  <- if (length(id_col) > 0) as.character(md[[id_col[1]]]) else rownames(md)

  base <- data.frame(
    cell_id         = cell_id,
    phf1_pos        = as.logical(as.character(md$PHF1)),   # boolean post-stain
    dist_to_phf1_um = as.numeric(md$dist_to_phf1_um),
    percent_neg     = as.numeric(md$percent.neg),
    nUMI_log        = log2(as.numeric(md$nCount_RNA) + 1),
    sample_id       = as.character(md$sample_id),
    # Braak is carried ONLY to colour the donor lines in the PHF1+/- paired figure; it is
    # not a covariate in any model here, and is tolerated as missing (the figure falls
    # back to colouring by donor) so a plotting nicety can never abort the analysis.
    Braak           = if ("Braak" %in% colnames(md)) as.character(md$Braak) else NA_character_,
    Sex             = as.character(md[[sex_col]]),
    Age             = as.numeric(md[[age_col]]),
    PMI             = as.numeric(md[[pmi_col]]),
    stringsAsFactors = FALSE
  )
  cbind(base, score_group(seu_ct, present_sets, "RXM_", seed))
}

# Observed LMM for one module.
fit_obs_lmer <- function(sub) {
  m <- tryCatch(
    lmerTest::lmer(score ~ dist_scaled + nUMI_log + percent_neg + Sex + Age_s + PMI_s +
                     (1 | sample_id), data = sub, REML = TRUE),
    error = function(e) { message("    obs lmer failed: ", conditionMessage(e)); NULL })
  if (is.null(m)) return(NULL)
  cf <- coef(summary(m))
  if (!"dist_scaled" %in% rownames(cf)) return(NULL)
  list(model = m,
       estimate = cf["dist_scaled", "Estimate"], se = cf["dist_scaled", "Std. Error"],
       t = cf["dist_scaled", "t value"], df = cf["dist_scaled", "df"],
       p = cf["dist_scaled", "Pr(>|t|)"], singular = isSingular(m))
}

# PHF1+ vs PHF1- contrast for one module -- what the strip's point SHAPE encodes.
#
# Fitted on ALL cells of the celltype with NO distance cap, deliberately: the cap bounds
# the distance axis, but applying it here would select the PHF1-negative comparison group
# on proximity to PHF1+ cells, which is the very gradient the rest of the figure is
# about, and would bias this contrast toward the null.
#
# Every donor contributes both groups, so this is a donor-PAIRED comparison: the log
# reports a paired Cohen's dz (mean paired difference / SD of differences) and the count
# of donors moving the same way alongside the cell-level estimate.
fit_phf1_contrast <- function(df_all, m) {
  sub <- df_all[, c("phf1_pos", "nUMI_log", "percent_neg", "Sex", "Age_s", "PMI_s",
                    "sample_id", "Braak")]
  sub$score <- df_all[[m]]
  sub <- sub[!is.na(sub$score), ]
  if (dplyr::n_distinct(sub$phf1_pos) < 2) {
    message(sprintf("    PHF1 contrast skipped for %s: only one PHF1 group present (%d+/%d-)",
                    m, sum(sub$phf1_pos), sum(!sub$phf1_pos)))
    return(NULL)
  }

  fit <- tryCatch(
    lmerTest::lmer(score ~ phf1_pos + nUMI_log + percent_neg + Sex + Age_s + PMI_s +
                     (1 | sample_id), data = sub, REML = TRUE),
    error = function(e) { message("    PHF1 contrast lmer failed: ", conditionMessage(e)); NULL })
  if (is.null(fit)) return(NULL)
  cf <- coef(summary(fit))
  rn <- grep("^phf1_pos", rownames(cf), value = TRUE)   # logical predictor -> 'phf1_posTRUE'
  if (!length(rn)) return(NULL)

  vc     <- as.data.frame(lme4::VarCorr(fit))
  sd_tot <- sqrt(sum(vc$vcov[is.na(vc$var2)]))
  est    <- cf[rn[1], "Estimate"]; se <- cf[rn[1], "Std. Error"]

  # donor-paired magnitude. per_donor is also returned: it IS the paired figure's data.
  per_donor <- sub %>%
    group_by(sample_id) %>%
    summarise(Braak = dplyr::first(Braak),
              n_pos = sum(phf1_pos), n_neg = sum(!phf1_pos),
              mean_pos = mean(score[phf1_pos]), mean_neg = mean(score[!phf1_pos]),
              .groups = "drop") %>%
    mutate(diff = mean_pos - mean_neg) %>%
    filter(n_pos >= 1, n_neg >= 1, is.finite(diff))
  dz <- if (nrow(per_donor) >= 2 && stats::sd(per_donor$diff) > 0)
          mean(per_donor$diff) / stats::sd(per_donor$diff) else NA_real_

  # ---- non-parametric cross-checks -------------------------------------------------
  # (a) CELL-level rank-sum (Mann-Whitney). Distribution-free. Effect size:
  #     AUC = P(a random PHF1+ cell > a random PHF1- cell), which is independent of n.
  #     rank-biserial r = 2*AUC - 1, in [-1, 1].
  xs <- sub$score[sub$phf1_pos]; ys <- sub$score[!sub$phf1_pos]
  wc <- suppressWarnings(stats::wilcox.test(xs, ys, exact = FALSE, correct = TRUE))
  auc <- as.numeric(wc$statistic) / (length(xs) * length(ys))

  # (b) DONOR-level signed-rank on the per-donor paired differences. n = donors, so this
  #     is the non-parametric analogue of dz.
  #     With 9 donors the smallest attainable two-sided p is 2/2^9 = 0.0039, so it CANNOT
  #     go below that however large the effect -- read it with the sign count, not alone.
  sr <- if (nrow(per_donor) >= 2)
          suppressWarnings(stats::wilcox.test(per_donor$diff)) else NULL
  nd <- nrow(per_donor); tot_rank <- nd * (nd + 1) / 2
  sr_rb <- if (!is.null(sr) && tot_rank > 0)
             (2 * as.numeric(sr$statistic) - tot_rank) / tot_rank else NA_real_

  list(estimate = est, se = se, ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
       wilcox_cell_W = as.numeric(wc$statistic), wilcox_cell_p = wc$p.value,
       wilcox_cell_auc = auc, wilcox_cell_rank_biserial = 2 * auc - 1,
       wilcox_donor_V = if (is.null(sr)) NA_real_ else as.numeric(sr$statistic),
       wilcox_donor_p = if (is.null(sr)) NA_real_ else sr$p.value,
       wilcox_donor_rank_biserial = sr_rb,
       t = cf[rn[1], "t value"], df = cf[rn[1], "df"], p = cf[rn[1], "Pr(>|t|)"],
       cohens_d = est / sd_tot, sd_total = sd_tot,
       dz = dz, mean_paired_diff = if (nrow(per_donor)) mean(per_donor$diff) else NA_real_,
       n_donors_paired = nrow(per_donor),
       n_donors_same_dir = if (nrow(per_donor))
         sum(sign(per_donor$diff) == sign(mean(per_donor$diff))) else NA_integer_,
       n_phf1_pos = sum(sub$phf1_pos), n_phf1_neg = sum(!sub$phf1_pos),
       singular = isSingular(fit), per_donor = per_donor)
}

# Model-free rolling (sliding-window) mean of y over x, evaluated on a grid.
roll_mean <- function(x, y, grid, half_window) {
  out <- lapply(grid, function(g) {
    idx <- which(x >= g - half_window & x <= g + half_window); n <- length(idx)
    if (n < 1) return(c(roll_mean = NA_real_, sem = NA_real_, n_window = 0))
    m <- mean(y[idx]); s <- if (n > 1) stats::sd(y[idx]) / sqrt(n) else NA_real_
    c(roll_mean = m, sem = s, n_window = n)
  })
  as.data.frame(do.call(rbind, out))
}

prep_celltype_df <- function(df) {
  df <- df[!is.na(df$dist_to_phf1_um) & !is.na(df$percent_neg) &
             !is.na(df$nUMI_log) & !is.na(df$Age) & !is.na(df$PMI) & !is.na(df$Sex), ]
  df$Sex       <- droplevels(factor(df$Sex))
  df$sample_id <- factor(df$sample_id)
  df
}

# Covariate filter for the PHF1+ vs PHF1- contrast.
#
# Deliberately NOT prep_celltype_df(). That one additionally requires a non-NA
# dist_to_phf1_um, which is exactly what a PHF1+ cell does NOT have -- it has no "nearest
# OTHER PHF1+ neuron" -- so reusing it silently deletes every PHF1+ cell and leaves the
# contrast with a single group. The contrast does not use distance at all, so it must not
# filter on it.
prep_contrast_df <- function(df) {
  df <- df[!is.na(df$percent_neg) & !is.na(df$nUMI_log) &
             !is.na(df$Age) & !is.na(df$PMI) & !is.na(df$Sex), ]
  df$Sex       <- droplevels(factor(df$Sex))
  df$sample_id <- factor(df$sample_id)
  df
}

# ---- panel geometry: identical axes across every distance figure --------------------
#
# Two things have to be forced, because ggplot gets both wrong for this figure family.
#
# HEIGHT. The panel is pinned to PANEL_H_MIN in every figure, so the y-axis is the same
# length everywhere and the canvas simply grows downwards for however many rows the
# legend wraps to.
#
# WIDTH. The neuron panels extend x to negative values to hold the PHF1+ strip, so their
# expanded x range is ~1254 data-units against ~1100 for the glia and _nostrip panels. At
# a fixed 1.6 in panel width that makes a micron 14% shorter on the strip figures than on
# the others -- the curves are not comparable by eye between them. Scaling the panel width
# by the x range instead pins MICRONS-PER-INCH, so 0->cap occupies the same physical
# length in every figure and the strip simply adds width on the left.
# TOTAL WIDTH is fixed at 10 cm for every distance figure. That forces the legend to the
# BOTTOM: measured on the real 10-module stress group, a right-hand legend needs 11.9 cm
# of canvas (10.7 cm for the 4-module autophagy group) before the panel gets any width at
# all, so it cannot fit and, worse, its width VARIES with the module count -- which would
# give each group a different axis length. A bottom legend costs 2.4 cm of furniture
# regardless of how many modules it holds, which is what makes one axis length possible.
TOTAL_W_CM  <- 10
PANEL_H_MIN <- 1.66   # floor, so small-legend panels keep the family's compact size

# PANEL_W_REF (inches given to the 0 -> x_cap range) is no longer a constant: it is
# solved for in flush_panels() so that the WIDEST style exactly fills TOTAL_W_CM. Set
# here only as the starting value used when a panel is drawn outside the queue.
PANEL_W_REF <- 1.6

#' Panel width in inches, scaled by the plot's own x range so that MICRONS-PER-INCH is
#' the same in every figure. A strip panel is wider than a distance-only panel by
#' exactly the strip it carries, so 0 -> cap occupies identical physical length in both.
panel_width_in <- function(p, x_span_ref, w_ref) {
  xr <- tryCatch(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x.range,
                 error = function(e) NULL)
  if (is.null(xr) || is.null(x_span_ref) || !is.finite(x_span_ref) || x_span_ref <= 0)
    w_ref else w_ref * diff(xr) / x_span_ref
}

#' Columns of the gtable holding the right-hand guide box.
guide_cols <- function(g) unique(g$layout$l[grepl("guide-box-right|^guide-box$", g$layout$name)])

#' Height of the right-hand guide box, in inches. It shares the panel ROW, so the panel
#' must be at least this tall or the legend is clipped without warning.
legend_height_in <- function(p) {
  g   <- ggplot2::ggplotGrob(p)
  idx <- which(grepl("guide-box-right|^guide-box$", g$layout$name))
  if (!length(idx)) return(0)
  suppressWarnings(max(c(0, vapply(idx, function(i) {
    h <- tryCatch(grid::convertHeight(grid::grobHeight(g$grobs[[i]]), "in", valueOnly = TRUE),
                  error = function(e) 0)
    if (is.finite(h)) h else 0
  }, numeric(1)))))
}

#' Width of the right-hand guide box, in inches. This is what the legend costs the axis.
guide_width_in <- function(p) {
  g  <- ggplot2::ggplotGrob(p)
  gc <- guide_cols(g)
  if (!length(gc)) return(0)
  w <- tryCatch(grid::convertWidth(sum(g$widths[gc]), "in", valueOnly = TRUE),
                error = function(e) 0)
  if (is.finite(w)) w else 0
}

#' Everything in the canvas that is NEITHER the panel NOR the legend: y-axis text, ticks
#' and margins. Both variable parts are zeroed so what remains is the fixed overhead, and
#' the panel and legend widths can then be allocated deliberately.
furniture_width_in <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  g$widths[pcol] <- grid::unit(0, "in")
  gc <- guide_cols(g)
  if (length(gc)) g$widths[gc] <- grid::unit(0, "in")
  tryCatch(grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
           error = function(e) NA_real_)
}

save_fixed_panel <- function(p, path, panel_h = PANEL_H_MIN, x_span_ref = NULL,
                             w_ref = PANEL_W_REF) {
  # w_ref is solved ONCE across the whole queue (flush_panels), so 0 -> x_cap is the
  # same physical length in every figure: the x axis is identical everywhere. A strip
  # panel is wider than a distance-only panel by exactly the strip it carries, which is
  # what makes the DRAWN AXIS match rather than merely the panel box.
  panel_w <- panel_width_in(p, x_span_ref, w_ref)

  g <- ggplot2::ggplotGrob(p)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  prow <- unique(g$layout$t[grepl("^panel", g$layout$name)])
  g$widths[pcol]  <- grid::unit(panel_w, "in")
  # Setting the panel row to an absolute unit is also what makes sum(g$heights)
  # convertible: by default it is unit(1, "null"), which has no inch value.
  g$heights[prow] <- grid::unit(panel_h, "in")

  # Absorb every spare millimetre into the GUTTER between panel and legend, so the
  # legend stays flush with the right edge and no figure ends in a void.
  #
  # There are two independent sources of slack. (1) The panel is solved from the widest
  # style, so a distance-only figure -- narrower, because it carries no strip -- leaves
  # ~0.5 cm over. (2) The panel solve budgets for the WIDEST legend in the run, so a
  # group with shorter labels (autophagy, 29 chars) leaves up to ~0.9 cm of legend
  # budget unspent. Neither can be given to the panel without breaking
  # the equal-axis guarantee, so it goes into the gutter instead, where it reads as
  # ordinary spacing.
  target_in <- grid::convertWidth(grid::unit(TOTAL_W_CM, "cm"), "in", valueOnly = TRUE)
  gc <- guide_cols(g)
  cur_in <- tryCatch(grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
                        error = function(e) NA_real_)
  if (is.finite(cur_in)) {
    slack <- target_in - cur_in
    if (slack >= 0) {
      # The column just left of the guide box is legend.box.spacing; widening it pushes
      # the legend right without touching the panel. Falls back to the last column only
      # if there is no legend to push against.
      pad_col <- if (length(gc) && min(gc) > 1L) min(gc) - 1L else length(g$widths)
      g$widths[pad_col] <- g$widths[pad_col] + grid::unit(slack, "in")
    } else {
      # Should not happen -- flush_panels() solves w_ref to prevent it -- but a silently
      # over-wide figure would break the axis-length guarantee, so say so.
      warning(sprintf("%s: content is %.2f cm wider than the %g cm target; the panel was ",
                      basename(path),
                      grid::convertWidth(grid::unit(-slack, "in"), "cm", valueOnly = TRUE),
                      TOTAL_W_CM),
              "shrunk to fit and its axis will NOT match the others.", call. = FALSE)
      g$widths[pcol] <- grid::unit(max(0.2, panel_w + slack), "in")
    }
  }

  total_h <- tryCatch(grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE),
                      error = function(e) NA_real_)
  if (!is.finite(total_h)) total_h <- panel_h + 0.44   # 0.44 = axes + margins
  ggsave(path, g, width = target_in, height = total_h, units = "in", device = "pdf")
  invisible(panel_w)
}

# Queue a distance panel; every queued panel is written at one common height by
# flush_panels() once all of them exist and the tallest legend is known.
PANEL_QUEUE <- list()
queue_fixed_panel <- function(p, path, x_span_ref) {
  PANEL_QUEUE[[length(PANEL_QUEUE) + 1L]] <<- list(p = p, path = path, x_span_ref = x_span_ref)
  invisible(NULL)
}

flush_panels <- function() {
  if (!length(PANEL_QUEUE)) return(invisible(NULL))
  # HEIGHT. A RIGHT-hand guide box shares the panel ROW, so a legend taller than the
  # panel is silently CLIPPED -- with 10 modules plus two significance guides the legend
  # needs ~1.98 in against a 1.66 in floor, and entries would simply vanish off the
  # bottom. Sizing each panel to its own legend fixes that but hands every group a
  # different y-axis length, so instead the panels are QUEUED, the tallest legend across
  # all of them is measured once, and every panel is written at that height: no
  # clipping, and one y-axis length everywhere.
  ph <- max(c(PANEL_H_MIN, vapply(PANEL_QUEUE, function(x) legend_height_in(x$p), numeric(1))))

  # Solve the reference panel width ONCE, across the whole queue, so that the widest
  # style exactly fills TOTAL_W_CM and no figure overflows. Each panel's width is then
  # w_ref scaled by its own x range, which keeps microns-per-inch identical across
  # styles: the strip figures are wider only by exactly the strip they carry, and every
  # figure of a given style therefore gets the same axis length.
  # WIDTH. One reference width for the whole run, so 0 -> x_cap is the same physical
  # length in every figure. It is solved from the TIGHTEST figure: the widest legend in
  # the run and the widest x range must both still fit 10 cm.
  #
  # The consequence is unavoidable and worth stating: a group whose labels are shorter
  # (autophagy's longest is 29 characters against stress's 41) does not need the legend
  # width the solve reserved for it, and that surplus becomes blank canvas. It is put in
  # the gutter between panel and legend by save_fixed_panel(), not at the page edge.
  # Spending it on the axis instead is possible, but only by letting each group have a
  # different axis length -- which is the thing this solve exists to prevent.
  target_in <- grid::convertWidth(grid::unit(TOTAL_W_CM, "cm"), "in", valueOnly = TRUE)
  scale_of  <- vapply(PANEL_QUEUE, function(x) panel_width_in(x$p, x$x_span_ref, 1), numeric(1))
  furn      <- vapply(PANEL_QUEUE, function(x) furniture_width_in(x$p), numeric(1))
  legend_w  <- suppressWarnings(max(c(0, vapply(PANEL_QUEUE,
                                                function(x) guide_width_in(x$p), numeric(1)))))
  if (!is.finite(legend_w)) legend_w <- 0

  ok    <- is.finite(scale_of) & is.finite(furn) & scale_of > 0
  w_ref <- if (any(ok)) min((target_in - furn[ok] - legend_w) / scale_of[ok]) else PANEL_W_REF
  if (!is.finite(w_ref) || w_ref <= 0.2) {
    warning("Could not solve a panel width that fits ", TOTAL_W_CM,
            " cm; falling back to ", PANEL_W_REF, " in.", call. = FALSE)
    w_ref <- PANEL_W_REF
  }
  widths <- vapply(PANEL_QUEUE, function(x)
    save_fixed_panel(x$p, x$path, panel_h = ph, x_span_ref = x$x_span_ref, w_ref = w_ref),
    numeric(1))

  cat(sprintf("\nWrote %d distance panel(s) at %g cm wide; panel height %.2f in; %.0f um per %.2f in;\n  legend budget %.2f in (the widest in the run).\n",
              length(PANEL_QUEUE), TOTAL_W_CM, ph,
              if (is.finite(MAX_DIST)) MAX_DIST else NA_real_, w_ref, legend_w))
  cat("Panel widths by style (in):\n")
  print(as.data.frame(table(`panel_width_in` = round(widths, 3))))
  cat("Two rows expected. 0 -> cap is the SAME physical length in every figure; a strip\n")
  cat("panel is wider only by the strip it carries, so the drawn x axis matches exactly.\n")
  invisible(NULL)
}

# PHF1+ reference (mean + 95% CI of the mean) for a given score column.
phf1_reference <- function(df, score_col) {
  s <- df[[score_col]][df$phf1_pos & !is.na(df[[score_col]])]
  n <- length(s)
  if (n < 1) return(NULL)
  m <- mean(s); sdv <- stats::sd(s); sem <- if (n > 1) sdv / sqrt(n) else NA_real_
  list(mean = m, sd = sdv, sem = sem, ci_lo = m - 1.96 * sem, ci_hi = m + 1.96 * sem, n = n)
}

##  ............................................................................
##  Load the object, derive the modules, score every celltype              ####
cat("\nLoading Seurat object:", args$seu, "\n")
seu <- readRDS(args$seu)
DefaultAssay(seu) <- "SCT"
stopifnot(
  "celltype column missing"        = "celltype"        %in% colnames(seu@meta.data),
  "PHF1 column missing"            = "PHF1"            %in% colnames(seu@meta.data),
  "dist_to_phf1_um column missing" = "dist_to_phf1_um" %in% colnames(seu@meta.data),
  "percent.neg column missing"     = "percent.neg"     %in% colnames(seu@meta.data),
  "nCount_RNA column missing"      = "nCount_RNA"      %in% colnames(seu@meta.data),
  "sample_id column missing"       = "sample_id"       %in% colnames(seu@meta.data)
)
sex_col <- pick_col(seu@meta.data, c("Sex"), "Sex")
age_col <- pick_col(seu@meta.data, c("Age"), "Age")
pmi_col <- pick_col(seu@meta.data, c("PMI", "PostMortemInterval", "PostMortem_Interval"), "PMI")

panel_genes <- rownames(seu[["SCT"]])
cat("CosMx panel genes (SCT assay):", length(panel_genes), "\n")

missing_ct <- setdiff(CELLTYPES, unique(as.character(seu$celltype)))
if (length(missing_ct))
  stop("Celltype(s) not present in the object: ", paste(missing_ct, collapse = ", "))

##  ............................................................................
##  Derive the module list from CATALOGUE_RULES                             ####
# Rebuilt in-process from the cached, provenance-stamped Reactome hierarchy + GMT, so
# this script runs standalone. It does NOT read results/reactome_programme_catalogue/.
# Because build_reactome_catalogue() is the same function the extraction script calls
# and the panel is the same 6175 genes, the gene lists here are IDENTICAL to that
# catalogue's genes_panel.
cat("\nDeriving Reactome modules from CATALOGUE_RULES...\n")
mods <- reactome_filtered_modules(panel_genes, gmt_path = args$gmt)

REACTOME_PRIMARY <- mods$REACTOME_PRIMARY
MODULE_GROUPS    <- mods$MODULE_GROUPS
MODULE_LABELS    <- mods$MODULE_LABELS
GROUP_LABELS     <- mods$GROUP_LABELS

stopifnot(
  "MODULE_GROUPS must partition REACTOME_PRIMARY" =
    setequal(unlist(MODULE_GROUPS, use.names = FALSE), names(REACTOME_PRIMARY)),
  "every module needs a label" = all(names(REACTOME_PRIMARY) %in% names(MODULE_LABELS)),
  "every group needs a label"  = all(names(MODULE_GROUPS)   %in% names(GROUP_LABELS)),
  "a group would exceed MODULE_COLOURS" =
    max(lengths(MODULE_GROUPS)) <= length(MODULE_COLOURS)
)
cat(sprintf("  %d modules: %s\n", length(REACTOME_PRIMARY),
            paste(sprintf("%s=%d", names(MODULE_GROUPS), lengths(MODULE_GROUPS)),
                  collapse = ", ")))

# GENES COME FROM THE CATALOGUE. Deliberately NOT re-resolved through
# reactome_genesets() + intersect(): the catalogue applies HGNC alias resolution
# (PARK2 -> PRKN, IL8 -> CXCL8, H2AFX -> H2AX, UFD1L -> UFD1) and a plain intersect
# does not, so a second pass would silently drop PRKN from PINK1-PRKN Mediated
# Mitophagy. These vectors are already panel-resolved and rule-filtered.
rx <- list(genes = mods$genes,
           provenance = tibble(module        = names(REACTOME_PRIMARY),
                               reactome_id   = unname(REACTOME_PRIMARY),
                               reactome_name = unname(mods$MODULE_NAMES),
                               enrichr_term  = mods$catalogue$enrichr_term,
                               n_genes       = as.integer(mods$catalogue$n_genes_reactome)))
print(as.data.frame(rx$provenance))

# prep_sets() is retained so every downstream structure is unchanged, but it is a
# no-op here by construction: the sets are already a subset of the panel and every
# one clears MIN_GENES_PRESENT (the rules impose a floor of 25). It still writes the
# overlap audit, which is the point.
prep <- prep_sets(rx$genes, panel_genes, MIN_GENES_PRESENT, "reactome", "reactome_primary")
SCORED_MODULES <- names(REACTOME_PRIMARY)[names(REACTOME_PRIMARY) %in% names(prep$present)]
dropped <- setdiff(names(REACTOME_PRIMARY), SCORED_MODULES)

overlap_tab <- rx$provenance %>%
  left_join(prep$overlap %>% dplyr::select(set, n_present, frac_present, scored),
            by = c("module" = "set")) %>%
  mutate(group = unname(vapply(module, function(m)
           names(MODULE_GROUPS)[vapply(MODULE_GROUPS, function(g) m %in% g, logical(1))][1],
           character(1))),
         module_label = unname(MODULE_LABELS[module]),
         min_genes_present = MIN_GENES_PRESENT) %>%
  dplyr::select(group, module, module_label, reactome_id, reactome_name, enrichr_term,
                n_total = n_genes, n_present, frac_present, min_genes_present, scored)
write.table(overlap_tab, file.path(args$output_dir, "stats_reactome_geneset_overlap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Gene-level audit: every gene of every set with an on_panel flag.
members <- bind_rows(lapply(names(rx$genes), function(m) tibble(
  module = m, module_label = unname(MODULE_LABELS[m]),
  reactome_id = unname(REACTOME_PRIMARY[m]),
  gene = rx$genes[[m]], on_panel = rx$genes[[m]] %in% panel_genes)))
write.table(members, file.path(args$output_dir, "stats_reactome_geneset_members.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nPanel overlap (floor =", MIN_GENES_PRESENT, "genes):\n")
print(as.data.frame(overlap_tab %>%
        dplyr::select(group, module, n_total, n_present, frac_present, scored)))
if (length(dropped))
  cat("\n[note] Dropped below the ", MIN_GENES_PRESENT, "-gene floor: ",
      paste(dropped, collapse = ", "), "\n", sep = "")
if (!length(SCORED_MODULES)) stop("No module cleared the panel-gene floor.")

MODULES_BY_GROUP <- lapply(MODULE_GROUPS, function(g) g[g %in% SCORED_MODULES])

# Colour is bound to the MODULE, not to its position in the legend. The legend is
# reordered per celltype by PHF1+ score (see run_group), so colouring by position would
# hand the same pathway a different colour in every panel and make the figures
# uncomparable at a glance -- the one thing a shared palette is for. Assigned once here,
# in the canonical (size-ordered) sequence, and looked up by module key thereafter.
MODULE_COLOUR_MAP <- unlist(lapply(MODULES_BY_GROUP, function(g)
  setNames(MODULE_COLOURS[seq_along(g)], g)), use.names = TRUE)
names(MODULE_COLOUR_MAP) <- sub("^[^.]*\\.", "", names(MODULE_COLOUR_MAP))

cat("\nScoring", length(SCORED_MODULES), "modules in", length(CELLTYPES), "celltypes...\n")
CT_DATA <- list()
for (ct in CELLTYPES) {
  CT_DATA[[ct]] <- extract_celltype_data(seu, ct, prep$present[SCORED_MODULES],
                                         seed = args$seed, sex_col = sex_col,
                                         age_col = age_col, pmi_col = pmi_col)
}
rm(seu); invisible(gc())

##  ............................................................................
##  Per-celltype model frame                                                ####
# PHF1-negative cells only, capped, with the CANONICAL distance transform. dist_sd is
# computed once per celltype (NOT per group), so a module's dist_scaled is identical in
# both groups. Age/PMI are scaled before the cap, matching the sibling scripts.
prep_model_df <- function(ct) {
  raw <- CT_DATA[[ct]]

  # Frame for the PHF1+ vs PHF1- contrast: ALL cells, no cap, and NO distance filter --
  # PHF1+ cells have no dist_to_phf1_um, so filtering on it would delete the very group
  # being contrasted. Age/PMI are scaled here and again on the distance frame below;
  # rescaling a covariate is a linear reparameterisation, so it moves neither model's
  # coefficient of interest.
  df_all <- prep_contrast_df(raw)
  df_all$Age_s <- as.numeric(scale(df_all$Age))
  df_all$PMI_s <- as.numeric(scale(df_all$PMI))
  n_pos <- sum(df_all$phf1_pos, na.rm = TRUE)
  if (is_neuron_celltype(ct) && any(raw$phf1_pos, na.rm = TRUE) && n_pos == 0)
    stop(sprintf(paste("%s: every PHF1+ cell was dropped from the contrast frame.",
                       "The covariate filter must not remove them -- see prep_contrast_df()."), ct))

  # Frame for the distance models: PHF1-negative cells only, capped, distance required.
  df_full <- prep_celltype_df(raw)
  df_full <- df_full[!df_full$phf1_pos, ]
  df_full$Age_s <- as.numeric(scale(df_full$Age))
  df_full$PMI_s <- as.numeric(scale(df_full$PMI))
  df <- if (is.finite(MAX_DIST)) df_full[df_full$dist_to_phf1_um <= MAX_DIST, ] else df_full

  # Canonical transform: plain natural log, scaled by SD only, NOT centred.
  # Distance is strictly > 0 for every PHF1-negative cell, so a non-positive value means
  # something upstream is wrong -- stop() rather than absorb it with an offset.
  if (any(!is.finite(df$dist_to_phf1_um) | df$dist_to_phf1_um <= 0))
    stop(sprintf("%s: non-positive or non-finite dist_to_phf1_um in %d modelled cell(s)",
                 ct, sum(!is.finite(df$dist_to_phf1_um) | df$dist_to_phf1_um <= 0)))
  dist_sd <- stats::sd(log(df$dist_to_phf1_um))
  df$dist_scaled <- log(df$dist_to_phf1_um) / dist_sd

  list(df = df, df_all = df_all, dist_sd = dist_sd, n_phf1_pos = n_pos,
       n_pre_cap = nrow(df_full), x_cap = if (is.finite(MAX_DIST)) MAX_DIST
                                          else as.numeric(quantile(df$dist_to_phf1_um, 0.99)))
}

# The PHF1+ strip is for NEURONS only. Gated on the celltype label (neuron_order from
# R/palettes.R), NOT on whether any cell is flagged PHF1+: PHF1 is a neuronal post-stain,
# so a handful of PHF1+ glia would be segmentation spillover, and drawing a strip off
# them would be an artefact.
is_neuron_celltype <- function(ct) ct %in% neuron_order

##  ............................................................................
##  Figure                                                                  ####
# Rolling means for every module of the group on one axis, RAW AddModuleScore scale
# (z-scoring flattens the modules onto each other, and the raw scale is what makes the
# PHF1+ strip readable). Line style = log-model BH significance, matching every other
# figure in this family. `refs` is NULL for celltypes with no PHF1+ cells (glia).
make_group_plot <- function(ct, group, modules, roll_df, ref_df, x_cap) {
  labs_lv <- unname(MODULE_LABELS[modules])
  cols    <- setNames(unname(MODULE_COLOUR_MAP[modules]), labs_lv)
  roll_df$module_label <- factor(roll_df$module_label, levels = labs_lv)

  # Complete the significance legend. See pad_levels() in R/palettes.R: with 20
  # rule-derived modules a whole group can come back with nothing significant, and a
  # "padj<0.05" key with no glyph reads as a rendering fault rather than as a real
  # (and interesting) null. Only the linetype/linewidth guide needs this -- the module
  # COLOUR guide is already complete from drop = FALSE alone.
  roll_df <- pad_levels(roll_df, "sig", c(SIG_BASIS, "ns"))

  p <- ggplot(roll_df, aes(dist_to_phf1_um, roll_mean,
                           colour = module_label, fill = module_label, group = module_label)) +
    geom_ribbon(aes(ymin = roll_mean - 1.96 * sem, ymax = roll_mean + 1.96 * sem),
                alpha = 0.15, colour = NA) +
    geom_line(aes(linetype = sig, linewidth = sig))

  # PHF1+ marginal strip at the LEFT edge: one point + 95% CI per module, dodged so the
  # error bars do not overlap, separated from the data by a thin rule. Colour comes from
  # the module scale (no extra legend entry); SHAPE encodes the PHF1+ vs PHF1- contrast
  # (filled = significant, hollow = n.s.) and gets its own legend, whose name must differ
  # from the linetype legend's or ggplot merges the two. clip="off" guards the label only.
  has_strip <- !is.null(ref_df) && nrow(ref_df) > 0
  if (has_strip) {
    ref_df$module_label <- factor(ref_df$module_label, levels = labs_lv)
    ref_df$con_sig      <- factor(ref_df$con_sig, levels = c(SIG_BASIS, "ns"))
    # Same reason as the linetype guide: if no module's PHF1+/- contrast is significant,
    # the filled-circle key would be labelled but empty.
    ref_df <- pad_levels(ref_df, "con_sig", c(SIG_BASIS, "ns"))
    p <- p +
      annotate("segment", x = -0.025 * x_cap, xend = -0.025 * x_cap,
               y = -Inf, yend = Inf, linewidth = 0.25, colour = "grey70") +
      geom_errorbar(data = ref_df,
                    aes(x = x_pos, y = ref_mean, ymin = ci_lo, ymax = ci_hi,
                        colour = module_label),
                    inherit.aes = FALSE, width = 0, linewidth = 0.4) +
      geom_point(data = ref_df,
                 aes(x = x_pos, y = ref_mean, colour = module_label, shape = con_sig),
                 inherit.aes = FALSE, size = 1.3, stroke = 0.5) +
      scale_shape_manual(values = setNames(CON_SHAPES, c(SIG_BASIS, "ns")),
                         name = CON_LEGEND, drop = FALSE, limits = c(SIG_BASIS, "ns")) +
      # Sits ABOVE the panel border (negative vjust + clip="off"), not on top of the data.
      # The extra top margin below is what stops it being cut off at the canvas edge.
      annotate("text", x = mean(ref_df$x_pos), y = Inf, vjust = -0.45,
               size = 2.6, fontface = "bold", colour = "grey15", label = "PHF1+")
  }

  brk <- pretty(c(0, x_cap)); brk <- brk[brk >= 0 & brk <= x_cap]
  xlo <- if (has_strip) -0.14 * x_cap else 0

  p +
    scale_colour_manual(values = cols, name = GROUP_LABELS[[group]], drop = FALSE) +
    scale_fill_manual(values = cols, guide = "none", drop = FALSE) +
    scale_linetype_manual(values = setNames(c("solid", "dashed"), c(SIG_BASIS, "ns")),
                          name = SIG_LEGEND, drop = FALSE, limits = c(SIG_BASIS, "ns")) +
    scale_linewidth_manual(values = setNames(c(0.8, 0.4), c(SIG_BASIS, "ns")),
                           name = SIG_LEGEND, drop = FALSE, limits = c(SIG_BASIS, "ns")) +
    scale_x_continuous(breaks = brk) +
    guides(colour = guide_legend(order = 1),
           linetype = guide_legend(order = 2), linewidth = guide_legend(order = 2),
           shape = guide_legend(order = 3,
                                override.aes = list(colour = "grey25", size = 1.3))) +
    labs(x = XLAB, y = YLAB) +
    coord_cartesian(xlim = c(xlo, x_cap), clip = "off") +
    theme_classic(base_size = 8) + fig_theme +
    # Legend RIGHT. It is given a COMMON width across every figure (solved in
    # flush_panels), and the panel takes what is left of the 10 cm -- so the x axis is
    # short, but it is the SAME short in every figure of a style, which is the property
    # that matters. Letting each legend take its natural width instead would make the
    # axis a different length in each group.
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "right", legend.box = "vertical",
          legend.text = element_text(size = 6), legend.title = element_text(size = 6),
          legend.key.width = grid::unit(14, "pt"), legend.key.height = grid::unit(8, "pt"),
          legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
          legend.spacing.y = grid::unit(1, "pt"),
          # Headroom for the "PHF1+" label above the panel. Applied to EVERY panel, not
          # just the strip ones: with a common panel height, a differing top margin would
          # still leave the strip and _nostrip/glia panels vertically offset from each
          # other on the page. Costs 8 pt of whitespace on the panels that don't need it.
          plot.margin = margin(t = 14, r = 6, b = 4, l = 5))
}

##  ............................................................................
##  PHF1+ vs PHF1- figures (neuron celltypes only)                          ####

# Forest: one row per module, the PHF1+ minus PHF1- LMM estimate with its 95% CI, with
# the per-donor differences as faint points behind (the same background-points idiom as
# decay_lambda_forest). Point shape = BH significance, matching the strip on the main
# figure. The donor points show the per-donor differences behind the estimate.
make_contrast_forest <- function(ct, group, con_stats, modules) {
  labs_lv <- rev(unname(MODULE_LABELS[modules]))    # rev so the first module sits at the top
  cols    <- setNames(unname(MODULE_COLOUR_MAP[modules]), unname(MODULE_LABELS[modules]))
  cs <- con_stats %>% mutate(module_label = factor(module_label, levels = labs_lv),
                             con_sig = factor(ifelse(significant %in% TRUE, SIG_BASIS, "ns"),
                                              levels = c(SIG_BASIS, "ns")))
  dp <- bind_rows(lapply(seq_len(nrow(cs)), function(i) {
          pd <- cs$per_donor[[i]]
          if (is.null(pd) || !nrow(pd)) return(NULL)
          data.frame(module_label = cs$module_label[i], diff = pd$diff)
        }))

  # Padded AFTER the donor points are extracted, deliberately. pad_levels() clones a
  # real row to build the filler, and `cs` carries a per_donor LIST column -- padding
  # first would clone that module's donor differences and draw them twice.
  cs <- pad_levels(cs, "con_sig", c(SIG_BASIS, "ns"))

  p <- ggplot(cs, aes(x = estimate, y = module_label, colour = module_label))
  if (!is.null(dp) && nrow(dp))
    p <- p + geom_point(data = dp, aes(x = diff, y = module_label, colour = module_label),
                        inherit.aes = FALSE, alpha = 0.35, size = 0.7, stroke = 0,
                        position = position_nudge(y = 0.22))
  p +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40", linetype = "dashed") +
    # geom_linerange, not geom_errorbarh: the latter is deprecated in current ggplot2 and
    # we want no end caps anyway.
    geom_linerange(aes(xmin = ci_lo, xmax = ci_hi), linewidth = 0.5) +
    geom_point(aes(shape = con_sig), size = 1.6, stroke = 0.5) +
    scale_colour_manual(values = cols, guide = "none", drop = FALSE) +
    scale_shape_manual(values = setNames(CON_SHAPES, c(SIG_BASIS, "ns")),
                       name = CON_LEGEND, drop = FALSE, limits = c(SIG_BASIS, "ns")) +
    labs(x = "PHF1+ minus PHF1- (module score)", y = NULL) +
    forest_theme + theme(legend.position = "bottom")
}

# Paired donor view: one facet per module, a line per donor joining its PHF1- mean to its
# PHF1+ mean, coloured by Braak. This is what dz measures, and it shows whether the
# donors agree.
make_contrast_bydonor <- function(ct, group, con_stats) {
  long <- bind_rows(lapply(seq_len(nrow(con_stats)), function(i) {
    pd <- con_stats$per_donor[[i]]
    if (is.null(pd) || !nrow(pd)) return(NULL)
    data.frame(module_label = con_stats$module_label[i], sample_id = pd$sample_id,
               Braak = pd$Braak,
               rbind(data.frame(status = "PHF1-", score = pd$mean_neg),
                     data.frame(status = "PHF1+", score = pd$mean_pos)))
  }))
  if (is.null(long) || !nrow(long)) return(NULL)
  long$status       <- factor(long$status, levels = c("PHF1-", "PHF1+"))
  long$module_label <- factor(long$module_label, levels = con_stats$module_label)

  use_braak <- any(!is.na(long$Braak)) && all(na.omit(long$Braak) %in% names(braak_line_palette))
  long$grp  <- if (use_braak) factor(long$Braak, levels = braak_levels) else factor(long$sample_id)

  p <- ggplot(long, aes(status, score, group = sample_id, colour = grp)) +
    geom_line(linewidth = 0.4, alpha = 0.85) +
    geom_point(size = 0.8, stroke = 0) +
    facet_wrap(~ module_label, scales = "free_y",
               ncol = min(4L, nlevels(long$module_label))) +
    labs(x = NULL, y = YLAB) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          strip.background = element_blank(), strip.text = element_text(size = 6),
          legend.text = element_text(size = 6), legend.title = element_text(size = 6),
          legend.position = "right", plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
  p <- if (use_braak)
         p + scale_colour_manual(values = braak_line_palette, name = "Braak", drop = FALSE)
       else p + scale_colour_discrete(name = "Donor")
  list(plot = p, long = long, n_facets = nlevels(long$module_label))
}

##  ............................................................................
##  One celltype x group unit: fit, plot, write the triple                  ####
lmm_coef_rows <- list()

run_group <- function(ct, group, mdf) {
  modules <- MODULES_BY_GROUP[[group]]
  if (!length(modules)) {
    cat("  [skip] ", ct, " / ", group, ": no module cleared the panel floor\n", sep = "")
    return(invisible(NULL))
  }
  df <- mdf$df; x_cap <- mdf$x_cap
  cat(sprintf("\n=== %s / %s (%d modules, %d cells, %d donors) ===\n",
              ct, group, length(modules), nrow(df), dplyr::n_distinct(df$sample_id)))

  # --- models -------------------------------------------------------------
  rows <- list()
  for (m in modules) {
    sub <- df[, c("cell_id", "dist_scaled", "nUMI_log", "percent_neg",
                  "Sex", "Age_s", "PMI_s", "sample_id")]
    sub$score <- df[[m]]
    obs <- fit_obs_lmer(sub)
    if (is.null(obs)) { cat("  [warn] LMM failed for ", m, "\n", sep = ""); next }
    # Cohen's d: estimate / sqrt(all random-effect variances + residual). Same idiom as
    # R/imc_phf1_morphology_exc.R (drop the covariance rows: is.na(var2)).
    vc     <- as.data.frame(lme4::VarCorr(obs$model))
    sd_tot <- sqrt(sum(vc$vcov[is.na(vc$var2)]))
    rows[[m]] <- tibble(
      group = group, celltype = ct, module = m,
      module_label   = unname(MODULE_LABELS[m]),
      reactome_id    = unname(REACTOME_PRIMARY[m]),
      distance_scale = "log_um",
      slope_per_sd = obs$estimate, se = obs$se,
      ci_lo_per_sd = obs$estimate - 1.96 * obs$se,
      ci_hi_per_sd = obs$estimate + 1.96 * obs$se,
      slope_per_unit = obs$estimate / mdf$dist_sd,
      cohens_d = obs$estimate / sd_tot, sd_total = sd_tot,
      t = obs$t, df = obs$df, p_LMM = obs$p, singular = obs$singular,
      dist_sd = mdf$dist_sd, n_cells = nrow(sub),
      n_donors = dplyr::n_distinct(sub$sample_id))
  }
  if (!length(rows)) { cat("  [skip] ", ct, " / ", group, ": every LMM failed\n", sep = ""); return(invisible(NULL)) }

  mod_stats <- bind_rows(rows)
  # BH across the modules of THIS celltype x group (the glia-script precedent).
  mod_stats$lmm_padj    <- p.adjust(mod_stats$p_LMM, method = "BH")
  mod_stats$significant <- mod_stats$lmm_padj < SIG_ALPHA
  mod_stats$verdict     <- ifelse(is.na(mod_stats$significant), "not testable",
                            ifelse(mod_stats$significant,
                                   "significant distance effect (model FDR<0.05)", "n.s. (model FDR)"))
  modules <- modules[modules %in% mod_stats$module]          # keep only modules that fitted
  sig_lut <- setNames(mod_stats$significant %in% TRUE, mod_stats$module)

  # --- legend order: PHF1+ module score, highest first ---------------------
  # The legend reads top-to-bottom in the order the modules score INSIDE tangle-bearing
  # neurons, so the entry at the top is the pathway most active in PHF1+ cells. This
  # sets the order for everything downstream in one place -- rolling-mean factor levels,
  # the strip's left-to-right x positions, and the contrast forest -- because all three
  # are built from `modules` after this point.
  #
  # Only celltypes with PHF1+ cells have a score to order by. Glia (and any neuron
  # subtype with none) keep the canonical size order rather than being given an
  # arbitrary one; their legends therefore need not match a neuron panel's, which is
  # stated in the stats log rather than hidden.
  phf1_order_src <- "canonical (n_genes_panel, descending) -- no PHF1+ cells to order by"
  if (is_neuron_celltype(ct) && any(CT_DATA[[ct]]$phf1_pos, na.rm = TRUE)) {
    ref_mean_of <- vapply(modules, function(m) {
      r <- phf1_reference(CT_DATA[[ct]], m)
      if (is.null(r)) NA_real_ else r$mean
    }, numeric(1))
    # na.last keeps any module without a reference at the bottom instead of dropping it.
    modules <- modules[order(-ref_mean_of, na.last = TRUE)]
    phf1_order_src <- "PHF1+ module score, descending"
  }
  cat("  Legend order: ", phf1_order_src, "\n", sep = "")

  # --- rolling means ------------------------------------------------------
  grid <- seq(0, x_cap, length.out = GRID_N)
  roll_df <- bind_rows(lapply(modules, function(m) {
    y <- df[[m]]; ok <- !is.na(y) & !is.na(df$dist_to_phf1_um)
    data.frame(module = m, module_label = unname(MODULE_LABELS[m]),
               dist_to_phf1_um = grid,
               roll_mean(df$dist_to_phf1_um[ok], y[ok], grid, HALF_WINDOW))
  }))
  roll_df$sig <- factor(ifelse(sig_lut[roll_df$module], SIG_BASIS, "ns"),
                        levels = c(SIG_BASIS, "ns"))

  # --- PHF1+ reference strip + PHF1+ vs PHF1- contrast (NEURONS only) ------
  ref_df <- NULL; con_stats <- NULL
  if (is_neuron_celltype(ct) && any(CT_DATA[[ct]]$phf1_pos, na.rm = TRUE)) {
    # Strip points spread across ~85% of the strip, from the panel's left edge
    # (xlo = -0.14 * x_cap) to just short of the separating rule (-0.025 * x_cap),
    # leaving a small gap at each end so no point sits on the rule or is clipped by the
    # edge. With up to 10 modules the spread is what keeps the error bars
    # distinguishable.
    n    <- length(modules)
    xpos <- if (n == 1) -0.08 * x_cap else seq(-0.130, -0.033, length.out = n) * x_cap
    ref_rows <- list(); con_rows <- list()
    for (i in seq_along(modules)) {
      m <- modules[i]
      r <- phf1_reference(CT_DATA[[ct]], m)
      if (is.null(r)) next
      ref_rows[[m]] <- data.frame(
        module = m, module_label = unname(MODULE_LABELS[m]),
        x_pos = xpos[i], ref_mean = r$mean, ci_lo = r$ci_lo, ci_hi = r$ci_hi,
        sd = r$sd, sem = r$sem, n_phf1_pos = r$n)
      cn <- fit_phf1_contrast(mdf$df_all, m)
      if (is.null(cn)) { cat("  [warn] PHF1+ vs PHF1- contrast failed for ", m, "\n", sep = ""); next }
      con_rows[[m]] <- tibble(
        group = group, celltype = ct, module = m, module_label = unname(MODULE_LABELS[m]),
        estimate = cn$estimate, se = cn$se, ci_lo = cn$ci_lo, ci_hi = cn$ci_hi,
        cohens_d = cn$cohens_d, sd_total = cn$sd_total,
        dz_donor_paired = cn$dz, mean_paired_diff = cn$mean_paired_diff,
        n_donors_paired = cn$n_donors_paired, n_donors_same_dir = cn$n_donors_same_dir,
        t = cn$t, df = cn$df, p_LMM = cn$p,
        wilcox_cell_W = cn$wilcox_cell_W, wilcox_cell_p = cn$wilcox_cell_p,
        wilcox_cell_auc = cn$wilcox_cell_auc,
        wilcox_cell_rank_biserial = cn$wilcox_cell_rank_biserial,
        wilcox_donor_V = cn$wilcox_donor_V, wilcox_donor_p = cn$wilcox_donor_p,
        wilcox_donor_rank_biserial = cn$wilcox_donor_rank_biserial,
        n_phf1_pos = cn$n_phf1_pos, n_phf1_neg = cn$n_phf1_neg, singular = cn$singular,
        per_donor = list(cn$per_donor))   # list-column: the paired figure's data
    }
    if (length(ref_rows)) ref_df <- bind_rows(ref_rows)
    if (length(con_rows)) {
      con_stats <- bind_rows(con_rows)
      # BH across the modules of this celltype x group -- its own family, parallel to
      # (and separate from) the distance test's.
      con_stats$padj        <- p.adjust(con_stats$p_LMM, method = "BH")
      con_stats$significant <- con_stats$padj < SIG_ALPHA
      # The Wilcoxons are cross-checks, BH-adjusted in their own families so each test
      # is corrected on the same footing. The FIGURE's shape encoding stays on the LMM.
      con_stats$wilcox_cell_padj  <- p.adjust(con_stats$wilcox_cell_p,  method = "BH")
      con_stats$wilcox_donor_padj <- p.adjust(con_stats$wilcox_donor_p, method = "BH")
    }
    # A module whose contrast could not be fitted is drawn hollow and flagged in the log.
    if (!is.null(ref_df)) {
      con_lut <- if (is.null(con_stats)) setNames(logical(0), character(0))
                 else setNames(con_stats$significant %in% TRUE, con_stats$module)
      ref_df$con_significant <- unname(con_lut[ref_df$module] %in% TRUE)
      ref_df$con_padj        <- if (is.null(con_stats)) NA_real_ else
        unname(setNames(con_stats$padj, con_stats$module)[ref_df$module])
      ref_df$con_sig         <- ifelse(ref_df$con_significant, SIG_BASIS, "ns")
    }
  }

  # --- figures ------------------------------------------------------------
  # Queued, not written yet: flush_panels() gives every distance panel the SAME panel
  # height (tallest legend across all of them) and the same microns-per-inch.
  x_span_ref <- 1.1 * x_cap    # the expanded x range of a plain 0 -> x_cap panel
  queue_fixed_panel(make_group_plot(ct, group, modules, roll_df, ref_df, x_cap),
                    file.path(args$output_dir,
                              sprintf("plot_reactome_%s_dist_%s.pdf", group, ct)),
                    x_span_ref)
  # Neuron panels are ALSO queued without the PHF1+ strip: the strip's y-range can compress
  # the curves, so the "_nostrip" copy is the one to use when the gradient itself is the
  # point. Same data, same models, same line encoding -- only the strip is dropped.
  if (!is.null(ref_df))
    queue_fixed_panel(make_group_plot(ct, group, modules, roll_df, NULL, x_cap),
                      file.path(args$output_dir,
                                sprintf("plot_reactome_%s_dist_%s_nostrip.pdf", group, ct)),
                      x_span_ref)

  # --- source data --------------------------------------------------------
  write.table(roll_df %>%
                mutate(group = group, celltype = ct, window_um = args$window_um,
                       log_significant = sig == SIG_BASIS) %>%
                dplyr::select(group, celltype, module, module_label, window_um,
                              log_significant, dist_to_phf1_um, roll_mean, sem, n_window),
              file.path(args$output_dir,
                        sprintf("source_data_reactome_%s_dist_%s_rollmean.tsv", group, ct)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  if (!is.null(ref_df))
    write.table(ref_df %>% mutate(group = group, celltype = ct) %>%
                  dplyr::select(group, celltype, module, module_label, x_pos,
                                ref_mean, ci_lo, ci_hi, sd, sem, n_phf1_pos,
                                con_significant, con_padj),
                file.path(args$output_dir,
                          sprintf("source_data_reactome_%s_dist_%s_phf1ref.tsv", group, ct)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  if (!is.null(con_stats)) {
    write.table(con_stats %>% dplyr::select(-per_donor),   # drop the list-column
                file.path(args$output_dir,
                sprintf("stats_reactome_%s_dist_%s_phf1contrast.tsv", group, ct)),
                sep = "\t", quote = FALSE, row.names = FALSE)

    # Forest of the PHF1+ minus PHF1- estimate, one row per module.
    fp <- make_contrast_forest(ct, group, con_stats, modules)
    ggsave(file.path(args$output_dir,
           sprintf("plot_reactome_%s_phf1contrast_forest_%s.pdf", group, ct)), fp,
           width = TOTAL_W_CM, height = 2.0 + 0.42 * nrow(con_stats), units = "cm", device = "pdf")

    # Paired donor view + its source data (one row per donor x module x PHF1 status).
    bd <- make_contrast_bydonor(ct, group, con_stats)
    if (!is.null(bd)) {
      g <- ggplot2::ggplotGrob(bd$plot)
      pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
      g$widths[pcol] <- grid::unit(1.05, "in")
      ggsave(file.path(args$output_dir,
             sprintf("plot_reactome_%s_phf1contrast_bydonor_%s.pdf", group, ct)), g,
             width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
             height = 0.6 + 1.25 * ceiling(bd$n_facets / min(4L, bd$n_facets)),
             units = "in", device = "pdf")
      write.table(bd$long %>% mutate(group = group, celltype = ct) %>%
                    dplyr::select(group, celltype, module_label, sample_id, Braak,
                                  status, score),
                  file.path(args$output_dir,
                  sprintf("source_data_reactome_%s_phf1contrast_bydonor_%s.tsv", group, ct)),
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }

  # --- structured stats tables -------------------------------------------
  write.table(mod_stats %>% dplyr::select(
                group, celltype, module, module_label, reactome_id, distance_scale,
                slope_per_sd, se, ci_lo_per_sd, ci_hi_per_sd, slope_per_unit,
                t, df, p_LMM, lmm_padj, significant, singular, dist_sd, n_cells, n_donors),
              file.path(args$output_dir,
                        sprintf("stats_reactome_%s_dist_%s_lmm.tsv", group, ct)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(mod_stats %>% dplyr::select(
                group, celltype, module, module_label,
                estimate_per_sd_log_dist = slope_per_sd,
                ci_lo_per_sd, ci_hi_per_sd, cohens_d, sd_total,
                n_cells, n_donors, p_LMM, lmm_padj, significant),
              file.path(args$output_dir,
                        sprintf("stats_reactome_%s_dist_%s_effectsize.tsv", group, ct)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # --- human-readable stats log ------------------------------------------
  sink(file.path(args$output_dir, sprintf("stats_reactome_%s_dist_%s.txt", group, ct)))
  cat("Reactome ", GROUP_LABELS[[group]], " module scores vs distance to PHF1+ neuron\n", sep = "")
  cat("Group:", group, " | celltype:", ct, " | distance scale: log_um\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

  # The module list is derived, not curated, so the rules that produced it are part
  # of the result and travel with every log. A reader cannot otherwise tell why a
  # given pathway is here and an obvious sibling is not.
  cat("LEGEND ORDER: ", phf1_order_src, "\n", sep = "")
  cat("  Top entry = highest mean module score in PHF1+ (tangle-bearing) neurons. The\n")
  cat("  same order drives the PHF1+ strip left-to-right and the contrast forest.\n")
  cat("  Colour is bound to the pathway, NOT to legend position, so a pathway keeps its\n")
  cat("  colour across celltypes even where the ordering differs.\n")
  if (!is.null(ref_df) && nrow(ref_df))
    print(as.data.frame(ref_df[order(-ref_df$ref_mean), c("module_label", "ref_mean")] |>
            transform(ref_mean = round(ref_mean, 4))), row.names = FALSE)
  cat("\n")

  cat("MODULE SELECTION (rule-derived, not hand-picked):\n")
  cat("  Reactome PCD / autophagy / stress subtrees, filtered by CATALOGUE_RULES:\n")
  for (nm in names(CATALOGUE_RULES))
    cat(sprintf("    %-22s %s\n", nm, format(CATALOGUE_RULES[[nm]], trim = TRUE)))
  cat("  Defaults left in force: exclude_disease = TRUE, require_gmt_match = TRUE,\n")
  cat("    min_genes_unique = 0 (NOT applied).\n")
  cat("  ", paste(strwrap(CATALOGUE_RULES_CAVEAT, width = 88, prefix = "  "),
                  collapse = "\n"), "\n", sep = "")
  cat("  Genes are the catalogue's panel-resolved sets (HGNC alias resolution applied),\n")
  cat("    not a fresh GMT intersect. Full member list: stats_reactome_geneset_members.tsv\n\n")

  cat("MODEL (one per module; LMM, log distance):\n")
  cat("  score ~ dist_to_phf1_um_scaled + nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)\n")
  cat("  Outcome = Seurat AddModuleScore of the Reactome geneset (SCT data slot; ctrl=", CTRL,
      " nbin=", NBIN, " seed=", args$seed, ").\n", sep = "")
  cat("  All modules of this celltype were scored in ONE AddModuleScore call.\n")
  cat("  Distance transform (canonical): log(dist_to_phf1_um) / sd(log(dist_to_phf1_um)),\n")
  cat("    natural log, scaled by SD only, NOT centred; dist_sd computed WITHIN this celltype\n")
  cat("    across both groups: dist_sd =", round(mdf$dist_sd, 6), "\n")
  cat("  Distance cap:", if (is.finite(MAX_DIST)) sprintf("%g um", MAX_DIST) else "none", "\n")
  cat("  Cells modelled:", nrow(df), "of", mdf$n_pre_cap, "PHF1-negative (pre-cap) | donors:",
      dplyr::n_distinct(df$sample_id), "\n")
  cat("  PHF1+ cells excluded from every model and rolling mean:", mdf$n_phf1_pos, "\n")
  cat("  Donor random intercept: (1|sample_id).\n\n")

  cat("GENE SETS (Reactome stable IDs -> cached Enrichr GMT; see stats_reactome_geneset_overlap.tsv):\n")
  print(as.data.frame(overlap_tab %>%
          dplyr::filter(module %in% modules) %>%
          dplyr::select(module, reactome_id, reactome_name, n_total, n_present, scored)))
  if (length(dropped))
    cat("\n  Dropped below the ", MIN_GENES_PRESENT, "-gene panel floor (all groups): ",
        paste(dropped, collapse = ", "), "\n", sep = "")
  cat("\n")

  cat("SIGNIFICANCE = BH-adjusted model p across the", nrow(mod_stats),
      "modules of this celltype x group; padj <", SIG_ALPHA, "-> significant.\n")
  cat("  Figure: solid + thick = significant, dashed + thin = n.s.\n\n")

  cat("Per-module results:\n")
  print(as.data.frame(mod_stats %>% dplyr::select(
    module, slope_per_sd, ci_lo_per_sd, ci_hi_per_sd, t, df, p_LMM, lmm_padj,
    significant, singular, verdict)))

  cat("\n=== Effect sizes ===\n")
  cat("slope_per_sd = change in module score per 1 SD of log distance (unstandardised),\n")
  cat("  with its 95% CI. cohens_d = slope_per_sd / sqrt(sum of the random-effect variances\n")
  cat("  + residual), i.e. the same slope in total-SD units.\n")
  cat("slope_per_unit = change per 1 log-um (slope_per_sd / dist_sd).\n\n")
  print(as.data.frame(mod_stats %>% dplyr::transmute(
    module,
    estimate_per_sd = round(slope_per_sd, 6),
    ci95 = sprintf("[%.6f, %.6f]", ci_lo_per_sd, ci_hi_per_sd),
    cohens_d = round(cohens_d, 4),
    per_log_um = round(slope_per_unit, 6),
    n_cells, n_donors)))

  if (!is.null(ref_df)) {
    cat("\nPHF1+ reference (drawn as the LEFT-EDGE marginal strip; NOT in the distance regression):\n")
    print(as.data.frame(ref_df %>%
            dplyr::select(module, ref_mean, ci_lo, ci_hi, sd, n_phf1_pos)))

    cat("\n=== PHF1+ vs PHF1- contrast (the strip's POINT SHAPE) ===\n")
    cat("  score ~ phf1_pos + nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)\n")
    cat("Fitted on ALL", nrow(mdf$df_all), "cells of this celltype, NO distance cap: the cap bounds\n")
    cat("  the distance axis, but applying it here would select the PHF1-negative comparison group\n")
    cat("  on proximity to PHF1+ cells -- the very gradient this figure is about -- and would bias\n")
    cat("  the contrast toward the null.\n")
    cat("Filled point = significant (BH across this celltype x group), hollow = n.s.\n")
    cat("Every donor contributes both groups, so this is a donor-PAIRED comparison: dz_donor_paired\n")
    cat("  = mean per-donor (PHF1+ minus PHF1-) difference / SD of those differences, and\n")
    cat("  n_donors_same_dir counts donors moving with the mean.\n\n")
    if (!is.null(con_stats)) {
      print(as.data.frame(con_stats %>% dplyr::select(
        module, estimate, ci_lo, ci_hi, cohens_d, dz_donor_paired,
        n_donors_same_dir, n_donors_paired, p_LMM, padj, significant)))

      cat("\n--- Wilcoxon cross-checks (non-parametric; the figure's shape stays on the LMM) ---\n")
      cat("(a) CELL-level rank-sum (Mann-Whitney), PHF1+ vs PHF1- per-cell scores.\n")
      cat("    Distribution-free (", con_stats$n_phf1_pos[1], " PHF1+ and ",
          con_stats$n_phf1_neg[1], " PHF1- cells from ", con_stats$n_donors_paired[1],
          " donors). Effect size:\n", sep = "")
      cat("    auc = P(a random PHF1+ cell > a random PHF1- cell),\n")
      cat("    which is independent of n; rank_biserial = 2*auc - 1. auc 0.5 = no separation.\n")
      print(as.data.frame(con_stats %>% dplyr::transmute(
        module, W = wilcox_cell_W, auc = round(wilcox_cell_auc, 4),
        rank_biserial = round(wilcox_cell_rank_biserial, 4),
        p = wilcox_cell_p, padj = wilcox_cell_padj,
        significant = wilcox_cell_padj < SIG_ALPHA)))

      cat("\n(b) DONOR-level signed-rank on the per-donor paired differences (n = ",
          con_stats$n_donors_paired[1], " donors).\n", sep = "")
      cat("    The non-parametric analogue of dz.\n")
      cat("    NOTE the floor: with n = ", con_stats$n_donors_paired[1],
          " the smallest attainable two-sided p is ",
          signif(2 / 2^con_stats$n_donors_paired[1], 3), ",\n", sep = "")
      cat("    so this test CANNOT go below it however large the effect. Read it together with\n")
      cat("    n_donors_same_dir above, not on its own.\n")
      print(as.data.frame(con_stats %>% dplyr::transmute(
        module, V = wilcox_donor_V,
        rank_biserial = round(wilcox_donor_rank_biserial, 4),
        p = wilcox_donor_p, padj = wilcox_donor_padj,
        significant = wilcox_donor_padj < SIG_ALPHA)))

      failed <- setdiff(ref_df$module, con_stats$module)
      if (length(failed))
        cat("\n  [warn] contrast could not be fitted (drawn hollow): ",
            paste(failed, collapse = ", "), "\n", sep = "")
    } else {
      cat("  [warn] no contrast could be fitted; every point is drawn hollow.\n")
    }
  } else if (!is_neuron_celltype(ct)) {
    cat("\nNo PHF1+ reference strip: ", ct, " is not a neuron celltype. PHF1 is a neuronal\n", sep = "")
    cat("  post-stain, so the strip is gated on the celltype label rather than on whether any\n")
    cat("  cell is flagged PHF1+ -- a few PHF1+ glia would be segmentation spillover, and a\n")
    cat("  strip drawn off them would be an artefact.\n")
  } else {
    cat("\nNo PHF1+ cells in this neuron celltype -- no reference strip.\n")
  }

  cat("\nNOTE (gene sets): these Reactome sets are EXTERNAL -- unlike the PHF1 marker set\n")
  cat("  they are not selected from the PHF1+/PHF1- contrast, so the PHF1+ strip is not\n")
  cat("  elevated by construction.\n")
  cat("NOTE (scale): AddModuleScore values from different gene sets are not on the same\n")
  cat("  absolute scale. Compare SHAPES across modules within a panel, not heights.\n")

  cat("\nsessionInfo():\n"); print(sessionInfo())
  sink()

  lmm_coef_rows[[paste(ct, group)]] <<- mod_stats

  cat(sprintf("  Wrote %s / %s: %d module(s), %d significant\n",
              ct, group, nrow(mod_stats), sum(mod_stats$significant %in% TRUE)))
  invisible(mod_stats)
}

##  ............................................................................
##  Run everything                                                          ####
for (ct in CELLTYPES) {
  mdf <- prep_model_df(ct)
  cat(sprintf("\n#### %s: %d PHF1-negative cells within cap (%d PHF1+ excluded) | dist_sd = %.4f\n",
              ct, nrow(mdf$df), mdf$n_phf1_pos, mdf$dist_sd))
  for (group in names(MODULE_GROUPS)) run_group(ct, group, mdf)
}

# Write every queued distance panel at one common panel height / x scale.
flush_panels()

if (length(lmm_coef_rows))
  write.table(bind_rows(lmm_coef_rows),
              file.path(args$output_dir, "stats_reactome_lmm_coeffs.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nAll outputs written to:", args$output_dir, "\n")
cat("Done.\n")
