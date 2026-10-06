#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_modulescore_vs_phf1_distance_modelp.R
#
# Figure panels: 6E
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# plot_modulescore_vs_phf1_distance_modelp.R
#
# Glial module scores vs distance to the nearest PHF1+ neuron. Significance / line
# annotations come from the MODEL's BH-adjusted p-value (lmm_padj for the LMMs; BH of
# the smooth p-value for the GAM). Modelled cells are capped at dist_to_phf1_um <=
# MAX_DIST (set by CAP_UM in run_modulescore_distance_modelp.sh).
#
# WHAT IT DOES
#   * Loads the gene-set collections (qs::qread; the Pandey oligodendrocyte sets are read from CSV):
#       Mancuso_2024 — 8 microglial states (HM, tCRM, CRM1, CRM2, RM, DAM, HLA, IRM)
#       Cameron_2024 — 8 reactive-astrocyte subclusters (2 neurotoxic, 6 neuroprotective)
#   * Maps identity-matched: Mancuso -> Micro only; Cameron -> Astro only.
#   * Mancuso: 8 states scored separately. Cameron: 2 POOLED meta-signatures
#     (Neurotoxic, Neuroprotective; full published lists) AND the 8 subclusters.
#   * Overlaps EVERY signature with the CosMx 6k panel before scoring; logs overlap.
#   * Scores with AddModuleScore on the SCT data slot, per cell-type subset.
#   * Models each module score vs distance, capped at MAX_DIST:
#       log     : score ~ log(dist)_scaled + nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)
#                 -- THE REPORTED FIT
#       decay   : score = P_j + A_j*exp(-dist/lambda_j) + nUMI_log + percent_neg,
#                 fitted PER DONOR by profile least squares, lambda meta-analysed
#     The "linear" and "gam" types are not run by default (only the log fit is used for
#     reporting). Their code paths in run_group_model() are intact; enable them by adding
#     them to the `for (mt in ...)` vector in the driver at the bottom.
#   * PER-DONOR VIEW: one faceted figure per group (a facet per module, a line per donor,
#     coloured by Braak), model-free, alongside the pooled rolling mean. See
#     R/rollmean_by_donor.R.
#   * Significance = BH-adjusted model p within each group (lmm_padj). Solid/bold
#     lines = padj<0.05, thin/dashed = ns.
#   * Distance: log (reported) AND decay; the linear and gam paths are available.
#
# LENGTH CONSTANT (model type "decay"). For every module the LOG model calls
# significant, the asymptotic decay above gives the length constant lambda in MICRONS
# with a CI, plus d95 = 2.996*lambda -- the distance at which 95% of the excess has
# resolved. MICROGLIA AND ASTROCYTES ONLY (oligodendrocytes are out of scope for this
# model type; they are still fitted by linear/log/gam). See R/decay_length_utils.R for
# the estimator, its gates (the decay must beat the log-linear model by a
# margin, lambda must be resolved inside the window, its CI must be narrow enough to
# quote, and it must not track the distance cap) and the deliberate use of RAW MICRONS
# rather than the canonical log/SD distance transform. The log model remains primary
# for "is there a gradient?"; decay answers "how far does it reach?".
#
# OUTPUT (per group x model), under the --output_dir:
#   plot_modulescore_dist_<group>_<CT>[_log|_decay].pdf             — coloured-line figure
#   plot_modulescore_bydonor_dist_<group>_<CT>.pdf                  — per-donor lines, Braak-coloured
#   source_data_modulescore_bydonor_dist_<group>_<CT>_rollmean.tsv  — exact drawn rows
#   stats_modulescore_bydonor_dist_<group>_<CT>.txt                 — window, donor table, caveats
#   source_data_modulescore_dist_<group>_<CT>[...]_fit.tsv          — EXACT drawn rows (fitted curves + CI)
#   source_data_modulescore_dist_<group>_<CT>[...]_persample.tsv    — per-donor x ring means
#   source_data_modulescore_dist_<group>_<CT>_decay_donors.tsv      — per-donor lambda (incl. exclusions)
#   stats_modulescore_dist_<group>_<CT>[...].txt                    — formulas, per-module summaries, sessionInfo
# Structured supplementary tables (written once):
#   stats_modulescore_geneset_overlap.tsv, stats_modulescore_lmm_coeffs.tsv,
#   stats_modulescore_decay_lambda.tsv                              — the headline lambda table
#   plot_decay_lambda_forest.pdf, source_data_decay_lambda_forest.tsv
#
# Use run_modulescore_distance_modelp.sh

##  ............................................................................
##  Packages + setup                                                        ####
suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(lme4)
  library(lmerTest)
  library(mgcv)
  library(argparse)
})

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")             # fig_theme, forest_theme, celltype_palette, ...
source("R/decay_length_utils.R")   # decay_analyse() and friends (model type "decay")
source("R/rollmean_by_donor.R")    # per-donor rolling means, coloured by Braak

set.seed(42)

##  ............................................................................
##  Arguments                                                               ####
parser <- ArgumentParser()
parser$add_argument("--seu",
  default = "seu_PHF1.rds",
  help = "Path to seu_PHF1.rds [default: seu_PHF1.rds]")
parser$add_argument("--geneset_dir",
  default = "<EXTERNAL_GENESET_DIR>",
  help = "Directory holding Mancuso_2024.qs and Cameron_2024.qs")
parser$add_argument("--output_dir",
  default = "plots/module_score_vs_phf1_distance_modelp_1000um",
  help = "Output directory")
parser$add_argument("--seed",    type = "integer", default = 42,
  help = "Seed (AddModuleScore + provenance) [default: 42]")
parser$add_argument("--max_dist_um", type = "double", default = 1000,
  help = "If > 0, restrict modelled cells to dist_to_phf1_um <= this (um); 0 = no cap")
parser$add_argument("--pandey_up_csv", default = "Pandey_hOligo2_UP_highPathology.csv",
  help = "CSV of UP genes in the disease-associated oligo (hOligo2) cluster; column 'gene'")
parser$add_argument("--pandey_dn_csv", default = "Pandey_hOligo2_DOWN_highPathology.csv",
  help = "CSV of DOWN genes in the disease-associated oligo (hOligo2) cluster; column 'gene'")
parser$add_argument("--window_um", type = "double", default = 75,
  help = "Rolling-mean window WIDTH (um) for the raw model-free overlay plots [default: 75]")
parser$add_argument("--decay_lambda_min", type = "double", default = DECAY_LAMBDA_MIN_UM,
  help = "Decay model: lower bound of the lambda search (um) [default: 5]")
parser$add_argument("--decay_lambda_max_mult", type = "double", default = DECAY_LAMBDA_MAX_MULT,
  help = "Decay model: lambda searched up to MULT * distance cap [default: 3]")
parser$add_argument("--decay_min_cells_donor", type = "integer", default = DECAY_MIN_CELLS_DONOR,
  help = "Decay model: min cells for a donor to contribute a lambda [default: 100]")
parser$add_argument("--decay_seg_baseline", default = "yes", choices = c("yes", "no"),
  help = paste("Decay model: use the stably-expressed-gene (SEG) null as the baseline,",
               "by entering the per-cell SEG module score as a covariate so the fitted",
               "baseline is P_j + b*SEG(d) rather than a flat plateau [default: yes]"))
parser$add_argument("--seg_qs", default = "seg_control/Reference_SEG_panel.qs",
  help = paste("SEG geneset .qs (named list, one element per celltype) used for the",
               "decay baseline; from R/prepare_seg_geneset_reference.R",
               "[default: seg_control/Reference_SEG_panel.qs]"))
args <- parser$parse_args()

MAX_DIST <- if (args$max_dist_um > 0) args$max_dist_um else Inf
if (is.finite(MAX_DIST))
  cat(sprintf("Distance cap: modelling cells with dist_to_phf1_um <= %g um\n", MAX_DIST))

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

# Constants -------------------------------------------------------------------
MIN_GENES_PRESENT <- 10                         # drop signatures with fewer panel genes
CTRL              <- 100
NBIN              <- 24
RING_BREAKS       <- c(0, 50, 100, 200, 300, 500, 700)
SIG_ALPHA         <- 0.05                        # significance threshold (model FDR)
GAM_K             <- 5                           # basis dim for the GAM distance smooth
HALF_WINDOW       <- args$window_um / 2          # rolling-mean half-window (um) for raw plots

# Identity-matched mapping
CONFIG <- list(Micro = "Mancuso_2024", Astro = "Cameron_2024")

# Distinct qualitative palette for submodules (NOT celltype_palette).
# Colourblind-friendly (Okabe-Ito, extended with Tol olive/indigo). Ordered so the
# 8-module glia panels (Mancuso Micro, Cameron subclusters) use positions 1-8, which
# avoid the low-contrast pale yellow (kept at position 9 for line legibility on white).
MODULE_COLOURS <- c(
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#CC79A7",  # reddish purple
  "#000000",  # black
  "#999933",  # olive (Tol)
  "#F0E442",  # yellow
  "#332288"   # indigo (Tol)
)
# Up/toxic = vermillion, down/protective = blue (Okabe-Ito; colourblind-safe red/blue).
MERGED_COLOURS <- c(Neurotoxic = "#D55E00", Neuroprotective = "#0072B2",
                    "hOligo2 Up" = "#D55E00", "hOligo2 Down" = "#0072B2")  # Pandey oligo modules

# Ring labels / midpoints
.ring_lower <- head(RING_BREAKS, -1)
.ring_upper <- tail(RING_BREAKS, -1)
RING_LAB    <- sprintf("%g-%g", .ring_lower, .ring_upper)
RING_MID    <- setNames((.ring_lower + .ring_upper) / 2, RING_LAB)

##  ............................................................................
##  Helpers                                                                 ####

# First present column name from a set of candidates (defensive against naming).
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
    collection   = collection,
    group        = group,
    set          = names(sets),
    n_total      = as.integer(n_total),
    n_present    = as.integer(n_pres),
    frac_present = round(n_pres / n_total, 4),
    scored       = scored
  )
  list(present = present[scored], overlap = overlap)
}

# AddModuleScore for one prefix; returns data.frame of scores renamed to set names.
# Retries with a smaller nbin if the control-gene binning fails.
score_group <- function(seu_ct, present_sets, prefix, seed) {
  set.seed(seed)
  seu_ct <- tryCatch(
    AddModuleScore(seu_ct, features = present_sets, name = prefix,
                   assay = "SCT", ctrl = CTRL, nbin = NBIN, seed = seed),
    error = function(e) {
      cat("    AddModuleScore failed at nbin=", NBIN, " (", conditionMessage(e),
          ") — retrying at nbin=15\n", sep = "")
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

# Subset Seurat to a celltype, score the supplied feature groups, return one
# data.frame with cell_id + covariates + all module-score columns.
extract_celltype_data <- function(seu, ct, feature_groups, seed,
                                  sex_col, age_col, pmi_col) {
  seu_ct <- subset(seu, celltype == ct)
  DefaultAssay(seu_ct) <- "SCT"
  md <- seu_ct@meta.data
  cat(sprintf("  %s cells: %d\n", ct, nrow(md)))

  id_cands <- c("cell_ID", "cell_id", "CellID")
  id_col   <- id_cands[id_cands %in% colnames(md)]
  cell_id  <- if (length(id_col) > 0) as.character(md[[id_col[1]]]) else rownames(md)

  base <- data.frame(
    cell_id         = cell_id,
    dist_to_phf1_um = as.numeric(md$dist_to_phf1_um),
    percent_neg     = as.numeric(md$percent.neg),
    nUMI_log        = log2(as.numeric(md$nCount_RNA) + 1),
    sample_id       = as.character(md$sample_id),
    # Braak is carried only to COLOUR the per-donor rolling-mean lines; it is not a
    # covariate in any model here. Tolerated as missing so a plotting nicety can
    # never abort the analysis -- the figure falls back to colouring by donor.
    Braak           = if ("Braak" %in% colnames(md)) as.character(md$Braak) else NA_character_,
    Sex             = as.character(md[[sex_col]]),
    Age             = as.numeric(md[[age_col]]),
    PMI             = as.numeric(md[[pmi_col]]),
    stringsAsFactors = FALSE
  )

  scores <- do.call(cbind, lapply(names(feature_groups), function(prefix)
    score_group(seu_ct, feature_groups[[prefix]], prefix, seed)))

  cbind(base, scores)   # row-aligned: all derived from seu_ct@meta.data order
}

# Observed LMM for one module. Returns fit stats + the fitted model.
fit_obs_lmer <- function(sub) {
  m <- tryCatch(
    lmerTest::lmer(score ~ dist_scaled + nUMI_log + percent_neg + Sex + Age_s + PMI_s +
                     (1 | sample_id), data = sub, REML = TRUE),
    error = function(e) { message("    obs lmer failed: ", conditionMessage(e)); NULL })
  if (is.null(m)) return(NULL)
  cf <- coef(summary(m))
  if (!"dist_scaled" %in% rownames(cf)) return(NULL)
  list(model = m,
       estimate = cf["dist_scaled", "Estimate"],
       se       = cf["dist_scaled", "Std. Error"],
       t        = cf["dist_scaled", "t value"],
       df       = cf["dist_scaled", "df"],
       p        = cf["dist_scaled", "Pr(>|t|)"],
       singular = isSingular(m))
}

# Fitted module-score trend across a distance grid, with Wald CI from the
# fixed-effect vcov (RE variance ignored). Covariates held at their means /
# the reference Sex level. grid_um is raw um; transformed to match the model.
predict_grid <- function(model, df, dist_sd, log_distance, x_cap) {
  lo      <- if (log_distance) min(df$dist_to_phf1_um, na.rm = TRUE) else 0
  grid_um <- seq(lo, x_cap, length.out = 200)
  dt      <- if (log_distance) log(grid_um) else grid_um
  nd <- data.frame(
    dist_scaled = dt / dist_sd,
    nUMI_log    = mean(df$nUMI_log),
    percent_neg = mean(df$percent_neg),
    Sex         = factor(levels(df$Sex)[1], levels = levels(df$Sex)),
    Age_s       = mean(df$Age_s),
    PMI_s       = mean(df$PMI_s)
  )
  fe_form <- ~ dist_scaled + nUMI_log + percent_neg + Sex + Age_s + PMI_s
  X    <- model.matrix(fe_form, nd)
  beta <- lme4::fixef(model)
  X    <- X[, names(beta), drop = FALSE]
  fit  <- as.numeric(X %*% beta)
  V    <- as.matrix(vcov(model))
  se   <- sqrt(rowSums((X %*% V) * X))
  data.frame(dist_to_phf1_um = grid_um, fitted_score = fit,
             ci_lo = fit - 1.96 * se, ci_hi = fit + 1.96 * se)
}

# ---- GAM path (mgcv::bam) -------------------------------------------------
# Penalized smooth on RAW distance + donor random-effect smooth: a data-driven
# nonlinear alternative to the linear/log LMM. Significance here is the smooth's
# parametric p (BH-adjusted).

# Pull (edf, F-statistic, p) for the distance smooth from a fitted bam.
.gam_dist_stat <- function(m) {
  st <- summary(m)$s.table
  rn <- grep("s\\(dist", rownames(st))
  if (length(rn) == 0) return(NULL)
  fcol <- if ("F" %in% colnames(st)) "F" else "Chi.sq"
  list(edf = st[rn[1], "edf"], stat = st[rn[1], fcol], p = st[rn[1], "p-value"])
}

fit_obs_gam <- function(sub) {
  m <- tryCatch(
    mgcv::bam(score ~ s(dist, k = GAM_K) + nUMI_log + percent_neg + Sex + Age_s + PMI_s +
                s(sample_id, bs = "re"),
              data = sub, method = "fREML", discrete = TRUE),
    error = function(e) { message("    obs gam failed: ", conditionMessage(e)); NULL })
  if (is.null(m)) return(NULL)
  s <- .gam_dist_stat(m)
  if (is.null(s)) return(NULL)
  list(model = m, edf = s$edf, stat = s$stat, p = s$p)
}

# Fitted smooth across the distance grid (population level: donor RE excluded),
# covariates at their means / reference Sex. Same output schema as predict_grid.
predict_grid_gam <- function(model, df, x_cap) {
  grid_um <- seq(0, x_cap, length.out = 200)
  nd <- data.frame(
    dist        = grid_um,
    nUMI_log    = mean(df$nUMI_log),
    percent_neg = mean(df$percent_neg),
    Sex         = factor(levels(df$Sex)[1], levels = levels(df$Sex)),
    Age_s       = mean(df$Age_s),
    PMI_s       = mean(df$PMI_s),
    sample_id   = factor(levels(df$sample_id)[1], levels = levels(df$sample_id))
  )
  pr  <- predict(model, newdata = nd, se.fit = TRUE, exclude = "s(sample_id)")
  fit <- as.numeric(pr$fit); se <- as.numeric(pr$se.fit)
  data.frame(dist_to_phf1_um = grid_um, fitted_score = fit,
             ci_lo = fit - 1.96 * se, ci_hi = fit + 1.96 * se)
}

# Model-free rolling (sliding-window) mean of y over x, evaluated on a grid.
# Returns mean, +/- 1 SEM, and window cell count at each grid point (95% CI = +/- 1.96 SEM).
roll_mean <- function(x, y, grid, half_window) {
  out <- lapply(grid, function(g) {
    idx <- which(x >= g - half_window & x <= g + half_window); n <- length(idx)
    if (n < 1) return(c(roll_mean = NA_real_, sem = NA_real_, n_window = 0))
    m <- mean(y[idx]); s <- if (n > 1) stats::sd(y[idx]) / sqrt(n) else NA_real_
    c(roll_mean = m, sem = s, n_window = n)
  })
  as.data.frame(do.call(rbind, out))
}

##  ............................................................................
##  Load gene sets, build pooled signatures                                 ####
cat("Loading gene sets...\n")
mancuso <- qs::qread(file.path(args$geneset_dir, "Mancuso_2024.qs"))
cameron <- qs::qread(file.path(args$geneset_dir, "Cameron_2024.qs"))
cat("  Mancuso sets:", paste(names(mancuso), collapse = ", "), "\n")
cat("  Cameron sets:", paste(names(cameron), collapse = ", "), "\n")

nt_names <- grep("^Neurotoxic",      names(cameron), value = TRUE)
np_names <- grep("^Neuroprotective", names(cameron), value = TRUE)
stopifnot("No Neurotoxic/Neuroprotective Cameron sets found" =
            length(nt_names) > 0 && length(np_names) > 0)
cameron_pooled <- list(
  Neurotoxic      = unique(unlist(cameron[nt_names], use.names = FALSE)),
  Neuroprotective = unique(unlist(cameron[np_names], use.names = FALSE))
)
cameron_sub <- cameron   # 8 subclusters
names(cameron_sub) <- gsub("_cluster", "", names(cameron_sub))   # compact labels: Neurotoxic_cluster0 -> Neurotoxic0

# Pandey disease-associated oligodendrocyte (hOligo2, high pathology) DEG sets.
pandey <- list(
  `hOligo2 Up`   = read.csv(args$pandey_up_csv, stringsAsFactors = FALSE)$gene,
  `hOligo2 Down` = read.csv(args$pandey_dn_csv, stringsAsFactors = FALSE)$gene
)
stopifnot("Pandey CSVs must have a 'gene' column" =
            length(pandey[["hOligo2 Up"]]) > 0 && length(pandey[["hOligo2 Down"]]) > 0)
cat("  Pandey hOligo2: Up", length(pandey[["hOligo2 Up"]]), "genes, Down",
    length(pandey[["hOligo2 Down"]]), "genes\n")

##  ............................................................................
##  Load Seurat object, determine panel, score both celltypes               ####
cat("Loading Seurat object (this is large)...\n")
seu <- readRDS(args$seu)
DefaultAssay(seu) <- "SCT"

stopifnot("celltype column missing"        = "celltype"        %in% colnames(seu@meta.data),
          "dist_to_phf1_um column missing" = "dist_to_phf1_um" %in% colnames(seu@meta.data),
          "percent.neg column missing"     = "percent.neg"     %in% colnames(seu@meta.data),
          "nCount_RNA column missing"      = "nCount_RNA"      %in% colnames(seu@meta.data),
          "sample_id column missing"       = "sample_id"       %in% colnames(seu@meta.data))

sex_col <- pick_col(seu@meta.data, c("Sex"), "Sex")
age_col <- pick_col(seu@meta.data, c("Age"), "Age")
pmi_col <- pick_col(seu@meta.data, c("PMI", "PostMortemInterval", "PostMortem_Interval"), "PMI")
cat(sprintf("Donor covariate columns: Sex=%s Age=%s PMI=%s\n", sex_col, age_col, pmi_col))

panel_genes <- rownames(seu[["SCT"]])
cat("CosMx panel genes (SCT assay):", length(panel_genes), "\n")

mancuso_prep <- prep_sets(mancuso,        panel_genes, MIN_GENES_PRESENT, "Mancuso_2024", "mancuso")
campool_prep <- prep_sets(cameron_pooled, panel_genes, MIN_GENES_PRESENT, "Cameron_2024", "cameron_pooled")
camsub_prep  <- prep_sets(cameron_sub,    panel_genes, MIN_GENES_PRESENT, "Cameron_2024", "cameron_subclusters")
pandey_prep  <- prep_sets(pandey,         panel_genes, MIN_GENES_PRESENT, "Pandey_hOligo2", "pandey_oligo")

overlap_tab <- bind_rows(mancuso_prep$overlap, campool_prep$overlap, camsub_prep$overlap, pandey_prep$overlap)
write.table(overlap_tab, file.path(args$output_dir, "stats_modulescore_geneset_overlap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Gene-set / panel overlap:\n"); print(as.data.frame(overlap_tab))

# Save EVERY gene list used (full sets, with an on_panel flag = the genes actually scored).
members_df <- function(sets, group) bind_rows(lapply(names(sets), function(m)
  tibble(group = group, module = m, gene = sets[[m]], on_panel = sets[[m]] %in% panel_genes)))
geneset_members <- bind_rows(
  members_df(mancuso,        "mancuso"),
  members_df(cameron_pooled, "cameron_pooled"),
  members_df(cameron_sub,    "cameron_subclusters"),
  members_df(pandey,         "pandey_oligo"))
write.table(geneset_members, file.path(args$output_dir, "stats_modulescore_geneset_members.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved gene-list members (on_panel = scored):",
    sum(geneset_members$on_panel), "of", nrow(geneset_members), "genes on panel\n")

##  ............................................................................
##  SEG null baseline for the decay model                                   ####
# Stably expressed genes (scSEGIndex, Lin et al. 2019), per celltype, from
# R/prepare_seg_geneset_reference.R -- the same sets as the negative-control script
# R/modulescore_seg_control_vs_phf1_distance.R, so the two agree by construction.
# Scored as its own feature_group (its own AddModuleScore call), so the module
# scores above stay bit-identical to a run without it. Used ONLY by the decay type.
USE_SEG     <- identical(args$decay_seg_baseline, "yes")
SEG_PRESENT <- list()
if (USE_SEG) {
  if (!file.exists(args$seg_qs))
    stop("SEG geneset not found: ", args$seg_qs,
         "\n  Run R/prepare_seg_geneset_reference.R, or pass --decay_seg_baseline no.")
  seg <- qs::qread(args$seg_qs)
  stopifnot("SEG geneset must be a named list of character vectors" =
              is.list(seg) && !is.null(names(seg)) && all(vapply(seg, is.character, logical(1))))
  cat("SEG baseline sets:", paste(names(seg), collapse = ", "), "\n")
  seg_overlap <- list()
  for (ct in c("Micro", "Astro")) {            # the two celltypes the decay type covers
    if (!ct %in% names(seg)) {
      cat("  [warn] no SEG set for ", ct, "; free-plateau fallback\n", sep = ""); next
    }
    pr <- prep_sets(setNames(list(seg[[ct]]), "seg_baseline"), panel_genes,
                    MIN_GENES_PRESENT, "Reference_scSEGIndex", ct)
    seg_overlap[[ct]] <- pr$overlap %>% mutate(celltype = ct)
    if (length(pr$present) > 0) {
      SEG_PRESENT[[ct]] <- pr$present[["seg_baseline"]]
      cat(sprintf("  %s: %d SEG genes on panel\n", ct, length(SEG_PRESENT[[ct]])))
    } else cat("  [warn] ", ct, ": too few SEG genes on panel; free-plateau fallback\n", sep = "")
  }
  if (length(seg_overlap) > 0)
    write.table(bind_rows(seg_overlap) %>%
                  dplyr::select(celltype, collection, set, n_total, n_present, frac_present, scored),
                file.path(args$output_dir, "stats_modulescore_decay_seg_overlap.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
}
seg_fg <- function(ct) if (is.null(SEG_PRESENT[[ct]])) list() else
  list(SEGBL_ = list(seg_baseline = SEG_PRESENT[[ct]]))

cat("Scoring Micro (Mancuso)...\n")
micro_df <- extract_celltype_data(
  seu, "Micro",
  feature_groups = c(list(Mancuso_ = mancuso_prep$present), seg_fg("Micro")),
  seed = args$seed, sex_col = sex_col, age_col = age_col, pmi_col = pmi_col)

cat("Scoring Astro (Cameron pooled + subclusters)...\n")
astro_df <- extract_celltype_data(
  seu, "Astro",
  feature_groups = c(list(CameronPooled_ = campool_prep$present,
                          CameronSub_    = camsub_prep$present), seg_fg("Astro")),
  seed = args$seed, sex_col = sex_col, age_col = age_col, pmi_col = pmi_col)

cat("Scoring Oligo (Pandey hOligo2 UP/DOWN)...\n")
oligo_df <- extract_celltype_data(
  seu, "Oligo",
  feature_groups = list(Pandey_ = pandey_prep$present),
  seed = args$seed, sex_col = sex_col, age_col = age_col, pmi_col = pmi_col)

rm(seu); invisible(gc())

##  ............................................................................
##  Prepare per-celltype model frames (PHF1-negative, scaled covariates)    ####
prep_celltype_df <- function(df) {
  df <- df[!is.na(df$dist_to_phf1_um) & !is.na(df$percent_neg) &
             !is.na(df$nUMI_log) & !is.na(df$Age) & !is.na(df$PMI) & !is.na(df$Sex), ]
  df$Sex       <- droplevels(factor(df$Sex))
  df$sample_id <- factor(df$sample_id)
  df
}
micro_df <- prep_celltype_df(micro_df)
astro_df <- prep_celltype_df(astro_df)
oligo_df <- prep_celltype_df(oligo_df)
cat(sprintf("Modelled PHF1-negative cells (pre-cap): Micro=%d, Astro=%d, Oligo=%d\n",
            nrow(micro_df), nrow(astro_df), nrow(oligo_df)))

##  ............................................................................
##  Group definitions                                                       ####
GROUPS <- list(
  list(group = "mancuso",             ct = "Micro", df_name = "micro",
       modules = names(mancuso_prep$present), merged = FALSE),
  list(group = "cameron_pooled",      ct = "Astro", df_name = "astro",
       modules = intersect(c("Neurotoxic", "Neuroprotective"), names(campool_prep$present)),
       merged = TRUE),
  list(group = "cameron_subclusters", ct = "Astro", df_name = "astro",
       modules = names(camsub_prep$present), merged = FALSE),
  list(group = "pandey_oligo",        ct = "Oligo", df_name = "oligo",
       modules = names(pandey_prep$present), merged = TRUE)
)
DF_BY_NAME <- list(micro = micro_df, astro = astro_df, oligo = oligo_df)

# Accumulators
lmm_coef_rows <- list()
LOG_PANELS    <- list()   # captured log-model ggplots for the combined 3-panel A4 figure
decay_rows    <- list()   # headline lambda table (one row per group x module)
decay_donors  <- list()   # per-donor lambdas, for the forest's background points
decay_colours <- c()      # module -> colour, for the forest

##  ............................................................................
##  Per-group x model analysis (model-FDR significance; NO permutation)     ####
run_group_model <- function(gdef, model_type) {
  # The decay model type shares nothing with the lmer/gam paths (per-donor profile
  # least squares + meta-analysis, raw-micron distance), so it dispatches out here.
  if (model_type == "decay") return(run_group_decay(gdef))
  group <- gdef$group; ct <- gdef$ct; modules <- gdef$modules; merged <- gdef$merged
  df_full <- DF_BY_NAME[[gdef$df_name]]
  is_gam       <- model_type == "gam"
  log_distance <- model_type == "log"
  suffix      <- switch(model_type, linear = "",   log = "_log",   gam = "_gam")
  scale_label <- switch(model_type, linear = "um", log = "log_um", gam = "gam")
  model_label <- switch(model_type,
                        linear = "LMM (linear distance)",
                        log    = "LMM (log distance)",
                        gam    = sprintf("GAM (penalized smooth, mgcv::bam, k=%d)", GAM_K))
  formula_str <- if (is_gam)
    "score ~ s(dist_to_phf1_um, k=K) + nUMI_log + percent_neg + Sex + Age + PMI + s(sample_id, bs='re')"
  else
    "score ~ dist_to_phf1_um_scaled + nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)"
  cat(sprintf("\n=== %s / %s  (model: %s%s) ===\n", group, ct, scale_label,
              if (is.finite(MAX_DIST)) sprintf(", <=%gum", MAX_DIST) else ""))

  # Scale donor covariates on the full PHF1-negative set
  df_full$Age_s <- as.numeric(scale(df_full$Age))
  df_full$PMI_s <- as.numeric(scale(df_full$PMI))

  # Distance cap: restrict the modelled set to cells within MAX_DIST.
  # (Module scores are computed on the full celltype population; the cap defines the
  #  modelled distance window, not the AddModuleScore control background.)
  df <- if (is.finite(MAX_DIST)) df_full[df_full$dist_to_phf1_um <= MAX_DIST, ] else df_full

  # Distance handling: lmer paths scale (optionally log); GAM smooths raw distance.
  if (is_gam) {
    dist_sd <- NA_real_
  } else {
    dt      <- if (log_distance) log(df$dist_to_phf1_um) else df$dist_to_phf1_um
    dist_sd <- sd(dt)
    df$dist_scaled <- dt / dist_sd
  }
  x_cap <- if (is.finite(MAX_DIST)) MAX_DIST else as.numeric(quantile(df$dist_to_phf1_um, 0.99))

  fit_list <- list(); mod_stats <- list()

  for (m in modules) {
    if (is_gam) {
      sub <- df[, c("cell_id", "nUMI_log", "percent_neg", "Sex", "Age_s", "PMI_s", "sample_id")]
      sub$dist  <- df$dist_to_phf1_um
      sub$score <- df[[m]]
      obs <- fit_obs_gam(sub)
      if (is.null(obs)) { cat("  ", m, ": GAM failed, skipped\n"); next }
      grid <- predict_grid_gam(obs$model, df, x_cap)
      nrm  <- list(slope_per_sd = NA_real_, se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                   slope_per_unit = NA_real_, stat = obs$stat, dfv = obs$edf,
                   p = obs$p, singular = NA)
    } else {
      sub <- df[, c("cell_id", "dist_scaled", "nUMI_log", "percent_neg",
                    "Sex", "Age_s", "PMI_s", "sample_id")]
      sub$score <- df[[m]]
      obs <- fit_obs_lmer(sub)
      if (is.null(obs)) { cat("  ", m, ": LMM failed, skipped\n"); next }
      grid <- predict_grid(obs$model, df, dist_sd, log_distance, x_cap)
      nrm  <- list(slope_per_sd = obs$estimate, se = obs$se,
                   ci_lo = obs$estimate - 1.96 * obs$se,
                   ci_hi = obs$estimate + 1.96 * obs$se,
                   slope_per_unit = obs$estimate / dist_sd, stat = obs$t,
                   dfv = obs$df, p = obs$p, singular = obs$singular)
    }

    grid$module <- m
    fit_list[[m]] <- grid

    mod_stats[[m]] <- tibble(
      group = group, celltype = ct, module = m, distance_scale = scale_label,
      slope_per_sd = nrm$slope_per_sd, se = nrm$se,
      ci_lo_per_sd = nrm$ci_lo, ci_hi_per_sd = nrm$ci_hi,
      slope_per_unit = nrm$slope_per_unit,
      t = nrm$stat, df = nrm$dfv, p_LMM = nrm$p,
      singular = nrm$singular, dist_sd = dist_sd
    )
  }

  if (length(mod_stats) == 0) { cat("  No modules modelled; skipping group.\n"); return(invisible(NULL)) }

  mod_stats <- bind_rows(mod_stats)
  # BH within this group's module set -> model FDR (the significance source).
  mod_stats$lmm_padj    <- p.adjust(mod_stats$p_LMM, method = "BH")
  mod_stats$significant <- mod_stats$lmm_padj < SIG_ALPHA
  mod_stats$verdict     <- ifelse(is.na(mod_stats$significant), "not testable",
                           ifelse(mod_stats$significant,
                                  "significant distance effect (model FDR<0.05)",
                                  "n.s. (model FDR)"))

  sig_basis <- "padj<0.05"

  # ---- assemble fit data frame for plotting / source data ----
  fit_df <- bind_rows(fit_list) %>%
    left_join(dplyr::select(mod_stats, module, significant, p_LMM, lmm_padj), by = "module") %>%
    mutate(module = factor(module, levels = modules),
           sig = factor(ifelse(significant, sig_basis, "ns"), levels = c(sig_basis, "ns")))

  # ---- coloured-line figure ----
  cols <- if (merged) MERGED_COLOURS[modules]
          else setNames(MODULE_COLOURS[seq_along(modules)], modules)
  xlab <- expression("Distance to PHF1+ neuron (" * mu * "m)")   # mu via plotmath (survives default pdf device)
  subt <- NULL   # no model-descriptor title on the plots
  # module-legend title: gene-set source + celltype
  legend_title <- switch(group,
                         mancuso             = "Mancuso Micro",
                         cameron_pooled      = "Cameron Astro",
                         cameron_subclusters = "Cameron Astro",
                         pandey_oligo        = "Pandey Oligo",
                         NULL)

  p <- ggplot(fit_df, aes(dist_to_phf1_um, fitted_score, colour = module, group = module))
  # 95% CI shading (fitted-curve Wald CI) on every plot, behind the lines.
  p <- p + geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = module),
                       alpha = 0.15, colour = NA)
  p <- p +
    geom_line(aes(linetype = sig, linewidth = sig)) +
    scale_colour_manual(values = cols, name = legend_title, drop = FALSE) +
    # linetype + linewidth share one legend (same name + breaks) so it always shows
    # BOTH a bold-solid key (significant) and a thin-dashed key (ns), even when only
    # one is present in the data (drop = FALSE keeps the absent level).
    scale_linetype_manual(values = setNames(c("solid", "dashed"), c(sig_basis, "ns")),
                          name = NULL, drop = FALSE, limits = c(sig_basis, "ns")) +
    scale_linewidth_manual(values = setNames(c(0.8, 0.4), c(sig_basis, "ns")),
                           name = NULL, drop = FALSE, limits = c(sig_basis, "ns")) +
    # fixed legend order (else ggplot flips it between plots): module legend on
    # top, the significance (padj) legend below, consistently.
    guides(colour = guide_legend(order = 1),
           linetype = guide_legend(order = 2), linewidth = guide_legend(order = 2)) +
    labs(x = xlab, y = "Module score (fitted)", subtitle = subt) +
    coord_cartesian(xlim = c(0, x_cap)) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.text = element_text(size = 6),
          # widen the legend key so the dashed (ns) linetype shows >=2 dashes
          legend.key.width = grid::unit(20, "pt"),
          # cut the panel-to-legend gap (~11pt default -> 3pt) and legend padding
          legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
          # trim fig_theme's large left margin (was for angled celltype x-labels)
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
  p <- p + scale_fill_manual(values = cols, guide = "none")

  # capture the log-model plot for the combined 3-panel A4 figure (Micro|Astro|Oligo)
  if (model_type == "log" && group %in% c("mancuso", "cameron_subclusters", "pandey_oligo"))
    LOG_PANELS[[group]] <<- p

  # Fix the panel (x-axis) to a constant physical width so the x-axes align across
  # ALL plots regardless of legend size; total figure width adapts to the legend.
  # Panel narrow + short so the individual plots are compact.
  g <- ggplot2::ggplotGrob(p)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  g$widths[pcol] <- grid::unit(1.6, "in")   # wide enough for the x-axis title (avoids left/right clip)
  ggsave(file.path(args$output_dir,
                   sprintf("plot_modulescore_dist_%s_%s%s.pdf", group, ct, suffix)),
         g, width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
         height = 2.1, units = "in", device = "pdf")

  # ---- source data: exact drawn rows (fitted curves + CI) ----
  fit_out <- fit_df %>%
    transmute(celltype = ct, group = group, distance_scale = scale_label,
              module = as.character(module), dist_to_phf1_um,
              fitted_score, ci_lo, ci_hi, significant)
  write.table(fit_out, file.path(args$output_dir,
              sprintf("source_data_modulescore_dist_%s_%s%s_fit.tsv", group, ct, suffix)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- supplementary source data: per-sample x ring means ----
  ring_long <- lapply(modules, function(m) {
    d <- df[, c("sample_id", "dist_to_phf1_um")]
    d$score    <- df[[m]]
    d$dist_ring <- cut(d$dist_to_phf1_um, breaks = RING_BREAKS, labels = RING_LAB,
                       include.lowest = TRUE, right = TRUE)
    d <- d[!is.na(d$dist_ring), ]
    d %>% group_by(sample_id, dist_ring) %>%
      summarise(mean_score = mean(score), sd_score = sd(score), n_cells = dplyr::n(),
                .groups = "drop") %>%
      mutate(celltype = ct, group = group, module = m,
             ring_mid_um = RING_MID[as.character(dist_ring)])
  }) %>% bind_rows() %>%
    dplyr::select(celltype, group, module, sample_id, dist_ring, ring_mid_um,
                  mean_score, sd_score, n_cells)
  write.table(ring_long, file.path(args$output_dir,
              sprintf("source_data_modulescore_dist_%s_%s%s_persample.tsv", group, ct, suffix)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- stats log (triple) ----
  sink(file.path(args$output_dir,
                 sprintf("stats_modulescore_dist_%s_%s%s.txt", group, ct, suffix)))
  cat("Module score vs distance to PHF1+ neuron (significance = model FDR)\n")
  cat("Group:", group, " | celltype:", ct, " | distance scale:", scale_label, "\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat(model_label, " (per module):\n", sep = "")
  cat("  ", formula_str, "\n", sep = "")
  cat("  Outcome = Seurat AddModuleScore (SCT data slot; ctrl=", CTRL, " nbin=", NBIN,
      " seed=", args$seed, ").\n", sep = "")
  if (is_gam) {
    cat("  GAM smooth on RAW distance. For GAM rows in the tables: t = s(dist) F-statistic,\n")
    cat("  df = s(dist) edf, p_LMM = smooth p-value; slope_* columns are NA (no single slope).\n")
  } else {
    cat("  Distance scaled to SD units (dist_sd =", round(dist_sd, 4), ").\n")
  }
  if (is.finite(MAX_DIST)) {
    cat("  Distance cap: cells with dist_to_phf1_um <=", MAX_DIST, "um.\n")
    cat("  Cells modelled (within cap):", nrow(df), " of ", nrow(df_full),
        " PHF1-negative | donors:", dplyr::n_distinct(df$sample_id), "\n")
  } else {
    cat("  Cells (PHF1-negative):", nrow(df), " | donors:", dplyr::n_distinct(df$sample_id), "\n")
  }
  cat("  Donor random intercept (1|sample_id).\n\n")
  cat("Significance = BH-adjusted model p (lmm_padj) within each group; padj <", SIG_ALPHA,
      " -> significant.\n\n")
  cat("Per-module results:\n")
  print(as.data.frame(mod_stats %>% dplyr::select(
    module, slope_per_sd, ci_lo_per_sd, ci_hi_per_sd, t, df, p_LMM, lmm_padj,
    significant, singular, verdict)))
  cat("\nsessionInfo():\n"); print(sessionInfo())
  sink()

  # accumulate structured coefficient table
  lmm_coef_rows[[paste(group, ct, scale_label)]] <<- mod_stats %>%
    dplyr::select(group, celltype, module, distance_scale, slope_per_sd, se,
                  ci_lo_per_sd, ci_hi_per_sd, slope_per_unit, t, df, p_LMM, lmm_padj,
                  significant, singular)

  cat(sprintf("  Done %s/%s%s: %d modules, %d significant (model FDR)\n",
              group, ct, suffix, nrow(mod_stats), sum(mod_stats$significant, na.rm = TRUE)))
}

##  ............................................................................
##  Per-donor rolling means (descriptive; one line per donor)                ####
# One faceted figure per group: a facet per module, a line per donor, coloured by
# Braak. Deliberately model-free -- see R/rollmean_by_donor.R.
run_group_bydonor <- function(gdef) {
  group <- gdef$group; ct <- gdef$ct; modules <- gdef$modules
  df_full <- DF_BY_NAME[[gdef$df_name]]
  df <- if (is.finite(MAX_DIST)) df_full[df_full$dist_to_phf1_um <= MAX_DIST, ] else df_full
  df <- df[is.finite(df$dist_to_phf1_um), ]
  x_cap <- if (is.finite(MAX_DIST)) MAX_DIST else as.numeric(quantile(df$dist_to_phf1_um, 0.99))

  long <- rollmean_by_donor(df, modules, x_cap = x_cap, half_window = HALF_WINDOW)
  if (is.null(long)) {
    cat("  [skip bydonor] ", group, "/", ct, ": no donor had enough cells\n", sep = "")
    return(invisible(NULL))
  }
  long$group <- group; long$celltype <- ct

  ncol_facets <- if (length(modules) > 4) 4 else length(modules)
  p <- plot_rollmean_by_donor(long, expression("Distance to PHF1+ neuron (" * mu * "m)"),
                              "Module score (rolling mean)", x_cap,
                              facet_col = "module", facet_levels = modules,
                              ncol = ncol_facets)
  if (!is.null(p)) {
    g <- ggplot2::ggplotGrob(p)
    pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
    g$widths[pcol] <- grid::unit(1.35, "in")
    nrow_facets <- ceiling(length(modules) / ncol_facets)
    ggsave(file.path(args$output_dir,
                     sprintf("plot_modulescore_bydonor_dist_%s_%s.pdf", group, ct)), g,
           width  = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
           height = 0.5 + 1.5 * nrow_facets, units = "in", device = "pdf")
  }

  write.table(long %>% dplyr::select(group, celltype, module, sample_id, Braak, window_um,
                                     dist_to_phf1_um, roll_mean, sem, n_window, n_cells_donor),
              file.path(args$output_dir,
                sprintf("source_data_modulescore_bydonor_dist_%s_%s_rollmean.tsv", group, ct)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  sink(file.path(args$output_dir,
                 sprintf("stats_modulescore_bydonor_dist_%s_%s.txt", group, ct)))
  cat("Per-donor rolling-mean module scores vs distance to PHF1+ neuron\n")
  cat("Group:", group, " | celltype:", ct, "\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("Outcome = Seurat AddModuleScore (SCT data slot; ctrl=", CTRL, " nbin=", NBIN,
      " seed=", args$seed, ").\n", sep = "")
  cat("Modules:", paste(modules, collapse = ", "), "\n\n")
  rollmean_by_donor_stats(long, x_cap, HALF_WINDOW,
                          n_donors_total = dplyr::n_distinct(df$sample_id))
  cat("\nColour = Braak stage, using a darkened line variant of braak_palette\n")
  cat("  (braak_line_palette in R/rollmean_by_donor.R): the canonical stage-1 fill colour is\n")
  cat("  too pale to read as a thin line on white. Stage ordering and hue family are unchanged.\n")
  cat("\nsessionInfo():\n"); print(sessionInfo())
  sink()

  cat(sprintf("  Wrote per-donor rolling-mean figure: %s / %s (%d modules, %d donors)\n",
              group, ct, length(modules), dplyr::n_distinct(long$sample_id)))
  invisible(long)
}

##  ............................................................................
##  Length constant (model type "decay")                                    ####
# Asymptotic decay of each module score away from the nearest PHF1+ neuron, fitted
# per donor and meta-analysed -> lambda (um) with a CI, and d95 = 2.996*lambda.
# See R/decay_length_utils.R for the estimator and its gates.
#
# GATED on the log model, per module: a length constant is only meaningful where the
# primary model found a gradient, so modules the log model calls n.s. are recorded in
# the headline table with a "not fitted" verdict and nothing is fitted for them.
#
# MICROGLIA AND ASTROCYTES ONLY. Oligodendrocytes are deliberately out of scope for
# the length-constant analysis (they are still fitted by the linear/log/gam types).
DECAY_GROUPS    <- c("mancuso", "cameron_pooled", "cameron_subclusters")
# Solid/dashed here encodes whether lambda passed the gates -- NOT
# significance (that is the log model's job), so the dashed level is labelled
# "flagged", not "ns".
DECAY_SIG_BASIS <- "lambda quotable"
DECAY_NS_BASIS  <- "flagged"

run_group_decay <- function(gdef) {
  group <- gdef$group; ct <- gdef$ct; modules <- gdef$modules; merged <- gdef$merged
  suffix <- "_decay"; scale_label <- "decay_um"
  cat(sprintf("\n=== %s / %s  (model: %s%s) ===\n", group, ct, scale_label,
              if (is.finite(MAX_DIST)) sprintf(", <=%gum", MAX_DIST) else ""))
  if (!group %in% DECAY_GROUPS) {
    cat("  [skip decay] out of scope for the length constant (microglia + astrocytes only)\n")
    return(invisible(NULL))
  }

  df_full <- DF_BY_NAME[[gdef$df_name]]
  df <- if (is.finite(MAX_DIST)) df_full[df_full$dist_to_phf1_um <= MAX_DIST, ] else df_full
  df <- df[is.finite(df$dist_to_phf1_um), ]
  # Canonical guard (as in R/nn3_over_phf1_distance.R): the log-linear competitor in
  # the model comparison needs dist > 0.
  if (any(df$dist_to_phf1_um <= 0))
    stop("Non-positive dist_to_phf1_um in the decay fit for ", group, "/", ct, ".")
  x_cap <- if (is.finite(MAX_DIST)) MAX_DIST else as.numeric(quantile(df$dist_to_phf1_um, 0.99))

  # ---- gate: per-module LOG-model significance ----
  log_coef <- lmm_coef_rows[[paste(group, ct, "log_um")]]
  sig_lut  <- setNames(rep(FALSE, length(modules)), modules)
  padj_lut <- setNames(rep(NA_real_, length(modules)), modules)
  if (!is.null(log_coef)) {
    sig_lut[as.character(log_coef$module)]  <- log_coef$significant %in% TRUE
    padj_lut[as.character(log_coef$module)] <- log_coef$lmm_padj
  }
  fit_mods <- modules[sig_lut[modules]]
  for (m in modules[!sig_lut[modules]]) {
    cat(sprintf("  [skip decay] %s: no log-model gradient (log padj = %s)\n",
                m, format.pval(padj_lut[[m]], digits = 3)))
    decay_rows[[paste(group, ct, m)]] <<- data.frame(
      group = group, celltype = ct, module = m, log_lmm_padj = padj_lut[[m]],
      n_donors_included = NA_integer_, n_donors_total = NA_integer_,
      lambda_um = NA_real_, lambda_lo = NA_real_, lambda_hi = NA_real_,
      quotable = FALSE, verdict = "not fitted: no log-model gradient",
      stringsAsFactors = FALSE)
  }
  if (length(fit_mods) == 0) {
    cat("  no module in this group has a log-model gradient; nothing to fit\n")
    return(invisible(NULL))
  }

  # ---- baseline: SEG null (see the stats log for the full rationale) ----
  # With the SEG null in play the fitted baseline is P_j + b*SEG(d): the SHAPE of the
  # stably-expressed-gene null, scaled by an estimated b, so only the EXCESS decay
  # above that null becomes A and lambda. b is estimated rather than fixed at 1
  # because the scores share a decline shape but sit on very different scales.
  seg_ok  <- USE_SEG && "seg_baseline" %in% names(df)
  covar_p <- c("nUMI_log", "percent_neg")
  covar   <- covar_p
  if (seg_ok) {
    # TWO covariates: seg_profile (per-donor distance-binned mean) carries the null's
    # SYSTEMATIC distance shape -- the baseline proper -- while the raw per-cell
    # seg_baseline absorbs cell-level technical variation. The raw column alone cannot
    # do the first job: regression dilution shrinks its coefficient towards 0 and it
    # removes almost none of the drift (see seg_distance_profile()).
    df$seg_profile <- seg_distance_profile(df$dist_to_phf1_um, df$seg_baseline, df$sample_id)
    covar <- c(covar_p, "seg_baseline", "seg_profile")
  }
  if (USE_SEG && !seg_ok)
    cat("  [warn] no SEG baseline available here; falling back to a free flat plateau\n")

  # ---- fit each surviving module ----
  res_list <- list(); res_nb_list <- list(); NULL_BENCH <- list()
  for (m in fit_mods) {
    if (seg_ok) NULL_BENCH[[m]] <- null_benchmark_gradient(df, m, "seg_baseline")
    res_list[[m]] <- decay_analyse(df, m, cap_um = x_cap, covar_cols = covar,
                                   lambda_min      = args$decay_lambda_min,
                                   lambda_max_mult = args$decay_lambda_max_mult,
                                   min_cells       = args$decay_min_cells_donor)
    # companion fit without the SEG baseline, so the log can show what it changed
    if (seg_ok)
      res_nb_list[[m]] <- decay_analyse(df, m, cap_um = x_cap, covar_cols = covar_p,
                                        lambda_min      = args$decay_lambda_min,
                                        lambda_max_mult = args$decay_lambda_max_mult,
                                        min_cells       = args$decay_min_cells_donor)
    r <- res_list[[m]]
    cat(sprintf("  %s: lambda = %s um [%s, %s] | %s\n", m,
                if (is.finite(r$meta$lambda))    sprintf("%.0f", r$meta$lambda)    else "NA",
                if (is.finite(r$meta$lambda_lo)) sprintf("%.0f", r$meta$lambda_lo) else "NA",
                if (is.finite(r$meta$lambda_hi)) sprintf("%.0f", r$meta$lambda_hi) else "NA",
                r$gates$verdict))
    rnb <- res_nb_list[[m]]
    decay_rows[[paste(group, ct, m)]] <<- cbind(
      decay_summary_row(r, group = group, celltype = ct, module = m,
                        log_lmm_padj = padj_lut[[m]],
                        seg_baseline = if (seg_ok) "SEG null" else "free plateau"),
      data.frame(
        lambda_no_seg    = if (is.null(rnb)) NA_real_ else rnb$meta$lambda,
        lambda_no_seg_lo = if (is.null(rnb)) NA_real_ else rnb$meta$lambda_lo,
        lambda_no_seg_hi = if (is.null(rnb)) NA_real_ else rnb$meta$lambda_hi,
        amplitude_no_seg = if (is.null(rnb)) NA_real_ else rnb$amplitude,
        quotable_no_seg  = if (is.null(rnb)) NA else rnb$gates$quotable,
        # magnitude benchmark: how many times steeper than the stably-expressed-gene null
        bench_sig_drop_sd = if (is.null(NULL_BENCH[[m]])) NA_real_ else NULL_BENCH[[m]]$sig[["mean"]],
        bench_seg_drop_sd = if (is.null(NULL_BENCH[[m]])) NA_real_ else NULL_BENCH[[m]]$seg[["mean"]],
        bench_excess_sd   = if (is.null(NULL_BENCH[[m]])) NA_real_ else NULL_BENCH[[m]]$excess[["mean"]],
        bench_excess_lo   = if (is.null(NULL_BENCH[[m]])) NA_real_ else NULL_BENCH[[m]]$excess[["lo"]],
        bench_excess_hi   = if (is.null(NULL_BENCH[[m]])) NA_real_ else NULL_BENCH[[m]]$excess[["hi"]],
        bench_ratio       = if (is.null(NULL_BENCH[[m]])) NA_real_ else NULL_BENCH[[m]]$ratio,
        bench_donors_sig_steeper = if (is.null(NULL_BENCH[[m]])) NA_integer_ else
                                     NULL_BENCH[[m]]$n_donors_same_dir))
    don <- cbind(group = group, celltype = ct, module = m, r$donors)
    decay_donors[[paste(group, ct, m)]] <<- don
  }

  cols <- if (merged) MERGED_COLOURS[modules]
          else setNames(MODULE_COLOURS[seq_along(modules)], modules)
  decay_colours[fit_mods] <<- cols[fit_mods]

  # ---- overlaid figure: one curve per module, lambda carried in the legend ----
  drawn <- names(Filter(function(r) !is.null(r$curve), res_list))
  if (length(drawn) > 0) {
    fit_df <- bind_rows(lapply(drawn, function(m) {
      r <- res_list[[m]]
      r$curve %>% mutate(module = m,
                         sig = factor(ifelse(r$gates$quotable, DECAY_SIG_BASIS, DECAY_NS_BASIS),
                                      levels = c(DECAY_SIG_BASIS, DECAY_NS_BASIS)))
    }))
    fit_df$module <- factor(fit_df$module, levels = drawn)
    # per-module lambda rug along the bottom: shows how much the length constants
    # differ between the states, which is the point of an overlaid panel. Only
    # gate-passing lambdas are marked -- a tick is a claim about a position on the
    # x-axis, and a lambda that failed the gates has no defensible position.
    quot <- vapply(drawn, function(m) isTRUE(res_list[[m]]$gates$quotable), logical(1))
    rug_df <- data.frame(module = factor(drawn, levels = drawn),
                         lambda_um = vapply(drawn, function(m) res_list[[m]]$meta$lambda, numeric(1)),
                         quotable = quot)
    rug_df <- rug_df[is.finite(rug_df$lambda_um) & rug_df$quotable, , drop = FALSE]
    # Legend labels carry lambda, but ONLY for modules that passed the gates. plotmath
    # (not a literal 'µ') so the glyph survives the default pdf device.
    lab_exprs <- lapply(drawn, function(m) {
      l <- res_list[[m]]$meta$lambda
      if (is.finite(l) && isTRUE(res_list[[m]]$gates$quotable))
        bquote(.(m) ~ "(" * lambda == .(round(l)) ~ mu * "m)") else bquote(.(m))
    })
    legend_title <- switch(group, mancuso = "Mancuso Micro", cameron_pooled = "Cameron Astro",
                           cameron_subclusters = "Cameron Astro", NULL)
    xlab <- expression("Distance to PHF1+ neuron (" * mu * "m)")

    p <- ggplot(fit_df, aes(dist_to_phf1_um, fitted_score, colour = module, group = module)) +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = module), alpha = 0.15, colour = NA) +
      geom_line(aes(linetype = sig, linewidth = sig)) +
      scale_colour_manual(values = cols[drawn], labels = lab_exprs, name = legend_title,
                          drop = FALSE) +
      scale_fill_manual(values = cols[drawn], guide = "none", drop = FALSE) +
      scale_linetype_manual(values = setNames(c("solid", "dashed"),
                                              c(DECAY_SIG_BASIS, DECAY_NS_BASIS)),
                            labels = list(expression(lambda ~ "quotable"), DECAY_NS_BASIS),
                            name = NULL, drop = FALSE,
                            limits = c(DECAY_SIG_BASIS, DECAY_NS_BASIS)) +
      scale_linewidth_manual(values = setNames(c(0.8, 0.4), c(DECAY_SIG_BASIS, DECAY_NS_BASIS)),
                             labels = list(expression(lambda ~ "quotable"), DECAY_NS_BASIS),
                             name = NULL, drop = FALSE,
                             limits = c(DECAY_SIG_BASIS, DECAY_NS_BASIS)) +
      guides(colour = guide_legend(order = 1),
             linetype = guide_legend(order = 2), linewidth = guide_legend(order = 2)) +
      labs(x = xlab, y = "Module score (decay fit)") +
      coord_cartesian(xlim = c(0, x_cap)) +
      theme_classic(base_size = 8) + fig_theme +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
            legend.text = element_text(size = 6), legend.key.width = grid::unit(20, "pt"),
            legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
            plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
    if (nrow(rug_df))
      p <- p + geom_rug(data = rug_df, aes(x = lambda_um, colour = module),
                        inherit.aes = FALSE, sides = "b", linewidth = 0.4,
                        length = grid::unit(0.06, "npc"), show.legend = FALSE)

    g <- ggplot2::ggplotGrob(p)
    pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
    g$widths[pcol] <- grid::unit(1.6, "in")
    ggsave(file.path(args$output_dir,
                     sprintf("plot_modulescore_dist_%s_%s%s.pdf", group, ct, suffix)),
           g, width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
           height = 2.1, units = "in", device = "pdf")

    fit_out <- fit_df %>%
      transmute(celltype = ct, group = group, distance_scale = scale_label,
                module = as.character(module), dist_to_phf1_um, fitted_score, ci_lo, ci_hi,
                plateau, amplitude, lambda_um, quotable = sig == DECAY_SIG_BASIS)
    write.table(fit_out, file.path(args$output_dir,
                sprintf("source_data_modulescore_dist_%s_%s%s_fit.tsv", group, ct, suffix)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    cat("  no fitted curve for any module (lambda not estimable); figure skipped\n")
  }

  # ---- per-donor source data (INCLUDING the excluded donors + reasons) ----
  write.table(bind_rows(decay_donors[paste(group, ct, fit_mods)]),
              file.path(args$output_dir,
              sprintf("source_data_modulescore_dist_%s_%s%s_donors.tsv", group, ct, suffix)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- stats log: one per group, one block per fitted module ----
  sink(file.path(args$output_dir,
                 sprintf("stats_modulescore_dist_%s_%s%s.txt", group, ct, suffix)))
  cat("Module score LENGTH CONSTANT: decay away from PHF1+ neurons\n")
  cat("Group:", group, " | celltype:", ct, " | model: asymptotic exponential decay\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("Outcome = Seurat AddModuleScore (SCT data slot; ctrl=", CTRL, " nbin=", NBIN,
      " seed=", args$seed, ").\n", sep = "")
  cat(sprintf("Cells modelled (within cap): %d of %d | donors: %d\n",
              nrow(df), nrow(df_full), dplyr::n_distinct(df$sample_id)))
  cat("\n=== BASELINE: stably expressed gene (SEG) null ===\n")
  if (seg_ok) {
    cat("The baseline is NOT a free flat plateau. The per-cell SEG module score (scSEGIndex,\n")
    cat("Lin et al. 2019 GigaScience 8(9):giz106; per-celltype sets from\n")
    cat("R/prepare_seg_geneset_reference.R, the same sets as the negative-control script\n")
    cat("R/modulescore_seg_control_vs_phf1_distance.R) enters the model as a covariate, so the\n")
    cat("fitted baseline is P_j + b*SEG(d): the SHAPE of the stably-expressed-gene null,\n")
    cat("scaled by an estimated b. Only the EXCESS decay above that null becomes A and lambda.\n")
    cat("b is ESTIMATED, not fixed at 1: the two scores can share a decline shape while\n")
    cat("differing several-fold in magnitude, so raw differencing would remove almost none of\n")
    cat("it.\n")
    cat("TWO SEG covariates, because the null does two jobs and one column cannot do both:\n")
    cat("  * seg_profile  = per-donor DISTANCE-BINNED mean SEG -- the baseline proper, i.e. the\n")
    cat("    null's systematic distance shape with per-cell noise averaged out.\n")
    cat("  * seg_baseline = RAW per-cell SEG -- absorbs cell-level technical variation, but is\n")
    cat("    far too noisy to carry the distance shape (regression dilution would shrink its\n")
    cat("    coefficient towards 0), so alone it would look like an adjustment while doing\n")
    cat("    nothing. See seg_distance_profile() in R/decay_length_utils.R.\n")
    cat(sprintf("SEG genes on panel for %s: %d\n", ct, length(SEG_PRESENT[[ct]])))
    sg_near <- mean(df$seg_baseline[df$dist_to_phf1_um <= 100], na.rm = TRUE)
    sg_far  <- mean(df$seg_baseline[df$dist_to_phf1_um >= 500], na.rm = TRUE)
    cat(sprintf("SEG null's own gradient here: mean %.5f at <=100 um vs %.5f at >=500 um (%+.1f%%); ",
                sg_near, sg_far, 100 * (sg_far - sg_near) / abs(sg_near)))
    cat(sprintf("Spearman rho with distance = %+.3f\n",
                suppressWarnings(stats::cor(df$seg_baseline, df$dist_to_phf1_um,
                                            method = "spearman", use = "complete.obs"))))
    cat("ASSUMPTION: a gradient in genes selected to be stably expressed is technical, not\n")
    cat("  biological.\n")
    cat("Where the null and a module share a decline SHAPE, exp(-d/lambda) and seg_profile are\n")
    cat("  near-collinear at the fitted lambda, so the SEG-adjusted fit is reported as a\n")
    cat("  sensitivity analysis, alongside the magnitude benchmark per module below.\n")
  } else {
    cat("NOT USED for this group -- free flat plateau P_j.")
    cat(if (USE_SEG) " (requested, but no SEG set was available.)\n" else " (--decay_seg_baseline no)\n")
  }

  cat("\nGATE: only modules whose LOG model showed a distance effect (lmm_padj < ", SIG_ALPHA,
      ") are\n  fitted here. Modules in this group and their log padj:\n", sep = "")
  for (m in modules)
    cat(sprintf("  %-20s log padj = %-10s -> %s\n", m, format.pval(padj_lut[[m]], digits = 3),
                if (sig_lut[[m]]) "FITTED" else "not fitted"))
  for (m in fit_mods) {
    cat("\n", strrep("=", 78), "\n", sep = "")
    cat("MODULE: ", m, "\n", sep = "")
    cat(strrep("=", 78), "\n\n", sep = "")
    decay_write_stats_block(res_list[[m]], m)
    rnb <- res_nb_list[[m]]
    if (!is.null(rnb)) {
      cat("\n=== What the SEG baseline changed (same cells, same estimator) ===\n")
      fmt <- function(r) if (is.finite(r$meta$lambda))
        sprintf("%.1f um [%.1f, %.1f]", r$meta$lambda, r$meta$lambda_lo, r$meta$lambda_hi)
        else "not estimable"
      cat(sprintf("  lambda WITH SEG baseline    : %-28s (donors %d, quotable %s)\n",
                  fmt(res_list[[m]]), res_list[[m]]$gates$n_donors_included,
                  res_list[[m]]$gates$quotable))
      cat(sprintf("  lambda WITHOUT SEG baseline : %-28s (donors %d, quotable %s)\n",
                  fmt(rnb), rnb$gates$n_donors_included, rnb$gates$quotable))
      cat(sprintf("  amplitude  with SEG %+.4f  |  without SEG %+.4f\n",
                  res_list[[m]]$amplitude, rnb$amplitude))
      cat(sprintf("  shared-lambda anchor  with SEG %.0f um  |  without SEG %.0f um\n",
                  res_list[[m]]$shared$lambda, rnb$shared$lambda))
    }
    if (seg_ok) {
      cat("\n=== MAGNITUDE BENCHMARK against the SEG null ===\n")
      decay_write_null_benchmark(NULL_BENCH[[m]], m)
      cat("\nPer-donor near-to-far drops (SD units):\n")
      print(NULL_BENCH[[m]]$per_donor, row.names = FALSE, digits = 4)
    }
  }
  cat("\nsessionInfo():\n"); print(sessionInfo())
  sink()

  cat(sprintf("  Done %s/%s%s: %d module(s) fitted, %d quotable\n", group, ct, suffix,
              length(fit_mods),
              sum(vapply(res_list, function(r) isTRUE(r$gates$quotable), logical(1)))))
  invisible(res_list)
}

##  ............................................................................
##  Raw (model-free) rolling-mean overlay, per group                        ####
# Same format as the PHF1 module-score-over-distance raw plot (model-free rolling
# mean + 95% CI), but with EVERY module of the group overlaid on one axis (colour =
# module) and NO PHF1+ reference line (glia). Each module's rolling-mean line is drawn
# SOLID+thick / DASHED+thin using that module's LOG-model significance (from the log
# analysis, lmm_coef_rows key '<group> <ct> log_um'), matching the log fitted plot.
# One plot per group: Mancuso Micro, Cameron Astro (subclusters), Pandey Oligo.
make_glia_raw_plot <- function(gdef) {
  group <- gdef$group; ct <- gdef$ct; modules <- gdef$modules; merged <- gdef$merged
  df_full <- DF_BY_NAME[[gdef$df_name]]
  df <- if (is.finite(MAX_DIST)) df_full[df_full$dist_to_phf1_um <= MAX_DIST, ] else df_full
  x_cap <- if (is.finite(MAX_DIST)) MAX_DIST else as.numeric(quantile(df$dist_to_phf1_um, 0.99))
  grid <- seq(0, x_cap, length.out = 200)
  cols <- if (merged) MERGED_COLOURS[modules]
          else setNames(MODULE_COLOURS[seq_along(modules)], modules)
  sig_basis <- "padj<0.05"

  # per-module LOG-model significance (thick/solid vs dashed), from the log analysis
  log_coef <- lmm_coef_rows[[paste(group, ct, "log_um")]]
  sig_lut  <- setNames(rep(FALSE, length(modules)), modules)
  if (!is.null(log_coef)) sig_lut[as.character(log_coef$module)] <- log_coef$significant %in% TRUE
  sig_fac <- function(m) factor(ifelse(sig_lut[as.character(m)], sig_basis, "ns"),
                                levels = c(sig_basis, "ns"))

  # RAW rolling mean (NOT z-scored): z-scoring flattens the several modules onto each
  # other, so plot the native AddModuleScore scale on a shared axis (the range copes).
  roll_df <- bind_rows(lapply(modules, function(m) {
    y <- df[[m]]; ok <- !is.na(y) & !is.na(df$dist_to_phf1_um)
    data.frame(module = m, roll_mean(df$dist_to_phf1_um[ok], y[ok], grid, HALF_WINDOW),
               dist_to_phf1_um = grid)
  }))
  roll_df$module <- factor(roll_df$module, levels = modules)
  roll_df$sig    <- sig_fac(roll_df$module)

  legend_title <- switch(group, mancuso = "Mancuso Micro", cameron_pooled = "Cameron Astro",
                         cameron_subclusters = "Cameron Astro", pandey_oligo = "Pandey Oligo", NULL)
  xlab <- expression("Distance to PHF1+ neuron (" * mu * "m)")   # mu via plotmath (survives default pdf device)

  p <- ggplot(roll_df, aes(dist_to_phf1_um, roll_mean, colour = module, fill = module, group = module)) +
    geom_ribbon(aes(ymin = roll_mean - 1.96 * sem, ymax = roll_mean + 1.96 * sem),
                alpha = 0.15, colour = NA) +
    geom_line(aes(linetype = sig, linewidth = sig)) +
    scale_colour_manual(values = cols, name = legend_title, drop = FALSE) +
    scale_fill_manual(values = cols, guide = "none", drop = FALSE) +
    scale_linetype_manual(values = setNames(c("solid", "dashed"), c(sig_basis, "ns")),
                          name = "Log-model FDR", drop = FALSE, limits = c(sig_basis, "ns")) +
    scale_linewidth_manual(values = setNames(c(0.8, 0.4), c(sig_basis, "ns")),
                           name = "Log-model FDR", drop = FALSE, limits = c(sig_basis, "ns")) +
    guides(colour = guide_legend(order = 1),
           linetype = guide_legend(order = 2), linewidth = guide_legend(order = 2)) +
    labs(x = xlab, y = "Module score") +
    coord_cartesian(xlim = c(0, x_cap)) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.text = element_text(size = 6),
          legend.key.width = grid::unit(18, "pt"), legend.key.height = grid::unit(8, "pt"),
          legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))

  g <- ggplot2::ggplotGrob(p)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  g$widths[pcol] <- grid::unit(1.6, "in")   # wide enough for the x-axis title (avoids left/right clip)
  ggsave(file.path(args$output_dir, sprintf("plot_modulescore_raw_dist_%s_%s.pdf", group, ct)),
         g, width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
         height = 2.1, units = "in", device = "pdf")

  out <- roll_df %>% mutate(group = group, celltype = ct, window_um = args$window_um,
                            log_significant = sig == sig_basis) %>%
    dplyr::select(group, celltype, module, window_um, log_significant,
                  dist_to_phf1_um, roll_mean, sem, n_window)
  write.table(out, file.path(args$output_dir,
              sprintf("source_data_modulescore_raw_dist_%s_%s_rollmean.tsv", group, ct)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  sink(file.path(args$output_dir, sprintf("stats_modulescore_raw_dist_%s_%s.txt", group, ct)))
  cat("Raw (model-free) rolling-mean module scores vs distance to PHF1+ neuron\n")
  cat("Group:", group, " | celltype:", ct, "\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("Rolling mean +/- 95% CI (mean +/- 1.96*SEM) of per-cell RAW AddModuleScore, per module\n")
  cat("(NOT z-scored; native scale, shared axis). Y axis = 'Module score (rolling mean, 95% CI)'.\n")
  cat("Window:", args$window_um, "um; distance cap:",
      if (is.finite(MAX_DIST)) MAX_DIST else "none", "um; grid: 200 points. NO PHF1+ reference (glia).\n")
  cat("Line style = LOG-model significance (solid+thick = model FDR<0.05, dashed+thin = n.s.),\n")
  cat("  taken from the log-distance LMM analysis (stats_modulescore_dist_", group, "_", ct,
      "_log.txt).\n", sep = "")
  cat("Cells:", nrow(df), " | donors:", dplyr::n_distinct(df$sample_id), "\n")
  cat("Modules (solid = log-significant):",
      paste(modules[sig_lut[modules]], collapse = ", "), "\n")
  cat("Modules:", paste(modules, collapse = ", "), "\n\n")
  cat("sessionInfo():\n"); print(sessionInfo())
  sink()

  cat(sprintf("  Wrote glia raw rolling-mean plot: %s / %s (%d modules)\n", group, ct, length(modules)))
  p   # return the ggplot for the combined figure
}

##  ............................................................................
##  Run everything                                                          ####
# Raw overlay plots: one per glial type -- Mancuso Micro, Cameron Astro (SUBCLUSTERS
# only, not the pooled Neurotoxic/Neuroprotective), Pandey Oligo.
RAW_GROUPS <- c("mancuso", "cameron_subclusters", "pandey_oligo")
RAW_PANELS <- list()
for (gdef in GROUPS) {
  # MODEL TYPES RUN. Only the log-distance fit is reported; "linear" and "gam" can be
  # added to this vector -- their code paths in run_group_model() are intact.
  # "decay" must come after "log": it is gated on the log model's per-module verdict.
  for (mt in c("log", "decay")) {
    run_group_model(gdef, mt)
  }
  if (gdef$group %in% RAW_GROUPS) RAW_PANELS[[gdef$group]] <- make_glia_raw_plot(gdef)
  run_group_bydonor(gdef)   # descriptive per-donor view, every group
}

# Structured coefficient table
write.table(bind_rows(lmm_coef_rows),
            file.path(args$output_dir, "stats_modulescore_lmm_coeffs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

##  ............................................................................
##  Length-constant headline table + forest                                 ####
if (length(decay_rows) > 0) {
  dtab <- bind_rows(decay_rows)
  write.table(dtab, file.path(args$output_dir, "stats_modulescore_decay_lambda.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n#### Length constants (lambda, um) ####\n")
  print(as.data.frame(dtab %>% dplyr::select(dplyr::any_of(
    c("group", "celltype", "module", "lambda_um", "lambda_lo", "lambda_hi", "d95_um",
      "n_donors_included", "quotable", "verdict")))), row.names = FALSE)

  ftab <- dtab[is.finite(dtab$lambda_um), , drop = FALSE]
  if (nrow(ftab) > 0) {
    ftab$row_label <- sprintf("%s | %s", ftab$celltype, ftab$module)
    ftab$series    <- ftab$module
    dpts <- bind_rows(decay_donors)
    dpts <- dpts[dpts$included, , drop = FALSE]
    if (nrow(dpts)) dpts$row_label <- sprintf("%s | %s", dpts$celltype, dpts$module)
    fp <- decay_lambda_forest(ftab, dpts, cap_um = if (is.finite(MAX_DIST)) MAX_DIST else NULL,
                              colours = decay_colours)
    if (!is.null(fp)) {
      fp <- fp + forest_theme + theme(legend.position = "bottom")
      ggsave(file.path(args$output_dir, "plot_decay_lambda_forest.pdf"), fp,
             width = 9.0, height = 1.2 + 0.42 * nrow(ftab), units = "cm", device = "pdf")
      write.table(ftab %>% dplyr::select(dplyr::any_of(
                    c("row_label", "group", "celltype", "module", "lambda_um",
                      "lambda_lo", "lambda_hi", "d95_um", "d95_lo", "d95_hi",
                      "n_donors_included", "n_donors_total", "quotable", "verdict"))),
                  file.path(args$output_dir, "source_data_decay_lambda_forest.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("Wrote lambda forest:", nrow(ftab), "row(s)\n")
    }
  }
}

# Combined 3-panel LOG figure (Micro Mancuso | Astro Cameron subclusters | Oligo Pandey),
# sized to fit within an A4-page width; legends on the right (compact Astro labels).
if (all(c("mancuso", "cameron_subclusters", "pandey_oligo") %in% names(LOG_PANELS))) {
  combo <- LOG_PANELS[["mancuso"]] + LOG_PANELS[["cameron_subclusters"]] +
           LOG_PANELS[["pandey_oligo"]] + patchwork::plot_layout(nrow = 1)
  ggsave(file.path(args$output_dir, "plot_modulescore_dist_COMBINED_log_Micro_Astro_Oligo.pdf"),
         combo, width = 19.05, height = 5, units = "cm", device = "pdf")
  cat("Wrote combined 3-panel log figure (A4 width): Micro | Astro subclusters | Oligo\n")
}

# Combined RAW rolling-mean (z-scored) figure: Micro (Mancuso) | Astro (Cameron subclusters).
# Oligo (Pandey) is deliberately excluded here.
if (all(c("mancuso", "cameron_subclusters") %in% names(RAW_PANELS))) {
  combo_raw <- RAW_PANELS[["mancuso"]] + RAW_PANELS[["cameron_subclusters"]] +
               patchwork::plot_layout(nrow = 1)
  ggsave(file.path(args$output_dir, "plot_modulescore_raw_COMBINED_Micro_Astro.pdf"),
         combo_raw, width = 12.7, height = 5, units = "cm", device = "pdf")
  cat("Wrote combined RAW rolling-mean figure (Micro | Astro), 12.7 cm\n")
}

cat("\nAll outputs written to:", args$output_dir, "\n")
cat("Done.\n")
