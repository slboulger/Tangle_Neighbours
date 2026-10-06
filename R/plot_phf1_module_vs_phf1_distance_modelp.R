#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_phf1_module_vs_phf1_distance_modelp.R
#
# Shared utility - sourced by the scripts above, no panel of its own
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# plot_phf1_module_vs_phf1_distance_modelp.R
#
# ptau-signature module scores over distance to nearest PHF1+ neuron, within the
# neuron subtypes eligible for PHF1 analysis. Direct sibling of
# plot_modulescore_vs_phf1_distance_modelp.R (SAME scoring, SAME LMM/GAM models,
# SAME model-FDR significance and figure style). Scores THREE signatures per
# eligible celltype, each analysed independently:
#
#   * phf1       -- the CELLTYPE-SPECIFIC PHF1 marker set from
#                   find_phf1_markers_by_celltype.R (phf1_markers/<ct>/<ct>_phf1_geneset.txt;
#                   upregulated in PHF1+ vs PHF1- within that celltype).
#   * otero_up   -- Otero-Garcia AT8 ptau signature, UP genes  (default set otero_L23_up).
#   * otero_down -- Otero-Garcia AT8 ptau signature, DOWN genes (default set otero_L23_down).
#     (Otero sets are GLOBAL: the same genes are scored in every celltype. Read from
#      otero_signatures/otero_at8_signatures.rds; --otero_up_set / --otero_down_set select which.)
#
# For EVERY signature x eligible celltype the FULL individual-plot set is produced
# (fitted linear/log/gam figures + an actual-data raw figure), exactly as for PHF1.
#
# INTEGRATED plot (one per eligible celltype): the integrated signatures' rolling-mean
# trends + 95% CI (NO raw points) on one axis, each Z-SCORED per signature (mean 0 / SD 1
# on that celltype's PHF1-negative cells within the cap) so the distance-trends are
# comparable in shape. Each signature also carries its PHF1+ neurons' reference
# (mean, 95% CI) transformed on the SAME z-scale. INTEGRATED_SIGS selects which
# signatures appear here (default phf1 + otero_up; Otero DOWN is EXCLUDED from the
# overlay but still gets its own full individual-plot set).
#
# The PHF1+ neurons define distance 0 and are EXCLUDED from every regression /
# rolling mean (fit/trend on PHF1-negative cells only); they appear only as the
# reference band / line.
#
# NOTES (also written into each stats log):
#   * The PHF1 geneset is SELECTED from the PHF1+/PHF1- contrast within the celltype.
#     The distance regression uses PHF1-NEGATIVE cells only (selection contrast !=
#     distance); the PHF1+ reference band is elevated by construction and is drawn as a
#     reference. (The Otero sets are external, so not subject to this selection.)
#   * Significance is the PARAMETRIC model FDR (BH within a celltype x signature).
#     With a single module per celltype x signature, model FDR == the raw model p.
#   * AddModuleScore values from different gene sets are not on the same absolute
#     scale; the integrated plot z-scores each signature to make the SHAPES
#     comparable (magnitudes are lost by design).
#
# LENGTH CONSTANT (model type "decay"). For every signature x celltype the LOG model
# calls significant, an asymptotic decay
#     score = P_j + A_j * exp(-dist/lambda_j) + b1*nUMI_log + b2*percent_neg
# is fitted PER DONOR by profile least squares and lambda meta-analysed across donors
# (DerSimonian-Laird + Knapp-Hartung), giving the length constant lambda in MICRONS with
# a CI, plus d95 = 2.996*lambda -- the distance at which 95% of the excess has resolved.
# See R/decay_length_utils.R for the estimator, its gates (the decay must beat the
# log-linear model by a margin, lambda must be resolved inside the window, its CI must
# be sufficiently narrow, and it must not track the distance cap) and the deliberate
# use of RAW MICRONS rather than the canonical log/SD distance transform.
# The log model remains primary for "is there a gradient?"; decay answers "how far?".
#
# MODEL TYPES RUN: log-distance LMM (the reported fit) and decay. The "linear" and
# "gam" types are not run (to cut run time) -- only the log fit is used for reporting.
# Their code paths in run_ct_model() are intact; enable them by adding them to the
# `for (mt in ...)` vector in the driver at the bottom.
#
# PER-DONOR VIEW: one faceted figure per celltype (a facet per signature, a line per
# donor, coloured by Braak), model-free; it shows each donor's trend alongside the
# pooled rolling mean. See R/rollmean_by_donor.R.
#
# OUTPUT (per eligible celltype <ct>, per signature <sig> in {phf1, oteroUP, oteroDOWN}),
# under --output_dir:
#   Fitted (per model type; suffix in {"_log", "_decay"}):
#     plot_<sig>modulescore_dist_<ct><suffix>.pdf            -- fitted curve + CI + PHF1+ band
#     source_data_<sig>modulescore_dist_<ct><suffix>_fit.tsv -- exact drawn rows (curve + CI)
#     source_data_<sig>modulescore_dist_<ct><suffix>_persample.tsv -- per-donor x ring means
#     stats_<sig>modulescore_dist_<ct><suffix>.txt           -- formula, model FDR, PHF1+ ref, sessionInfo
#   Actual data (model-free):
#     plot_<sig>modulescore_dist_<ct>_raw.pdf                -- points + rolling mean + 95% CI + PHF1+ band
#     source_data_<sig>modulescore_dist_<ct>_raw.tsv         -- one row per drawn cell
#     source_data_<sig>modulescore_dist_<ct>_raw_rollmean.tsv-- rolling-mean grid
#   Integrated (per celltype; TWO variants -- <tag> in {zscore, raw}; INTEGRATED_SIGS):
#     plot_integrated_<tag>_dist_<ct>.pdf                    -- rolling means + 95% CI + PHF1+ refs
#       (zscore = each signature z-scored so shapes compare; raw = native AddModuleScore scale)
#     source_data_integrated_<tag>_dist_<ct>_rollmean.tsv, ..._phf1ref.tsv
#     stats_integrated_<tag>_dist_<ct>.txt
#     plot_integrated_<tag>_COMBINED.pdf                     -- per-celltype integrated panels, side by side
#   Per-donor (descriptive, one figure per celltype, faceted by signature):
#     plot_modulescore_bydonor_dist_<ct>.pdf                   -- one line per donor, Braak-coloured
#     source_data_modulescore_bydonor_dist_<ct>_rollmean.tsv   -- exact drawn rows
#     stats_modulescore_bydonor_dist_<ct>.txt                  -- window, donor table, caveats
#   Length constant (only where the log model is significant):
#     plot_<sig>modulescore_dist_<ct>_decay.pdf                -- decay curve + lambda CI + plateau
#     source_data_<sig>modulescore_dist_<ct>_decay_fit.tsv     -- exact drawn curve rows
#     source_data_<sig>modulescore_dist_<ct>_decay_donors.tsv  -- per-donor lambda (incl. exclusions)
#     stats_<sig>modulescore_dist_<ct>_decay.txt               -- meta-analysis, effect sizes, gates
#   Written once:
#     stats_modulescore_geneset_overlap.tsv, stats_modulescore_lmm_coeffs.tsv,
#     source_data_modulescore_phf1ref.tsv,
#     stats_modulescore_decay_lambda.tsv                       -- the headline lambda table
#     plot_decay_lambda_forest.pdf, source_data_decay_lambda_forest.tsv
#
# Use run_phf1_module_distance_modelp.sh

##  ............................................................................
##  Packages + setup                                                        ####
suppressPackageStartupMessages({
  library(Seurat)
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
source("R/palettes.R")                       # fig_theme, forest_theme, neuron_order, ...
source("R/eligible_neuron_subtypes.R")       # eligible_neuron_subtypes()
source("R/decay_length_utils.R")             # decay_analyse() and friends (model type "decay")
source("R/rollmean_by_donor.R")              # per-donor rolling means, coloured by Braak

set.seed(42)

##  ............................................................................
##  Arguments                                                               ####
parser <- ArgumentParser()
parser$add_argument("--seu", default = "seu_PHF1.rds",
  help = "Path to seu_PHF1.rds [default: seu_PHF1.rds]")
parser$add_argument("--markers_dir", default = "phf1_markers",
  help = "Directory holding <ct>/<ct>_phf1_geneset.txt [default: phf1_markers]")
parser$add_argument("--otero_rds", default = "otero_signatures/otero_at8_signatures.rds",
  help = "RDS: named list of Otero-Garcia AT8 signature gene vectors")
parser$add_argument("--otero_up_set", default = "otero_L23_up",
  help = "Name of the Otero UP set to score [default: otero_L23_up]")
parser$add_argument("--otero_down_set", default = "otero_L23_down",
  help = "Name of the Otero DOWN set to score [default: otero_L23_down]")
parser$add_argument("--output_dir",
  default = "plots/phf1_module_vs_phf1_distance_modelp_1000um",
  help = "Output directory")
parser$add_argument("--seed", type = "integer", default = 42,
  help = "Seed (AddModuleScore + provenance) [default: 42]")
parser$add_argument("--max_dist_um", type = "double", default = 1000,
  help = "Restrict modelled cells to dist_to_phf1_um <= this (um); 0 = no cap [default: 1000]")
parser$add_argument("--window_um", type = "double", default = 75,
  help = "Rolling-mean window WIDTH (um) for the actual-data / integrated plots [default: 75]")
parser$add_argument("--min_phf1_cells", type = "integer", default = 6,
  help = "Eligibility: min PHF1+ cells per sample [default: 6]")
parser$add_argument("--min_samples", type = "integer", default = 5,
  help = "Eligibility: min samples meeting min_phf1_cells [default: 5]")
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

# Constants (mirror plot_modulescore_vs_phf1_distance_modelp.R) ---------------
MIN_GENES_PRESENT <- 10                          # drop genesets with fewer panel genes
CTRL              <- 100
NBIN              <- 24
RING_BREAKS       <- c(0, 50, 100, 200, 300, 500, 700, 1000)
SIG_ALPHA         <- 0.05                         # significance threshold (model FDR)
GAM_K             <- 5                            # basis dim for the GAM distance smooth
HALF_WINDOW       <- args$window_um / 2           # rolling-mean half-window (um)

PHF1_REF_COL <- "grey30"                          # PHF1+ neurons reference band/line
# internal factor levels for the (legend-less) individual plots
SERIES_FIT <- "signature (non-PHF1+ neurons)"
SERIES_RAW <- "non-PHF1+ neurons"
SERIES_REF <- "PHF1+ neurons (mean, 95% CI)"

# Ring labels / midpoints
.ring_lower <- head(RING_BREAKS, -1)
.ring_upper <- tail(RING_BREAKS, -1)
RING_LAB    <- sprintf("%g-%g", .ring_lower, .ring_upper)
RING_MID    <- setNames((.ring_lower + .ring_upper) / 2, RING_LAB)

XLAB <- expression("Distance to PHF1+ neuron (" * mu * "m)")   # mu via plotmath (survives default pdf device)

##  ............................................................................
##  Helpers (as in plot_modulescore_vs_phf1_distance_modelp.R unless noted) ####

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

# Subset Seurat to a celltype, score the supplied signature list (one AddModuleScore
# call, one column per signature key), return one data.frame with cell_id + PHF1
# status + covariates + module-score columns.
extract_celltype_data <- function(seu, ct, present_sets, seed,
                                  sex_col, age_col, pmi_col, seg_set = NULL) {
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
    # Braak is carried only to COLOUR the per-donor rolling-mean lines; it is not a
    # covariate in any model here. Tolerated as missing so a plotting nicety can
    # never abort the analysis -- the figure falls back to colouring by donor.
    Braak           = if ("Braak" %in% colnames(md)) as.character(md$Braak) else NA_character_,
    Sex             = as.character(md[[sex_col]]),
    Age             = as.numeric(md[[age_col]]),
    PMI             = as.numeric(md[[pmi_col]]),
    stringsAsFactors = FALSE
  )
  scores <- score_group(seu_ct, present_sets, "SIG_", seed)   # cols = signature keys
  # Stably-expressed-gene (SEG) null baseline for the decay model. Scored in its OWN
  # AddModuleScore call so the signature scores above are bit-identical to a run
  # without it: score_group() re-seeds every call and each call builds its own
  # expression-matched control bins, so adding this cannot perturb linear/log/gam.
  if (!is.null(seg_set) && length(seg_set) > 0)
    scores <- cbind(scores, score_group(seu_ct, list(seg_baseline = seg_set), "SEGBL_", seed))
  cbind(base, scores)
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

# Fitted trend across a distance grid, Wald CI from fixed-effect vcov.
predict_grid <- function(model, df, dist_sd, log_distance, x_cap) {
  lo      <- if (log_distance) min(df$dist_to_phf1_um[df$dist_to_phf1_um > 0], na.rm = TRUE) else 0
  grid_um <- seq(lo, x_cap, length.out = 200)
  dt      <- if (log_distance) log(grid_um) else grid_um
  nd <- data.frame(
    dist_scaled = dt / dist_sd, nUMI_log = mean(df$nUMI_log),
    percent_neg = mean(df$percent_neg),
    Sex = factor(levels(df$Sex)[1], levels = levels(df$Sex)),
    Age_s = mean(df$Age_s), PMI_s = mean(df$PMI_s))
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

predict_grid_gam <- function(model, df, x_cap) {
  grid_um <- seq(0, x_cap, length.out = 200)
  nd <- data.frame(
    dist = grid_um, nUMI_log = mean(df$nUMI_log), percent_neg = mean(df$percent_neg),
    Sex = factor(levels(df$Sex)[1], levels = levels(df$Sex)),
    Age_s = mean(df$Age_s), PMI_s = mean(df$PMI_s),
    sample_id = factor(levels(df$sample_id)[1], levels = levels(df$sample_id)))
  pr  <- predict(model, newdata = nd, se.fit = TRUE, exclude = "s(sample_id)")
  fit <- as.numeric(pr$fit); se <- as.numeric(pr$se.fit)
  data.frame(dist_to_phf1_um = grid_um, fitted_score = fit,
             ci_lo = fit - 1.96 * se, ci_hi = fit + 1.96 * se)
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

# Fix panel width so x-axes align; save at a compact height.
save_fixed_panel <- function(p, path) {
  g <- ggplot2::ggplotGrob(p)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  g$widths[pcol] <- grid::unit(1.6, "in")   # wide enough for the x-axis title (avoids left/right clip)
  ggsave(path, g, width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
         height = 2.1, units = "in", device = "pdf")
}

# PHF1+ reference band (mean + 95% CI of the mean) for a given score column.
phf1_reference <- function(df, score_col) {
  s <- df[[score_col]][df$phf1_pos & !is.na(df[[score_col]])]
  n <- length(s)
  if (n < 1) return(NULL)
  m <- mean(s); sdv <- stats::sd(s); sem <- if (n > 1) sdv / sqrt(n) else NA_real_
  list(mean = m, sd = sdv, sem = sem, ci_lo = m - 1.96 * sem, ci_hi = m + 1.96 * sem, n = n)
}

##  ............................................................................
##  Signatures + eligible celltypes                                         ####
otero <- readRDS(args$otero_rds)
stopifnot("Otero UP set not found"   = args$otero_up_set   %in% names(otero),
          "Otero DOWN set not found" = args$otero_down_set %in% names(otero))
cat("Otero sets: UP =", args$otero_up_set, "(", length(otero[[args$otero_up_set]]),
    "genes) | DOWN =", args$otero_down_set, "(", length(otero[[args$otero_down_set]]), "genes)\n")

read_geneset <- function(markers_dir, ct) {
  ct_safe <- gsub("[^A-Za-z0-9_-]", "_", ct)
  f <- file.path(markers_dir, ct_safe, paste0(ct_safe, "_phf1_geneset.txt"))
  if (!file.exists(f)) return(NULL)
  g <- trimws(readLines(f)); unique(g[nzchar(g)])
}

# Signature registry. genes(ct) returns the raw gene vector for that celltype.
SIGNATURES <- list(
  phf1 = list(key = "phf1", label = "PHF1 signature", file_tag = "phf1modulescore",
              colour = "#D55E00", collection = "phf1_markers",             # Okabe-Ito vermillion (CB-safe)
              genes = function(ct) read_geneset(args$markers_dir, ct)),
  otero_up = list(key = "otero_up", label = "Otero-Garcia Exc1/2 UP", file_tag = "oteroUPmodulescore",
              colour = "#0072B2", collection = "otero_at8",                # Okabe-Ito blue
              genes = function(ct) otero[[args$otero_up_set]]),
  otero_down = list(key = "otero_down", label = "Otero-Garcia DOWN", file_tag = "oteroDOWNmodulescore",
              colour = "#CC79A7", collection = "otero_at8",                # Okabe-Ito reddish-purple
              genes = function(ct) otero[[args$otero_down_set]])
)
SIG_ORDER        <- c("phf1", "otero_up", "otero_down")   # individual-plot order (all signatures)
INTEGRATED_SIGS  <- c("phf1", "otero_up")                 # integrated overlay: Otero DOWN excluded

project_root <- getwd()
elig <- eligible_neuron_subtypes(project_root,
                                 min_phf1_cells = args$min_phf1_cells,
                                 min_samples    = args$min_samples,
                                 neuron_levels  = neuron_order)
if (length(elig) == 0) stop("No eligible neuron subtypes found.")
cat("Eligible neuron subtypes:\n  ", paste(elig, collapse = ", "), "\n", sep = "")

##  ............................................................................
##  Load Seurat object, determine panel                                     ####
cat("Loading Seurat object (this is large)...\n")
seu <- readRDS(args$seu)
DefaultAssay(seu) <- "SCT"

stopifnot("celltype column missing"        = "celltype"        %in% colnames(seu@meta.data),
          "PHF1 column missing"            = "PHF1"            %in% colnames(seu@meta.data),
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

##  ............................................................................
##  SEG null baseline for the decay model                                   ####
# Stably expressed genes (scSEGIndex, Lin et al. 2019), per celltype, from
# R/prepare_seg_geneset_reference.R -- the same sets used by the negative-control
# script R/modulescore_seg_control_vs_phf1_distance.R, so the two agree by
# construction. Used ONLY by the decay model type (see run_ct_decay).
USE_SEG      <- identical(args$decay_seg_baseline, "yes")
SEG_PRESENT  <- list()
seg_overlap  <- list()
if (USE_SEG) {
  if (!requireNamespace("qs", quietly = TRUE))
    stop("Package 'qs' is required to read the SEG geneset (--decay_seg_baseline no to skip).")
  if (!file.exists(args$seg_qs))
    stop("SEG geneset not found: ", args$seg_qs,
         "\n  Run R/prepare_seg_geneset_reference.R, or pass --decay_seg_baseline no.")
  seg <- qs::qread(args$seg_qs)
  stopifnot("SEG geneset must be a named list of character vectors" =
              is.list(seg) && !is.null(names(seg)) && all(vapply(seg, is.character, logical(1))))
  cat("SEG baseline sets:", paste(names(seg), collapse = ", "), "\n")
  for (ct in elig) {
    if (!ct %in% names(seg)) {
      cat("  [warn] no SEG set for ", ct, "; its decay fit falls back to a free plateau\n", sep = "")
      next
    }
    pr <- prep_sets(setNames(list(seg[[ct]]), "seg_baseline"), panel_genes,
                    MIN_GENES_PRESENT, "Reference_scSEGIndex", ct)
    seg_overlap[[ct]] <- pr$overlap %>% mutate(celltype = ct)
    if (length(pr$present) > 0) {
      SEG_PRESENT[[ct]] <- pr$present[["seg_baseline"]]
      cat(sprintf("  %s: %d SEG genes on panel\n", ct, length(SEG_PRESENT[[ct]])))
    } else {
      cat("  [warn] ", ct, ": too few SEG genes on the panel; free-plateau fallback\n", sep = "")
    }
  }
  if (length(seg_overlap) > 0)
    write.table(bind_rows(seg_overlap) %>%
                  dplyr::select(celltype, collection, set, n_total, n_present, frac_present, scored),
                file.path(args$output_dir, "stats_modulescore_decay_seg_overlap.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
}

##  ............................................................................
##  Score each eligible celltype (all signatures in one AddModuleScore pass) ####
CT_DATA      <- list()   # per-ct scored data frame (full population)
CT_SIGS      <- list()   # per-ct: signature keys actually scored
overlap_rows <- list()

for (ct in elig) {
  raw_sets <- lapply(SIGNATURES, function(s) { g <- s$genes(ct); if (is.null(g)) character(0) else unique(g) })
  # panel-filter each signature; record overlap; keep those with >= MIN_GENES_PRESENT.
  present <- list()
  for (k in names(SIGNATURES)) {
    prep <- prep_sets(setNames(list(raw_sets[[k]]), k), panel_genes, MIN_GENES_PRESENT,
                      SIGNATURES[[k]]$collection, ct)
    overlap_rows[[paste(ct, k)]] <- prep$overlap %>% mutate(celltype = ct, signature = k)
    if (length(prep$present) > 0) present[[k]] <- prep$present[[k]]
  }
  if (length(present) == 0) { cat("  [skip] ", ct, ": no usable signature\n", sep = ""); next }
  cat("Scoring ", ct, " (", paste(sprintf("%s=%d", names(present), lengths(present)),
                                   collapse = ", "), ")...\n", sep = "")
  CT_DATA[[ct]] <- extract_celltype_data(seu, ct, present, seed = args$seed,
                                         sex_col = sex_col, age_col = age_col, pmi_col = pmi_col,
                                         seg_set = SEG_PRESENT[[ct]])
  CT_SIGS[[ct]] <- names(present)
}
rm(seu); invisible(gc())

if (length(overlap_rows) > 0)
  write.table(bind_rows(overlap_rows) %>%
                dplyr::select(celltype, signature, collection, set, n_total, n_present, frac_present, scored),
              file.path(args$output_dir, "stats_modulescore_geneset_overlap.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
if (length(CT_DATA) == 0) stop("No eligible celltype had a usable signature.")

##  ............................................................................
##  Plot builders                                                           ####

# Fitted-line figure with the PHF1+ reference band (legend-less; encoding:
# <colour> solid/dashed = signature fitted trend in non-PHF1+ neurons (solid = model
# FDR<0.05, dashed = n.s.) with 95% CI; grey dotted line + band = PHF1+ neurons mean/95% CI).
make_fitted_plot <- function(fit_df, ref, x_cap, sig_basis, module_colour, ylab) {
  ref_df <- data.frame(dist_to_phf1_um = c(0, x_cap), fitted_score = ref$mean,
                       ci_lo = ref$ci_lo, ci_hi = ref$ci_hi,
                       series = factor(SERIES_REF, levels = c(SERIES_FIT, SERIES_REF)))
  fit_df$series <- factor(SERIES_FIT, levels = c(SERIES_FIT, SERIES_REF))
  series_cols <- setNames(c(module_colour, PHF1_REF_COL), c(SERIES_FIT, SERIES_REF))
  ggplot() +
    geom_ribbon(data = ref_df, aes(dist_to_phf1_um, ymin = ci_lo, ymax = ci_hi, fill = series),
                alpha = 0.13, colour = NA) +
    geom_line(data = ref_df, aes(dist_to_phf1_um, fitted_score, colour = series),
              linetype = "dotted", linewidth = 0.4) +
    geom_ribbon(data = fit_df, aes(dist_to_phf1_um, ymin = ci_lo, ymax = ci_hi, fill = series),
                alpha = 0.15, colour = NA) +
    geom_line(data = fit_df, aes(dist_to_phf1_um, fitted_score, colour = series,
                                 linetype = sig, linewidth = sig)) +
    scale_colour_manual(values = series_cols, guide = "none", drop = FALSE) +
    scale_fill_manual(values = series_cols, guide = "none", drop = FALSE) +
    scale_linetype_manual(values = setNames(c("solid", "dashed"), c(sig_basis, "ns")),
                          guide = "none", drop = FALSE, limits = c(sig_basis, "ns")) +
    scale_linewidth_manual(values = setNames(c(0.8, 0.4), c(sig_basis, "ns")),
                           guide = "none", drop = FALSE, limits = c(sig_basis, "ns")) +
    labs(x = XLAB, y = ylab) +
    coord_cartesian(xlim = c(0, x_cap)) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none",
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
}

# Actual-data figure: raw points + model-free rolling mean + PHF1+ reference band (legend-less).
# The rolling-mean line is SOLID+thick when the LOG-model distance effect is significant
# (sig_flag = model FDR<0.05), DASHED+thin otherwise -- same annotation as the log fitted plot.
make_raw_plot <- function(raw_df, roll_df, ref, x_cap, module_colour, ylab, sig_flag) {
  ref_df <- data.frame(dist_to_phf1_um = c(0, x_cap), ci_lo = ref$ci_lo, ci_hi = ref$ci_hi,
                       ymean = ref$mean, series = factor(SERIES_REF, levels = c(SERIES_RAW, SERIES_REF)))
  raw_df$series  <- factor(SERIES_RAW, levels = c(SERIES_RAW, SERIES_REF))
  roll_df$series <- factor(SERIES_RAW, levels = c(SERIES_RAW, SERIES_REF))
  series_cols <- setNames(c(module_colour, PHF1_REF_COL), c(SERIES_RAW, SERIES_REF))
  roll_lty <- if (isTRUE(sig_flag)) "solid" else "dashed"   # log-model FDR<0.05 -> solid
  roll_lwd <- if (isTRUE(sig_flag)) 0.8 else 0.4
  ggplot() +
    geom_ribbon(data = ref_df, aes(dist_to_phf1_um, ymin = ci_lo, ymax = ci_hi, fill = series),
                alpha = 0.13, colour = NA) +
    geom_line(data = ref_df, aes(dist_to_phf1_um, ymean, colour = series),
              linetype = "dotted", linewidth = 0.4) +
    geom_point(data = raw_df, aes(dist_to_phf1_um, score, colour = series),
               alpha = 0.15, size = 0.3, stroke = 0) +
    geom_ribbon(data = roll_df, aes(dist_to_phf1_um, ymin = roll_mean - 1.96 * sem,
                                    ymax = roll_mean + 1.96 * sem, fill = series),
                alpha = 0.25, colour = NA) +
    geom_line(data = roll_df, aes(dist_to_phf1_um, roll_mean, colour = series),
              linetype = roll_lty, linewidth = roll_lwd) +
    scale_colour_manual(values = series_cols, guide = "none", drop = FALSE) +
    scale_fill_manual(values = series_cols, guide = "none", drop = FALSE) +
    labs(x = XLAB, y = ylab) +
    coord_cartesian(xlim = c(0, x_cap)) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none",
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
}

##  ............................................................................
##  Analysis functions                                                      ####
lmm_coef_rows <- list()
phf1ref_rows  <- list()
decay_rows    <- list()   # headline lambda table (one row per signature x celltype)
decay_donors  <- list()   # per-donor lambdas, for the forest's background points

# One fitted model type, one celltype, one signature.
run_ct_model <- function(ct, df_full_all, ref, model_type, sg) {
  # The decay model type shares nothing with the lmer/gam paths (per-donor profile
  # least squares + meta-analysis, raw-micron distance), so it dispatches out here.
  if (model_type == "decay") return(run_ct_decay(ct, df_full_all, ref, sg))
  m <- sg$key; is_gam <- model_type == "gam"; log_distance <- model_type == "log"
  suffix      <- switch(model_type, linear = "",   log = "_log",   gam = "_gam")
  scale_label <- switch(model_type, linear = "um", log = "log_um", gam = "gam")
  model_label <- switch(model_type, linear = "LMM (linear distance)", log = "LMM (log distance)",
                        gam = sprintf("GAM (penalized smooth, mgcv::bam, k=%d)", GAM_K))
  formula_str <- if (is_gam)
    "score ~ s(dist_to_phf1_um, k=K) + nUMI_log + percent_neg + Sex + Age + PMI + s(sample_id, bs='re')"
  else
    "score ~ dist_to_phf1_um_scaled + nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)"
  cat(sprintf("\n=== %s / %s  (model: %s%s) ===\n", sg$key, ct, scale_label,
              if (is.finite(MAX_DIST)) sprintf(", <=%gum", MAX_DIST) else ""))

  df_full <- df_full_all
  df_full$Age_s <- as.numeric(scale(df_full$Age))
  df_full$PMI_s <- as.numeric(scale(df_full$PMI))
  df <- if (is.finite(MAX_DIST)) df_full[df_full$dist_to_phf1_um <= MAX_DIST, ] else df_full
  if (log_distance) df <- df[df$dist_to_phf1_um > 0, ]

  if (is_gam) { dist_sd <- NA_real_ } else {
    dt <- if (log_distance) log(df$dist_to_phf1_um) else df$dist_to_phf1_um
    dist_sd <- stats::sd(dt); df$dist_scaled <- dt / dist_sd
  }
  x_cap <- if (is.finite(MAX_DIST)) MAX_DIST else as.numeric(quantile(df$dist_to_phf1_um, 0.99))

  if (is_gam) {
    sub <- df[, c("cell_id", "nUMI_log", "percent_neg", "Sex", "Age_s", "PMI_s", "sample_id")]
    sub$dist <- df$dist_to_phf1_um; sub$score <- df[[m]]
    obs <- fit_obs_gam(sub)
    if (is.null(obs)) { cat("  GAM failed, skipped\n"); return(invisible(NULL)) }
    grid <- predict_grid_gam(obs$model, df, x_cap)
    nrm  <- list(slope_per_sd = NA_real_, se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                 slope_per_unit = NA_real_, stat = obs$stat, dfv = obs$edf, p = obs$p, singular = NA)
  } else {
    sub <- df[, c("cell_id", "dist_scaled", "nUMI_log", "percent_neg", "Sex", "Age_s", "PMI_s", "sample_id")]
    sub$score <- df[[m]]
    obs <- fit_obs_lmer(sub)
    if (is.null(obs)) { cat("  LMM failed, skipped\n"); return(invisible(NULL)) }
    grid <- predict_grid(obs$model, df, dist_sd, log_distance, x_cap)
    nrm  <- list(slope_per_sd = obs$estimate, se = obs$se,
                 ci_lo = obs$estimate - 1.96 * obs$se, ci_hi = obs$estimate + 1.96 * obs$se,
                 slope_per_unit = obs$estimate / dist_sd, stat = obs$t, dfv = obs$df,
                 p = obs$p, singular = obs$singular)
  }

  mod_stats <- tibble(
    signature = sg$key, celltype = ct, module = m, distance_scale = scale_label,
    slope_per_sd = nrm$slope_per_sd, se = nrm$se, ci_lo_per_sd = nrm$ci_lo, ci_hi_per_sd = nrm$ci_hi,
    slope_per_unit = nrm$slope_per_unit, t = nrm$stat, df = nrm$dfv, p_LMM = nrm$p,
    singular = nrm$singular, dist_sd = dist_sd)
  mod_stats$lmm_padj    <- p.adjust(mod_stats$p_LMM, method = "BH")   # single test -> == p
  mod_stats$significant <- mod_stats$lmm_padj < SIG_ALPHA
  mod_stats$verdict     <- ifelse(is.na(mod_stats$significant), "not testable",
                           ifelse(mod_stats$significant,
                                  "significant distance effect (model FDR<0.05)", "n.s. (model FDR)"))
  sig_basis <- "padj<0.05"

  fit_df <- grid %>%
    mutate(sig = factor(ifelse(mod_stats$significant, sig_basis, "ns"), levels = c(sig_basis, "ns")))

  p <- make_fitted_plot(fit_df, ref, x_cap, sig_basis, sg$colour,
                        paste0(sg$label, " score (fitted)"))
  save_fixed_panel(p, file.path(args$output_dir,
                   sprintf("plot_%s_dist_%s%s.pdf", sg$file_tag, ct, suffix)))

  fit_out <- fit_df %>%
    transmute(signature = sg$key, celltype = ct, module = m, distance_scale = scale_label,
              dist_to_phf1_um, fitted_score, ci_lo, ci_hi, significant = mod_stats$significant,
              phf1ref_mean = ref$mean, phf1ref_ci_lo = ref$ci_lo, phf1ref_ci_hi = ref$ci_hi)
  write.table(fit_out, file.path(args$output_dir,
              sprintf("source_data_%s_dist_%s%s_fit.tsv", sg$file_tag, ct, suffix)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  d <- df[, c("sample_id", "dist_to_phf1_um")]; d$score <- df[[m]]
  d$dist_ring <- cut(d$dist_to_phf1_um, breaks = RING_BREAKS, labels = RING_LAB,
                     include.lowest = TRUE, right = TRUE)
  d <- d[!is.na(d$dist_ring), ]
  ring_long <- d %>% group_by(sample_id, dist_ring) %>%
    summarise(mean_score = mean(score), sd_score = stats::sd(score), n_cells = dplyr::n(), .groups = "drop") %>%
    mutate(signature = sg$key, celltype = ct, module = m, ring_mid_um = RING_MID[as.character(dist_ring)]) %>%
    dplyr::select(signature, celltype, module, sample_id, dist_ring, ring_mid_um, mean_score, sd_score, n_cells)
  write.table(ring_long, file.path(args$output_dir,
              sprintf("source_data_%s_dist_%s%s_persample.tsv", sg$file_tag, ct, suffix)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  sink(file.path(args$output_dir, sprintf("stats_%s_dist_%s%s.txt", sg$file_tag, ct, suffix)))
  cat(sg$label, "score vs distance to PHF1+ neuron (MODEL-FDR version, no permutation)\n")
  cat("Signature:", sg$key, " | celltype:", ct, " | distance scale:", scale_label, "\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat(model_label, ":\n  ", formula_str, "\n", sep = "")
  cat("  Outcome = Seurat AddModuleScore of the ", sg$label, " geneset (SCT data slot; ctrl=", CTRL,
      " nbin=", NBIN, " seed=", args$seed, ").\n", sep = "")
  if (is_gam) {
    cat("  GAM smooth on RAW distance. t = s(dist) F-statistic, df = s(dist) edf,\n")
    cat("  p_LMM = smooth p-value; slope_* columns are NA.\n")
  } else cat("  Distance scaled to SD units (dist_sd =", round(dist_sd, 4), ").\n")
  if (is.finite(MAX_DIST)) {
    cat("  Distance cap: cells with dist_to_phf1_um <=", MAX_DIST, "um.\n")
    cat("  Cells modelled (within cap):", nrow(df), " of ", nrow(df_full),
        " PHF1-negative | donors:", dplyr::n_distinct(df$sample_id), "\n")
  } else cat("  Cells (PHF1-negative):", nrow(df), " | donors:", dplyr::n_distinct(df$sample_id), "\n")
  cat("  Inferential unit = sample (n donors); cell-level n is NOT the biological n.\n\n")
  cat("PHF1+ reference (drawn as the horizontal band; NOT in the regression):\n")
  cat("  n PHF1+ cells =", ref$n, " | mean =", round(ref$mean, 4),
      " | 95% CI = [", round(ref$ci_lo, 4), ",", round(ref$ci_hi, 4), "]\n\n")
  cat("Significance = BH-adjusted model p (lmm_padj); padj <", SIG_ALPHA, "-> significant.\n")
  if (sg$key == "phf1") {
    cat("NOTE (selection): the geneset is SELECTED from the PHF1+/PHF1- contrast within this\n")
    cat("  celltype. The regression is on PHF1-NEGATIVE cells only; the PHF1+ reference band\n")
    cat("  is elevated by construction and is drawn as a reference.\n")
  }
  cat("Single test per celltype x signature -> padj == p.\n\n")
  cat("Per-module results:\n")
  print(as.data.frame(mod_stats %>% dplyr::select(
    signature, module, slope_per_sd, ci_lo_per_sd, ci_hi_per_sd, t, df, p_LMM, lmm_padj,
    significant, singular, verdict)))
  cat("\nsessionInfo():\n"); print(sessionInfo())
  sink()

  lmm_coef_rows[[paste(sg$key, ct, scale_label)]] <<- mod_stats %>%
    dplyr::select(signature, celltype, module, distance_scale, slope_per_sd, se,
                  ci_lo_per_sd, ci_hi_per_sd, slope_per_unit, t, df, p_LMM, lmm_padj, significant, singular)
  cat(sprintf("  Done %s/%s%s: padj=%.3g, %s\n", sg$key, ct, suffix,
              mod_stats$lmm_padj[1], mod_stats$verdict[1]))
  invisible(list(significant = isTRUE(mod_stats$significant[1]),
                 lmm_padj = mod_stats$lmm_padj[1], p_LMM = mod_stats$p_LMM[1]))
}

##  ............................................................................
##  Per-donor rolling means (descriptive; one line per donor)                ####
# One faceted figure per celltype: a facet per signature, a line per donor,
# coloured by Braak. Model-free -- see R/rollmean_by_donor.R.
run_ct_bydonor <- function(ct, df_full_all, keys) {
  if (!length(keys)) return(invisible(NULL))
  df <- if (is.finite(MAX_DIST)) df_full_all[df_full_all$dist_to_phf1_um <= MAX_DIST, ]
        else df_full_all
  df <- df[is.finite(df$dist_to_phf1_um), ]
  x_cap <- if (is.finite(MAX_DIST)) MAX_DIST else as.numeric(quantile(df$dist_to_phf1_um, 0.99))

  long <- rollmean_by_donor(df, keys, x_cap = x_cap, half_window = HALF_WINDOW)
  if (is.null(long)) { cat("  [skip bydonor] ", ct, ": no donor had enough cells\n", sep = ""); return(invisible(NULL)) }

  labs_map <- vapply(SIGNATURES[keys], function(s) s$label, character(1))
  long$signature_label <- unname(labs_map[long$module])
  long$celltype <- ct

  p <- plot_rollmean_by_donor(long, XLAB, "Module score (rolling mean)", x_cap,
                              facet_col = "signature_label",
                              facet_levels = unname(labs_map[keys]), ncol = length(keys))
  if (!is.null(p)) {
    g <- ggplot2::ggplotGrob(p)
    pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
    g$widths[pcol] <- grid::unit(1.35, "in")   # per-facet panel width; x-axes align across figures
    ggsave(file.path(args$output_dir, sprintf("plot_modulescore_bydonor_dist_%s.pdf", ct)), g,
           width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
           height = 2.2, units = "in", device = "pdf")
  }

  write.table(long %>% dplyr::select(celltype, module, signature_label, sample_id, Braak,
                                     window_um, dist_to_phf1_um, roll_mean, sem,
                                     n_window, n_cells_donor),
              file.path(args$output_dir,
                        sprintf("source_data_modulescore_bydonor_dist_%s_rollmean.tsv", ct)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  sink(file.path(args$output_dir, sprintf("stats_modulescore_bydonor_dist_%s.txt", ct)))
  cat("Per-donor rolling-mean module scores vs distance to PHF1+ neuron\n")
  cat("Celltype:", ct, " | signatures:", paste(keys, collapse = ", "), "\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("Outcome = Seurat AddModuleScore (SCT data slot; ctrl=", CTRL, " nbin=", NBIN,
      " seed=", args$seed, "). PHF1+ cells are EXCLUDED (they define distance 0).\n\n", sep = "")
  rollmean_by_donor_stats(long, x_cap, HALF_WINDOW,
                          n_donors_total = dplyr::n_distinct(df$sample_id))
  cat("\nColour = Braak stage, using a darkened line variant of braak_palette\n")
  cat("  (braak_line_palette in R/rollmean_by_donor.R): the canonical stage-1 fill colour is\n")
  cat("  too pale to read as a thin line on white. Stage ordering and hue family are unchanged.\n")
  cat("\nsessionInfo():\n"); print(sessionInfo())
  sink()

  cat(sprintf("  Wrote per-donor rolling-mean figure: %s (%d signatures, %d donors)\n",
              ct, length(keys), dplyr::n_distinct(long$sample_id)))
  invisible(long)
}

##  ............................................................................
##  Length constant (model type "decay")                                    ####
# Asymptotic decay of the module score away from the nearest PHF1+ neuron, fitted
# per donor and meta-analysed -> lambda (um) with a CI, and d95 = 2.996*lambda.
# GATED on the log model: a length constant is only meaningful where the primary
# model found a gradient, so signatures the log model calls n.s. are recorded in
# the headline table with a "not fitted" verdict and nothing is fitted.

DECAY_SIG_BASIS <- "lambda quotable"   # reuses the fitted plots' solid/dashed encoding

# Overlay the numbers the analysis exists to produce on top of make_fitted_plot():
# plateau line, lambda and d95 markers, a rug of the per-donor lambdas, and the
# lambda / d95 annotation.
add_decay_annotations <- function(p, res, x_cap, colour, ref = NULL) {
  m <- res$meta
  if (!is.finite(m$lambda)) return(p)
  # Park the label top-right, but DUCK UNDER the PHF1+ reference band when that band
  # is the topmost feature -- otherwise the two overprint. The right-hand half of
  # the panel is empty once the curve has decayed, so there is room there.
  y_lab <- Inf
  if (!is.null(ref) && is.finite(ref$ci_lo) &&
      ref$ci_lo > max(res$curve$ci_hi, na.rm = TRUE)) y_lab <- ref$ci_lo

  # A lambda that fails the gates is not printed on the panel. These panels are
  # legend-less, so a dashed line alone is too subtle a cue; instead the panel states
  # that lambda is not identified and shows no number, no lambda or d95 marker, and no
  # plateau line.
  if (!isTRUE(res$gates$quotable))
    return(p + annotate("text", x = x_cap, y = y_lab, hjust = 1, vjust = 1.15,
                        size = 2.0, colour = "grey25",
                        label = "lambda not identified\n(see stats log)"))

  inc <- res$donors[res$donors$included & is.finite(res$donors$lambda_um), , drop = FALSE]
  lab <- sprintf("atop(lambda == %.0f ~ mu * 'm' ~ '[' * %.0f * '-' * %.0f * ']', d[95] == %.0f ~ mu * 'm')",
                 m$lambda, m$lambda_lo, m$lambda_hi, res$d95_um)
  if (is.finite(res$plateau))
    p <- p + geom_hline(yintercept = res$plateau, linetype = "dotted",
                        linewidth = 0.3, colour = colour)
  p <- p +
    geom_vline(xintercept = m$lambda, linetype = "dashed", linewidth = 0.3, colour = colour) +
    geom_vline(xintercept = min(res$d95_um, x_cap), linetype = "dotdash",
               linewidth = 0.3, colour = "grey45")
  if (nrow(inc))
    p <- p + geom_rug(data = inc, aes(x = lambda_um), inherit.aes = FALSE,
                      sides = "b", linewidth = 0.3, colour = "grey45",
                      length = grid::unit(0.05, "npc"))
  p + annotate("text", x = x_cap, y = y_lab, label = lab, parse = TRUE,
               hjust = 1, vjust = 1.15, size = 2.0, colour = "black")
}

run_ct_decay <- function(ct, df_full_all, ref, sg) {
  m <- sg$key; suffix <- "_decay"; scale_label <- "decay_um"
  key <- paste(sg$key, ct)
  cat(sprintf("\n=== %s / %s  (model: %s%s) ===\n", sg$key, ct, scale_label,
              if (is.finite(MAX_DIST)) sprintf(", <=%gum", MAX_DIST) else ""))

  # ---- gate on the LOG model (the project's primary model) ----
  log_row <- lmm_coef_rows[[paste(sg$key, ct, "log_um")]]
  log_padj <- if (is.null(log_row)) NA_real_ else log_row$lmm_padj[1]
  if (is.null(log_row) || !isTRUE(log_row$significant[1])) {
    cat(sprintf("  [skip decay] no log-model gradient (log padj = %s); lambda not fitted\n",
                format.pval(log_padj, digits = 3)))
    decay_rows[[key]] <<- data.frame(
      signature = sg$key, celltype = ct, module = m, log_lmm_padj = log_padj,
      n_donors_included = NA_integer_, n_donors_total = NA_integer_,
      lambda_um = NA_real_, lambda_lo = NA_real_, lambda_hi = NA_real_,
      quotable = FALSE, verdict = "not fitted: no log-model gradient",
      stringsAsFactors = FALSE)
    return(invisible(NULL))
  }

  # ---- data: raw microns, within the cap, PHF1-negative cells only ----
  df <- if (is.finite(MAX_DIST)) df_full_all[df_full_all$dist_to_phf1_um <= MAX_DIST, ]
        else df_full_all
  df <- df[is.finite(df$dist_to_phf1_um) & !is.na(df[[m]]), ]
  # Canonical guard (as in R/nn3_over_phf1_distance.R): every modelled PHF1-negative
  # cell must sit at a strictly positive distance -- the log-linear competitor in the
  # model comparison needs it, and a zero would mean a PHF1+ cell leaked in.
  if (any(df$dist_to_phf1_um <= 0))
    stop("Non-positive dist_to_phf1_um in the decay fit for ", sg$key, "/", ct,
         "; a PHF1-negative cell cannot be at distance 0.")
  x_cap <- if (is.finite(MAX_DIST)) MAX_DIST else as.numeric(quantile(df$dist_to_phf1_um, 0.99))

  # BASELINE. With the SEG null in play the fitted baseline is P_j + b*SEG(d): the
  # SHAPE of the stably-expressed-gene null, scaled by an estimated b, so only the
  # EXCESS decay above that null is attributed to the signature. b is estimated
  # rather than fixed at 1 because the two module scores sit on very different
  # absolute scales (SEG is ~15x smaller here) while sharing a decline shape;
  # differencing them raw would be wrong, scaling them is not.
  seg_ok  <- USE_SEG && "seg_baseline" %in% names(df)
  covar_p <- c("nUMI_log", "percent_neg")
  covar   <- covar_p
  if (seg_ok) {
    # TWO covariates, because the null does two different jobs and one column cannot
    # do both. seg_profile (the per-donor distance-binned mean) carries the null's
    # SYSTEMATIC distance shape -- the baseline proper; the raw per-cell seg_baseline
    # is far too noisy for that job (regression dilution would shrink its coefficient
    # to ~0 and remove none of the drift, see seg_distance_profile()) but it does
    # usefully absorb per-CELL technical variation.
    df$seg_profile <- seg_distance_profile(df$dist_to_phf1_um, df$seg_baseline, df$sample_id)
    covar <- c(covar_p, "seg_baseline", "seg_profile")
  }
  if (USE_SEG && !seg_ok)
    cat("  [warn] no SEG baseline available here; falling back to a free flat plateau\n")

  res <- decay_analyse(df, m, cap_um = x_cap, covar_cols = covar,
                       lambda_min      = args$decay_lambda_min,
                       lambda_max_mult = args$decay_lambda_max_mult,
                       min_cells       = args$decay_min_cells_donor)
  # Companion fit WITHOUT the SEG baseline, so the log can show what the baseline did.
  res_nb <- if (seg_ok) decay_analyse(df, m, cap_um = x_cap, covar_cols = covar_p,
                                      lambda_min      = args$decay_lambda_min,
                                      lambda_max_mult = args$decay_lambda_max_mult,
                                      min_cells       = args$decay_min_cells_donor) else NULL
  mm <- res$meta

  # ---- figure ----
  if (!is.null(res$curve)) {
    fit_df <- res$curve %>%
      mutate(sig = factor(ifelse(res$gates$quotable, DECAY_SIG_BASIS, "ns"),
                          levels = c(DECAY_SIG_BASIS, "ns")))
    # Grey out the curve when lambda failed the gates, so the panel cannot be read as
    # a length-constant result at a glance (the signature colour means "result").
    curve_col <- if (isTRUE(res$gates$quotable)) sg$colour else "grey60"
    p <- make_fitted_plot(fit_df, ref, x_cap, DECAY_SIG_BASIS, curve_col,
                          paste0(sg$label, " score (decay fit)"))
    p <- add_decay_annotations(p, res, x_cap, curve_col, ref)
    save_fixed_panel(p, file.path(args$output_dir,
                     sprintf("plot_%s_dist_%s%s.pdf", sg$file_tag, ct, suffix)))

    fit_out <- res$curve %>%
      transmute(signature = sg$key, celltype = ct, module = m, distance_scale = scale_label,
                dist_to_phf1_um, fitted_score, ci_lo, ci_hi, plateau, amplitude, lambda_um,
                quotable = res$gates$quotable,
                phf1ref_mean = ref$mean, phf1ref_ci_lo = ref$ci_lo, phf1ref_ci_hi = ref$ci_hi)
    write.table(fit_out, file.path(args$output_dir,
                sprintf("source_data_%s_dist_%s%s_fit.tsv", sg$file_tag, ct, suffix)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    cat("  no fitted curve (lambda not estimable); figure skipped\n")
  }

  # ---- per-donor source data (INCLUDING the excluded donors + reasons) ----
  don <- cbind(signature = sg$key, celltype = ct, module = m, res$donors)
  write.table(don, file.path(args$output_dir,
              sprintf("source_data_%s_dist_%s%s_donors.tsv", sg$file_tag, ct, suffix)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  decay_donors[[key]] <<- don

  # ---- stats log ----
  # Fraction of the PHF1+ / background gap spanned by the decaying microenvironment:
  # available only here, because only this script has a PHF1+ reference.
  gap      <- ref$mean - res$plateau
  gap_frac <- if (is.finite(gap) && abs(gap) > 1e-12) res$amplitude / gap else NA_real_

  nb <- NULL   # magnitude benchmark vs the SEG null; filled inside the log below
  sink(file.path(args$output_dir, sprintf("stats_%s_dist_%s%s.txt", sg$file_tag, ct, suffix)))
  cat(sg$label, "score: LENGTH CONSTANT of the decay away from PHF1+ neurons\n")
  cat("Signature:", sg$key, " | celltype:", ct, " | model: asymptotic exponential decay\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("Outcome = Seurat AddModuleScore of the ", sg$label, " geneset (SCT data slot; ctrl=",
      CTRL, " nbin=", NBIN, " seed=", args$seed, ").\n", sep = "")
  cat(sprintf("Fitted because the LOG model found a gradient here (log lmm_padj = %s).\n\n",
              format.pval(log_padj, digits = 3)))
  cat(sprintf("Cells modelled (PHF1-negative, within cap): %d | donors: %d\n\n",
              nrow(df), dplyr::n_distinct(df$sample_id)))

  cat("=== BASELINE: stably expressed gene (SEG) null ===\n")
  if (seg_ok) {
    cat("The baseline is NOT a free flat plateau. The per-cell SEG module score (scSEGIndex,\n")
    cat("Lin et al. 2019 GigaScience 8(9):giz106; per-celltype sets from\n")
    cat("R/prepare_seg_geneset_reference.R, the same sets as the negative-control script\n")
    cat("R/modulescore_seg_control_vs_phf1_distance.R) enters the model as a covariate, so the\n")
    cat("fitted baseline is P_j + b*SEG(d): the SHAPE of the stably-expressed-gene null,\n")
    cat("scaled by an estimated b. Only the EXCESS decay above that null becomes A and lambda.\n")
    cat("TWO SEG covariates are used, because the null does two jobs and one column cannot do\n")
    cat("  both:\n")
    cat("   * seg_profile  = the per-donor DISTANCE-BINNED mean SEG score. This is the baseline\n")
    cat("     proper -- the null's systematic distance shape, with per-cell noise averaged out.\n")
    cat("   * seg_baseline = the RAW per-cell SEG score, which absorbs cell-level technical\n")
    cat("     variation but is far too noisy to carry the distance shape: the systematic part\n")
    cat("     is a small fraction of the per-cell scatter, so regression dilution would shrink\n")
    cat("     its coefficient towards 0 and remove almost none of the drift. Using the raw\n")
    cat("     column alone would look like an adjustment while doing nothing.\n")
    cat("b is ESTIMATED, not fixed at 1, because the SEG and signature module scores sit on\n")
    cat("  different absolute scales; raw differencing of the two scores would not align them.\n")
    cat(sprintf("SEG genes on panel: %d\n", length(SEG_PRESENT[[ct]])))
    sg_near <- mean(df$seg_baseline[df$dist_to_phf1_um <= 100], na.rm = TRUE)
    sg_far  <- mean(df$seg_baseline[df$dist_to_phf1_um >= 500], na.rm = TRUE)
    cat(sprintf("SEG null's own gradient here: mean %.5f at <=100 um vs %.5f at >=500 um (%+.1f%%);",
                sg_near, sg_far, 100 * (sg_far - sg_near) / abs(sg_near)))
    cat(sprintf(" Spearman rho with distance = %+.3f\n",
                suppressWarnings(stats::cor(df$seg_baseline, df$dist_to_phf1_um,
                                            method = "spearman", use = "complete.obs"))))
    cat("ASSUMPTION: a gradient in genes selected to be stably expressed is technical, not\n")
    cat("  biological. If cells near tangles were genuinely more transcriptionally active as a\n")
    cat("  whole, this adjustment would remove real signal (over-adjustment).\n")
    cat("NOTE on the adjusted lambda: where the null and the signature share a decline SHAPE,\n")
    cat("  exp(-d/lambda) and seg_profile are near-COLLINEAR at the fitted lambda, so the split\n")
    cat("  of a shared gradient between them is poorly determined. The adjusted fit is a\n")
    cat("  sensitivity analysis; the magnitude benchmark below is the direct comparison\n")
    cat("  against the null.\n")
  } else {
    cat("NOT USED for this fit -- free flat plateau P_j.")
    cat(if (USE_SEG) " (requested, but no SEG set was available here.)\n" else
        " (--decay_seg_baseline no)\n")
  }
  cat("\n")
  decay_write_stats_block(res, paste(sg$label, "score"))

  if (!is.null(res_nb)) {
    cat("\n=== What the SEG baseline changed (same cells, same estimator) ===\n")
    fmt <- function(r) if (is.finite(r$meta$lambda))
      sprintf("%.1f um [%.1f, %.1f]", r$meta$lambda, r$meta$lambda_lo, r$meta$lambda_hi) else "not estimable"
    cat(sprintf("  lambda WITH SEG baseline    : %-28s (donors %d, quotable %s)\n",
                fmt(res), res$gates$n_donors_included, res$gates$quotable))
    cat(sprintf("  lambda WITHOUT SEG baseline : %-28s (donors %d, quotable %s)\n",
                fmt(res_nb), res_nb$gates$n_donors_included, res_nb$gates$quotable))
    cat(sprintf("  amplitude  with SEG %+.4f  |  without SEG %+.4f\n",
                res$amplitude, res_nb$amplitude))
    cat(sprintf("  shared-lambda anchor  with SEG %.0f um  |  without SEG %.0f um\n",
                res$shared$lambda, res_nb$shared$lambda))
    cat("  With collinear shapes, little change in lambda is expected (see the note above).\n")
  }

  if (seg_ok) {
    cat("\n=== MAGNITUDE BENCHMARK against the SEG null ===\n")
    nb <- null_benchmark_gradient(df, m, "seg_baseline")
    decay_write_null_benchmark(nb, sg$label)
    cat("\nPer-donor near-to-far drops (SD units):\n")
    print(nb$per_donor, row.names = FALSE, digits = 4)
    write.table(cbind(signature = sg$key, celltype = ct, module = m, nb$per_donor),
                file.path(args$output_dir,
                sprintf("source_data_%s_dist_%s%s_nullbenchmark.tsv", sg$file_tag, ct, suffix)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  cat("\n=== Position relative to the PHF1+ neurons themselves ===\n")
  cat(sprintf("PHF1+ mean = %.4f [%.4f, %.4f] (n = %d) | modelled background plateau = %.4f\n",
              ref$mean, ref$ci_lo, ref$ci_hi, ref$n, res$plateau))
  cat(sprintf("Fraction of the PHF1+ / background gap spanned by the decaying field = %s\n",
              if (is.finite(gap_frac)) sprintf("%.2f", gap_frac) else "NA"))
  cat("  i.e. how much of the way from background up to a tangle-bearing neuron the\n")
  cat("  nearest PHF1-NEGATIVE neurons already are. > 1 would mean the field overshoots\n")
  cat("  the PHF1+ cells themselves; << 1 means a small fraction of the gap.\n")
  if (sg$key == "phf1") {
    cat("\nNOTE (selection): the geneset is SELECTED from the PHF1+/PHF1- contrast within\n")
    cat("  this celltype. The decay is fitted on PHF1-NEGATIVE cells only; the PHF1+ level is\n")
    cat("  elevated by construction and the gap fraction above inherits that. The Otero sets\n")
    cat("  are external and do not.\n")
  }
  cat("\nsessionInfo():\n"); print(sessionInfo())
  sink()

  decay_rows[[key]] <<- cbind(
    decay_summary_row(res, signature = sg$key, celltype = ct, module = m,
                      log_lmm_padj = log_padj,
                      seg_baseline = if (seg_ok) "SEG null" else "free plateau"),
    data.frame(
      lambda_no_seg    = if (is.null(res_nb)) NA_real_ else res_nb$meta$lambda,
      lambda_no_seg_lo = if (is.null(res_nb)) NA_real_ else res_nb$meta$lambda_lo,
      lambda_no_seg_hi = if (is.null(res_nb)) NA_real_ else res_nb$meta$lambda_hi,
      amplitude_no_seg = if (is.null(res_nb)) NA_real_ else res_nb$amplitude,
      quotable_no_seg  = if (is.null(res_nb)) NA else res_nb$gates$quotable,
      # magnitude benchmark: how many times steeper than the stably-expressed-gene null
      bench_sig_drop_sd = if (is.null(nb)) NA_real_ else nb$sig[["mean"]],
      bench_seg_drop_sd = if (is.null(nb)) NA_real_ else nb$seg[["mean"]],
      bench_excess_sd    = if (is.null(nb)) NA_real_ else nb$excess[["mean"]],
      bench_excess_lo    = if (is.null(nb)) NA_real_ else nb$excess[["lo"]],
      bench_excess_hi    = if (is.null(nb)) NA_real_ else nb$excess[["hi"]],
      bench_ratio        = if (is.null(nb)) NA_real_ else nb$ratio,
      bench_donors_sig_steeper = if (is.null(nb)) NA_integer_ else nb$n_donors_same_dir))
  cat(sprintf("  Done %s/%s%s: lambda = %s um [%s, %s] | %s\n", sg$key, ct, suffix,
              if (is.finite(mm$lambda)) sprintf("%.0f", mm$lambda) else "NA",
              if (is.finite(mm$lambda_lo)) sprintf("%.0f", mm$lambda_lo) else "NA",
              if (is.finite(mm$lambda_hi)) sprintf("%.0f", mm$lambda_hi) else "NA",
              res$gates$verdict))
  invisible(res)
}

# Model-free actual-data plot for one celltype x signature. The rolling-mean line is
# drawn solid+thick / dashed+thin using the LOG-model significance (log_sig), matching
# the log fitted plot's annotation. log_sig may be NULL (log model failed) -> dashed.
run_ct_raw <- function(ct, df_full_all, ref, sg, log_sig) {
  m <- sg$key
  df <- if (is.finite(MAX_DIST)) df_full_all[df_full_all$dist_to_phf1_um <= MAX_DIST, ] else df_full_all
  df <- df[!is.na(df[[m]]) & !is.na(df$dist_to_phf1_um), ]
  x_cap <- if (is.finite(MAX_DIST)) MAX_DIST else as.numeric(quantile(df$dist_to_phf1_um, 0.99))
  sig_flag <- isTRUE(log_sig$significant)

  raw_df <- data.frame(cell_id = df$cell_id, dist_to_phf1_um = df$dist_to_phf1_um, score = df[[m]])
  grid   <- seq(0, x_cap, length.out = 200)
  roll   <- roll_mean(raw_df$dist_to_phf1_um, raw_df$score, grid, HALF_WINDOW)
  roll_df <- data.frame(dist_to_phf1_um = grid, roll)

  p <- make_raw_plot(raw_df, roll_df, ref, x_cap, sg$colour, "Module score", sig_flag)
  save_fixed_panel(p, file.path(args$output_dir,
                   sprintf("plot_%s_dist_%s_raw.pdf", sg$file_tag, ct)))

  pos <- CT_DATA[[ct]]; pos <- pos[pos$phf1_pos & !is.na(pos[[m]]), ]
  raw_out <- rbind(
    data.frame(signature = sg$key, celltype = ct, cell_id = raw_df$cell_id,
               dist_to_phf1_um = raw_df$dist_to_phf1_um, score = raw_df$score, phf1_pos = FALSE),
    data.frame(signature = sg$key, celltype = ct, cell_id = pos$cell_id,
               dist_to_phf1_um = pos$dist_to_phf1_um, score = pos[[m]], phf1_pos = TRUE))
  write.table(raw_out, file.path(args$output_dir,
              sprintf("source_data_%s_dist_%s_raw.tsv", sg$file_tag, ct)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(signature = sg$key, celltype = ct, window_um = args$window_um,
                         log_significant = sig_flag, log_lmm_padj = log_sig$lmm_padj %||% NA_real_, roll_df),
              file.path(args$output_dir, sprintf("source_data_%s_dist_%s_raw_rollmean.tsv", sg$file_tag, ct)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Done %s/%s raw plot: %d PHF1-negative cells, window=%g um\n",
              sg$key, ct, nrow(raw_df), args$window_um))
}

# Integrated per-celltype overlay: rolling means + 95% CI (no raw points) for every
# integrated signature + each signature's PHF1+ reference (mean, 95% CI), colour-coded.
#   mode = "z"  : z-score each signature (mean 0/SD 1 on PHF1-negative cells within cap)
#                 so the trends are shape-comparable regardless of AddModuleScore scale.
#   mode = "raw": plot RAW AddModuleScore values on a shared native-scale axis.
build_integrated <- function(ct, neg_df, keys, mode = "z") {
  tag  <- if (mode == "z") "zscore" else "raw"
  ylab <- if (mode == "z") "Module score (z, per signature)" else "Module score"   # raw: shortened axis label
  x_cap <- if (is.finite(MAX_DIST)) MAX_DIST else as.numeric(quantile(neg_df$dist_to_phf1_um, 0.99))
  neg_cap <- if (is.finite(MAX_DIST)) neg_df[neg_df$dist_to_phf1_um <= MAX_DIST, ] else neg_df
  grid <- seq(0, x_cap, length.out = 200)
  pos  <- CT_DATA[[ct]]; pos <- pos[pos$phf1_pos, ]

  roll_rows <- list(); ref_rows <- list(); par_rows <- list()
  for (k in keys) {
    lab <- SIGNATURES[[k]]$label
    y <- neg_cap[[k]]; ok <- !is.na(y) & !is.na(neg_cap$dist_to_phf1_um)
    mu <- mean(y[ok]); sdv <- stats::sd(y[ok])
    if (!is.finite(sdv) || sdv == 0) next
    tf <- function(v) if (mode == "z") (v - mu) / sdv else v      # z-score or identity
    r  <- roll_mean(neg_cap$dist_to_phf1_um[ok], tf(y[ok]), grid, HALF_WINDOW)
    roll_rows[[k]] <- data.frame(signature = lab, dist_to_phf1_um = grid, r)   # roll_mean, sem, n_window
    pv <- pos[[k]][!is.na(pos[[k]])]; n <- length(pv)             # PHF1+ reference (same transform)
    if (n >= 1) {
      zp <- tf(pv); mp <- mean(zp); sem <- if (n > 1) stats::sd(zp) / sqrt(n) else NA_real_
      ref_rows[[k]] <- data.frame(signature = lab, ref_mean = mp,
                                  ci_lo = mp - 1.96 * sem, ci_hi = mp + 1.96 * sem, n_phf1_pos = n)
    }
    par_rows[[k]] <- data.frame(signature = lab, key = k, neg_mean = mu, neg_sd = sdv, n_neg_cells = sum(ok))
  }
  if (length(roll_rows) == 0) return(invisible(NULL))

  sig_levels <- vapply(keys, function(k) SIGNATURES[[k]]$label, character(1))
  sig_cols   <- setNames(vapply(keys, function(k) SIGNATURES[[k]]$colour, character(1)), sig_levels)
  roll_df <- bind_rows(roll_rows); roll_df$signature <- factor(roll_df$signature, levels = sig_levels)
  ref_df  <- bind_rows(ref_rows);  ref_df$signature  <- factor(ref_df$signature,  levels = sig_levels)
  ref_band <- ref_df %>% tidyr::crossing(dist_to_phf1_um = c(0, x_cap))

  p <- ggplot() +
    geom_ribbon(data = ref_band, aes(dist_to_phf1_um, ymin = ci_lo, ymax = ci_hi, fill = signature),
                alpha = 0.10, colour = NA) +
    geom_line(data = ref_band, aes(dist_to_phf1_um, ref_mean, colour = signature),
              linetype = "dotted", linewidth = 0.4) +
    # 95% CI of the rolling mean (mean +/- 1.96*SEM), behind the trend lines
    geom_ribbon(data = roll_df, aes(dist_to_phf1_um, ymin = roll_mean - 1.96 * sem,
                                    ymax = roll_mean + 1.96 * sem, fill = signature),
                alpha = 0.15, colour = NA) +
    geom_line(data = roll_df, aes(dist_to_phf1_um, roll_mean, colour = signature), linewidth = 0.7) +
    scale_colour_manual(values = sig_cols, name = NULL, drop = FALSE) +
    scale_fill_manual(values = sig_cols, guide = "none", drop = FALSE) +
    labs(x = XLAB, y = ylab) +
    coord_cartesian(xlim = c(0, x_cap)) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.text = element_text(size = 6), legend.key.width = grid::unit(14, "pt"),
          legend.key.height = grid::unit(9, "pt"), legend.position = "right",
          legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))

  save_fixed_panel(p, file.path(args$output_dir, sprintf("plot_integrated_%s_dist_%s.pdf", tag, ct)))

  write.table(roll_df %>% mutate(celltype = ct, mode = mode, window_um = args$window_um) %>%
                dplyr::select(celltype, mode, signature, window_um, dist_to_phf1_um, roll_mean, sem, n_window),
              file.path(args$output_dir, sprintf("source_data_integrated_%s_dist_%s_rollmean.tsv", tag, ct)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(ref_df %>% mutate(celltype = ct, mode = mode) %>%
                dplyr::select(celltype, mode, signature, ref_mean, ci_lo, ci_hi, n_phf1_pos),
              file.path(args$output_dir, sprintf("source_data_integrated_%s_dist_%s_phf1ref.tsv", tag, ct)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  sink(file.path(args$output_dir, sprintf("stats_integrated_%s_dist_%s.txt", tag, ct)))
  cat("Integrated ptau-signature trends vs distance to PHF1+ neuron -- celltype:", ct,
      " | mode:", mode, "\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  if (mode == "z") {
    cat("Each signature Z-SCORED (mean 0 / SD 1) on this celltype's PHF1-negative cells within the ",
        MAX_DIST, " um cap; then a ", args$window_um, " um rolling mean + 95% CI over distance (no raw\n", sep = "")
    cat("points). PHF1+ scores transformed with the SAME per-signature mean/SD and drawn as a horizontal\n")
    cat("reference (mean, 95% CI). Rolling means hover near 0; the PHF1+ reference shows the SD gap.\n")
    cat("AddModuleScore magnitudes are not comparable across signatures -> z-scoring for this overlay.\n\n")
  } else {
    cat("RAW AddModuleScore values (NOT z-scored): ", args$window_um, " um rolling mean + 95% CI over\n", sep = "")
    cat("distance for PHF1-negative cells, with each signature's PHF1+ reference (raw mean, 95% CI),\n")
    cat("on a shared native-scale axis.\n\n")
  }
  cat("Per-signature parameters (PHF1-negative cells within cap):\n")
  print(as.data.frame(bind_rows(par_rows)))
  cat("\nPHF1+ reference:\n"); print(as.data.frame(ref_df))
  cat("\nInferential stats are in the per-signature stats_*_dist_", ct, "*.txt files.\n", sep = "")
  cat("\nsessionInfo():\n"); print(sessionInfo())
  sink()

  cat(sprintf("  Done integrated %s plot %s: %d signatures\n", mode, ct, length(roll_rows)))
  p
}

##  ............................................................................
##  Run everything                                                          ####
INTEGRATED_PANELS_Z   <- list()   # z-scored overlay panels
INTEGRATED_PANELS_RAW <- list()   # raw (non-z) overlay panels

for (ct in names(CT_DATA)) {
  all_keys <- SIG_ORDER[SIG_ORDER %in% CT_SIGS[[ct]]]              # every scored signature
  int_keys <- INTEGRATED_SIGS[INTEGRATED_SIGS %in% CT_SIGS[[ct]]]  # integrated overlay subset
  df_all <- prep_celltype_df(CT_DATA[[ct]])
  neg_df <- df_all[!df_all$phf1_pos, ]
  cat(sprintf("\n#### %s: %d PHF1-negative cells (pre-cap) | signatures: %s\n",
              ct, nrow(neg_df), paste(all_keys, collapse = ", ")))

  for (k in all_keys) {
    sig <- SIGNATURES[[k]]
    ref <- phf1_reference(CT_DATA[[ct]], k)
    if (is.null(ref)) { cat("  [skip] ", ct, "/", k, ": no PHF1+ cells for reference\n", sep = ""); next }
    phf1ref_rows[[paste(ct, k)]] <- tibble(signature = k, celltype = ct, n_phf1_pos = ref$n,
                                           mean = ref$mean, sd = ref$sd, sem = ref$sem,
                                           ci_lo = ref$ci_lo, ci_hi = ref$ci_hi)
    log_res <- NULL
    # MODEL TYPES RUN. Only the log-distance fit is reported, so "linear" and "gam"
    # are not run (to cut run time) -- add them to this vector if needed; their code
    # paths in run_ct_model() are intact.
    # "decay" must come after "log": it is gated on the log model's verdict.
    for (mt in c("log", "decay")) {
      res <- run_ct_model(ct, neg_df, ref, mt, sig)
      if (mt == "log") log_res <- res
    }
    run_ct_raw(ct, neg_df, ref, sig, log_res)
  }

  # Descriptive per-donor view, once per celltype across all its signatures.
  run_ct_bydonor(ct, neg_df, all_keys)

  INTEGRATED_PANELS_Z[[ct]]   <- build_integrated(ct, neg_df, int_keys, mode = "z")
  INTEGRATED_PANELS_RAW[[ct]] <- build_integrated(ct, neg_df, int_keys, mode = "raw")
}

# Written-once tables.
if (length(lmm_coef_rows) > 0)
  write.table(bind_rows(lmm_coef_rows),
              file.path(args$output_dir, "stats_modulescore_lmm_coeffs.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
if (length(phf1ref_rows) > 0)
  write.table(bind_rows(phf1ref_rows),
              file.path(args$output_dir, "source_data_modulescore_phf1ref.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

##  ............................................................................
##  Length-constant headline table + forest                                 ####
if (length(decay_rows) > 0) {
  dtab <- bind_rows(decay_rows)
  write.table(dtab, file.path(args$output_dir, "stats_modulescore_decay_lambda.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n#### Length constants (lambda, um) ####\n")
  print(as.data.frame(dtab %>% dplyr::select(dplyr::any_of(
    c("signature", "celltype", "lambda_um", "lambda_lo", "lambda_hi", "d95_um",
      "n_donors_included", "quotable", "verdict")))), row.names = FALSE)

  ftab <- dtab[is.finite(dtab$lambda_um), , drop = FALSE]
  if (nrow(ftab) > 0) {
    ftab$row_label <- sprintf("%s | %s", ftab$signature, ftab$celltype)
    ftab$series    <- ftab$signature
    dpts <- bind_rows(decay_donors)
    dpts <- dpts[dpts$included, , drop = FALSE]
    if (nrow(dpts)) dpts$row_label <- sprintf("%s | %s", dpts$signature, dpts$celltype)
    sig_cols <- vapply(SIGNATURES, function(s) s$colour, character(1))
    fp <- decay_lambda_forest(ftab, dpts, cap_um = if (is.finite(MAX_DIST)) MAX_DIST else NULL,
                              colours = sig_cols)
    if (!is.null(fp)) {
      fp <- fp + forest_theme + theme(legend.position = "bottom")
      ggsave(file.path(args$output_dir, "plot_decay_lambda_forest.pdf"), fp,
             width = 9.0, height = 1.2 + 0.42 * nrow(ftab), units = "cm", device = "pdf")
      write.table(ftab %>% dplyr::select(dplyr::any_of(
                    c("row_label", "signature", "celltype", "module", "lambda_um",
                      "lambda_lo", "lambda_hi", "d95_um", "d95_lo", "d95_hi",
                      "n_donors_included", "n_donors_total", "quotable", "verdict"))),
                  file.path(args$output_dir, "source_data_decay_lambda_forest.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("Wrote lambda forest:", nrow(ftab), "row(s)\n")
    }
  }
}

# Combined integrated figures (per-celltype panels side by side), one per mode.
save_integrated_combined <- function(panel_list, path) {
  panels <- Filter(Negate(is.null), panel_list)
  if (length(panels) < 1) return(invisible())
  combo <- Reduce(`+`, panels) + patchwork::plot_layout(nrow = 1, guides = "collect")
  ggsave(path, combo, width = 9.0 * length(panels), height = 6, units = "cm", device = "pdf")
  cat("Wrote combined integrated figure:", basename(path), "(", length(panels), "panel(s) )\n")
}
save_integrated_combined(INTEGRATED_PANELS_Z,
                         file.path(args$output_dir, "plot_integrated_zscore_COMBINED.pdf"))
save_integrated_combined(INTEGRATED_PANELS_RAW,
                         file.path(args$output_dir, "plot_integrated_raw_COMBINED.pdf"))

cat("\nAll outputs written to:", args$output_dir, "\n")
cat("Done.\n")
