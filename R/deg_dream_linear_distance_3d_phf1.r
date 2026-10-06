#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# deg_dream_linear_distance_3d_phf1.r
#
# Figure panels: S1C
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# deg_dream_linear_distance_3d_phf1.r
#
# ============================================================================
# OUT-OF-PLANE TANGLE CORRECTION -- SENSITIVITY FORK
# ============================================================================
# A fork of deg_dream_linear_distance_phf1.r that differs in EXACTLY TWO ways:
#
#   --dist_table  a TSV of cell_id -> distance that REPLACES dist_to_phf1_um
#   --dist_col    which column of it to use (one augmented draw per fit)
#
# The model, the covariates, the gene filter, the transform, the cell filters
# and the output contract are untouched. Anything else that differs is a bug.
#
# WHY. dist_to_phf1_um is an IN-PLANE distance from a single 2D section, so the
# nearest tangle may lie above or below the imaged tissue and some neurons
# labelled PHF1-negative are genuinely tangle-bearing. Stage 1
# (R/phf1_3d_distance.R) adds H = n_obs (1 - p_detect)/p_detect hidden anchors
# per donor, drawn from the eligible PHF1-negative pool and weighted by
# kernel-smoothed observed anchor intensity, and recomputes the distance
# against the augmented set through the SAME compute_dist_to_phf1_um().
# The unseen-plane tables (t3_unseen_plane_anchors() in R/tangle_3d_utils.R)
# are passed in the same way; run_deg_linear_distance_3d.sh uses those for
# Fig. S1C.
#
# ----------------------------------------------------------------------------
# HOW TO READ THE OUTPUT
# ----------------------------------------------------------------------------
# Hidden anchors are a POSITION-WEIGHTED RANDOM DRAW, not an identification of
# which cells are mislabelled (~131 anchors are added). The augmentation
# therefore attenuates the distance slope through anchor-set dilution on its
# own, independently of the underlying effect. The size of that attenuation is
# measured by simulation:
#
#     ARM B DILUTION CONSTANT = slope_aug / slope_obs = 0.906 [0.896, 0.915]
#
# from R/test_tangle_3d_utils.R section 4, stored in
# results/tangle_3d/armB_dilution_constant.tsv and echoed into every log below.
# It applies to the hidden-anchor tables, not to the unseen-plane tables.
#
# ----------------------------------------------------------------------------
# THE IDENTITY CHECK IS THE POINT OF dist3d/aug_p1p000_identity.tsv
# ----------------------------------------------------------------------------
# At p_detect = 1 no hidden anchors are drawn, so that table holds the canonical
# distances -- sourced from the SCE's own stored column, agreeing to 2e-13 um
# (about one ULP; not bit-identical because the value round-trips through a
# text file). This script must therefore reproduce
#   deg/de_linear_distance/<ct>/<ct>_dist_to_phf1_um_scaled.tsv
# to within floating-point noise: expect logFC to agree to ~1e-12 and the gene
# ranking and significance calls to be unchanged. That is the proof the fork
# changed nothing but the distance. RUN IT FIRST (array index 1).
#
# Original header follows.
# ----------------------------------------------------------------------------
# deg_dream_linear_distance_phf1.r
#
# Cell-level mixed model DEG using a linear term on continuous distance
# from nearest PHF1+ neuron, among PHF1-NEGATIVE cells only.
#
# Model (per gene):
#   ~ dist_to_phf1_um_scaled + nUMI_log + percent_neg + Sex + Age + PMI
#     + (1 | sample_id)
#
# With --log_distance the distance predictor becomes log(dist_to_phf1_um)
# scaled to SD units. This is biologically motivated when the effect is
# expected to decay non-linearly (e.g. diffusion-based signals): a unit
# change near the source is more meaningful than the same unit change far
# away. Natural log is used (not log1p): distance to the nearest PHF1+ neuron
# is strictly > 0 for the modelled PHF1-negative cells, so no offset is needed.
#
# Within-sample variation in distance is preserved by operating at the
# cell level. The random effect (1|sample_id) accounts for donor clustering.
#
# Normalisation: CPM (no TMM). Weights via voomWithDreamWeights.
# Structure, covariates, output format and file naming follow the spline-distance
# DEG variant, with the spline replaced by a single linear term.
#
# Prediction curves and per-gene plateau distances are written in the same
# format as the spline variant for downstream compatibility.
#
# Expects annotated SCE files produced by label_phf1_neighbours.r
# (colData column: dist_to_phf1_um)
#
# Gene filter: --min_pct_cells sets the detection-rate cut (default 5, i.e.
# 0.05 * n_cells). The default writes to the unsuffixed output directory; any
# other value writes to a _pct<X> directory instead, so threshold experiments
# never overwrite the default outputs.
#
# Exclusion radius: --min_distance_um (default 0) additionally drops PHF1-negative
# cells CLOSER than the given radius to the nearest PHF1+ neuron, on top of the
# --max_dist_um upper cap. It is a sensitivity analysis: the cells immediately
# adjacent to a tangle are the ones most exposed to segmentation bleed-through and
# to the optical PHF1 spillover floor, so refitting at 10/20/30/50 um asks whether
# the distance slope survives once they are removed. Distance is a cell-inclusion
# filter here, exactly like --max_dist_um; the model formula is unchanged.
# 0 reproduces the canonical run exactly and writes to deg/de_linear_distance/;
# any other radius writes to deg/de_linear_distance_min<R>um/.
#
# Two things are needed for the coefficients to be comparable across radii:
#   (a) The detection-rate gene filter is computed on the FULL PHF1-negative,
#       within-cap set BEFORE the radius filter, and that background gene set is
#       reused at every radius. Computed after the radius filter it would shrink
#       with the radius and the comparison would be confounded by a moving gene
#       set. A fingerprint of the background is printed so identity across logs is
#       checkable by eye.
#   (b) dist_sd is recomputed from the retained cells, so logFC (per SD of
#       transformed distance) is NOT comparable across radii. The back-divided
#       logFC_per_unit / CI.L_per_unit / CI.R_per_unit columns ARE, and are what
#       any cross-radius comparison must use. dist_sd is echoed per run and
#       carried in the results table alongside min_distance_um and n_cells.
#
# Cortical layer is not modelled here, and there is no laminar annotation column
# in the CosMx colData; the radius diagnostics probe for one, report its absence,
# and skip.
#
# Use run_deg_linear_distance.sh (canonical 5%, celltype x radius array)
#  or run_deg_linear_distance_pct.sh (alternative threshold)

##  ............................................................................
##  Load packages                                                           ####
library(tidyverse)
library(argparse)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(Matrix)
library(edgeR)
library(limma)
library(reformulas)
library(variancePartition)
library(BiocParallel)

.muffleReformulasMigrationWarning <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("findbars|nobars", conditionMessage(w)) &&
        grepl("reformulas", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
}

##  ............................................................................
##  Parse command-line arguments                                            ####
parser <- ArgumentParser()
required <- parser$add_argument_group("Required", "required arguments")

required$add_argument(
  "--sce",
  help     = "Path to annotated _sce_neighbours.qs file for one celltype",
  required = TRUE
)
required$add_argument(
  "--confounding_vars",
  help     = "Comma-separated confounding variables e.g. Sex,Age,PMI",
  metavar  = "Sex,Age,PMI",
  required = TRUE
)
required$add_argument(
  "--output_dir",
  help     = "Base output directory",
  required = TRUE
)
required$add_argument(
  "--ncores",
  help    = "Cores for dream [default: 4]",
  default = 4,
  type    = "integer"
)
required$add_argument(
  "--max_dist_um",
  help    = "Maximum distance from PHF1+ neuron to include in model (um) [default: 1000]",
  default = 1000,
  type    = "double"
)
parser$add_argument(
  "--log_distance",
  help   = "Apply natural-log transformation to distance before scaling [default: FALSE]",
  action = "store_true",
  default = FALSE
)
parser$add_argument(
  "--min_pct_cells",
  help    = paste("Gene filter: minimum percent of modelled cells in which a gene",
                  "must be detected (counts > 0) to be tested [default: 5]"),
  default = 5,
  type    = "double"
)
parser$add_argument(
  "--dist_table",
  help    = paste("TSV with a cell_id column and one or more distance columns.",
                  "Replaces dist_to_phf1_um for the cells it names. Produced by",
                  "R/phf1_3d_distance.R. Required."),
  required = TRUE
)
parser$add_argument(
  "--dist_col",
  help    = "Column of --dist_table to use as the distance [default: rep_01]",
  default = "rep_01"
)
parser$add_argument(
  "--config_tag",
  help    = "Short tag for the output dir, e.g. p0p750_intensity. Required.",
  required = TRUE
)
parser$add_argument(
  "--dilution_tsv",
  help    = "Arm B dilution constant from R/test_tangle_3d_utils.R",
  default = "results/tangle_3d/armB_dilution_constant.tsv"
)
parser$add_argument(
  "--min_distance_um",
  help    = paste("Exclusion radius: drop PHF1-negative cells closer than this to the",
                  "nearest PHF1+ neuron (um). 0 = no exclusion, reproducing the",
                  "canonical run exactly [default: 0]"),
  default = 0,
  type    = "double"
)

args <- parser$parse_args()

celltype <- gsub("_sce_neighbours\\.qs$", "", basename(args$sce))
args$confounding_vars <- trimws(strsplit(args$confounding_vars, ",")[[1]])

if (args$min_pct_cells < 0 || args$min_pct_cells >= 100) {
  stop("--min_pct_cells must be in [0, 100).")
}

if (args$min_distance_um < 0) {
  stop("--min_distance_um must be >= 0.")
}
if (args$min_distance_um >= args$max_dist_um) {
  stop("--min_distance_um must be < --max_dist_um (", args$max_dist_um,
       " um) or no cells are retained.")
}

# The canonical 5% run keeps writing to de_linear_distance/ so existing outputs
# and downstream scripts are untouched. Any other threshold is routed to a
# suffixed directory (e.g. de_linear_distance_pct3p5/), matching the variant
# convention already used by de_linear_distance_raw/.
# The variant marker is part of the BASE name and applied unconditionally --
# a replaced distance column is a different analysis, not a dial of this one,
# so it belongs in the base name and not in the _pct / _min suffixes. The guard
# below refuses to write outside that namespace, so a
# mis-set --output_dir can never overwrite deg/de_linear_distance/.
DEFAULT_MIN_PCT <- 5
dir_tag <- if (isTRUE(all.equal(args$min_pct_cells, DEFAULT_MIN_PCT))) {
  sprintf("de_linear_distance_3d_%s", args$config_tag)
} else {
  sprintf("de_linear_distance_3d_%s_pct%s", args$config_tag,
          gsub("\\.", "p", format(args$min_pct_cells, trim = TRUE, scientific = FALSE)))
}

# Same principle for the exclusion radius: 0 leaves the tag (and therefore the
# canonical output path) untouched, any other radius gets its own directory so a
# sweep can never overwrite the reference run. A decimal radius renders as e.g.
# min10p5um, matching the pct3p5 convention above.
DEFAULT_MIN_DIST_UM <- 0
if (args$min_distance_um > DEFAULT_MIN_DIST_UM) {
  dir_tag <- paste0(dir_tag, sprintf("_min%sum",
    gsub("\\.", "p", format(args$min_distance_um, trim = TRUE, scientific = FALSE))))
}

stopifnot(
  "dir_tag escaped the 3D variant namespace - refusing to write" =
    grepl("^de_linear_distance_3d_", dir_tag)
)
outdir <- file.path(args$output_dir, dir_tag, celltype)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cli::cli_text("Output dir: {.path {outdir}}")

##  ............................................................................
##  Load SCE and subset to PHF1-negative cells only                        ####
cli::cli_text("Reading {.strong {celltype}}")

sce <- qs::qread(args$sce)
rownames(sce) <- rowData(sce)$gene
sce$sample_id <- droplevels(as.factor(sce$sample_id))

stopifnot(
  "dist_to_phf1_um column missing — run label_phf1_neighbours.r first" =
    "dist_to_phf1_um" %in% colnames(colData(sce))
)

##  ............................................................................
##  THE ONE SUBSTANTIVE CHANGE: swap in the augmented-anchor distance         ####
# Done HERE -- after the SCE is read, before ANY cell filter -- so that every
# downstream step (the max_dist cap, the gene background, dist_sd, the model)
# sees the corrected column and nothing has to be told about it twice.
dil <- tryCatch({
  tb <- read.delim(args$dilution_tsv, stringsAsFactors = FALSE)
  setNames(tb$value, tb$quantity)
}, error = function(e) NULL)
if (is.null(dil)) {
  stop("Arm B dilution constant not found at ", args$dilution_tsv,
       ". Run: Rscript R/test_tangle_3d_utils.R --reps 40\n",
       "The constant is required: the ",
       "augmentation attenuates the slope on its own and this records by how much.")
}

dt <- read.delim(args$dist_table, stringsAsFactors = FALSE, check.names = FALSE)
if (!"cell_id" %in% colnames(dt))
  stop("--dist_table has no cell_id column: ", args$dist_table)
if (!args$dist_col %in% colnames(dt))
  stop("--dist_col '", args$dist_col, "' not in ", args$dist_table,
       ". Available: ", paste(setdiff(colnames(dt), "cell_id"), collapse = ", "))

cid <- if ("cell_id" %in% colnames(colData(sce))) sce$cell_id else colnames(sce)
new_d <- dt[[args$dist_col]][match(cid, dt$cell_id)]
old_d <- sce$dist_to_phf1_um

# Every modelled cell must be covered. A silent NA here would drop cells from
# this arm but not from the canonical one, and the two would stop being
# comparable for reasons that never appear in the log.
was_modelled <- !is.na(old_d) & old_d <= args$max_dist_um
n_uncov <- sum(was_modelled & is.na(new_d))
if (n_uncov > 0)
  stop(n_uncov, " cells that the canonical run models are absent from ",
       args$dist_table, ". The distance table and the SCE disagree on the cell set.")

cat("\n=== Augmented-anchor distance substitution ===\n")
cat("Distance table: ", args$dist_table, "\n")
cat("Column:         ", args$dist_col, "\n")
cat("Config tag:     ", args$config_tag, "\n")
cat(sprintf("Cells matched:   %d / %d\n", sum(!is.na(new_d)), length(cid)))
cat(sprintf("Median distance: %.2f -> %.2f um (%.2f%% shorter) over %d modelled cells\n",
            median(old_d[was_modelled]), median(new_d[was_modelled]),
            100 * (1 - median(new_d[was_modelled]) / median(old_d[was_modelled])),
            sum(was_modelled)))
if (any(new_d[was_modelled] <= 0, na.rm = TRUE))
  stop("Non-positive augmented distance; log() requires dist > 0. ",
       "Fix upstream in Stage 1, do not absorb it with an offset.")

cat(sprintf("\nARM B DILUTION CONSTANT: %.3f [%.3f, %.3f]\n",
            dil[["dilution_slope_aug_over_obs"]], dil[["dilution_ci_lo"]],
            dil[["dilution_ci_hi"]]))
cat("  Augmentation attenuates the distance slope by about this much ON ITS OWN,\n")
cat("  whether or not the field effect is real, because hidden anchors are a\n")
cat("  position-weighted random draw and not an identification of mislabelled\n")
cat("  cells. It applies to hidden-anchor tables, not to unseen-plane tables.\n\n")

# THE GENE BACKGROUND IS FROZEN ON THE CANONICAL CELL SET.
# The augmented distances are shorter, so more cells fall inside --max_dist_um:
# 24,807 against the canonical 24,672 at p_detect = 0.75. That change in the
# MODELLED SET is correct and intended -- under a corrected anchor set those
# cells genuinely are within 1000 um of a tangle. But letting the DETECTION-RATE
# GENE FILTER move with it is not: the denominator shifts, the tested gene set
# shifts with it (2,462 against the canonical 2,468 here), and any cross-arm
# comparison of coefficients is then confounded by a moving gene set.
#
# This is exactly the trap the canonical script's own header describes for the
# exclusion radius, and it takes the same remedy: compute the background ONCE,
# on the canonical (observed-distance) within-cap set, and reuse it at every
# p_detect. The fingerprint printed later must therefore be IDENTICAL across
# every arm including the identity run -- that is checkable by eye from the logs.
canon_in_cap <- !is.na(old_d) & old_d <= args$max_dist_um
frozen_counts <- counts(sce)[, canon_in_cap, drop = FALSE]
frozen_n      <- sum(canon_in_cap)
frozen_keep   <- Matrix::rowSums(frozen_counts > 0) >=
                   (args$min_pct_cells / 100) * frozen_n
FROZEN_BKG    <- rownames(sce)[frozen_keep]
rm(frozen_counts); gc(verbose = FALSE)
cat(sprintf("Gene background FROZEN on the canonical set: %d genes from %d cells\n",
            length(FROZEN_BKG), frozen_n))
cat("  (identical across every arm by construction -- compare the fingerprint below)\n\n")

# PHF1+ cells keep their NA so they stay excluded, exactly as canonically.
sce$dist_to_phf1_um <- ifelse(is.na(old_d), NA_real_, new_d)

# Stage 1: exclude PHF1+ cells (dist_to_phf1_um is NA for them) and restrict to
# cells within --max_dist_um of a tangle. This is the reference set — invariant
# across exclusion radii, and the set the gene background is computed on.
MAX_DIST_UM <- args$max_dist_um
MIN_DIST_UM <- args$min_distance_um
sce_full <- sce[, !is.na(sce$dist_to_phf1_um) & sce$dist_to_phf1_um <= MAX_DIST_UM]

cat("Cells after excluding PHF1+ and > ", MAX_DIST_UM, " um:", ncol(sce_full), "\n")
cat("Distance summary (um):\n")
print(summary(sce_full$dist_to_phf1_um))
cat("Cells per sample:\n")
print(table(sce_full$sample_id))

##  ............................................................................
##  Gene filtering — computed PRE-radius, on the full within-cap set          ####
# Detection-rate filter: a gene must be detected in at least --min_pct_cells
# percent of the modelled (PHF1-negative, within-cap) cells. Note the panel is
# targeted, so at low thresholds essentially the whole panel is retained; the
# per-gene detection rate is tightly packed around 0.04-0.05, which makes the
# gene count strongly sensitive to the exact cut in the 3-5% range.
#
# This block sits BEFORE the exclusion-radius filter deliberately. Computed after
# it, the denominator would shrink with the radius and the tested gene set would
# differ between runs, confounding any cross-radius comparison of coefficients.
# The fingerprint below is printed so identity of the background across radii can
# be confirmed from the logs alone.
raw_counts <- counts(sce_full)
n_cells    <- ncol(sce_full)
expr_cells <- rowSums(raw_counts > 0)
keep       <- expr_cells >= (args$min_pct_cells / 100) * n_cells
# NOT rownames(sce_full)[keep] -- see the frozen-background block above. The
# recomputed set is reported alongside so any drift is visible, but the FROZEN
# one is what is tested.
recomputed_bkg <- rownames(sce_full)[keep]
bkg_genes      <- FROZEN_BKG
if (!identical(sort(recomputed_bkg), sort(bkg_genes)))
  cat(sprintf("NOTE: background recomputed on the augmented set would be %d genes (%d differ); using the FROZEN %d.\n",
              length(recomputed_bkg),
              length(union(setdiff(recomputed_bkg, bkg_genes),
                           setdiff(bkg_genes, recomputed_bkg))),
              length(bkg_genes)))

cat("Total genes:", nrow(sce_full), "\n")
cat(sprintf("Gene filter: detected in >= %s%% of %d cells (>= %.1f cells)\n",
            format(args$min_pct_cells, trim = TRUE), n_cells,
            (args$min_pct_cells / 100) * n_cells))
cat("Genes kept:", length(bkg_genes), "\n")
cat(sprintf("Background genes (computed PRE-radius, invariant across radii): %d\n",
            length(bkg_genes)))
cat(sprintf("Background fingerprint: n=%d first=%s last=%s sum_nchar=%d\n",
            length(bkg_genes), bkg_genes[1], bkg_genes[length(bkg_genes)],
            sum(nchar(bkg_genes))))

if (length(bkg_genes) == 0) {
  stop("No genes passed the detection-rate filter. Exiting.")
}

##  ............................................................................
##  Stage 2: exclusion radius                                                ####
# Drop cells CLOSER than --min_distance_um to the nearest PHF1+ neuron. At the
# default 0 nothing is dropped (distance is strictly > 0 for these cells), so the
# canonical run is reproduced exactly.
donor_levels <- levels(droplevels(as.factor(sce_full$sample_id)))
n_before     <- ncol(sce_full)
dist_before  <- sce_full$dist_to_phf1_um
donor_before <- table(factor(as.character(sce_full$sample_id), levels = donor_levels))

keep_cell <- sce_full$dist_to_phf1_um >= MIN_DIST_UM
sce_neg   <- sce_full[, keep_cell]
sce_neg$sample_id <- droplevels(as.factor(sce_neg$sample_id))

n_after     <- ncol(sce_neg)
donor_after <- table(factor(as.character(sce_neg$sample_id), levels = donor_levels))

cat("\n=== Exclusion radius diagnostics ===\n")
cat(sprintf("Exclusion radius (min_distance_um): %s um%s\n",
            format(MIN_DIST_UM, trim = TRUE),
            if (MIN_DIST_UM == 0) "  (no exclusion — canonical run)" else ""))
cat(sprintf("Cells before radius filter: %d\n", n_before))
cat(sprintf("Cells after  radius filter: %d\n", n_after))
cat(sprintf("Cells dropped:              %d (%.2f%%)\n",
            n_before - n_after,
            if (n_before > 0) 100 * (n_before - n_after) / n_before else NA_real_))

cat("\nCells per sample_id, before vs after:\n")
donor_tab <- data.frame(
  sample_id  = donor_levels,
  before     = as.integer(donor_before[donor_levels]),
  after      = as.integer(donor_after[donor_levels]),
  row.names  = NULL,
  stringsAsFactors = FALSE
)
donor_tab$dropped     <- donor_tab$before - donor_tab$after
donor_tab$pct_dropped <- ifelse(donor_tab$before > 0,
                                round(100 * donor_tab$dropped / donor_tab$before, 2), NA_real_)
print(donor_tab, row.names = FALSE)

cat("\nDistance summary (um) BEFORE radius filter:\n")
print(summary(dist_before))
cat("Distance summary (um) AFTER radius filter:\n")
print(summary(sce_neg$dist_to_phf1_um))

# Braak composition before vs after. Braak is donor-level metadata and is NOT a
# covariate in this model — this is a composition diagnostic only, so a missing
# column is a warning, never an abort.
if ("Braak" %in% colnames(colData(sce_full))) {
  braak_before <- as.character(colData(sce_full)$Braak)
  braak_after  <- as.character(colData(sce_neg)$Braak)
  braak_lv     <- sort(unique(braak_before[!is.na(braak_before)]))
  braak_tab <- data.frame(
    Braak      = braak_lv,
    before     = as.integer(table(factor(braak_before, levels = braak_lv))[braak_lv]),
    after      = as.integer(table(factor(braak_after,  levels = braak_lv))[braak_lv]),
    row.names  = NULL,
    stringsAsFactors = FALSE
  )
  braak_tab$prop_before <- round(braak_tab$before / max(1, sum(braak_tab$before)), 4)
  braak_tab$prop_after  <- round(braak_tab$after  / max(1, sum(braak_tab$after)),  4)
  cat("\nBraak composition before vs after radius filter:\n")
  print(braak_tab, row.names = FALSE)
} else {
  cat("\nBraak column absent from colData — skipping Braak composition diagnostic.\n")
}

# Laminar / layer composition. There is no laminar annotation column in the CosMx
# colData (layer appears only inside the celltype label, which is constant within
# a run), so this normally reports absence and skips. Layer is deliberately never
# added to the model formula in this project — if a future object does carry such
# a column it is tabulated here as a diagnostic only.
LAYER_CANDIDATES <- c("layer", "Layer", "laminar", "cortical_layer",
                      "cortical_depth", "depth", "depth_rel", "scaled_Y")
layer_cols <- intersect(LAYER_CANDIDATES, colnames(colData(sce_full)))
if (length(layer_cols) == 0) {
  cat("\nNo laminar/layer annotation column in colData - skipping laminar",
      "composition diagnostic\n(layer is deliberately not modelled in this project).\n")
  cat("  probed:", paste(LAYER_CANDIDATES, collapse = ", "), "\n")
} else {
  for (lc in layer_cols) {
    cat(sprintf("\n%s composition before vs after radius filter:\n", lc))
    lv <- sort(unique(as.character(colData(sce_full)[[lc]])))
    print(data.frame(
      level  = lv,
      before = as.integer(table(factor(as.character(colData(sce_full)[[lc]]), levels = lv))[lv]),
      after  = as.integer(table(factor(as.character(colData(sce_neg)[[lc]]),  levels = lv))[lv]),
      row.names = NULL
    ), row.names = FALSE)
  }
  cat("(diagnostic only — NOT added to the model formula)\n")
}
cat("\n")

if (ncol(sce_neg) < 50) {
  stop("Fewer than 50 PHF1-negative cells — insufficient for dream. Exiting.")
}

n_samples <- length(unique(sce_neg$sample_id))
if (n_samples < 3) {
  stop("Fewer than 3 donors — random effect not estimable. Exiting.")
}

# Per-donor sufficiency. At larger radii a sparse celltype legitimately runs out
# of cells in some donors; aborting here with the donors named is the intended
# outcome, not a bug. Fitting on 2 donors or on donors contributing a handful of
# cells would produce a coefficient that is not comparable with the other radii.
MIN_CELLS_PER_DONOR <- 20
thin <- donor_levels[as.integer(donor_after[donor_levels]) < MIN_CELLS_PER_DONOR]
if (length(thin) > 0) {
  stop("Exclusion radius ", MIN_DIST_UM, " um leaves fewer than ", MIN_CELLS_PER_DONOR,
       " cells in ", length(thin), " donor(s): ",
       paste(sprintf("%s (n=%d)", thin, as.integer(donor_after[thin])), collapse = ", "),
       ". Expected for sparse celltypes at large radii — this radius is not",
       " estimable for ", celltype, ". Exiting.")
}

n_donors_ok <- sum(as.integer(donor_after[donor_levels]) >= MIN_CELLS_PER_DONOR)
if (n_donors_ok < 3) {
  stop("Exclusion radius ", MIN_DIST_UM, " um leaves only ", n_donors_ok,
       " donor(s) with >= ", MIN_CELLS_PER_DONOR, " cells (",
       paste(donor_levels[as.integer(donor_after[donor_levels]) >= MIN_CELLS_PER_DONOR],
             collapse = ", "),
       "); the random effect is not estimable on fewer than 3. Exiting.")
}

##  ............................................................................
##  Prepare colData covariates                                              ####
cd <- as.data.frame(colData(sce_neg))

cd$nUMI_log    <- log2(colSums(counts(sce_neg)) + 1)
cd$nGene       <- colSums(counts(sce_neg) > 0)
cd$percent_neg <- as.numeric(cd$percent.neg)

cd$Age <- as.numeric(scale(as.numeric(cd$Age)))
cd$PMI <- as.numeric(scale(as.numeric(cd$PMI)))

# Scale distance to SD units so the coefficient is interpretable as
# logCPM change per SD of distance, and to aid numerical convergence.
# If --log_distance, apply natural log first (biologically motivated for
# diffusion-decay signals; also compresses the long right tail). distance to the
# nearest PHF1+ neuron is strictly > 0 for the modelled PHF1-negative cells, so
# natural log is used (no log1p offset needed); a non-positive value is a data
# error and is caught here.
if (args$log_distance) {
  if (any(cd$dist_to_phf1_um <= 0, na.rm = TRUE))
    stop("Non-positive dist_to_phf1_um encountered; log() requires dist > 0.")
  cd$dist_to_phf1_um_transformed <- log(cd$dist_to_phf1_um)
  cat("Distance transformation: log(dist_to_phf1_um)\n")
  cat("log(distance) summary:\n")
  print(summary(cd$dist_to_phf1_um_transformed))
} else {
  cd$dist_to_phf1_um_transformed <- cd$dist_to_phf1_um
  cat("Distance transformation: none (raw um)\n")
}
dist_sd <- sd(cd$dist_to_phf1_um_transformed, na.rm = TRUE)
cd$dist_to_phf1_um_scaled <- cd$dist_to_phf1_um_transformed / dist_sd
cat(sprintf("Distance SD (after transformation): %.4f (used for scaling)\n", dist_sd))

colData(sce_neg) <- DataFrame(cd, row.names = rownames(cd))

##  ............................................................................
##  Apply the pre-radius gene background                                     ####
# bkg_genes was computed above on sce_full (pre-radius), so the tested gene set is
# identical at every exclusion radius.
sce_filt <- sce_neg[bkg_genes, ]

##  ............................................................................
##  Build DGEList                                                           ####
# CPM normalisation: library-size scaling only, no TMM adjustment.
dge <- DGEList(counts = counts(sce_filt))

##  ............................................................................
##  Build formula                                                           ####
conf_str <- paste(args$confounding_vars, collapse = " + ")
qc_str   <- "nUMI_log + percent_neg"

form_str <- sprintf(
  "~ dist_to_phf1_um_scaled + %s + %s + (1|sample_id)",
  qc_str, conf_str
)

cat("Dream formula:", form_str, "\n")
form <- as.formula(form_str)

form_vars <- setdiff(all.vars(form), "sample_id")
missing   <- setdiff(form_vars, colnames(colData(sce_filt)))
if (length(missing) > 0) {
  stop("Missing colData columns for formula: ", paste(missing, collapse = ", "))
}

pd <- as.data.frame(colData(sce_filt))

##  ............................................................................
##  Diagnostic: collinearity between distance and QC covariates             ####
cat("\nDiagnostic | Pearson r between dist_to_phf1_um and QC/confounders:\n")
diag_vars <- intersect(c("nUMI_log", "percent_neg", args$confounding_vars), colnames(pd))
for (v in diag_vars) {
  if (is.numeric(pd[[v]])) {
    r <- cor(pd$dist_to_phf1_um, pd[[v]], use = "complete.obs")
    cat(sprintf("  %-20s r = %+.3f%s\n", v, r,
                ifelse(abs(r) > 0.3, "  <-- WARNING: collinear", "")))
  }
}
cat("\n")

##  ............................................................................
##  voom + dream                                                            ####
param <- MulticoreParam(args$ncores, progressbar = TRUE)

cat("Running voomWithDreamWeights...\n")
vobj <- .muffleReformulasMigrationWarning(voomWithDreamWeights(
  counts  = dge,
  formula = form,
  data    = pd,
  BPPARAM = param
))

cat("Running dream...\n")
fit_dream <- .muffleReformulasMigrationWarning(dream(
  exprObj = vobj,
  formula = form,
  data    = pd,
  BPPARAM = param
))
fit_dream <- eBayes(fit_dream)

##  ............................................................................
##  Extract results: linear distance coefficient                            ####
cat("Coefficient names in fit:\n")
print(colnames(coef(fit_dream)))

dist_coef <- grep("dist_to_phf1_um_scaled", colnames(coef(fit_dream)), value = TRUE)
if (length(dist_coef) == 0) {
  stop("dist_to_phf1_um_scaled not found in fit. Available: ",
       paste(colnames(coef(fit_dream)), collapse = ", "))
}

res <- topTable(
  fit_dream,
  coef          = dist_coef,
  number        = nrow(fit_dream),
  adjust.method = "BH",
  confint       = TRUE
)

res <- res %>%
  mutate(
    gene         = rownames(.),
    contrast     = dist_coef,
    model        = form_str,
    # Back-transform: divide by dist_sd to undo scaling.
    # If --log_distance: units are logCPM per unit log(um).
    # If raw distance:   units are logCPM per um.
    logFC_per_unit = logFC / dist_sd,
    CI.L_per_unit  = CI.L  / dist_sd,
    CI.R_per_unit  = CI.R  / dist_sd,
    distance_scale = ifelse(args$log_distance, "log_um", "um"),
    # Provenance for the exclusion-radius sweep. Appended at the END of the
    # column order so every existing name-based reader is unaffected.
    # dist_sd is carried because logFC is expressed per SD of transformed
    # distance and dist_sd changes with the radius: logFC is therefore NOT
    # comparable across radii, logFC_per_unit is.
    min_distance_um = args$min_distance_um,
    n_cells         = ncol(sce_filt),
    dist_sd         = dist_sd,
    # 3D-variant provenance, appended at the END so every existing name-based
    # reader is unaffected.
    config_tag      = args$config_tag,
    dist_table      = basename(args$dist_table),
    dist_col        = args$dist_col,
    arm             = "B_augmented_anchors",
    dilution_const  = unname(dil[["dilution_slope_aug_over_obs"]])
  ) %>%
  dplyr::rename(pval = P.Value, padj = adj.P.Val) %>%
  dplyr::select(gene, logFC, logFC_per_unit, CI.L, CI.R,
                CI.L_per_unit, CI.R_per_unit,
                pval, padj, AveExpr, t, B, contrast, model, distance_scale,
                min_distance_um, n_cells, dist_sd,
                config_tag, dist_table, dist_col, arm, dilution_const)

file_name <- file.path(outdir,
  sprintf("%s_%s.tsv", celltype, gsub("[^A-Za-z0-9]", "_", dist_coef)))
write.table(res, file_name,
            col.names = TRUE, row.names = FALSE, sep = "\t")
cat("Written:", file_name, "\n")

cat("Significant genes (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")
cat("Significant genes (padj < 0.1):",  sum(res$padj < 0.1,  na.rm = TRUE), "\n")

##  ............................................................................
##  Predict fitted expression curves across a distance grid                ####
# Write _spline_joint_F.tsv and _spline_predictions.tsv in the same format
# as deg_dream_spline_distance_phf1.r for downstream compatibility.

SIG_THRESHOLD <- 0.1

top_genes <- res %>%
  filter(padj < SIG_THRESHOLD) %>%
  arrange(pval) %>%
  pull(gene)

cat("Genes with padj <", SIG_THRESHOLD, ":", length(top_genes), "\n")

# _spline_joint_F.tsv — F = t^2 for a single-coefficient test
res_joint <- res %>%
  transmute(
    gene    = gene,
    F       = t^2,
    pval    = pval,
    padj    = padj,
    AveExpr = AveExpr,
    test    = "linear_distance_t",
    model   = model
  )

joint_file <- file.path(outdir, paste0(celltype, "_spline_joint_F.tsv"))
write.table(res_joint, joint_file,
            col.names = TRUE, row.names = FALSE, sep = "\t")
cat("Written joint (linear) results:", joint_file, "\n")

if (length(top_genes) > 0) {

  # Grid spans the observed distance range. Start at the minimum observed
  # distance (not 0) so log(dist) is defined under --log_distance.
  dist_min  <- min(pd$dist_to_phf1_um, na.rm = TRUE)
  dist_max  <- max(pd$dist_to_phf1_um, na.rm = TRUE)
  dist_grid <- seq(dist_min, dist_max, length.out = 200)

  # Convert grid to the same transformed scale used in the model
  if (args$log_distance) {
    dist_grid_transformed <- log(dist_grid)
  } else {
    dist_grid_transformed <- dist_grid
  }
  dist_grid_scaled <- dist_grid_transformed / dist_sd

  # For a linear model: fitted = intercept + slope_per_scaled_unit * dist_grid_scaled
  # slope_per_scaled_unit = logFC (from topTable, in scaled-distance units)
  intercept_vec <- coef(fit_dream)[top_genes, "(Intercept)", drop = TRUE]
  slope_vec     <- setNames(res$logFC, res$gene)[top_genes]

  pred_long <- lapply(top_genes, function(g) {
    data.frame(
      gene            = g,
      dist_to_phf1_um = dist_grid,           # always raw µm for plotting
      fitted_logCPM   = intercept_vec[g] + slope_vec[g] * dist_grid_scaled,
      stringsAsFactors = FALSE
    )
  })
  pred_long <- bind_rows(pred_long)

  pred_file <- file.path(outdir, paste0(celltype, "_spline_predictions.tsv"))
  write.table(pred_long, pred_file,
              col.names = TRUE, row.names = FALSE, sep = "\t")
  cat("Written linear prediction curves:", pred_file, "\n")

  ##  ..........................................................................
  ##  Plateau distances                                                     ####
  # For a linear model the derivative is constant so the 10%-of-max-derivative
  # criterion used in the spline script always returns 0. Instead report the
  # distance where cumulative absolute change reaches 90% of total change.

  plateau_df <- pred_long %>%
    group_by(gene) %>%
    arrange(dist_to_phf1_um) %>%
    mutate(
      abs_change_cum = cumsum(abs(c(0, diff(fitted_logCPM)))),
      total_change   = max(abs_change_cum)
    ) %>%
    filter(total_change > 0, abs_change_cum >= 0.9 * total_change) %>%
    summarise(
      # dplyr::first must be qualified: SingleCellExperiment is attached after
      # tidyverse, so S4Vectors::first masks it and a bare first() aborts with
      # "unable to find an inherited method ... signature x = numeric".
      plateau_dist_um = dplyr::first(dist_to_phf1_um),
      .groups = "drop"
    ) %>%
    left_join(res_joint %>% dplyr::select(gene, F, padj), by = "gene") %>%
    arrange(padj)

  plateau_file <- file.path(outdir, paste0(celltype, "_plateau_distances.tsv"))
  write.table(plateau_df, plateau_file,
              col.names = TRUE, row.names = FALSE, sep = "\t")
  cat("Written per-gene plateau distances:", plateau_file, "\n")

  cat("\nPlateau distance summary across significant genes (um):\n")
  print(summary(plateau_df$plateau_dist_um))

  suggested_threshold <- median(plateau_df$plateau_dist_um, na.rm = TRUE)
  cat(sprintf("\nSuggested binary threshold (median 90%% cumulative change): %.1f um\n",
              suggested_threshold))

} else {
  cat("No significant genes at padj <", SIG_THRESHOLD, "— prediction step skipped.\n")
}

##  ............................................................................
##  Session summary                                                         ####
cat("\n--- Session summary ---\n")
cat("Celltype:              ", celltype, "\n")
cat("Cells modelled:        ", ncol(sce_filt), "\n")
cat("Max distance (um):     ", MAX_DIST_UM, "\n")
cat("Min distance (exclusion radius, um):", MIN_DIST_UM, "\n")
cat("Cells before radius filter:", n_before, "\n")
cat("Cells after radius filter: ", n_after, "\n")
cat("Donors after radius filter:", n_samples, "\n")
cat("Gene filter (pct cells):", args$min_pct_cells, "\n")
cat("Background genes (pre-radius):", length(bkg_genes), "\n")
cat("Output dir tag:        ", dir_tag, "\n")
cat("Genes tested:          ", nrow(fit_dream), "\n")
cat("Formula:               ", form_str, "\n")
cat("Donors:                ", n_samples, "\n")
cat("log_distance:          ", args$log_distance, "\n")
cat("dist SD (transformed): ", round(dist_sd, 4), "\n")
cat("  NOTE: dist_sd is recomputed from the retained cells, so logFC (per SD) is\n")
cat("  NOT comparable across exclusion radii; logFC_per_unit is.\n")
cat("Significant (padj<0.1):", sum(res$padj < 0.1, na.rm = TRUE), "\n")
cat("Config tag:            ", args$config_tag, "\n")
cat("Distance table:        ", basename(args$dist_table), "[", args$dist_col, "]\n")
cat(sprintf("Arm B dilution constant: %.3f -- the attenuation expected from the\n",
            dil[["dilution_slope_aug_over_obs"]]))
cat("  hidden-anchor augmentation alone; not applicable to unseen-plane tables.\n")
