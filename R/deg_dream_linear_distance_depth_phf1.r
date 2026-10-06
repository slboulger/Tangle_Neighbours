#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# deg_dream_linear_distance_depth_phf1.r
#
# Figure panels: S1B
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# deg_dream_linear_distance_depth_phf1.r
#
# DEPTH-ADJUSTED variant of deg_dream_linear_distance_phf1.r. Identical in every
# respect except that one fixed effect is added: depth_s, the standardised
# relative position of each cell along its sample's manual pia->white-matter
# compass axis.
#
# Cell-level mixed model DEG using a linear term on continuous distance
# from nearest PHF1+ neuron, among PHF1-NEGATIVE cells only.
#
# Model (per gene):
#   ~ dist_to_phf1_um_scaled + depth_s + nUMI_log + percent_neg + Sex + Age + PMI
#     + (1 | sample_id)
#
# WHY
# ---
# A sensitivity analysis asking whether the distance gradient survives conditioning
# on large-scale position within the section. It is the CosMx counterpart of the IMC
# depth ladder in imc_covariate_sensitivity.R (Figure S8A), which standardises
# depth within the modelled cell set and reports it as a labelled sensitivity.
# docs/MODELS.md lists depth_rel among the image-derived covariates "used as
# labelled sensitivities", so no model specification changes.
#
# WHAT depth_s IS
# ---------------
# depth_rel comes from a manual per-sample compass direction for L1 (`l1_dir` in
# cortical_depth_utils.R). depth_s is a LARGE-SCALE WITHIN-SECTION POSITION
# covariate, not a cortical-layer annotation. It is correlated with the distance
# predictor (|rho| 0.02 to 0.57 across samples), so it tests whether the distance
# gradient survives conditioning on where in the section a cell sits.
# The laminar validation table (Spearman rho between the laminar excitatory ordinal
# and depth_rel) is printed into every log by cat_depth_orientation().
#
# depth_rel is computed against a FROZEN all-cell bounding box, not from the
# modelled subset -- see the header of cortical_depth_utils.R for why that is
# load-bearing. Only the z-scoring is set-dependent, and that is a linear
# rescaling which cannot move the distance coefficient or its p-value.
#
# RELATION TO docs/MODELS.md
# --------------------------
# Set 1 carries no depth term, and docs/MODELS.md (Set 3) lists scaled_Y/depth_s as a
# sensitivity term, never in a headline. This script is a labelled sensitivity
# (Fig. S1B), as depth is in IMC Set 3.
#
# It writes ONLY to deg/de_linear_distance_depth*/ (enforced by a stopifnot on
# dir_tag below) and never to deg/de_linear_distance/, which backs Fig. 3A, 3B
# and 4C-E.
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
# 0 uses the canonical cell set and writes to deg/de_linear_distance_depth/;
# any other radius writes to deg/de_linear_distance_depth_min<R>um/.
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
# There is no laminar annotation column in the CosMx colData; the radius
# diagnostics probe for one and report its absence, and depth_s is derived from
# coordinates rather than from any such column.
#
# Use run_deg_linear_distance_depth.sh (canonical 5%, one job per celltype).
# The contrast against the unadjusted parent run is collate_depth_covariate.r.

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

# Locate this script's own directory so cortical_depth_utils.R is found whether
# the job is launched from the project root or from R/. Same shim idiom as
# palettes.R's loader in the collate scripts.
.get_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  m  <- grep("^--file=", ca, value = TRUE)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
  getwd()
}
source(file.path(.get_script_dir(), "cortical_depth_utils.R"))

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
  "--min_distance_um",
  help    = paste("Exclusion radius: drop PHF1-negative cells closer than this to the",
                  "nearest PHF1+ neuron (um). 0 = no exclusion, reproducing the",
                  "canonical run exactly [default: 0]"),
  default = 0,
  type    = "double"
)
parser$add_argument(
  "--coords_csv",
  help    = paste("All-cell coordinate export used to freeze the per-sample depth",
                  "bounding box [default: phf1_v2/PHF1/seu_coords.csv]"),
  default = DEPTH_COORDS_DEFAULT
)
parser$add_argument(
  "--depth_bbox",
  help    = paste("Cache path for the frozen depth bounding box",
                  "[default: phf1_v2/cortical_depth_bbox.tsv]"),
  default = DEPTH_BBOX_DEFAULT
)
parser$add_argument(
  "--parent_results",
  help    = paste("Optional: directory of the unadjusted parent run for this celltype",
                  "(e.g. phf1_v2/deg/de_linear_distance). When given, the cell set and",
                  "dist_sd are asserted identical and a paired comparison TSV is written."),
  default = ""
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

# The default 5% run writes to de_linear_distance_depth/. Any other threshold is
# routed to a suffixed directory (e.g. de_linear_distance_depth_pct3p5/).
DEFAULT_MIN_PCT <- 5
dir_tag <- if (isTRUE(all.equal(args$min_pct_cells, DEFAULT_MIN_PCT))) {
  "de_linear_distance_depth"
} else {
  sprintf("de_linear_distance_depth_pct%s",
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

# _depth is part of the BASE tag, applied unconditionally, so no combination of
# --min_pct_cells and --min_distance_um can route this variant's output into the
# canonical deg/de_linear_distance/ tree.
stopifnot(
  "dir_tag escaped the depth variant namespace - refusing to write" =
    grepl("^de_linear_distance_depth", dir_tag)
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
bkg_genes  <- rownames(sce_full)[keep]

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

##  ............................................................................
##  Cortical-depth covariate                                                ####
# depth_rel is projected against the FROZEN all-cell bounding box, never against
# this celltype's PHF1-negative subset -- otherwise depth_rel would silently mean
# something different in every celltype and at every exclusion radius. See the
# header of cortical_depth_utils.R.
cat("\n--- Cortical depth ---\n")
depth_bbox <- load_cortical_depth_bbox(coords_csv = args$coords_csv,
                                       cache      = args$depth_bbox)
cat(sprintf("Frozen depth bbox: %s (%d samples)\n", args$depth_bbox, nrow(depth_bbox)))

DEPTH_COORD_X <- "x_slide_mm"
DEPTH_COORD_Y <- "y_slide_mm"
missing_coord <- setdiff(c(DEPTH_COORD_X, DEPTH_COORD_Y), colnames(cd))
if (length(missing_coord))
  stop("Coordinate column(s) missing from colData, cannot compute depth: ",
       paste(missing_coord, collapse = ", "))
if (anyNA(cd[[DEPTH_COORD_X]]) || anyNA(cd[[DEPTH_COORD_Y]]))
  stop("NA in ", DEPTH_COORD_X, "/", DEPTH_COORD_Y, " for modelled cells.")

cd <- add_cortical_depth_from_bbox(cd, depth_bbox,
                                   coord_x = DEPTH_COORD_X, coord_y = DEPTH_COORD_Y)
stopifnot(
  "depth_rel has NA after the bbox projection"  = !anyNA(cd$depth_rel),
  "depth_rel outside [0,1] after the projection" =
    min(cd$depth_rel) >= 0 && max(cd$depth_rel) <= 1
)

# Standardised within the modelled cell set, matching imc_covariate_sensitivity.R
# (Figure S8A). A linear rescaling of a frozen quantity, so it cannot move the
# distance coefficient or its p-value -- only depth_s's own units.
depth_rel_sd <- sd(cd$depth_rel, na.rm = TRUE)
if (!isTRUE(depth_rel_sd > 0))
  stop("depth_rel has zero variance in the modelled cells; depth_s is undefined.")
cd$depth_s <- as.numeric(scale(cd$depth_rel))
cat(sprintf("depth_rel: mean %.4f, sd %.4f, range [%.4f, %.4f]\n",
            mean(cd$depth_rel), depth_rel_sd, min(cd$depth_rel), max(cd$depth_rel)))

# Orientation validation, printed on every run. Warns, never stops: the axis is a
# usable position covariate whether or not it carries laminar signal.
#
# Computed on the ALL-CELL frame, not on `cd`. This SCE holds exactly one
# celltype, so the laminar ordinal would have a single level here and every rho
# would come back NA. The axis is a property of the cohort, so it is measured
# once on the same frame the bbox came from and cached beside it.
depth_valid <- cortical_depth_validation(coords_csv = args$coords_csv,
                                         bbox       = depth_bbox)
cat_depth_orientation(depth_valid)

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
  "~ dist_to_phf1_um_scaled + depth_s + %s + %s + (1|sample_id)",
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
diag_vars <- intersect(c("nUMI_log", "percent_neg", "depth_s", args$confounding_vars),
                       colnames(pd))
for (v in diag_vars) {
  if (is.numeric(pd[[v]])) {
    r <- cor(pd$dist_to_phf1_um, pd[[v]], use = "complete.obs")
    cat(sprintf("  %-20s r = %+.3f%s\n", v, r,
                ifelse(abs(r) > 0.3, "  <-- WARNING: collinear", "")))
  }
}
cat("\n")

##  ............................................................................
##  Diagnostic: how much of the distance predictor does depth absorb?       ####
# This is the depth analogue of the FOV clone's between_r2() decomposition, and
# it is the number to quote alongside any attenuation of the distance
# coefficient. depth_s and the distance predictor are BOTH cell-level fixed
# effects here, so the relevant quantity is simply their shared variance: R^2 of
# dist ~ depth_s is the fraction of the predictor the adjustment can remove.
cat("Diagnostic | overlap between depth_s and the distance predictor:\n")
depth_rho_dist <- suppressWarnings(
  cor(pd$dist_to_phf1_um_scaled, pd$depth_s, method = "spearman", use = "complete.obs"))
depth_r_dist   <- cor(pd$dist_to_phf1_um_scaled, pd$depth_s, use = "complete.obs")
depth_r2_dist  <- depth_r_dist^2

# WITHIN-DONOR is the quantity that matters, and it is the HEADLINE here.
# The model carries (1|sample_id), so it conditions on the donor-demeaned
# predictors; the pooled correlation mixes in between-donor differences the
# random intercept has already absorbed. Measured on this cohort the pooled
# figure understates the within-donor R^2 by 2x to 48x and for Oligo it even
# flips SIGN (+0.052 pooled vs -0.078 within). Quoting the pooled number as
# "how much of the predictor depth can remove" is therefore wrong.
.demean <- function(v, g) { g <- as.character(g); v - tapply(v, g, mean)[g] }
dz_w <- .demean(pd$dist_to_phf1_um_scaled, pd$sample_id)
ds_w <- .demean(pd$depth_s,                pd$sample_id)
depth_r_within   <- cor(dz_w, ds_w)
depth_r2_within  <- depth_r_within^2
depth_rho_within <- suppressWarnings(cor(dz_w, ds_w, method = "spearman"))

cat(sprintf("  WITHIN-DONOR Pearson r   = %+.4f   <-- the overlap the model sees\n",
            depth_r_within))
cat(sprintf("  WITHIN-DONOR R^2         = %.5f  (%.2f%% of the predictor is shared with\n",
            depth_r2_within, 100 * depth_r2_within))
cat(sprintf("                             depth; %.2f%% is independent of it)\n",
            100 * (1 - depth_r2_within)))
cat(sprintf("  WITHIN-DONOR Spearman rho= %+.4f\n", depth_rho_within))
cat(sprintf("  [pooled, for reference only: r = %+.4f, R^2 = %.5f, rho = %+.4f]\n",
            depth_r_dist, depth_r2_dist, depth_rho_dist))
if (depth_r2_within > 2 * depth_r2_dist)
  cat("  NOTE: the pooled figure understates the within-donor overlap; use the\n",
      "        within-donor value when reading any attenuation.\n", sep = "")
if (sign(depth_r_within) != sign(depth_r_dist))
  cat("  NOTE: pooled and within-donor overlap have OPPOSITE SIGN (Simpson's\n",
      "        reversal across donors). Only the within-donor value is meaningful\n",
      "        for a model carrying (1|sample_id).\n", sep = "")
cat("  per-sample Spearman rho:\n")
for (s in sort(unique(as.character(pd$sample_id)))) {
  i <- which(as.character(pd$sample_id) == s)
  rr <- if (length(i) > 2 && sd(pd$depth_s[i]) > 0)
          suppressWarnings(cor(pd$dist_to_phf1_um_scaled[i], pd$depth_s[i],
                               method = "spearman")) else NA_real_
  cat(sprintf("    %-10s n = %6d   rho = %+.3f\n", s, length(i), rr))
}
if (abs(depth_r_within) > 0.7)
  cat("  *** WARNING: depth_s and distance are strongly collinear; the split of the\n",
      "      shared gradient between them is not identified. ***\n", sep = "")
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
    # Depth provenance, appended at the END so every name-based reader of the
    # canonical table is unaffected and the two arms stay column-compatible.
    depth_rho_dist  = depth_rho_dist,
    depth_r2_dist   = depth_r2_dist,
    depth_r_within  = depth_r_within,
    depth_r2_within = depth_r2_within
  ) %>%
  dplyr::rename(pval = P.Value, padj = adj.P.Val) %>%
  dplyr::select(gene, logFC, logFC_per_unit, CI.L, CI.R,
                CI.L_per_unit, CI.R_per_unit,
                pval, padj, AveExpr, t, B, contrast, model, distance_scale,
                min_distance_um, n_cells, dist_sd,
                depth_rho_dist, depth_r2_dist, depth_r_within, depth_r2_within)

file_name <- file.path(outdir,
  sprintf("%s_%s.tsv", celltype, gsub("[^A-Za-z0-9]", "_", dist_coef)))
write.table(res, file_name,
            col.names = TRUE, row.names = FALSE, sep = "\t")
cat("Written:", file_name, "\n")

##  ............................................................................
##  The depth coefficient itself                                            ####
# Written separately rather than merged into the distance table, which must keep
# the parent's column contract. This is what says whether depth_s did anything at
# all, and it feeds the depth-specific columns of collate_depth_covariate.r.
depth_coef <- grep("^depth_s$", colnames(coef(fit_dream)), value = TRUE)
if (length(depth_coef) != 1)
  stop("Expected exactly one depth_s coefficient, found: ",
       paste(colnames(coef(fit_dream)), collapse = ", "))

res_depth <- topTable(fit_dream, coef = depth_coef, number = nrow(fit_dream),
                      adjust.method = "BH", confint = TRUE) %>%
  mutate(gene = rownames(.), contrast = depth_coef, model = form_str,
         n_cells = ncol(sce_filt),
         depth_rho_dist = depth_rho_dist, depth_r2_dist = depth_r2_dist,
         depth_r_within = depth_r_within, depth_r2_within = depth_r2_within) %>%
  dplyr::rename(pval = P.Value, padj = adj.P.Val) %>%
  dplyr::select(gene, logFC, CI.L, CI.R, pval, padj, AveExpr, t, B,
                contrast, model, n_cells, depth_rho_dist, depth_r2_dist,
                depth_r_within, depth_r2_within)

depth_file <- file.path(outdir, sprintf("%s_depth_s.tsv", celltype))
write.table(res_depth, depth_file, col.names = TRUE, row.names = FALSE, sep = "\t")
cat("Written:", depth_file, "\n")

# One-row summary consumed by collate_depth_covariate.r, mirroring the role of
# the FOV clone's <ct>_re_variance_summary.tsv.
depth_summary <- data.frame(
  celltype             = celltype,
  n_cells              = ncol(sce_filt),
  n_genes              = nrow(res_depth),
  depth_rho_dist       = depth_rho_dist,
  depth_r_dist         = depth_r_dist,
  depth_r2_dist        = depth_r2_dist,
  depth_rho_within     = depth_rho_within,
  depth_r_within       = depth_r_within,
  depth_r2_within      = depth_r2_within,
  depth_rel_sd         = depth_rel_sd,
  median_abs_t_depth   = median(abs(res_depth$t), na.rm = TRUE),
  median_abs_t_dist    = median(abs(res$t), na.rm = TRUE),
  n_sig_depth_0.1      = sum(res_depth$padj < 0.1,  na.rm = TRUE),
  n_sig_depth_0.05     = sum(res_depth$padj < 0.05, na.rm = TRUE),
  pct_sig_depth_0.1    = round(100 * mean(res_depth$padj < 0.1, na.rm = TRUE), 2),
  laminar_rho_pooled   = depth_valid$pooled_rho,
  laminar_axis_warn    = isTRUE(depth_valid$warn),
  dist_sd              = dist_sd,
  min_distance_um      = args$min_distance_um,
  stringsAsFactors     = FALSE
)
depth_summary_file <- file.path(outdir, sprintf("%s_depth_summary.tsv", celltype))
write.table(depth_summary, depth_summary_file,
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
cat("Written:", depth_summary_file, "\n")
cat(sprintf("Genes significant for depth_s (padj < 0.1): %d of %d (%.1f%%)\n",
            depth_summary$n_sig_depth_0.1, nrow(res_depth),
            depth_summary$pct_sig_depth_0.1))

cat("Significant genes (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")
cat("Significant genes (padj < 0.1):",  sum(res$padj < 0.1,  na.rm = TRUE), "\n")

##  ............................................................................
##  Integrity + paired comparison against the unadjusted parent run         ####
# depth_rel is defined for every cell, so adding depth_s cannot drop a cell: the
# modelled set must be IDENTICAL to the parent's. That identity is what lets
# collate_depth_covariate.r contrast the two arms without refitting the parent
# here, so it is asserted rather than assumed. n_cells and dist_sd are the two
# quantities that would move if the cell set had changed.
if (nzchar(args$parent_results)) {
  # Accept either the parent results DIRECTORY or the celltype's TSV directly --
  # run_deg_linear_distance_fov.sh passes a file, so both conventions are in use.
  parent_file <- if (grepl("\\.tsv$", args$parent_results)) {
    args$parent_results
  } else {
    file.path(args$parent_results, celltype,
              sprintf("%s_%s.tsv", celltype, gsub("[^A-Za-z0-9]", "_", dist_coef)))
  }
  if (!file.exists(parent_file)) {
    warning("Parent results not found, skipping the paired comparison: ", parent_file,
            call. = FALSE)
  } else {
    par_res <- read.delim(parent_file, stringsAsFactors = FALSE)
    cat("\n--- Integrity vs the unadjusted parent run ---\n")
    cat("Parent: ", parent_file, "\n")

    p_cells <- unique(par_res$n_cells); p_sd <- unique(par_res$dist_sd)
    cat(sprintf("  n_cells  parent %s | depth %d\n",
                paste(p_cells, collapse = "/"), ncol(sce_filt)))
    cat(sprintf("  dist_sd  parent %s | depth %.10f\n",
                paste(sprintf("%.10f", p_sd), collapse = "/"), dist_sd))
    if (length(p_cells) != 1 || p_cells != ncol(sce_filt))
      stop("Cell count differs from the parent run (", paste(p_cells, collapse = "/"),
           " vs ", ncol(sce_filt), "). Adding depth_s must not change the cell set; ",
           "the two arms are not comparable and the contrast would be invalid.")
    if (length(p_sd) != 1 || !isTRUE(abs(p_sd / dist_sd - 1) < 1e-9))
      stop("dist_sd differs from the parent run (", paste(p_sd, collapse = "/"),
           " vs ", dist_sd, "). The predictor is on a different scale, so logFC ",
           "is not comparable between the arms.")
    cat("  Integrity assertions PASSED: same cells, same distance scaling.\n")

    j <- res %>%
      dplyr::select(gene, y = logFC_per_unit, padj_d = padj, logFC_d = logFC) %>%
      dplyr::inner_join(
        par_res %>% dplyr::select(gene, x = logFC_per_unit,
                                  padj_p = padj, logFC_p = logFC),
        by = "gene")
    cat(sprintf("  Genes paired: %d (parent %d, depth %d)\n",
                nrow(j), nrow(par_res), nrow(res)))

    n_sig_p <- sum(j$padj_p < 0.1, na.rm = TRUE)
    n_sig_d <- sum(j$padj_d < 0.1, na.rm = TRUE)
    n_both  <- sum(j$padj_p < 0.1 & j$padj_d < 0.1, na.rm = TRUE)
    cmp <- data.frame(
      celltype              = celltype,
      n_paired              = nrow(j),
      pearson_logFC         = suppressWarnings(cor(j$x, j$y)),
      spearman_logFC        = suppressWarnings(cor(j$x, j$y, method = "spearman")),
      ols_slope_depth_on_par = unname(coef(stats::lm(y ~ x, data = j))[2]),
      median_abs_logFC_ratio = median(abs(j$y), na.rm = TRUE) /
                               median(abs(j$x), na.rm = TRUE),
      n_sig_parent          = n_sig_p,
      n_sig_depth           = n_sig_d,
      n_sig_both            = n_both,
      pct_parentsig_retained = ifelse(n_sig_p > 0, round(100 * n_both / n_sig_p, 2), NA_real_),
      jaccard_sig           = ifelse((n_sig_p + n_sig_d - n_both) > 0,
                                     round(n_both / (n_sig_p + n_sig_d - n_both), 4), NA_real_),
      n_sign_flips_among_sig = sum(j$padj_p < 0.1 & sign(j$x) != sign(j$y), na.rm = TRUE),
      depth_rho_dist        = depth_rho_dist,
      depth_r2_dist         = depth_r2_dist,
      depth_r_within        = depth_r_within,
      depth_r2_within       = depth_r2_within,
      stringsAsFactors      = FALSE
    )
    print(cmp, row.names = FALSE, digits = 4)
    cmp_file <- file.path(outdir, sprintf("%s_depth_vs_parent.tsv", celltype))
    write.table(cmp, cmp_file, col.names = TRUE, row.names = FALSE,
                sep = "\t", quote = FALSE)
    cat("Written:", cmp_file, "\n")
    cat("\n  Read the attenuation together with depth_r2_dist above.\n")
  }
} else {
  cat("\n--parent_results not given; skipping the integrity check and paired comparison.\n")
}
cat("\n")

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
  #
  # setNames(as.numeric(...)) rather than [..., drop = TRUE]: with EXACTLY ONE
  # gene in top_genes, matrix single-row-single-column subsetting drops the
  # dimnames, so intercept_vec[g] would return NA and every fitted value would
  # become NA. Forcing the names is correct for any length >= 1.
  #
  # Note the curve is drawn at depth_s = 0 (its mean, since depth_s is z-scored),
  # so these predictions are the depth-adjusted expectation for an average-depth
  # cell and remain directly comparable to the parent's.
  intercept_vec <- setNames(as.numeric(coef(fit_dream)[top_genes, "(Intercept)"]),
                            top_genes)
  slope_vec     <- setNames(res$logFC, res$gene)[top_genes]

  if (anyNA(intercept_vec) || anyNA(slope_vec)) {
    stop("NA in the prediction intercepts/slopes for ",
         sum(is.na(intercept_vec)) + sum(is.na(slope_vec)), " of ",
         2 * length(top_genes), " lookups - the fitted curves would be all-NA. ",
         "Refusing to write a silently corrupt _spline_predictions.tsv.")
  }

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
cat("--- Depth adjustment ---\n")
cat("depth_s r with distance (within-donor):", sprintf("%+.4f", depth_r_within), "\n")
cat("R^2 within-donor:         ", sprintf("%.5f", depth_r2_within),
    sprintf("(%.2f%% of the predictor is independent of depth)\n",
            100 * (1 - depth_r2_within)))
cat("  [pooled, reference only:  ", sprintf("%.5f", depth_r2_dist), "]\n")
cat("Genes sig for depth_s (padj<0.1):", depth_summary$n_sig_depth_0.1, "\n")
cat("Laminar validation rho:", sprintf("%+.3f", depth_valid$pooled_rho),
    if (isTRUE(depth_valid$warn))
      " <-- positional covariate, not a laminar annotation\n" else "\n")
