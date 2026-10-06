#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# build_cosmx_deg_modules.R
#
# Produces: Table S9 input - builds the DEG-derived gene modules
# (genesets/cosmx_deg_modules/) that phf1_module_otero_overlap_stats.R tests.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# build_cosmx_deg_modules.R
#
# Writes the CosMx DEG-derived module gene sets to disk. Thin CLI over
# R/cosmx_deg_module_utils.R, which holds the rule and documents it in full.
#
# WHAT IT BUILDS, per size N in --n_sweep:
#   phf1cat        top N by recentred t from the CATEGORICAL DEG (PHF1+ vs PHF1-)
#   distance       top N by recentred t from the LINEAR-DISTANCE DEG (up near tangles)
#   field_strict   distance minus the FULL padj-significant PHF1 set   <- headline field set
#   field_matched  distance minus the size-matched PHF1 set            <- sensitivity
#   seg draws      --seg_draws random N-gene subsets of the celltype's SEG set
#
# The selection unit is the PROBE, which is the unit the DEG actually tested, so
# N is exact in selection space. Multi-gene CosMx probes ("MZT2A/B") are expanded
# into individual symbols only when the flat gene list is written, because that
# is the form the snRNA reference needs. Both counts are recorded (n_probes and
# n_genes) - for CBLN2 exactly one multi-gene probe falls in the top 150
# of either set, so the drift is at most one gene, but it is reported not assumed.
#
# WHY field_strict IS THE HEADLINE. Subtracting only the size-matched PHF1 set
# leaves genes that ARE significantly up in tangle-bearing cells and merely
# missed the top N - calling that "field only" would overstate it. field_strict
# subtracts the full padj<0.1 PHF1 list, so nothing in it is cell-autonomous at
# the project threshold. It is smaller and that is the price.
#
# RUN (cheap - two TSVs and a 3.8 KB qs; no HPC job)
#   Rscript R/build_cosmx_deg_modules.R
#   Rscript R/build_cosmx_deg_modules.R --celltype Exc-IT-L3-5-CHGA-IL1RAPL2

suppressPackageStartupMessages({
  library(qs)
  library(argparse)
})

PROJECT_ROOTS <- c(
  "<PROJECT_ROOT>/phf1_v2",
  "<PROJECT_ROOT>/phf1_v2"
)
.root <- PROJECT_ROOTS[dir.exists(PROJECT_ROOTS)]
if (length(.root) == 0) stop("No project root found.", call. = FALSE)
setwd(.root[1])
source("R/cosmx_deg_module_utils.R")

##  ............................................................................
##  Arguments                                                               ####
parser <- ArgumentParser()
parser$add_argument("--celltype", default = "Exc-IT-L2-3-CBLN2-HOPX",
  help = "CosMx celltype to build modules for [default: Exc-IT-L2-3-CBLN2-HOPX]")
parser$add_argument("--cat_dir", default = "deg/de_cat/de_PHF1",
  help = "Categorical PHF1+ vs PHF1- DEG root [default: deg/de_cat/de_PHF1]")
parser$add_argument("--dist_dir", default = "deg/de_linear_distance",
  help = "Canonical linear-distance DEG root [default: deg/de_linear_distance]")
parser$add_argument("--padj", type = "double", default = 0.1,
  help = "FDR cut defining eligibility and the full PHF1 set [default: 0.1]")
parser$add_argument("--n_sweep", default = "25,53,100,150",
  help = paste("Comma-separated set sizes. The headline is --n_headline; the",
               "rest are the sensitivity sweep. Sizes above the number of",
               "padj-significant genes are RANK-EXTENDED and flagged as such.",
               "[default: 25,53,100,150]"))
parser$add_argument("--n_headline", type = "integer", default = 53,
  help = paste("The headline set size. Must be in --n_sweep. 53 is the number of",
               "padj<0.1 up-near-tangle genes for CBLN2, i.e. the largest fully",
               "significant matched pair. [default: 53]"))
parser$add_argument("--seg_qs", default = "seg_control/Reference_SEG_panel.qs",
  help = "SEG null gene sets, named list keyed by CosMx celltype")
parser$add_argument("--seg_draws", type = "integer", default = 100,
  help = "Number of size-matched SEG draws per sweep level [default: 100]")
parser$add_argument("--out_dir", default = "genesets/cosmx_deg_modules",
  help = "Output root [default: genesets/cosmx_deg_modules]")
parser$add_argument("--seed", type = "integer", default = 42, help = "Seed")
args <- parser$parse_args()

set.seed(args$seed)

N_SWEEP <- as.integer(trimws(strsplit(args$n_sweep, ",", fixed = TRUE)[[1]]))
if (anyNA(N_SWEEP) || any(N_SWEEP <= 0))
  stop("--n_sweep must be positive integers, got: ", args$n_sweep, call. = FALSE)
N_SWEEP <- sort(unique(N_SWEEP))
if (!args$n_headline %in% N_SWEEP)
  stop("--n_headline (", args$n_headline, ") must be one of --n_sweep (",
       paste(N_SWEEP, collapse = ", "), ")", call. = FALSE)
if (args$seg_draws < 1) stop("--seg_draws must be >= 1", call. = FALSE)

CT   <- args$celltype
ODIR <- module_dir(args$out_dir, CT)
dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)

cat("\n=== build_cosmx_deg_modules ===\n")
cat("  project root :", getwd(), "\n")
cat("  celltype     :", CT, "\n")
cat("  padj cut     :", args$padj, "\n")
cat("  size sweep   :", paste(N_SWEEP, collapse = ", "),
    " (headline", args$n_headline, ")\n")
cat("  output       :", ODIR, "\n\n")

##  ............................................................................
##  Read + recentre                                                         ####

cat_df  <- read_deg_table(cat_deg_path(args$cat_dir,   CT), "categorical")
dist_df <- read_deg_table(dist_deg_path(args$dist_dir, CT), "linear-distance")

cat_df  <- recentre_deg(cat_df,  flip = FALSE)   # positive t = up in tangle-bearing
dist_df <- recentre_deg(dist_df, flip = TRUE)    # negative t = up near tangles

diag_df <- rbind(
  diagnose_recentring(cat_df,  FALSE, args$padj, "categorical (PHF1+ vs PHF1-)"),
  diagnose_recentring(dist_df, TRUE,  args$padj, "linear distance (up near tangles)"))

cat("Recentring diagnostic:\n")
print(diag_df[, c("table", "n_tested", "median_t", "n_padj_sig",
                  "n_sig_raw_pos", "n_sig_recentred_pos", "recentring_is_noop")],
      row.names = FALSE)
if (!all(diag_df$recentring_is_noop))
  cat("\n  NOTE: recentring CHANGED the significant split for at least one table.\n",
      "  That is a departure from the CBLN2 baseline, where it is a no-op.\n",
      "  Read the provenance file before using these sets.\n", sep = "")
cat("\n")

##  ............................................................................
##  Build                                                                   ####

mods <- build_deg_modules(cat_df, dist_df, args$padj, N_SWEEP)

seg_all <- qs::qread(args$seg_qs)
if (!is.list(seg_all) || is.null(names(seg_all)))
  stop("SEG geneset must be a named list of character vectors: ", args$seg_qs,
       call. = FALSE)
if (!CT %in% names(seg_all))
  stop("SEG geneset has no entry for ", CT, "; has: ",
       paste(names(seg_all), collapse = ", "), call. = FALSE)
seg_ct <- unique(as.character(seg_all[[CT]]))
cat("SEG null pool for", CT, ":", length(seg_ct), "genes\n")

seg_draws <- seg_matched_draws(seg_ct, N_SWEEP, args$seg_draws, args$seed)

##  ............................................................................
##  Write                                                                   ####

## Expand multi-gene probes only at write time - the flat list is what the
## snRNA scorer consumes, and it needs individual symbols.
mods$meta$n_genes <- NA_integer_
for (i in seq_len(nrow(mods$meta))) {
  key <- sprintf("%s_N%d", mods$meta$set[i], mods$meta$N[i])
  g   <- expand_probe_names(mods$sets[[key]])
  mods$meta$n_genes[i] <- length(g)
  writeLines(g, module_geneset_path(args$out_dir, CT, mods$meta$set[i], mods$meta$N[i]))
}
## n_genes_probe = probes selected (exactly N for the two matched sets).
## n_genes       = symbols written, after multi-gene probe expansion.
mods$meta <- mods$meta[, c("set", "N", "n_genes_probe", "n_genes", "n_sig",
                           "all_sig", "s_floor", "s_top", "note")]

utils::write.table(mods$membership, membership_path(args$out_dir, CT),
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(seg_draws, seg_draws_path(args$out_dir, CT),
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(mods$meta,
                   file.path(ODIR, sprintf("%s_modules_summary.tsv", ct_tag(CT))),
                   sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nSet summary:\n")
print(mods$meta[, c("set", "N", "n_genes_probe", "n_genes", "n_sig", "all_sig",
                    "s_floor")], row.names = FALSE)

## Sets that will be refused downstream for being too small. MIN_GENES_PRESENT
## is 10 in every scoring script in this project, and the intersection with the
## snRNA universe can only shrink these further.
MIN_GENES_PRESENT <- 10
thin <- mods$meta[mods$meta$n_genes < MIN_GENES_PRESENT, , drop = FALSE]
if (nrow(thin)) {
  cat("\n  WARNING: below MIN_GENES_PRESENT =", MIN_GENES_PRESENT,
      "and will be SKIPPED by the scorer:\n")
  print(thin[, c("set", "N", "n_genes")], row.names = FALSE)
}

##  ............................................................................
##  Provenance                                                              ####

sink(provenance_path(args$out_dir, CT))
cat("CosMx DEG-derived module gene sets\n")
cat("==================================\n\n")
cat("Celltype   :", CT, "\n")
cat("Built      :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Project    :", getwd(), "\n")
cat("Seed       :", args$seed, "\n\n")

cat("SOURCE TABLES\n")
cat("  categorical    :", cat_deg_path(args$cat_dir,   CT), "\n")
cat("  linear distance:", dist_deg_path(args$dist_dir, CT), "\n")
cat("  SEG null pool  :", args$seg_qs, " (", length(seg_ct), " genes)\n\n", sep = "")

cat("SELECTION RULE\n")
cat("  1. Direction relative to the PANEL-WIDE MEDIAN t:\n")
cat("       categorical  s =  (t - median(t))   positive = up in tangle-bearing\n")
cat("       distance     s = -(t - median(t))   positive = up near tangles\n")
cat("  2. Rank by s descending; take the top N, N common to both sets.\n")
cat("  3. padj <", args$padj, "defines eligibility for n_sig and for the full\n")
cat("     PHF1 set used as the field_strict subtrahend.\n\n")

cat("RECENTRING DIAGNOSTIC\n")
cat("  Both tables carry a large panel-wide POSITIVE median t. The recentring is\n")
cat("  nonetheless a NO-OP FOR SET MEMBERSHIP: every gene surviving FDR has |t|\n")
cat("  far above the median, and subtracting a constant cannot reorder a\n")
cat("  one-sided ranking. It is retained as the DIRECTION rule and as the\n")
cat("  reported MAGNITUDE, which is what makes the two sets' effect floors\n")
cat("  comparable across tables with different medians. It is NOT a re-test.\n\n")
print(diag_df, row.names = FALSE)
cat("\n  recentring_is_noop == TRUE means the raw-sign and recentred-sign splits\n")
cat("  are identical, i.e. the counts above prove the claim rather than asserting it.\n\n")

cat("DISTANCE TABLE\n")
cat("  The linear model is fitted on PHF1-NEGATIVE cells within the distance cap,\n")
cat("  and its median t is +", round(diag_df$median_t[2], 4), ".\n\n", sep = "")

cat("SET SIZES\n")
print(mods$meta, row.names = FALSE)
cat("\n  n_genes_probe = probes selected (the DEG test unit; exactly N for the two\n")
cat("    matched sets). n_genes = symbols written, after expanding multi-gene\n")
cat("    CosMx probes such as MZT2A/B.\n")
cat("  all_sig = FALSE means the set is RANK-EXTENDED past padj <", args$padj, ";\n")
cat("    only n_sig of its members survive FDR, so the set as a whole is not\n")
cat("    significant.\n")
cat("  field_strict and field_matched are DERIVED and are NOT size-matched - the\n")
cat("    equal-size property holds for phf1cat and distance only.\n\n")

cat("HEADLINE\n")
hd <- mods$meta[mods$meta$N == args$n_headline, , drop = FALSE]
cat("  N =", args$n_headline, "- the largest fully significant matched pair.\n")
print(hd[, c("set", "n_genes_probe", "n_genes", "n_sig", "all_sig", "s_floor")],
      row.names = FALSE)
cat("\n  Effect-size floors at the headline N are the two s_floor values for\n")
cat("  phf1cat and distance. Their closeness is the evidence that the two sets\n")
cat("  are matched in magnitude as well as in size.\n\n")

cat("SEG NULL\n")
cat("  ", args$seg_draws, " random N-gene subsets of the ", length(seg_ct),
    "-gene ", CT, " SEG set,\n", sep = "")
cat("  per sweep level, seeded per (N, draw). The full 305-gene set is NOT used:\n")
cat("  module-score variance scales with set size, so a 305-gene null is smoother\n")
cat("  than the signatures it is compared with and is not a fair reference. The\n")
cat("  draws instead give an EMPIRICAL NULL DISTRIBUTION for the predictor\n")
cat("  coefficient once scored and modelled.\n")
cat("  Table:", seg_draws_path(args$out_dir, CT), "\n\n")

cat("FILES WRITTEN\n")
for (f in sort(list.files(ODIR))) cat("  ", f, "\n", sep = "")

cat("\nsessionInfo():\n")
print(sessionInfo())
sink()

cat("\nWrote provenance:", provenance_path(args$out_dir, CT), "\n")
cat("Done.\n\n")
