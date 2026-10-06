#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# prepare_seg_geneset.R
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
# prepare_seg_geneset.R
#
# ONE-OFF. Extracts the published human stably expressed gene (SEG) list and writes
# it as a geneset .qs in the same shape as Mancuso_2024.qs / Cameron_2024.qs, so that
# prep_sets() in the module-score scripts consumes it unchanged.
#
# SOURCE
#   Lin Y, Ghazanfar S, Strbenac D, Wang A, Patrick E, Lin DM, Speed T, Yang JYH,
#   Yang P. Evaluating stably expressed genes in single cells.
#   GigaScience 2019;8(9):giz106.  (scSEGIndex)
#   Distributed as the `segList` data object in the Bioconductor package scMerge.
#
# WHY
#   Negative control for the module-score-over-distance analyses: a gene set that is
#   by construction NOT expected to vary with pathology or with distance to the
#   nearest PHF1+ neuron. If the SEG module is flat where the reported signatures are
#   sloped, the reported gradients are not an artefact of scoring any gene set on
#   these cells.
#
# WHAT IT DOES
#   * Requires scMerge. If absent it STOPS with an install instruction rather than
#     substituting a different gene list (the control has to be the published one).
#   * INSPECTS segList before subsetting: element naming differs across scMerge
#     versions and there is a parallel Ensembl-gene-ID version of each list. The
#     human gene-SYMBOL vector is resolved by name and then validated as symbols.
#   * Prints length() and the first 20 symbols.
#   * Writes list(SEG = <symbols>) to <geneset_dir>/Lin_SEG_2019.qs (qs::qsave), plus
#     a plain-text provenance sidecar.
#
#   Rscript R/prepare_seg_geneset.R

suppressPackageStartupMessages({
  library(qs)
  library(argparse)
})

# Project root: HPC path first (the project convention), local Mac mount as fallback,
# so this can be run either on the HPC or locally on <RDS_ROOT> without editing.
PROJECT_ROOTS <- c(
  "<PROJECT_ROOT>/phf1_v2",
  "<PROJECT_ROOT>/phf1_v2"
)
.root <- PROJECT_ROOTS[dir.exists(PROJECT_ROOTS)]
if (length(.root) == 0)
  stop("No project root found. Tried:\n  ", paste(PROJECT_ROOTS, collapse = "\n  "), call. = FALSE)
setwd(.root[1])
cat("Project root:", getwd(), "\n")

set.seed(42)

##  ............................................................................
##  Arguments                                                               ####
parser <- ArgumentParser()
# Project-local, alongside the other project gene-set dirs (otero_signatures/,
# phf1_markers/).
parser$add_argument("--out_dir", default = "seg_control",
  help = "Directory for the SEG geneset, relative to the project root [default: seg_control]")
parser$add_argument("--out_name", default = "Lin_SEG_2019.qs",
  help = "Output .qs filename [default: Lin_SEG_2019.qs]")
parser$add_argument("--set_name", default = "SEG",
  help = "Name of the single set inside the output list [default: SEG]")
args <- parser$parse_args()

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(args$out_dir))
  stop("could not create out_dir: ", args$out_dir, call. = FALSE)
# Probe writability now, so a read-only target fails legibly here rather than after
# the whole extraction has run.
.probe <- file.path(args$out_dir, ".write_probe")
if (!file.create(.probe, showWarnings = FALSE))
  stop("out_dir is not writable: ", normalizePath(args$out_dir, mustWork = FALSE),
       "\n  (pass --out_dir to somewhere writable.)", call. = FALSE)
unlink(.probe)
cat("Output dir:", normalizePath(args$out_dir), "\n")

out_qs  <- file.path(args$out_dir, args$out_name)
out_txt <- file.path(args$out_dir, sub("\\.qs$", "_provenance.txt", args$out_name))

##  ............................................................................
##  Require scMerge - do NOT substitute another list                        ####
if (!requireNamespace("scMerge", quietly = TRUE)) {
  stop(paste0(
    "\n",
    "scMerge is not installed in this R environment, so the published SEG list\n",
    "(Lin et al. 2019 GigaScience, scSEGIndex) cannot be extracted.\n",
    "\n",
    "This script deliberately does NOT substitute another housekeeping/stable gene\n",
    "list: the point of the control is that it is the published one.\n",
    "\n",
    "Install into the analysis environment and re-run:\n",
    "\n",
    "  R -e 'if (!requireNamespace(\"BiocManager\", quietly=TRUE)) ",
    "install.packages(\"BiocManager\", repos=\"https://cloud.r-project.org\"); ",
    "BiocManager::install(c(\"scMerge\", \"ruv\", \"M3Drop\"), ask=FALSE, update=FALSE)'\n",
    "\n",
    "  (`ruv` and `M3Drop` are scMerge dependencies.)\n",
    "\n",
    "Then:  Rscript R/prepare_seg_geneset.R\n"),
    call. = FALSE)
}

cat("scMerge version:", as.character(utils::packageVersion("scMerge")), "\n\n")

##  ............................................................................
##  Load and INSPECT segList before subsetting                              ####
# Element naming differs across scMerge versions, and each SEG list is shipped in
# both a gene-symbol and an Ensembl-gene-ID flavour. Print the structure, then
# resolve by name and validate - never index positionally.
segList <- NULL
utils::data("segList", package = "scMerge", envir = environment())
if (is.null(segList)) stop("data('segList', package='scMerge') did not bind `segList`")

cat("==== str(segList, max.level = 2) ====\n")
print(utils::str(segList, max.level = 2, list.len = 50))
cat("\nTop-level names(segList):\n")
print(names(segList))
for (nm in names(segList)) {
  el <- segList[[nm]]
  if (is.list(el)) {
    cat(sprintf("\n  segList[['%s']] is a list of %d; names:\n", nm, length(el)))
    print(names(el))
    cat("  lengths:\n"); print(vapply(el, length, integer(1)))
  } else {
    cat(sprintf("\n  segList[['%s']] is a %s of length %d\n", nm, class(el)[1], length(el)))
  }
}
cat("\n")

##  ............................................................................
##  Resolve the HUMAN gene-SYMBOL vector                                    ####
# Flatten to depth 2 so both shapes are handled: segList$human$human_scSEG (a nested
# list, the usual layout) and a flat segList$human_scSEG.
flat <- list()
for (nm in names(segList)) {
  el <- segList[[nm]]
  if (is.list(el)) {
    for (nm2 in names(el)) flat[[paste(nm, nm2, sep = "$")]] <- el[[nm2]]
  } else {
    flat[[nm]] <- el
  }
}
flat <- flat[vapply(flat, is.character, logical(1))]
stopifnot("segList contains no character vectors" = length(flat) > 0)

cat("Character-vector candidates in segList (path -> length):\n")
print(vapply(flat, length, integer(1)))

is_ensembl <- function(x) mean(grepl("^ENS[A-Z]*G[0-9]{6,}", x)) > 0.5
paths      <- names(flat)
# human, single-cell SEG (scSEGIndex), symbol flavour: require human + scSEG in the
# path, exclude the explicit Ensembl-ID variants, then exclude anything that still
# looks like Ensembl IDs by content.
cand <- paths[grepl("human", paths, ignore.case = TRUE) &
              grepl("scSEG", paths, ignore.case = TRUE) &
              !grepl("ensembl|geneid|entrez", paths, ignore.case = TRUE)]
cand <- cand[!vapply(flat[cand], is_ensembl, logical(1))]

if (length(cand) == 0) {
  # Fall back to any human SEG list that is symbol-like, still never positional.
  cand <- paths[grepl("human", paths, ignore.case = TRUE) &
                grepl("seg",   paths, ignore.case = TRUE) &
                !grepl("ensembl|geneid|entrez", paths, ignore.case = TRUE)]
  cand <- cand[!vapply(flat[cand], is_ensembl, logical(1))]
  if (length(cand) > 0)
    cat("\nNOTE: no 'human*scSEG' symbol element; fell back to human SEG candidates.\n")
}
if (length(cand) == 0)
  stop(paste0("Could not resolve a human gene-symbol SEG vector from segList. ",
              "Candidates were:\n  ", paste(paths, collapse = "\n  "),
              "\nInspect the printed structure above and set the element explicitly."),
       call. = FALSE)
if (length(cand) > 1) {
  cat("\nMultiple human symbol SEG candidates; taking the LONGEST:\n")
  print(vapply(flat[cand], length, integer(1)))
  cand <- cand[which.max(vapply(flat[cand], length, integer(1)))]
}

seg_path <- cand[1]
seg      <- sort(unique(trimws(flat[[seg_path]])))
seg      <- seg[nzchar(seg) & !is.na(seg)]

cat(sprintf("\nResolved element: segList[['%s']]\n", gsub("\\$", "']][['", seg_path)))

##  ............................................................................
##  Validate                                                               ####
stopifnot(
  "resolved SEG element is not a character vector" = is.character(seg),
  "resolved SEG element is implausibly short (<200 genes)" = length(seg) > 200,
  "resolved SEG element looks like Ensembl IDs, not gene symbols" =
    mean(grepl("^ENS[A-Z]*G[0-9]{6,}", seg)) < 0.01
)

cat("\nlength(SEG) =", length(seg), "\n")
cat("First 20 symbols:\n")
print(head(seg, 20))

##  ............................................................................
##  Write geneset .qs (same shape as Mancuso_2024.qs) + provenance          ####
seg_list <- list(seg)
names(seg_list) <- args$set_name        # named list of character vectors

qs::qsave(seg_list, out_qs)
cat("\nWrote:", out_qs, "\n")

# Round-trip check: prep_sets() will call intersect() on each element, so the object
# must read back as a named list of character vectors.
chk <- qs::qread(out_qs)
stopifnot("round-trip failed: not a named list" = is.list(chk) && !is.null(names(chk)),
          "round-trip failed: element is not character" = is.character(chk[[args$set_name]]),
          "round-trip failed: length changed" = length(chk[[args$set_name]]) == length(seg))
cat("Round-trip verified: list(", args$set_name, " = <character[", length(seg), "]>)\n", sep = "")

writeLines(c(
  "Lin_SEG_2019.qs - provenance",
  "",
  "Content : list(SEG = <human stably expressed gene symbols>)",
  sprintf("n_genes : %d", length(seg)),
  sprintf("Source  : scMerge %s, data('segList'), element segList[['%s']]",
          as.character(utils::packageVersion("scMerge")),
          gsub("\\$", "']][['", seg_path)),
  "Citation: Lin Y, Ghazanfar S, Strbenac D, Wang A, Patrick E, Lin DM, Speed T,",
  "          Yang JYH, Yang P. Evaluating stably expressed genes in single cells.",
  "          GigaScience 2019;8(9):giz106.  (scSEGIndex)",
  sprintf("Extracted: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("Script   : R/prepare_seg_geneset.R"),
  sprintf("R        : %s", R.version.string),
  "",
  "Used as the published negative-control module in",
  "R/modulescore_seg_control_vs_phf1_distance.R",
  "",
  "First 20 symbols:",
  paste(head(seg, 20), collapse = ", ")
), out_txt)
cat("Wrote:", out_txt, "\n")

cat("\nsessionInfo():\n"); print(sessionInfo())
cat("\nDone.\n")
