#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# build_group_sce_neighbours.r
#
# Upstream pipeline - builds the objects every panel reads
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# build_group_sce_neighbours.r
#
# Build the two POOLED neuron objects:
#
#   Excitatory  <- cbind of all Exc-* celltype_sce_neighbours objects
#   Inhibitory  <- cbind of all Inh-* celltype_sce_neighbours objects
#
# The six glia / vascular groups (Astro, Oligo, OPC, Micro, Endo, VLMC) are
# single celltypes and are read straight from celltype_sce_neighbours/ —
# nothing is duplicated for them. Unassigned and Unassigned Neuron are not pooled.
#
# Group membership is derived from `neuron_order` in R/palettes.R (^Exc- / ^Inh-)
# rather than hardcoded, so a future subtype rename cannot silently drop a
# population.
#
# To keep the pooled excitatory object manageable, `logcounts` and all
# reducedDims are dropped — downstream DEG uses counts() only.
# The fine `celltype` column is retained for subtype-composition diagnostics.
#
# Run ONCE, interactively (OpenOnDemand, ~64 GB). Not a PBS job.
#
#   Rscript R/build_group_sce_neighbours.r

##  ............................................................................
##  Load packages                                                           ####
suppressPackageStartupMessages({
  library(argparse)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

PROJECT_DIR <- "<PROJECT_ROOT>/phf1_v2"

##  ............................................................................
##  Parse command-line arguments                                            ####
parser <- ArgumentParser()

parser$add_argument(
  "--input_dir",
  help    = "Directory of per-celltype *_sce_neighbours.qs files",
  default = file.path(PROJECT_DIR, "celltype_sce_neighbours")
)
parser$add_argument(
  "--output_dir",
  help    = "Directory to write the pooled group objects into",
  default = file.path(PROJECT_DIR, "celltype_group_sce_neighbours")
)
parser$add_argument(
  "--palettes",
  help    = "Path to palettes.R (source of neuron_order)",
  default = file.path(PROJECT_DIR, "R", "palettes.R")
)

args <- parser$parse_args()

source(args$palettes)  # exports neuron_order, glia_order, ...

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

##  ............................................................................
##  Group definitions                                                       ####
GROUPS <- list(
  Excitatory = grep("^Exc-", neuron_order, value = TRUE),
  Inhibitory = grep("^Inh-", neuron_order, value = TRUE)
)

if (any(lengths(GROUPS) == 0)) {
  stop("No Exc-* or Inh-* entries found in neuron_order — has palettes.R been renamed?")
}

cat("Group membership (from neuron_order in palettes.R):\n")
for (g in names(GROUPS)) {
  cat(sprintf("  %-12s %d subtypes: %s\n", g, length(GROUPS[[g]]),
              paste(GROUPS[[g]], collapse = ", ")))
}
cat("\n")

##  ............................................................................
##  Build each group                                                        ####
# Strip everything dream does not use. logcounts and the reducedDims (PCA /
# HARMONY / UMAP) are the bulk of the object and are never read downstream.
.slim <- function(sce) {
  assays(sce)      <- assays(sce)["counts"]
  reducedDims(sce) <- list()
  altExps(sce)     <- list()
  sce
}

for (group in names(GROUPS)) {

  members <- GROUPS[[group]]
  paths   <- file.path(args$input_dir, paste0(members, "_sce_neighbours.qs"))

  missing <- members[!file.exists(paths)]
  if (length(missing) > 0) {
    stop("Missing input objects for group ", group, ": ",
         paste(missing, collapse = ", "))
  }

  cat("=== Building", group, "===\n")

  parts <- vector("list", length(members))
  for (i in seq_along(members)) {
    cat(sprintf("  reading %-30s ", members[i]))
    s <- .slim(qs::qread(paths[i]))
    rownames(s) <- rowData(s)$gene
    parts[[i]] <- s
    cat(sprintf("%d genes x %d cells\n", nrow(s), ncol(s)))
  }
  names(parts) <- members

  # cbind() requires identical rows. All members are column subsets of the same
  # sce.qs so this should hold; assert rather than assume.
  ref_rows <- rownames(parts[[1]])
  for (i in seq_along(parts)) {
    if (!identical(rownames(parts[[i]]), ref_rows)) {
      stop("rownames mismatch between ", members[1], " and ", members[i],
           " — the per-celltype objects are not row-aligned.")
    }
  }

  # colData columns can differ if a later annotation pass touched only some
  # files. Reduce to the intersection, preserving the first object's order.
  common_cols <- Reduce(intersect, lapply(parts, function(s) colnames(colData(s))))
  dropped     <- setdiff(colnames(colData(parts[[1]])), common_cols)
  if (length(dropped) > 0) {
    cat("  colData columns not shared by all members, dropped: ",
        paste(dropped, collapse = ", "), "\n", sep = "")
  }
  parts <- lapply(parts, function(s) { colData(s) <- colData(s)[, common_cols, drop = FALSE]; s })

  merged <- do.call(cbind, unname(parts))

  if (anyDuplicated(colnames(merged)) > 0) {
    stop("Duplicate cell IDs after cbind for group ", group,
         " — cell barcodes are not unique across subtypes.")
  }

  expected_cells <- sum(vapply(parts, ncol, integer(1)))
  stopifnot("cell count changed during cbind" = ncol(merged) == expected_cells)

  merged$celltype_group <- group

  out_file <- file.path(args$output_dir, paste0(group, "_sce_neighbours.qs"))
  qs::qsave(merged, out_file)

  cat(sprintf("  -> %d genes x %d cells written to %s\n", nrow(merged), ncol(merged), out_file))
  cat("  subtype composition:\n")
  print(table(droplevels(as.factor(merged$celltype))))
  cat("  cells per donor:\n")
  print(table(droplevels(as.factor(merged$sample_id))))
  cat("  assays:", paste(assayNames(merged), collapse = ", "),
      "| reducedDims:", ifelse(length(reducedDimNames(merged)) == 0, "(none)",
                               paste(reducedDimNames(merged), collapse = ", ")), "\n\n")

  rm(parts, merged); gc(verbose = FALSE)
}

cat("--- Session info ---\n")
print(sessionInfo())
