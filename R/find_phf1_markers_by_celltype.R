#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# find_phf1_markers_by_celltype.R
#
# Produces: Table S9
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# find_phf1_markers_by_celltype.R
#
# For each celltype, finds markers of PHF1+ vs PHF1- cells using
# Seurat FindMarkers (Wilcoxon, SCT assay).
#
# This is the gene-selection step that defines the PHF1 module: within each cell type,
# tangle-bearing vs tangle-free cells, keeping genes higher in tangle-bearing cells. Inference on
# the module is the Set 1 model of module score over distance (docs/MODELS.md).
#
# Outputs per celltype:
#   phf1_markers/<celltype>/<celltype>_phf1_markers.tsv  — full marker table
#   phf1_markers/<celltype>/<celltype>_phf1_geneset.txt  — upregulated genes only
#                                                          (p_val_adj < 0.05,
#                                                           avg_log2FC > 0.25)
#                                                          for use with AddModuleScore
#
# Global output:
#   phf1_markers/all_celltypes_phf1_markers.tsv — all celltypes merged

library(Seurat)
library(tidyverse)

# ---- Paths -------------------------------------------------------------------

seu_path <- "<PROJECT_ROOT>/phf1_v2/seu_PHF1.rds"

out_base <- "<PROJECT_ROOT>/phf1_v2/phf1_markers"

# ---- Thresholds --------------------------------------------------------------

PADJ_THRESHOLD <- 0.05
LOGFC_THRESHOLD <- 0.25   # Seurat default; also minimum for geneset inclusion
MIN_PCT         <- 0.2    # Seurat default
MIN_CELLS       <- 3      # minimum cells per group to attempt FindMarkers

# ---- Load Seurat object ------------------------------------------------------

cat("Loading Seurat object...\n")
seu <- readRDS(seu_path)

DefaultAssay(seu) <- "SCT"

cat("Cells:", ncol(seu), "\n")
cat("Genes:", nrow(seu), "\n")

stopifnot("celltype column missing from metadata" = "celltype" %in% colnames(seu@meta.data))
stopifnot("PHF1 column missing from metadata"     = "PHF1"     %in% colnames(seu@meta.data))

celltypes <- sort(unique(seu$celltype))
cat("Celltypes found:", length(celltypes), "\n")
print(celltypes)

# ---- Loop over celltypes -----------------------------------------------------

all_results <- list()

for (ct in celltypes) {

  cat("\n=== Celltype:", ct, "===\n")

  ct_safe <- gsub("[^A-Za-z0-9_-]", "_", ct)
  out_dir  <- file.path(out_base, ct_safe)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  seu_ct <- subset(seu, celltype == ct)

  n_true  <- sum(seu_ct$PHF1 == "TRUE",  na.rm = TRUE)
  n_false <- sum(seu_ct$PHF1 == "FALSE", na.rm = TRUE)
  cat(sprintf("  PHF1+: %d  PHF1-: %d\n", n_true, n_false))

  if (n_true < MIN_CELLS || n_false < MIN_CELLS) {
    cat(sprintf("  Skipping — fewer than %d cells in one group\n", MIN_CELLS))
    next
  }

  Idents(seu_ct) <- "PHF1"

  markers <- tryCatch(
    FindMarkers(
      seu_ct,
      ident.1         = "TRUE",
      ident.2         = "FALSE",
      assay           = "SCT",
      test.use        = "wilcox",
      min.pct         = MIN_PCT,
      logfc.threshold = LOGFC_THRESHOLD,
      verbose         = FALSE
    ),
    error = function(e) {
      cat("  ERROR in FindMarkers:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(markers) || nrow(markers) == 0) {
    cat("  No markers returned\n")
    next
  }

  markers <- markers %>%
    mutate(
      gene     = rownames(.),
      celltype = ct
    ) %>%
    dplyr::rename(padj = p_val_adj) %>%
    dplyr::select(gene, celltype, avg_log2FC, pct.1, pct.2, p_val, padj)

  cat(sprintf("  Markers returned: %d\n", nrow(markers)))
  cat(sprintf("  Significant (padj < %g, logFC > %g): %d\n",
              PADJ_THRESHOLD, LOGFC_THRESHOLD,
              sum(markers$padj < PADJ_THRESHOLD & markers$avg_log2FC > LOGFC_THRESHOLD,
                  na.rm = TRUE)))

  # Full marker table
  markers_file <- file.path(out_dir, paste0(ct_safe, "_phf1_markers.tsv"))
  write.table(markers, markers_file,
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  cat("  Written:", markers_file, "\n")

  # Upregulated geneset (for AddModuleScore)
  geneset <- markers %>%
    dplyr::filter(padj < PADJ_THRESHOLD, avg_log2FC > LOGFC_THRESHOLD) %>%
    arrange(desc(avg_log2FC)) %>%
    pull(gene)

  cat(sprintf("  Geneset size (upregulated only): %d genes\n", length(geneset)))

  if (length(geneset) > 0) {
    geneset_file <- file.path(out_dir, paste0(ct_safe, "_phf1_geneset.txt"))
    writeLines(geneset, geneset_file)
    cat("  Written:", geneset_file, "\n")
  } else {
    cat("  No significant upregulated genes — geneset file not written\n")
  }

  all_results[[ct]] <- markers
}

# ---- Merged output -----------------------------------------------------------

if (length(all_results) > 0) {
  merged <- bind_rows(all_results)
  merged_file <- file.path(out_base, "all_celltypes_phf1_markers.tsv")
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
  write.table(merged, merged_file,
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  cat("\nWritten merged results:", merged_file, "\n")
  cat("Total markers across all celltypes:", nrow(merged), "\n")
} else {
  cat("\nNo results to merge.\n")
}

cat("\nDone.\n")
