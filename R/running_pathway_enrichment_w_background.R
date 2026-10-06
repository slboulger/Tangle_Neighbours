# ---------------------------------------------------------------------------
# running_pathway_enrichment_w_background.R
#
# Produces: Table S5
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# Load enrichment script (background-aware function)
source("<PROJECT_ROOT>/phf1_v2/R/pathway_enrichment_enrichR_with_background.r")

library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)

# Root directory with DE results
base_dir <- "<PROJECT_ROOT>/phf1_v2/deg/de_cat/de_PHF1"
# base_dir <- "<PROJECT_ROOT>/phf1_v2/deg/de_cat_distal300/de_PHF1"
# base_dir <- "<PROJECT_ROOT>/phf1_v2/deg/de_cat/de_phf1_p95_pos_5mad"
# base_dir <- "<PROJECT_ROOT>/phf1_v2/deg/de_cat/de_phf1_p95_pos_3mad"

# ---- STEP 1: Enumerate celltypes ----
# Enrichment background is PER-CELLTYPE (each celltype's own tested-gene universe),
# computed inside the loop below. This matches the linear-distance pathway pipeline
# (pathway_enrichment_enrichR_linear.r, background = that celltype's tested genes).
# A cross-celltype union would add genes not expressed in a given celltype and inflate
# the background; per-celltype is the expressed universe for each celltype.
celltypes <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

# ---- STEP 2: Loop through celltypes ----
min_genes <- 3  # minimum number of genes to run enrichment

for (ct in celltypes) {
  deg_files <- list.files(ct, pattern = "\\.tsv$", full.names = TRUE, recursive = FALSE)
  
  for (deg in deg_files) {
    message("\nRunning pathway enrichment for: ", deg)
    
    deg_file <- read.delim(deg, header = TRUE)

    # Per-celltype tested-gene background (this celltype's expressed universe).
    # Matches the linear pipeline.
    background <- unique(deg_file$gene)
    cat("Per-celltype background size:", length(background), "\n")

    up_genes <- deg_file %>% dplyr::filter(padj < 0.1, logFC > 0.25) %>% pull(gene)
    down_genes <- deg_file %>% dplyr::filter(padj < 0.1, logFC < -0.25) %>% pull(gene)
    
    # Skip too-small gene sets
    if (length(up_genes) < min_genes) {
      message(" → Skipping UP genes (too few genes: ", length(up_genes), ")")
      up_genes <- NULL
    }
    if (length(down_genes) < min_genes) {
      message(" → Skipping DOWN genes (too few genes: ", length(down_genes), ")")
      down_genes <- NULL
    }
    if (is.null(up_genes) && is.null(down_genes)) next
    
    # Output directories
    outdir_up <- file.path(ct, "enrichr", "UP")
    outdir_down <- file.path(ct, "enrichr", "DOWN")
    dir.create(outdir_up, showWarnings = FALSE, recursive = TRUE)
    dir.create(outdir_down, showWarnings = FALSE, recursive = TRUE)
    
    # ---- Helper to run enrichment and save results ----
    run_and_save <- function(gene_list, outdir, direction) {
      if (is.null(gene_list) || length(gene_list) == 0) return(NULL)
      
      message(" → Enrichment for ", direction, " genes: ", length(gene_list))
      
      res <- tryCatch({
        pathway_analysis_enrichr(
          interest_gene = gene_list,
          enrichment_database = c(
            "GO_Molecular_Function_2025",
            "GO_Cellular_Component_2025",
            "GO_Biological_Process_2025"#,
            # "Reactome_Pathways_2024"
          ),
          min_size = 15,
          max_size = 250,
          min_overlap = 3,
          background = background
        )
      }, error = function(e) {
        message("   → Enrichment failed: ", e$message)
        return(NULL)
      })
      
      if (is.null(res)) return(NULL)
      
      base_prefix <- tools::file_path_sans_ext(basename(deg))
      
      # Save each database separately
      for (db in setdiff(names(res), "plot")) {
        db_res <- res[[db]]
        message("   → Columns in ", db, ": ", paste(colnames(db_res), collapse = ", "))
        
        if (is.null(db_res) || nrow(db_res) == 0) {
          message("   → Skipping database ", db, " (no results)")
          next
        }
        
        # Save TSV
        tsv_file <- file.path(outdir, paste0(base_prefix, "_", direction, "_", db, ".tsv"))
        write.table(db_res, file = tsv_file, sep = "\t", quote = FALSE, row.names = FALSE)
        
        # Save PNG
        png_file <- file.path(outdir, paste0(base_prefix, "_", direction, "_", db, ".png"))
        png(png_file, width = 10, height = 8, res = 300, units = "in")
        print(dotplot_enrichr(db_res))
        dev.off()
      }
      
      # Save merged results if any valid database exists
      valid_dbs <- setdiff(names(res), "plot")[sapply(res[setdiff(names(res), "plot")], function(x) !is.null(x) && nrow(x) > 0)]
      if (length(valid_dbs) == 0) return(NULL)
      
      merged_res <- bind_rows(res[valid_dbs], .id = "Database")
      merged_file <- file.path(outdir, paste0(base_prefix, "_", direction, "_merged.tsv"))
      write.table(merged_res, file = merged_file, sep = "\t", quote = FALSE, row.names = FALSE)
      
      merged_png <- file.path(outdir, paste0(base_prefix, "_", direction, "_merged.png"))
      png(merged_png, width = 10, height = 8, res = 300, units = "in")
      print(dotplot_enrichr(merged_res))
      dev.off()
    }
    
    run_and_save(up_genes, outdir_up, "UP")
    run_and_save(down_genes, outdir_down, "DOWN")
  }
}
