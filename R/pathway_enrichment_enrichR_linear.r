# ---------------------------------------------------------------------------
# pathway_enrichment_enrichR_linear.r
#
# Produces: Table S7
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
################################################################################
#' Functional enrichment analysis using enrichR — linear distance DEG version
#'
#' Reads per-celltype *_dist_to_phf1_um_scaled.tsv results from
#' deg_dream_linear_distance_phf1.r and runs enrichR ORA separately for:
#'   - UP:   genes with positive logFC (higher expression further from PHF1+)
#'   - DOWN: genes with negative logFC (lower expression further from PHF1+,
#'           i.e. upregulated near PHF1+ cells)
#'
#' Background is all genes tested in each celltype (not a fixed universe).
#'
#' Output structure mirrors the categorical enrichR analysis:
#'   <ct_dir>/enrichr/UP/
#'     <db>.tsv          — per-database significant pathways
#'     <db>.png          — per-database dotplot
#'     merged_enrichr.tsv  — all databases combined
#'     merged_enrichr.png  — merged dotplot
#'   <ct_dir>/enrichr/DOWN/
#'     (same)
#'
#' Adapted minimally from pathway_enrichment_enrichR_with_background.r
################################################################################

library(tidyverse)
library(enrichR)
library(stringr)
library(ggplot2)
library(cowplot)
library(purrr)
library(cli)

# ---- Configuration -----------------------------------------------------------

base_dir <- "<PROJECT_ROOT>/phf1_v2/deg/de_linear_distance"
# base_dir <- "<PROJECT_ROOT>/phf1_v2/deg/de_linear_distance_fov"
# base_dir <- "<PROJECT_ROOT>/phf1_v2/deg/de_linear_distance_raw"

PADJ_THRESHOLD <- 0.1

# logFC threshold: set to 0 to use all significant genes regardless of slope
# magnitude. Raise (e.g. 0.05) to require a minimum effect size.
LOGFC_THRESHOLD <- 0

ENRICHMENT_DATABASES <- c(
  "GO_Molecular_Function_2025",
  "GO_Cellular_Component_2025",
  "GO_Biological_Process_2025"#,
  # "WikiPathways_2024_Human",
  # "KEGG_2021_Human"
)

# ---- Helper functions (as in pathway_enrichment_enrichR_with_background.r) ---

format_res_table_enrichr <- function(res, min_size = 15, max_size = 250, min_overlap = 3) {
  res_table <- res %>% as.data.frame() %>%
    dplyr::transmute(
      geneset     = .get_geneset(Term),
      description = gsub("\\(GO:.*|Homo sapiens.R-HSA.*|WP.*", "", Term),
      database    = database,
      size        = as.numeric(gsub(".*\\/", "", Overlap)),
      overlap     = as.numeric(gsub("\\/.*", "", Overlap)),
      odds_ratio  = round(Odds.Ratio, 2),
      pval        = as.numeric(format(P.value, format = "e", digits = 2)),
      FDR         = as.numeric(format(Adjusted.P.value, format = "e", digits = 2))
    )

  res_table$genes <- res$Genes

  text_out <- "cardiac|lens|Mesenchymal|myotube|glioma|melanoma|thyroid|glomerular|renal|retina|nigra|vessel|estrogen|steroid|androgen|artery|bone|skeletal|muscle|aorta|cartilage|pancreatic|myoblast|embryonic|amyotrophic|neural tube|circadian|ectoderm|stem cell|vitamin|chylomicron|coronary|osteoclast|addiction|tumor|myometrial|prolactin|glioblastoma|sensory|cancer|carcinoma|hepatitis|oocyte|cardiomyocyte|heart|cardiac|eye|kidney|ear|auditory|Allograft|lupus|graft|rett|nose"

  # Non-CNS tissue/organ and developmental-morphology GO terms that cannot be
  # represented in a cortical dataset. Excluded (pre-specified, on biology, not
  # tuned to results) to restrict the tested universe to CNS-plausible processes.
  # Immune terms are deliberately NOT excluded (relevant for microglia).
  text_out <- paste(text_out,
    "endocardial|genitalia|mammary|pharyngeal|mesenchym|chondrocyte|ossification|biomineral|angiogenesis|vasoconstriction|bilateral symmetry|left/right|left / right|pattern formation|epithelial|tight junction",
    sep = "|")

  res_table <- res_table %>%
    dplyr::filter(!grepl(text_out, description, ignore.case = TRUE)) %>%
    dplyr::filter(size >= min_size, size < max_size, overlap >= min_overlap, odds_ratio > 2) %>%
    dplyr::mutate(FDR = p.adjust(pval, method = "BH"))

  res_table <- res_table[res_table$FDR <= 0.05, ]
  return(res_table)
}

.get_geneset <- function(term) {
  geneset <- purrr::map_chr(as.character(term),
                            ~ str_extract(., "GO:.*|R-HSA.*|WP.*"))
  geneset <- gsub("\\)|Homo sapiens", "", geneset)
  as.character(geneset)
}

dotplot_enrichr <- function(dt) {
  dt <- dt %>%
    dplyr::filter(!is.na(pval)) %>%
    dplyr::top_n(., min(10, nrow(.)), -pval) %>%
    dplyr::group_by(clusters) %>%
    dplyr::top_n(., min(2, nrow(.)), -pval)
  dt$description <- stringr::str_wrap(dt$description, 40)

  ggplot2::ggplot(dt, aes(
    x = odds_ratio,
    y = stats::reorder(description, odds_ratio)
  )) +
    geom_point(aes(fill = FDR, size = overlap),
               shape = 21, alpha = 0.7, color = "black") +
    scale_size(name = "Overlap", range = c(3, 8)) +
    xlab("Total Odds Ratio") +
    ylab("") +
    scale_fill_gradient(
      low = "navy", high = "gold", name = "FDR",
      guide = guide_colorbar(reverse = TRUE),
      limits = c(0, 0.1),
      aesthetics = c("fill")
    ) +
    guides(size = guide_legend(override.aes = list(fill = "gold", color = "gold"))) +
    cowplot::theme_cowplot() +
    cowplot::background_grid()
}

cluster_pathway <- function(enrichment_res,
                            cut_height     = 0.75,
                            cluster_method = "complete",
                            plot_dendo     = FALSE) {
  dt <- enrichment_res %>% dplyr::select(c(description, genes))

  if (nrow(dt) == 0) {
    cli::cli_alert("Number of enriched terms is zero")
    enrichment_res <- NULL
  } else if (nrow(dt) < 5) {
    cli::cli_alert("Number of enriched terms is too small for clustering")
    enrichment_res$clusters <- 1
  } else {
    dt_list <- split(as.character(dt$genes), as.character(dt$description))
    dt_list <- lapply(dt_list, function(x) unlist(strsplit(x, ";")[[1]]))

    mat <- t(splitstackshape:::charMat(listOfValues = dt_list, fill = 0L))
    colnames(mat) <- names(dt_list)
    mat <- as.data.frame(mat)

    idx <- which(colSums(mat) == nrow(mat))
    if (length(idx) > 0) mat <- mat[, -idx]

    if (nrow(mat) < 5) {
      enrichment_res$clusters <- 1
    } else {
      kappa_mat <- colpair_map(mat, cohen.kappa.pair)
      kappa_mat <- as.data.frame(kappa_mat)
      rownames(kappa_mat) <- kappa_mat$term
      kappa_mat$term <- NULL
      kappa_mat <- as.matrix(kappa_mat)
      kappa_tree <- hclust(as.dist(1 - kappa_mat), method = cluster_method)

      geneset_cluster <- stats::cutree(kappa_tree, h = quantile(kappa_tree$height, cut_height))
      geneset_cluster <- data.frame(description = names(geneset_cluster), clusters = geneset_cluster)

      cli::cli_alert("Total {max(geneset_cluster$clusters)} geneset clusters found.")

      if (plot_dendo) {
        plot(kappa_tree, labels = FALSE)
        abline(h = stats::quantile(kappa_tree$height, cut_height))
      }

      enrichment_res <- dplyr::left_join(enrichment_res, geneset_cluster, by = "description")
      enrichment_res$clusters[is.na(enrichment_res$clusters)] <- max(enrichment_res$clusters, na.rm = TRUE) + 1
      attr(enrichment_res, "kappa_tree") <- kappa_tree
    }
  }
  return(enrichment_res)
}

cohen.kappa.pair <- function(x, y) {
  tab <- table(x, y)
  k   <- rhoR::kappa_ct(tab)
  round(k, 3)
}

colpair_map <- function(.data, .f, ..., .diagonal = NA) {
  out <- purrr::map_dfr(.data, summarise_col, .f, .data, ...)
  as_cordf(out, diagonal = .diagonal)
}

summarise_col <- function(x, f, data) {
  dplyr::summarise(data, dplyr::across(.cols = dplyr::everything(), .fns = f, x))
}

as_cordf <- function(x, diagonal = NA) {
  if (inherits(x, "cor_df")) { warning("x is already a correlation data frame."); return(x) }
  x <- as.data.frame(x)
  if (ncol(x) != nrow(x)) stop("Input object x is not square.")
  if (ncol(x) > 1) diag(x) <- diagonal
  new_cordf(x, names(x))
}

new_cordf <- function(x, term = NULL) {
  if (!is.null(term)) x <- first_col(x, term)
  class(x) <- c("cor_df", class(x))
  x
}

first_col <- function(df, ..., var = "term") {
  stopifnot(is.data.frame(df))
  if (tibble::has_name(df, var)) stop("Column named ", var, " already exists!")
  new_col <- tibble::tibble(...)
  names(new_col) <- var
  dplyr::as_tibble(c(new_col, df))
}

# ---- Main: loop over celltypes -----------------------------------------------

celltypes <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
celltypes <- celltypes[celltypes != ""]

if (length(celltypes) == 0) stop("No celltype subdirectories found in: ", base_dir)
cat("Found", length(celltypes), "celltypes\n")

for (ct_dir in celltypes) {

  ct <- basename(ct_dir)
  cat("\n=== Celltype:", ct, "===\n")

  deg_file <- list.files(ct_dir,
                         pattern   = "_dist_to_phf1_um_scaled\\.tsv$",
                         full.names = TRUE)

  if (length(deg_file) == 0) {
    cat("  No *_dist_to_phf1_um_scaled.tsv found — skipping\n")
    next
  }

  res <- read_tsv(deg_file[1], show_col_types = FALSE)

  # Background: all genes tested in this celltype
  background_genes <- res$gene

  # Significant genes split by direction
  sig <- res %>% dplyr::filter(!is.na(padj), padj < PADJ_THRESHOLD, abs(logFC) > LOGFC_THRESHOLD)

  gene_lists <- list(
    UP   = sig %>% dplyr::filter(logFC >  LOGFC_THRESHOLD) %>% pull(gene),
    DOWN = sig %>% dplyr::filter(logFC < -LOGFC_THRESHOLD) %>% pull(gene)
  )

  cat(sprintf("  Significant: %d UP (pos logFC), %d DOWN (neg logFC)\n",
              length(gene_lists$UP), length(gene_lists$DOWN)))

  for (direction in c("UP", "DOWN")) {

    gene_list <- gene_lists[[direction]]
    out_dir   <- file.path(ct_dir, "enrichr", direction)

    if (length(gene_list) < 3) {
      cat(sprintf("  %s: fewer than 5 genes — skipping enrichment\n", direction))
      next
    }

    cat(sprintf("  Running enrichR for %s (%d genes)...\n", direction, length(gene_list)))

    raw <- enrichR::enrichr(
      genes           = gene_list,
      databases       = ENRICHMENT_DATABASES,
      background      = background_genes,
      include_overlap = TRUE
    )

    raw <- purrr::discard(raw, ~ nrow(.) == 0)

    if (length(raw) == 0) {
      cat(sprintf("  %s: no enrichment results returned\n", direction))
      next
    }

    enrichr_res <- lapply(names(raw), function(db) {
      raw[[db]] %>% mutate(database = db)
    })
    names(enrichr_res) <- names(raw)

    enrichr_res <- purrr::map(enrichr_res, format_res_table_enrichr)
    enrichr_res <- purrr::discard(enrichr_res, ~ nrow(.) == 0)

    if (length(enrichr_res) == 0) {
      cat(sprintf("  %s: no pathways significant at FDR <= 0.05\n", direction))
      next
    }

    enrichr_res <- lapply(enrichr_res, cluster_pathway,
                          cut_height = 0.75, cluster_method = "complete", plot_dendo = FALSE)
    enrichr_res <- purrr::discard(enrichr_res, is.null)

    if (length(enrichr_res) == 0) next

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Save per-database TSV and PNG
    for (db in names(enrichr_res)) {
      dt      <- enrichr_res[[db]]
      db_safe <- gsub("[^A-Za-z0-9_]", "_", db)

      out_tsv <- file.path(out_dir, paste0(db_safe, ".tsv"))
      write.table(dt, out_tsv,
                  col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
      cat(sprintf("    Written: %s\n", basename(out_tsv)))

      if (!is.null(dt) && nrow(dt) > 0) {
        p <- dotplot_enrichr(dt) +
          labs(title = sprintf("%s | %s | %s", ct, direction, db))
        out_png <- file.path(out_dir, paste0(db_safe, ".png"))
        ggsave(out_png, p, width = 10, height = 6, dpi = 300, device = "png")
        cat(sprintf("    Written: %s\n", basename(out_png)))
      }
    }

    # Merged TSV across all databases
    merged <- bind_rows(enrichr_res)
    merged_tsv <- file.path(out_dir, "merged_enrichr.tsv")
    write.table(merged, merged_tsv,
                col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    cat(sprintf("    Written: merged_enrichr.tsv (%d rows)\n", nrow(merged)))

    # Merged dotplot
    if (nrow(merged) > 0) {
      p_merged <- dotplot_enrichr(merged) +
        labs(title = sprintf("%s | %s | all databases", ct, direction))
      merged_png <- file.path(out_dir, "merged_enrichr.png")
      ggsave(merged_png, p_merged, width = 10, height = 6, dpi = 300, device = "png")
      cat(sprintf("    Written: merged_enrichr.png\n"))
    }
  }
}

cat("\nDone.\n")
