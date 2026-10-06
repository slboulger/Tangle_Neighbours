# ---------------------------------------------------------------------------
# pathway_enrichment_enrichR_with_background.r
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
################################################################################
#' Functional enrichment analysis using enrichR
#'
#' Performs impacted pathway analysis with a list of genes.
#'
#' @param gene_file A data frame or the path of a .tsv file containing
#' a list of genes, their fold-change, p-value and adjusted p-value.
#' Column names should be gene, logFC, pval and padj respectively.
#' @param enrichment_database Name of the database for enrichment. User can
#' specify one or more database names from [enrichR::listEnrichrDbs()].
#' @param is_output If TRUE a folder will be created and results of enrichment
#' analysis will be saved otherwise a R list will be returned. Default FALSE
#' @param output_dir Path for the output directory. Default is current dir.
#'
#' @return enrichment_result a list of data.frames containing enrichment output
#' and a list of plots of top 10 significant genesets.
#'
#' @family Impacted pathway analysis
#'
#' @importFrom enrichR enrichr listEnrichrDbs
#' @importFrom cli cli_alert_info cli_text
#' @importFrom ggplot2 ggplot ggsave
#' @importFrom cowplot theme_cowplot background_grid
#' @importFrom stringr str_wrap
#' @importFrom assertthat assert_that
#' @importFrom dplyr %>% mutate
#' @importFrom purrr map map_chr discard
#'
#' @export
pathway_analysis_enrichr <- function(interest_gene = NULL,
                                     enrichment_database = c(
                                       "GO_Molecular_Function_2025",
                                       "GO_Cellular_Component_2025",
                                       "GO_Biological_Process_2025",
                                       "WikiPathways_2024_Human",
                                       "KEGG_2021_Human"
                                     ),
                                     min_size = 10,
                                     max_size = 250,
                                     min_overlap = 3,
                                     background = NULL) {
  library(enrichR)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  
  res <- enrichR::enrichr(
    genes = as.character(interest_gene),
    databases = enrichment_database,
    background = background,
    include_overlap = TRUE
  )
  
  res <- purrr::discard(res, function(x) {
    nrow(x) == 0
  })
  
  enrichr_res <- lapply(names(res), function(x){dt <- res[[x]] %>% 
    mutate(database = x)})
  
  names(enrichr_res) <- names(res)
  
  enrichr_res <- purrr::map(
    enrichr_res,
    ~ format_res_table_enrichr(., min_size, max_size, min_overlap)
  )
  
  enrichr_res <- purrr::discard(enrichr_res, function(x) {
    nrow(x) == 0
  })
  
  if (length(enrichr_res) == 0) {
    cli::cli_text(
      "{.strong No significant impacted pathways found at FDR <= 0.05! }"
    )
    enrichr_res <- NULL
  } else {
    
    enrichr_res <- lapply(enrichr_res, cluster_pathway,
                  cut_height = 0.75,
                  cluster_method = "complete",
                  plot_dendo = FALSE)
    
    enrichr_res$plot <- lapply(
      enrichr_res,
      function(dt) dotplot_enrichr(dt)
    )
  }
  
  return(enrichr_res)
}


#' Format result table
#' @keywords internal

format_res_table_enrichr <- function(res, min_size, max_size, min_overlap) {
  res_table <- res %>% as.data.frame() %>%
    dplyr::transmute(
      #geneset = res,
      geneset = .get_geneset(Term),
      description = gsub("\\(GO:.*|Homo sapiens.R-HSA.*|WP.*" , "", Term),
      database = database,
      size = as.numeric(gsub(".*\\/", "", Overlap)),
      overlap = as.numeric(gsub("\\/.*", "", Overlap)),
      odds_ratio = round(Odds.Ratio, 2),
      pval = as.numeric(format(P.value, format = "e", digits = 4)),
      FDR = as.numeric(format(Adjusted.P.value, format = "e", digits = 4))
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
    dplyr::filter(!grepl(text_out, description, ignore.case = T)) %>%
    dplyr::filter(size >= min_size,  size < max_size, overlap >= min_overlap, odds_ratio > 2) %>%
    dplyr::mutate(FDR = p.adjust(pval, method = "BH"))
  
  res_table <- res_table[res_table$FDR <= 0.05, ]
  return(res_table)
}

.get_geneset <- function(term) {
  # geneset <- purrr::map_chr(
  #   as.character(term),
  #   ~ strsplit(., "(", fixed = TRUE)[[1]][2]
  # )
  
  geneset <- purrr::map_chr(as.character(term),
                            ~ str_extract(. , "GO:.*|R-HSA.*|WP.*"))
  geneset <- gsub("\\)|Homo sapiens", "", geneset)
  geneset <- as.character(geneset)
  return(geneset)
}

#' dotplot for ORA. x axis perturbation, y axis description
#' @importFrom stats reorder
#' @keywords internal


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
               shape = 21, alpha = 0.7, color = "black"
    ) +
    scale_size(name = "Overlap", range = c(3, 8)) +
    xlab("Total Odds Ratio") +
    ylab("") +
    scale_fill_gradient(
      low = "navy", high = "gold", name = "FDR",
      guide = guide_colorbar(reverse = TRUE),
      limits = c(0, 0.1),
      aesthetics = c("fill")
    ) +
    guides(size = guide_legend(
      override.aes = list(fill = "gold", color = "gold")
    )) +
    cowplot::theme_cowplot() +
    cowplot::background_grid()
}


####
cluster_pathway <- function(enrichment_res,
                            cut_height = 0.75,
                            cluster_method = "complete",
                            plot_dendo = TRUE){
  
  dt <- enrichment_res %>%
    dplyr::select(c(description, genes))
  
  if(nrow(dt) == 0){
    cli::cli_alert("Number of enriched term is zero")
    enrichment_res <- NULL
  } else if (nrow(dt) < 5) {
    cli::cli_alert("Number of enriched term is too small")
    enrichment_res$clusters <- 1
  } else {
    
    ##creating a list of all genesets with gene names
    dt_list <- split(as.character(dt$genes), as.character(dt$description))
    dt_list <- lapply(dt_list, function(x){ x <- unlist(strsplit(x, ";")[[1]])})
    
    ##preparing matrix for cohen's kappa calculation
    mat <- t(splitstackshape:::charMat(listOfValues = dt_list, fill = 0L))
    colnames(mat) <- names(dt_list)
    mat <- as.data.frame(mat)
    
    idx <- which(colSums(mat) == nrow(mat))
    
    if (length(idx) > 0) {
      mat <- mat[ , -idx]
    }
    
    if (nrow(mat) < 5) {
      enrichment_res$clusters <- 1
    } else {
      
      kappa_mat <- colpair_map(mat, cohen.kappa.pair)
      kappa_mat <- as.data.frame(kappa_mat)
      rownames(kappa_mat) <- kappa_mat$term
      kappa_mat$term <- NULL
      kappa_mat <- as.matrix(kappa_mat)
      kappa_tree <- hclust(as.dist(1-kappa_mat), method=cluster_method)
      
      geneset_cluster <- stats::cutree(kappa_tree, h = quantile(kappa_tree$height, cut_height))
      geneset_cluster <- data.frame(description = names(geneset_cluster), clusters = geneset_cluster)
      
      cli::cli_alert("Total {max(geneset_cluster$cluster)} geneset clusters found.")
      
      if (plot_dendo){
        plot(kappa_tree, labels=F)
        abline(h = stats::quantile(kappa_tree$height, cut_height))
      }
      
      enrichment_res <- dplyr::left_join(enrichment_res, geneset_cluster)
      
      enrichment_res$clusters[is.na(enrichment_res$clusters)] <- max(enrichment_res$clusters,na.rm = T) + 1
      
      attr(enrichment_res, "kappa_tree") <- kappa_tree
    }
  }
  
  return(enrichment_res)
  
}


cohen.kappa.pair <- function(x, y){
  tab <- table(x, y)
  k <- rhoR::kappa_ct(tab)
  k <- round(k, 3)
  return(k)
}

colpair_map <- function(.data, .f, ..., .diagonal = NA){
  
  out <- purrr::map_dfr(.data, summarise_col, .f, .data, ...)
  
  as_cordf(out, diagonal = .diagonal)
  
}

summarise_col <- function(x, f, data){
  
  dplyr::summarise(data, dplyr::across(.cols = dplyr::everything(),
                                       .fns = f,
                                       x))
}

as_cordf <- function(x, diagonal = NA) {
  if(inherits(x, "cor_df")){
    warning("x is already a correlation data frame.")
    return(x)
  }
  x <- as.data.frame(x)
  if(ncol(x) != nrow(x)) {
    stop("Input object x is not a square. ",
         "The number of columns must be equal to the number of rows.")
  }
  if (ncol(x) > 1) diag(x) <- diagonal
  new_cordf(x, names(x))
}

new_cordf <- function(x, term = NULL) {
  if (!is.null(term)) {
    x <- first_col(x, term)
  }
  class(x) <- c("cor_df", class(x))
  x
}

first_col <- function(df, ..., var = "term") {
  stopifnot(is.data.frame(df))
  
  if (tibble::has_name(df, var))
    stop("There is a column named ", var, " already!")
  
  new_col <- tibble::tibble(...)
  names(new_col) <- var
  new_df <- c(new_col, df)
  dplyr::as_tibble(new_df)
}