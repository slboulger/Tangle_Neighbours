#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# reactome_module_genesets.R
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
# reactome_module_genesets.R
#
# Resolve Reactome stable IDs (R-HSA-...) to gene vectors, for use as module-score
# gene sets. Sourceable; pure definitions, no side effects.
#
# WHY THIS EXISTS. The Enrichr Reactome GMT we already cache is keyed by pathway NAME
# only -- it carries no R-HSA IDs at all. A curated, hypothesis-driven pathway list is
# naturally written as stable IDs (names drift between Reactome releases; IDs do not),
# so something has to bridge the two. That bridge already exists in pieces:
#
#   R-HSA id  --[ genesets/reactome/ReactomePathways.txt ]-->  pathway name
#   name      --[ norm_name() in R/reactome_pcd_family.R  ]-->  Enrichr GMT term
#   term      --[ genesets/enrichr_gmt/*.gmt              ]-->  gene vector
#
# get_reactome_pcd_family() in R/reactome_pcd_family.R does exactly this join, but it
# also walks the pathway DAG and stop()s below 15 matched sets, so it is not reusable
# for a small hand-picked ID list. This file is the same two-hop join with no subtree
# walk and no size floor.
#
# Both cached files are provenance-stamped (.meta.txt: url, date, md5), so the gene
# sets behind a figure are reproducible from the repository alone -- no network call
# at analysis time.
#
# Used by: R/plot_reactome_stress_death_vs_phf1_distance.R

suppressPackageStartupMessages({
  library(tibble)
})

source("R/reactome_pcd_family.R")   # get_reactome_file(), norm_name(), REACTOME_CACHE

REACTOME_GMT_DEFAULT <- "genesets/enrichr_gmt/Reactome_Pathways_2024.gmt"

#' Resolve Reactome stable IDs to gene vectors via the cached Enrichr GMT.
#'
#' @param ids Named character vector: names are the module keys you want back,
#'   values are Reactome stable IDs, e.g. c(intrinsic_apop = "R-HSA-109606").
#' @param gmt_path Cached Enrichr Reactome GMT.
#' @param species Species column value in ReactomePathways.txt.
#' @param refresh TRUE to re-download the Reactome hierarchy file.
#' @return list(genes = <named list of character>, provenance = <tibble>)
#'   provenance: module, reactome_id, reactome_name, enrichr_term, n_genes
#'
#' Fails loudly on any unresolved ID. A silently dropped pathway would become a
#' silently missing module in the figure, which is worse than not running.
reactome_genesets <- function(ids,
                              gmt_path = REACTOME_GMT_DEFAULT,
                              species  = "Homo sapiens",
                              refresh  = FALSE) {
  stopifnot("ids must be a named character vector" =
              is.character(ids) && length(ids) > 0 &&
              !is.null(names(ids)) && all(nzchar(names(ids))))
  if (anyDuplicated(names(ids)))
    stop("Duplicate module keys in ids: ",
         paste(unique(names(ids)[duplicated(names(ids))]), collapse = ", "))
  if (!file.exists(gmt_path))
    stop("Cached Reactome GMT not found: ", gmt_path,
         "\n  Download it once from",
         "\n  https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_Pathways_2024",
         "\n  (the published run used the copy fetched 2026-08-10, md5 c30ddaf5ce1a1e8fd7f7850b39ebe4f6).")

  # quote="" / comment.char="": Reactome pathway names contain apostrophes and '#'.
  # The file has no header; columns are stId / name / species.
  pth <- utils::read.delim(get_reactome_file("ReactomePathways.txt", refresh = refresh),
                           header = FALSE, quote = "", comment.char = "",
                           col.names = c("reactome_id", "reactome_name", "species"),
                           stringsAsFactors = FALSE)
  pth <- pth[pth$species == species, , drop = FALSE]
  if (!nrow(pth)) stop("No rows for species '", species, "' in ReactomePathways.txt")

  gmt <- suppressWarnings(fgsea::gmtPathways(gmt_path))
  gmt <- lapply(gmt, function(g) unique(g[nzchar(g)]))
  if (!length(gmt)) stop("Cached GMT parsed to zero sets: ", gmt_path)
  gmt_key <- setNames(names(gmt), norm_name(names(gmt)))

  idx  <- match(unname(ids), pth$reactome_id)
  nm   <- pth$reactome_name[idx]
  term <- unname(gmt_key[norm_name(nm)])

  bad_id <- names(ids)[is.na(idx)]
  if (length(bad_id))
    stop("Reactome ID not found for species '", species, "': ",
         paste(sprintf("%s (%s)", bad_id, ids[bad_id]), collapse = ", "))
  bad_term <- which(is.na(term))
  if (length(bad_term))
    stop("Reactome pathway has no matching term in ", basename(gmt_path), ": ",
         paste(sprintf("%s (%s = '%s')", names(ids)[bad_term], ids[bad_term], nm[bad_term]),
               collapse = "; "))

  genes <- setNames(lapply(term, function(t) gmt[[t]]), names(ids))
  list(genes      = genes,
       provenance = tibble(module        = names(ids),
                           reactome_id   = unname(ids),
                           reactome_name = nm,
                           enrichr_term  = term,
                           n_genes       = as.integer(lengths(genes))))
}
