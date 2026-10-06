# ---------------------------------------------------------------------------
# reactome_catalogue_build.R
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
#' Build the Reactome programme catalogue. Sourceable; pure definitions.
#'
#' Kept as a sourceable library so the catalogue can be built from inside a
#' plotting script without shelling out or reading a saved file. There is ONE
#' implementation, used by:
#'   R/reactome_filtered_modules.R             (derives plotting modules from it)
#'
#' Everything here is a property of the Reactome release and the gene panel. No
#' expression, no celltype, no DEG run. Given the same cached hierarchy, the same
#' GMT and the same panel vector, it returns the same table every time -- which is
#' what lets a figure script reproduce the catalogue's gene lists exactly rather
#' than approximately.
#'
#' Requires dplyr / tibble, and R/reactome_pcd_family.R for get_reactome_file()
#' and norm_name().
################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble)
})

REACTOME_ROOTS_DEFAULT <- c(pcd       = "R-HSA-5357801",   # Programmed Cell Death
                            autophagy = "R-HSA-9612973",   # Autophagy
                            stress    = "R-HSA-8953897")   # Cellular responses to stimuli
REACTOME_DISEASE_ROOT  <- "R-HSA-1643685"

#' Crude HGNC family root. Stage 1 strips a trailing numeric run and any letters
#' after it (CASP3 -> CASP, MAPK1 -> MAPK). Stage 2, only if stage 1 changed
#' nothing and the symbol is long enough to survive it, strips one trailing
#' letter (YWHAB/YWHAE/YWHAZ -> YWHA).
#'
#' Deliberately crude and WILL both over- and under-merge (BCL2 -> BCL but
#' BCL2L1 -> BCL2L, so they do not group). It exists to make a set that is really
#' one gene family repeated -- YWHA* at 5 of 7 -- obvious at a glance, not to be a
#' nomenclature authority. The winning root and its members are both returned so
#' every call it makes is auditable.
reactome_family_root <- function(x) {
  r <- sub("[0-9]+[A-Za-z]*$", "", x)
  r <- ifelse(r == x & nchar(x) >= 5L, substr(x, 1L, nchar(x) - 1L), r)
  ifelse(nchar(r) >= 3L, r, x)
}

#' Build the catalogue: one row per pathway in the union of the root subtrees.
#'
#' @param panel Character vector of panel gene symbols. REQUIRED -- the catalogue
#'   is defined relative to a panel and there is no sensible default.
#' @param gmt_path Cached Enrichr Reactome GMT.
#' @param roots Named character vector of programme roots (key = programme name).
#' @param disease_root Reactome Disease root; descendants are FLAGGED, not dropped.
#' @param refresh TRUE to re-download the Reactome hierarchy files.
#' @param verbose Print progress.
#' @return tibble, one row per pathway, with identity / hierarchy / gene content /
#'   family columns. `retained` and the redundancy columns are NOT set here -- run
#'   filter_catalogue() from R/reactome_catalogue_utils.R for those.
build_reactome_catalogue <- function(panel,
                                     gmt_path     = "genesets/enrichr_gmt/Reactome_Pathways_2024.gmt",
                                     roots        = REACTOME_ROOTS_DEFAULT,
                                     disease_root = REACTOME_DISEASE_ROOT,
                                     refresh      = FALSE,
                                     verbose      = TRUE) {

  stopifnot("panel must be a non-empty character vector" =
              is.character(panel) && length(panel) > 0)
  panel <- unique(panel[nzchar(panel) & !is.na(panel)])
  say <- function(...) if (verbose) message(...)

  ##  -- 1. Hierarchy ---------------------------------------------------------
  say("== Reading Reactome hierarchy ==")

  # quote="" / comment.char="": pathway names contain apostrophes and '#'.
  pth <- utils::read.delim(get_reactome_file("ReactomePathways.txt", refresh = refresh),
                           header = FALSE, quote = "", comment.char = "",
                           col.names = c("stId", "pathway_name", "species"),
                           stringsAsFactors = FALSE)
  rel <- utils::read.delim(get_reactome_file("ReactomePathwaysRelation.txt", refresh = refresh),
                           header = FALSE, quote = "", comment.char = "",
                           col.names = c("parent", "child"), stringsAsFactors = FALSE)

  pth <- pth[pth$species == "Homo sapiens", , drop = FALSE]
  rel <- rel[startsWith(rel$parent, "R-HSA-") & startsWith(rel$child, "R-HSA-"), , drop = FALSE]
  rel <- unique(rel)
  pth$pathway_name <- trimws(pth$pathway_name)   # Reactome carries stray trailing spaces

  # A release that retires or re-IDs a root must fail loudly, not silently return
  # a one-term "programme".
  miss <- c(roots, disease = disease_root)
  miss <- miss[!miss %in% pth$stId]
  if (length(miss))
    stop("Reactome root(s) not found among Homo sapiens pathways: ",
         paste(sprintf("%s (%s)", names(miss), miss), collapse = ", "),
         "\n  The release may have re-IDed them -- check ",
         file.path(REACTOME_CACHE, "ReactomePathways.txt"), " before changing the roots.")

  kids    <- split(rel$child,  rel$parent)
  parents <- split(rel$parent, rel$child)

  # The hierarchy is a DAG, not a tree, so the visited set is what stops the walk.
  descendants <- function(root) {
    ids <- character(0); queue <- root
    while (length(queue)) {
      node <- queue[1]; queue <- queue[-1]
      if (node %in% ids) next
      ids <- c(ids, node); queue <- c(queue, kids[[node]])
    }
    ids
  }

  # Shortest (BFS) and longest (DAG DP) path from root. A single "depth" is
  # undefined in a DAG -- a node with two parents at different levels genuinely
  # sits at two depths -- so both are reported and neither is averaged. The
  # longest path uses Kahn topological order, which also gives a free cycle check.
  depths_from <- function(root) {
    nodes <- descendants(root)
    inset <- setNames(rep(TRUE, length(nodes)), nodes)

    dmin <- setNames(rep(NA_integer_, length(nodes)), nodes); dmin[root] <- 0L
    frontier <- root; lev <- 0L
    while (length(frontier)) {
      lev <- lev + 1L
      nxt <- unique(unlist(kids[frontier], use.names = FALSE))
      nxt <- nxt[!is.na(nxt)]
      nxt <- nxt[!is.na(inset[nxt]) & is.na(dmin[nxt])]
      if (!length(nxt)) break
      dmin[nxt] <- lev; frontier <- nxt
    }

    ch <- lapply(nodes, function(n) intersect(kids[[n]], nodes)); names(ch) <- nodes
    indeg <- setNames(integer(length(nodes)), nodes)
    for (n in nodes) for (c in ch[[n]]) indeg[c] <- indeg[c] + 1L

    dmax <- setNames(rep(-Inf, length(nodes)), nodes); dmax[root] <- 0
    queue <- names(indeg)[indeg == 0L]; ordered <- 0L
    while (length(queue)) {
      n <- queue[1]; queue <- queue[-1]; ordered <- ordered + 1L
      for (c in ch[[n]]) {
        if (is.finite(dmax[n])) dmax[c] <- max(dmax[c], dmax[n] + 1)
        indeg[c] <- indeg[c] - 1L
        if (indeg[c] == 0L) queue <- c(queue, c)
      }
    }
    if (ordered < length(nodes))
      stop("Cycle detected in the Reactome subtree below ", root,
           " (", length(nodes) - ordered, " nodes unordered). Longest-path depth ",
           "is undefined -- inspect ReactomePathwaysRelation.txt.")

    dmax[!is.finite(dmax)] <- NA_real_
    tibble(stId = nodes, depth_min = unname(dmin[nodes]),
           depth_max = as.integer(unname(dmax[nodes])))
  }

  prog <- lapply(names(roots), function(k) depths_from(roots[[k]]) %>% mutate(programme = k))
  names(prog) <- names(roots)
  prog_ids <- lapply(prog, `[[`, "stId")

  # A node can in principle sit under more than one programme (DAG): depth_min is
  # the min over programmes, depth_max the max, and root_programme lists all.
  depth_tbl <- bind_rows(prog) %>%
    group_by(stId) %>%
    summarise(depth_min = min(depth_min, na.rm = TRUE),
              depth_max = max(depth_max, na.rm = TRUE),
              root_programme = paste(sort(unique(programme)), collapse = "|"),
              .groups = "drop")

  disease_ids <- descendants(disease_root)

  hier <- tibble(stId = sort(unique(unlist(prog_ids, use.names = FALSE)))) %>%
    left_join(pth[, c("stId", "pathway_name", "species")], by = "stId") %>%
    filter(!is.na(pathway_name)) %>%          # non-Hs ids reachable via the DAG
    left_join(depth_tbl, by = "stId") %>%
    mutate(
      # parents / children come from the FULL relation graph, not the extracted
      # subtree: a set can be a leaf inside this subtree while having children
      # elsewhere, and calling that a leaf would be wrong.
      parent_ids = lapply(stId, function(i) sort(unique(parents[[i]]))),
      n_parents  = lengths(parent_ids),
      n_children = vapply(stId, function(i) length(unique(kids[[i]])), integer(1)),
      node_type  = ifelse(n_children == 0L, "leaf", "internal"),
      is_disease = stId %in% disease_ids)

  for (k in names(roots)) hier[[paste0("in_", k)]] <- hier$stId %in% prog_ids[[k]]
  in_cols <- paste0("in_", names(roots))
  hier <- hier %>% arrange(across(all_of(in_cols), ~ !.x), depth_min, pathway_name)

  say(sprintf("  %d Homo sapiens pathways across %d roots (%s); %d disease-flagged",
              nrow(hier), length(roots), paste(names(roots), collapse = ", "),
              sum(hier$is_disease)))

  ##  -- 2. Gene content ------------------------------------------------------
  if (!file.exists(gmt_path))
    stop("Cached Reactome GMT not found: ", gmt_path,
         "\n  Download it once from",
         "\n  https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_Pathways_2024",
         "\n  (the published run used the copy fetched 2026-08-10, md5 c30ddaf5ce1a1e8fd7f7850b39ebe4f6).")
  gmt <- suppressWarnings(fgsea::gmtPathways(gmt_path))
  gmt <- lapply(gmt, function(g) unique(g[nzchar(g)]))
  if (!length(gmt)) stop("Cached GMT parsed to zero sets: ", gmt_path)

  # Enrichr title-cases Reactome names and strips punctuation; norm_name()
  # collapses both sides to the same key. See R/reactome_pcd_family.R.
  gmt_key <- setNames(names(gmt), norm_name(names(gmt)))
  hier$enrichr_term <- unname(gmt_key[norm_name(hier$pathway_name)])
  hier$gmt_matched  <- !is.na(hier$enrichr_term)
  hier$genes_reactome <- lapply(hier$enrichr_term,
                                function(t) if (is.na(t)) character(0) else gmt[[t]])
  hier$n_genes_reactome <- lengths(hier$genes_reactome)
  say(sprintf("  %d/%d pathways matched a term in %s",
              sum(hier$gmt_matched), nrow(hier), basename(gmt_path)))

  ##  -- 3. Panel + HGNC alias resolution -------------------------------------
  # ALIAS2EG carries current official symbols as well as withdrawn/previous ones,
  # so mapping BOTH sides to Entrez and matching there resolves a GMT symbol that
  # is simply an older name for a panel gene (PARK2 -> PRKN, IL8 -> CXCL8,
  # H2AFX -> H2AX, UFD1L -> UFD1). Aliases mapping to more than one Entrez are
  # ambiguous and are NOT used -- a wrong rescue is worse than a miss.
  alias_tbl <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egALIAS2EG)
  alias_n   <- table(alias_tbl$alias_symbol)
  alias_uni <- alias_tbl[alias_tbl$alias_symbol %in% names(alias_n)[alias_n == 1L], ]
  sym2eg    <- setNames(alias_uni$gene_id, alias_uni$alias_symbol)

  panel_eg <- sym2eg[panel]; panel_eg <- panel_eg[!is.na(panel_eg)]
  eg2panel <- setNames(names(panel_eg), unname(panel_eg))
  eg2panel <- eg2panel[!duplicated(names(eg2panel))]

  resolve_panel <- function(g) {
    direct <- intersect(g, panel)
    rest   <- setdiff(g, panel)
    eg     <- sym2eg[rest]
    hit    <- !is.na(eg) & eg %in% names(eg2panel)
    from   <- rest[hit]
    to     <- unname(eg2panel[eg[hit]])
    keep   <- !is.na(to) & nzchar(to)
    list(genes    = sort(unique(c(direct, to[keep]))),
         rescued  = unique(from[keep]),
         unmapped = setdiff(rest, from[keep]))
  }

  pmap <- lapply(hier$genes_reactome, resolve_panel)
  hier$genes_panel     <- lapply(pmap, `[[`, "genes")
  hier$n_genes_panel   <- lengths(hier$genes_panel)
  hier$n_alias_rescued <- vapply(pmap, function(x) length(x$rescued), integer(1))
  hier$n_unmapped      <- vapply(pmap, function(x) length(x$unmapped), integer(1))
  hier$genes_unmapped  <- lapply(pmap, `[[`, "unmapped")
  hier$pct_panel <- ifelse(hier$n_genes_reactome > 0,
                           hier$n_genes_panel / hier$n_genes_reactome, NA_real_)
  say(sprintf("  panel %d genes; alias-rescued %d, unmapped %d (summed over sets)",
              length(panel), sum(hier$n_alias_rescued), sum(hier$n_unmapped)))

  ##  -- 4. Family dominance --------------------------------------------------
  fam <- lapply(hier$genes_panel, function(g) {
    if (!length(g)) return(list(frac = NA_real_, root = NA_character_, genes = NA_character_))
    r <- reactome_family_root(g)
    t <- sort(table(r), decreasing = TRUE)
    list(frac = as.numeric(t[1]) / length(g), root = names(t)[1],
         genes = paste(sort(g[r == names(t)[1]]), collapse = ";"))
  })
  hier$max_family_frac  <- vapply(fam, `[[`, numeric(1), "frac")
  hier$family_top       <- vapply(fam, `[[`, character(1), "root")
  hier$family_top_genes <- vapply(fam, `[[`, character(1), "genes")

  attr(hier, "roots")    <- roots
  attr(hier, "gmt_path") <- gmt_path
  attr(hier, "n_panel")  <- length(panel)
  hier
}
