# ---------------------------------------------------------------------------
# reactome_catalogue_utils.R
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
#' Reactome programme catalogue -- filtering and redundancy utilities.
#'
#' Sourced by both the catalogue builder and the filtering step so a set cannot be
#' scored one way at extraction and another way at filtering:
#'   R/reactome_catalogue_build.R       (build the catalogue)
#'   R/reactome_filtered_modules.R      (apply a filter and derive modules)
#'
#' THE OPERATIVE GENE LIST IS `genes_panel` -- the Reactome set intersected with
#' the CosMx panel (after HGNC alias resolution). Nothing here is conditioned on
#' expression, detection rate, celltype or any DEG run: a set's size, redundancy
#' and family composition are properties of the panel, which is fixed. That is
#' what makes the catalogue one row per pathway instead of one row per pathway
#' per celltype.
#'
#' WHY FILTERING IS A FUNCTION AND NOT A SPREADSHEET SUBSET. Four columns are NOT
#' properties of a pathway: jaccard_max, jaccard_partner, cluster_id and
#' n_genes_unique are all defined RELATIVE TO THE OTHER SETS YOU KEPT. Raise the
#' size floor and a set that had zero unique genes can acquire twenty, because
#' the sets that were covering them dropped out. Subsetting the table by hand
#' leaves those columns describing the filter that built the file. Every filter
#' change recomputes them here.
#'
#' Defines only functions and constants -- sourcing runs no analysis.
#' Requires dplyr / tibble / igraph.
################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble)
})

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

isTRUE_vec <- function(x) !is.na(x) & x

##  ............................................................................
##  The filter axes                                                         ####

#' Every tunable filter, in one place: name, default (= no-op), and what it means.
#'
#' Anything listed here is settable three ways and behaves identically in all
#' three: as an argument to filter_catalogue(), as a --<name> CLI flag on a
#' filter driver, and as a column of a sweep grid. Add a filter
#' here and it appears in all three at once.
#'
#' `default` is a list column because the defaults are deliberately of mixed type
#' (numeric / character / logical) -- one filter per row, each keeping its own.
CATALOGUE_FILTERS <- tibble::tibble(
  name = c("min_genes_panel", "max_genes_panel", "min_pct_panel", "min_genes_unique",
           "max_jaccard", "min_depth", "max_depth", "node_type", "programmes",
           "exclude_disease", "require_gmt_match", "max_family_frac",
           "representatives_only", "expr"),
  default = list(5, Inf, 0, 0,
                 1, 0, Inf, "any", "any",
                 TRUE, TRUE, 1,
                 FALSE, ""),
  type = c("numeric", "numeric", "numeric", "numeric",
           "numeric", "numeric", "numeric", "character", "character",
           "logical", "logical", "numeric",
           "logical", "character"),
  description = c(
    "n_genes_panel >= x    set members present on the panel",
    "n_genes_panel <= x    drop huge, unspecific sets",
    "pct_panel >= x        fraction of the FULL Reactome set that is on the panel",
    "n_genes_unique >= x   panel members in no other surviving set  [PEELED]",
    "jaccard_max <= x      no surviving pair overlaps more than this  [PEELED]",
    "depth_min >= x        shortest path from the programme root",
    "depth_min <= x",
    "'any' | 'leaf' | 'internal'   node_type, from the FULL relation graph",
    "'any' or a |-separated subset of programme keys, e.g. 'pcd|autophagy'",
    "drop descendants of the Reactome Disease root (circular in an AD cohort)",
    "drop sets with no term in the gene-set library (they carry zero genes)",
    "max_family_frac <= x  drop sets dominated by one gene family, e.g. YWHA*",
    "after clustering, keep only one set per redundancy cluster",
    "free-form R expression over any catalogue column, ANDed with the rest")
)

#' Defaults as a plain named list, ready to be modified and passed on.
catalogue_filter_defaults <- function() {
  setNames(CATALOGUE_FILTERS$default, CATALOGUE_FILTERS$name)
}

#' Print the filter menu -- the "what can I turn?" reference.
print_catalogue_filters <- function() {
  d <- CATALOGUE_FILTERS
  cat(sprintf("  %-22s %-10s %s\n", "FILTER", "DEFAULT", "MEANING"))
  for (i in seq_len(nrow(d)))
    cat(sprintf("  %-22s %-10s %s\n", d$name[i],
                format(d$default[[i]], trim = TRUE), d$description[i]))
  invisible(d)
}

##  ............................................................................
##  Redundancy, recomputed against whatever survived                        ####

jaccard2 <- function(a, b) {
  u <- length(union(a, b))
  if (u == 0L) return(NA_real_)
  length(intersect(a, b)) / u
}

#' Recompute the filter-relative columns against the currently retained sets.
#'
#' Operates on rows where `retained` is TRUE and leaves the rest NA -- a Jaccard
#' against sets you excluded is not a meaningful number, and filling it in would
#' invite exactly the mistake this function exists to prevent.
#'
#' Clusters are connected components of the graph whose edges are pairs at or
#' above `jaccard_threshold`. Components are transitive by construction: A~B and
#' B~C puts A, B and C together even if A and C fall below the threshold. Unlike
#' greedy pairwise collapsing, this does not depend on the order in which the sets
#' are visited.
annotate_redundancy <- function(tbl,
                                jaccard_threshold = 0.7,
                                representative    = c("largest", "smallest")) {
  representative <- match.arg(representative)
  stopifnot(all(c("stId", "retained", "genes_panel") %in% names(tbl)))

  tbl$jaccard_max     <- NA_real_
  tbl$jaccard_partner <- NA_character_
  tbl$cluster_id      <- NA_character_
  tbl$cluster_size    <- NA_integer_
  tbl$is_cluster_representative <- NA
  tbl$n_genes_unique  <- NA_integer_

  keep <- which(isTRUE_vec(tbl$retained))
  if (!length(keep)) return(tbl)

  gk  <- tbl$genes_panel[keep]
  ids <- tbl$stId[keep]
  k   <- length(keep)

  if (k >= 2L) {
    J <- matrix(0, k, k, dimnames = list(ids, ids))
    for (a in seq_len(k - 1L)) for (b in seq(a + 1L, k)) {
      J[a, b] <- J[b, a] <- jaccard2(gk[[a]], gk[[b]])
    }
    diag(J) <- NA_real_
    tbl$jaccard_max[keep] <- apply(J, 1, function(r) if (all(is.na(r))) NA_real_ else max(r, na.rm = TRUE))
    tbl$jaccard_partner[keep] <- ids[apply(J, 1, function(r) if (all(is.na(r))) NA_integer_ else which.max(r))]
    A <- !is.na(J) & J >= jaccard_threshold
  } else {
    A <- matrix(FALSE, 1, 1, dimnames = list(ids, ids))
  }
  diag(A) <- FALSE

  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  comp <- igraph::components(g)$membership
  tbl$cluster_id[keep]   <- sprintf("C%03d", comp)
  tbl$cluster_size[keep] <- as.integer(table(comp)[as.character(comp)])

  sz <- lengths(gk)
  rep_flag <- logical(k)
  for (cc in unique(comp)) {
    m <- which(comp == cc)
    # Ties broken by stId so the pick is deterministic across runs and machines.
    pick <- if (representative == "largest") m[order(-sz[m], ids[m])][1] else m[order(sz[m], ids[m])][1]
    rep_flag[pick] <- TRUE
  }
  tbl$is_cluster_representative[keep] <- rep_flag

  # n_genes_unique: panel members appearing in NO other surviving set. The single
  # most useful diagnostic in the table -- a set whose unique count is 0
  # contributed no genes its neighbours did not already contribute.
  tab <- table(unlist(gk, use.names = FALSE))
  tbl$n_genes_unique[keep] <- vapply(gk, function(g) sum(tab[g] == 1L), integer(1))
  tbl
}

#' Long table of every surviving pair at or above `report_floor`.
jaccard_pairs <- function(tbl, jaccard_threshold = 0.7, report_floor = 0.1) {
  keep <- which(isTRUE_vec(tbl$retained))
  if (length(keep) < 2L) return(NULL)
  gk <- tbl$genes_panel[keep]; ids <- tbl$stId[keep]; k <- length(keep)
  J <- matrix(0, k, k)
  for (a in seq_len(k - 1L)) for (b in seq(a + 1L, k)) J[a, b] <- J[b, a] <- jaccard2(gk[[a]], gk[[b]])
  pr <- which(upper.tri(J) & J >= report_floor, arr.ind = TRUE)
  if (!nrow(pr)) return(NULL)
  tibble(stId_a = ids[pr[, 1]], stId_b = ids[pr[, 2]],
         name_a = tbl$pathway_name[keep][pr[, 1]],
         name_b = tbl$pathway_name[keep][pr[, 2]],
         n_a = lengths(gk)[pr[, 1]], n_b = lengths(gk)[pr[, 2]],
         n_shared = mapply(function(i, j) length(intersect(gk[[i]], gk[[j]])), pr[, 1], pr[, 2]),
         jaccard = J[pr],
         collapsed = J[pr] >= jaccard_threshold) %>%
    arrange(desc(jaccard))
}

##  ............................................................................
##  The one entry point                                                     ####

#' Apply a filter set to the catalogue and rebuild the filter-relative columns.
#'
#' @param cat Catalogue tibble (from reactome_programme_catalogue.rds). List
#'   columns must be intact -- read the .rds, not the .tsv.
#' @param ... Any filter from CATALOGUE_FILTERS, by name. Unnamed/unknown
#'   arguments are an error, not a silent no-op: a typo'd filter name that did
#'   nothing would look like "the filter had no effect".
#' @param jaccard_threshold,representative Passed to annotate_redundancy().
#' @return The catalogue with `retained` set and jaccard_max / jaccard_partner /
#'   cluster_id / cluster_size / is_cluster_representative / n_genes_unique
#'   recomputed against the surviving sets.
#'
#' @examples
#'   cat <- readRDS("results/reactome_programme_catalogue/reactome_programme_catalogue.rds")
#'   a <- filter_catalogue(cat, min_genes_panel = 10)
#'   b <- filter_catalogue(cat, min_genes_panel = 10, node_type = "leaf", min_depth = 2)
#'   d <- filter_catalogue(cat, expr = "n_genes_panel >= 8 & max_family_frac < 0.5")
filter_catalogue <- function(cat, ...,
                             jaccard_threshold = 0.7,
                             representative    = c("largest", "smallest")) {
  representative <- match.arg(representative)
  user <- list(...)
  if (length(user) && (is.null(names(user)) || any(!nzchar(names(user)))))
    stop("All filter arguments must be named. See print_catalogue_filters().")
  bad <- setdiff(names(user), CATALOGUE_FILTERS$name)
  if (length(bad))
    stop("Unknown filter(s): ", paste(bad, collapse = ", "),
         "\n  Available: ", paste(CATALOGUE_FILTERS$name, collapse = ", "))

  f <- modifyList(catalogue_filter_defaults(), user)
  in_cols <- grep("^in_", names(cat), value = TRUE)

  keep <-
    cat$n_genes_panel >= f$min_genes_panel &
    cat$n_genes_panel <= f$max_genes_panel &
    coalesce(cat$pct_panel, 0) >= f$min_pct_panel &
    coalesce(cat$max_family_frac, 0) <= f$max_family_frac &
    coalesce(cat$depth_min, 0) >= f$min_depth &
    coalesce(cat$depth_min, 0) <= f$max_depth

  if (!identical(f$node_type, "any")) keep <- keep & cat$node_type == f$node_type
  if (isTRUE(f$exclude_disease))      keep <- keep & !cat$is_disease
  if (isTRUE(f$require_gmt_match) && "gmt_matched" %in% names(cat))
    keep <- keep & (is.na(cat$gmt_matched) | cat$gmt_matched)

  if (!identical(f$programmes, "any")) {
    want <- paste0("in_", trimws(strsplit(f$programmes, "|", fixed = TRUE)[[1]]))
    miss <- setdiff(want, in_cols)
    if (length(miss))
      stop("Unknown programme(s): ", paste(sub("^in_", "", miss), collapse = ", "),
           "\n  Available: ", paste(sub("^in_", "", in_cols), collapse = ", "))
    keep <- keep & Reduce(`|`, lapply(want, function(w) cat[[w]]))
  }

  if (nzchar(f$expr %||% "")) {
    e <- eval(parse(text = f$expr), envir = cat, enclos = parent.frame())
    if (!is.logical(e) || length(e) != nrow(cat))
      stop("--expr must evaluate to a logical vector of length nrow(cat); got ",
           class(e)[1], " of length ", length(e))
    keep <- keep & e
  }

  cat$retained <- isTRUE_vec(keep)
  cat <- annotate_redundancy(cat, jaccard_threshold = jaccard_threshold,
                             representative = representative)

  # representatives_only is applied AFTER clustering (it is defined by the
  # clustering) and then everything is recomputed once more, so n_genes_unique
  # describes the final set rather than the pre-collapse one.
  if (isTRUE(f$representatives_only)) {
    cat$retained <- isTRUE_vec(cat$retained) & isTRUE_vec(cat$is_cluster_representative)
    cat <- annotate_redundancy(cat, jaccard_threshold = jaccard_threshold,
                               representative = representative)
  }

  # max_jaccard CANNOT be applied as a one-shot predicate either, and for the
  # sharper reason: jaccard_max is a property of a PAIR, so both members of an
  # over-lapping pair violate it. Dropping every violator at once deletes both
  # halves of every redundant pair and keeps neither. Peeled instead: remove one
  # member of the worst-overlapping pair, recompute, repeat -- which leaves one
  # representative of each pair standing, the intended behaviour.
  #
  # Which member goes is tied to `representative`, so this filter and
  # representatives_only cannot disagree about who survives: keeping the largest
  # means dropping the smaller of the pair, and vice versa.
  if (f$max_jaccard < 1) {
    repeat {
      k <- which(isTRUE_vec(cat$retained))
      if (length(k) < 2L) break
      viol <- k[coalesce(cat$jaccard_max[k], 0) > f$max_jaccard]
      if (!length(viol)) break
      ord <- if (representative == "largest") {
        order(-cat$jaccard_max[viol], cat$n_genes_panel[viol], cat$stId[viol])
      } else {
        order(-cat$jaccard_max[viol], -cat$n_genes_panel[viol], cat$stId[viol])
      }
      cat$retained[viol[ord][1]] <- FALSE
      cat <- annotate_redundancy(cat, jaccard_threshold = jaccard_threshold,
                                 representative = representative)
    }
  }

  # min_genes_unique CANNOT be applied as a one-shot predicate. n_genes_unique is
  # a function of the surviving set and can only RISE as sets are removed, so
  # testing it against the values you walked in with drops sets that would have
  # passed once their coverers were gone. It is peeled instead: remove the single
  # worst violator, recompute, repeat. One at a time rather than all at once,
  # because dropping a set can rescue its neighbours -- removing them together
  # discards sets that never actually failed.
  #
  # Runs AFTER the max_jaccard peel: removing redundant sets can only raise
  # n_genes_unique, so doing it in this order drops the fewest sets overall.
  if (f$min_genes_unique > 0) {
    repeat {
      k <- which(isTRUE_vec(cat$retained))
      if (!length(k)) break
      viol <- k[coalesce(cat$n_genes_unique[k], 0L) < f$min_genes_unique]
      if (!length(viol)) break
      # Worst first: fewest unique, then smallest set, then stId for determinism.
      worst <- viol[order(cat$n_genes_unique[viol], cat$n_genes_panel[viol], cat$stId[viol])][1]
      cat$retained[worst] <- FALSE
      cat <- annotate_redundancy(cat, jaccard_threshold = jaccard_threshold,
                                 representative = representative)
    }
  }

  attr(cat, "filters") <- f
  attr(cat, "jaccard_threshold") <- jaccard_threshold
  attr(cat, "representative") <- representative
  cat
}

#' One-line summary of what a filter did. The thing you actually read when
#' comparing filters.
summarise_filter <- function(cat, label = NA_character_) {
  r <- isTRUE_vec(cat$retained)
  tibble(
    label       = label,
    n_sets      = nrow(cat),
    n_retained  = sum(r),
    # n_clusters is how many INDEPENDENT tests the retained sets amount to.
    n_clusters  = n_distinct(cat$cluster_id[r]),
    n_redundant = sum(r & coalesce(cat$cluster_size, 1L) > 1L),
    median_genes = stats::median(cat$n_genes_panel[r]),
    min_genes    = suppressWarnings(min(cat$n_genes_panel[r])),
    max_genes    = suppressWarnings(max(cat$n_genes_panel[r])),
    # A set with no unique member cannot contribute genes its neighbours do not
    # already contribute. This is the number to watch.
    n_zero_unique = sum(r & coalesce(cat$n_genes_unique, 0L) == 0L),
    median_unique = stats::median(cat$n_genes_unique[r]),
    max_jaccard   = suppressWarnings(max(cat$jaccard_max[r], na.rm = TRUE)),
    max_fam_frac  = suppressWarnings(max(cat$max_family_frac[r], na.rm = TRUE))
  ) %>% mutate(across(where(is.numeric), ~ ifelse(is.finite(.x), .x, NA_real_)))
}

#' Run a grid of filter settings and return one summary row per setting.
#'
#' @param cat Catalogue tibble.
#' @param grid data.frame whose columns are filter names (from CATALOGUE_FILTERS)
#'   and whose rows are settings to try. Build it with expand.grid().
#' @param ... Fixed filters applied to every row of the grid.
#'
#' @examples
#'   g <- expand.grid(min_genes_panel = c(5, 10, 15), node_type = c("any", "leaf"),
#'                    stringsAsFactors = FALSE)
#'   sweep_filters(cat, g)
sweep_filters <- function(cat, grid, ..., jaccard_threshold = 0.7,
                          representative = "largest", verbose = TRUE) {
  bad <- setdiff(names(grid), CATALOGUE_FILTERS$name)
  if (length(bad)) stop("Unknown filter column(s) in grid: ", paste(bad, collapse = ", "))
  fixed <- list(...)

  out <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    row <- as.list(grid[i, , drop = FALSE])
    row <- lapply(row, function(x) if (is.factor(x)) as.character(x) else x)
    lab <- paste(sprintf("%s=%s", names(row), vapply(row, function(x) format(x, trim = TRUE), "")),
                 collapse = ", ")
    if (verbose) message(sprintf("  [%d/%d] %s", i, nrow(grid), lab))
    f <- do.call(filter_catalogue,
                 c(list(cat), row, fixed,
                   list(jaccard_threshold = jaccard_threshold, representative = representative)))
    out[[i]] <- bind_cols(summarise_filter(f, label = lab), as_tibble(row))
  }
  bind_rows(out)
}

#' Collapse list columns so a tibble can be written as TSV.
flatten_list_cols <- function(df) {
  for (cn in names(df)) if (is.list(df[[cn]]))
    df[[cn]] <- vapply(df[[cn]], function(x) paste(x, collapse = ";"), character(1))
  as.data.frame(df)
}
