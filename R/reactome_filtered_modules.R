# ---------------------------------------------------------------------------
# reactome_filtered_modules.R
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
#' Reactome plotting modules, derived by RULE rather than hand-picked.
#'
#' Sourced by both figure scripts so they cannot disagree about the pathway list:
#'   R/plot_reactome_stress_death_vs_phf1_distance.R  -> plots/reactome_stress_death_vs_phf1_distance_1000um/
#'   R/plot_reactome_gene_contribution.R              -> plots/reactome_gene_contribution_1000um/
#'
#' WHY RULE-DERIVED. The module list is DERIVED: walk the three programme subtrees,
#' then apply CATALOGUE_RULES. The selection is stated as rules a reader can check
#' rather than a list a reader must take on trust, and it cannot silently drift from
#' the catalogue. Consequences of the rules are noted under NOTES below.
#'
#' SELF-CONTAINED AND EXACT. This rebuilds the catalogue in-process from the cached,
#' provenance-stamped Reactome hierarchy + GMT; it does NOT read the saved
#' results/reactome_programme_catalogue/*.rds. The figure scripts therefore run on
#' their own. Because build_reactome_catalogue() is the same function the extraction
#' script calls, and the CosMx panel is identical between seu_PHF1's SCT assay and
#' the celltype SCEs (verified: 6175 genes, setequal), the gene lists it
#' returns are IDENTICAL to the on-disk catalogue's `genes_panel`, not merely similar.
#'
#' GENES COME FROM THE CATALOGUE, NOT FROM A SECOND RESOLUTION PASS. Callers must
#' use `$genes` rather than re-resolving through reactome_genesets() + intersect().
#' The catalogue applies HGNC alias resolution (PARK2 -> PRKN, IL8 -> CXCL8,
#' H2AFX -> H2AX, UFD1L -> UFD1); a plain intersect() does not, and would silently
#' drop PRKN from PINK1-PRKN Mediated Mitophagy -- the gene the set is named after.
#'
#' NOTES on the rule-derived list:
#'   - There is no UBIQUITIN-PROTEASOME module. R-HSA-983168 ("Antigen
#'     processing: Ubiquitination & Proteasome degradation") is filed by Reactome
#'     under Immune System / Metabolism of proteins -- outside all three roots, so
#'     no rule over these subtrees reaches it.
#'   - The individual UPR arms are not separate modules: ATF6 (10 genes) and
#'     IRE1alpha fall below the size floor, so UPR appears as the parent term.
#'   - Sets are chosen for size, panel coverage and non-redundancy, not for
#'     mechanistic interest.
#'
#' Requires R/reactome_pcd_family.R, R/reactome_catalogue_build.R,
#' R/reactome_catalogue_utils.R (sourced here).
################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble)
})

source("R/reactome_pcd_family.R")        # get_reactome_file(), norm_name()
source("R/reactome_catalogue_build.R")   # build_reactome_catalogue()
source("R/reactome_catalogue_utils.R")   # filter_catalogue(), CATALOGUE_FILTERS

##  ............................................................................
##  THE RULES                                                               ####

#' The selection rules, in ONE place.
#'
#' Yields 20 sets: 11 stress, 6 PCD, 3 autophagy, every one its own redundancy
#' cluster.
#'
#'   min_genes_panel 25   enough panel genes for a stable module score
#'   max_genes_panel 80   excludes broad parents that are not a "programme"
#'   min_pct_panel  0.70  the panel must represent most of the real pathway,
#'                        so the score is the pathway and not a fragment of it
#'   max_jaccard    0.80  no surviving pair overlaps more than this (PEELED)
#'   representatives_only one set per redundancy cluster
#'   n_parents > 0        excludes the three programme roots themselves
#'
#' exclude_disease and require_gmt_match are left at their TRUE defaults.
#' min_genes_unique is deliberately NOT applied -- see the note below.
CATALOGUE_RULES <- list(
  min_genes_panel      = 25,
  max_genes_panel      = 80,
  min_pct_panel        = 0.70,
  max_jaccard          = 0.80,
  representatives_only = TRUE,
  expr                 = "n_parents > 0"
)

#' Note carried with the rules, printed into every stats log that uses them.
#' 5 of the 20 retained sets have NO panel gene absent from every other retained
#' set (Aggrephagy, RIPK1-mediated regulated necrosis, Mitophagy, Activation of
#' BH3-only proteins, HSF1-dependent transactivation), because clustering only
#' fires at Jaccard >= 0.7 and a set can be fully covered by several sub-threshold
#' partners, so their scores overlap with their neighbours'. `min_genes_unique`
#' would remove them; it is not applied because it would also drop named modalities.
CATALOGUE_RULES_CAVEAT <- paste(
  "Sets were selected by rule (size, panel coverage, non-redundancy), NOT for",
  "mechanistic interest. 5 of 20 retained sets carry no gene unique to them among",
  "the retained sets, so their scores overlap with their neighbours'.",
  "The ubiquitin-proteasome module is absent: Reactome files it outside all three",
  "programme roots, so no rule over these subtrees can reach it.")

##  ............................................................................
##  Palette                                                                 ####

#' Colourblind-safe qualitative palette, 12 colours: Okabe-Ito (8) extended with
#' Paul Tol muted (olive, indigo, teal, wine). The pale yellow #F0E442 is ordered
#' LAST-but-three because it is illegible as a thin line on white.
#' 12 rather than the sibling scripts' 10 because the rule-derived stress group has
#' 11 modules.
MODULE_PALETTE <- c(
  "#E69F00",  # orange           (Okabe-Ito)
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#CC79A7",  # reddish purple
  "#000000",  # black
  "#999933",  # olive            (Tol)
  "#332288",  # indigo           (Tol)
  "#44AA99",  # teal             (Tol)
  "#882255",  # wine             (Tol)
  "#F0E442"   # yellow           -- last: poor line contrast on white
)

##  ............................................................................
##  Helpers                                                                 ####

#' Readable, unique, filesystem-safe module key from a pathway name.
#' Keys end up in filenames (plot_genecontrib_<key>_<celltype>.pdf), so they must
#' be stable and free of punctuation. Collisions get a numeric suffix rather than
#' being silently merged.
module_key <- function(names_vec, max_chars = 30L) {
  k <- tolower(names_vec)
  k <- gsub("\\(.*?\\)", " ", k)                 # drop parenthetical asides
  k <- gsub("[^a-z0-9]+", "_", k)
  k <- gsub("^_+|_+$", "", k)
  k <- substr(k, 1L, max_chars)
  k <- gsub("_+$", "", k)
  dup <- duplicated(k)
  while (any(dup)) {                             # never merge two pathways
    k[dup] <- paste0(k[dup], "_", ave(k, k, FUN = seq_along)[dup])
    dup <- duplicated(k)
  }
  k
}

#' Legend label: the Reactome name, shortened where it runs long. Full names go to
#' the stats log and the overlap table, so nothing is lost.
#'
#' 41 characters is chosen, not rounded to: it is the exact length of "Detoxification of
#' Reactive Oxygen Species", the longest name worth keeping whole. One character fewer and
#' that label starts losing "Species". It also clips "HSP90 chaperone cycle for steroid
#' hormone receptors (SHR) in the presence of ligand" (82 chars), which is the intent --
#' left whole, that one name would set the width of every right-hand legend by itself and
#' steal the space from the x axis.
MODULE_LABEL_CHARS <- 41L
module_label <- function(x, max_chars = MODULE_LABEL_CHARS)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1L, max_chars - 3L), "..."), x)

#' Legend title for each programme group: the Reactome ROOT's own pathway name.
#'
#' Read from the catalogue rather than hard-coded, so the legend cannot drift from the
#' subtree it names and a Reactome rename propagates automatically. Reactome's exact
#' strings are title-cased for some roots and sentence-cased for others ("Programmed
#' Cell Death" but "Cellular responses to stimuli"); that inconsistency is Reactome's
#' and is reproduced verbatim, because a label that matches the source exactly is
#' checkable and one that has been tidied is not.
#'
#' @param cat_tbl Unfiltered catalogue (must still contain the root rows -- the rules
#'   drop them from the module list, which is why this is computed before filtering).
#' @param roots Named character vector, programme key -> root stable ID.
group_labels_from_roots <- function(cat_tbl, roots) {
  nm <- setNames(cat_tbl$pathway_name, cat_tbl$stId)
  lab <- unname(nm[roots])
  miss <- names(roots)[is.na(lab)]
  if (length(miss))
    stop("Root pathway name not found for programme(s): ", paste(miss, collapse = ", "),
         ". The catalogue must be built before the roots are filtered out.")
  setNames(lab, names(roots))
}

##  ............................................................................
##  The entry point                                                         ####

#' Build the catalogue, apply CATALOGUE_RULES, return a plotting module config.
#'
#' @param panel Character vector of panel gene symbols (e.g. rownames(seu[["SCT"]])).
#' @param gmt_path Cached Enrichr Reactome GMT.
#' @param rules Filter list; defaults to CATALOGUE_RULES.
#' @param verbose Print progress.
#' @return list(
#'   REACTOME_PRIMARY = named chr, module key -> Reactome stable ID
#'   MODULE_GROUPS    = list, group -> module keys (plotting/legend order)
#'   MODULE_LABELS    = named chr, module key -> legend label
#'   MODULE_NAMES     = named chr, module key -> FULL Reactome name
#'   genes            = named list, module key -> panel-resolved gene vector
#'                      (USE THIS; do not re-resolve through the GMT)
#'   catalogue        = the filtered catalogue tibble, all columns
#'   rules            = the rules actually applied)
reactome_filtered_modules <- function(panel,
                                      gmt_path = "genesets/enrichr_gmt/Reactome_Pathways_2024.gmt",
                                      rules    = CATALOGUE_RULES,
                                      roots    = REACTOME_ROOTS_DEFAULT,
                                      verbose  = TRUE) {

  cat_tbl <- build_reactome_catalogue(panel, gmt_path = gmt_path, roots = roots,
                                      verbose = verbose)
  # Taken from the UNFILTERED catalogue: the rules exclude the roots themselves
  # (n_parents > 0), so by the time `kept` exists their names are gone.
  root_labels <- group_labels_from_roots(cat_tbl, roots)

  filt <- do.call(filter_catalogue, c(list(cat_tbl), rules))
  kept <- filt[isTRUE_vec(filt$retained), ]

  if (!nrow(kept))
    stop("No pathway survived CATALOGUE_RULES. The rules or the panel have changed; ",
         "inspect before proceeding -- an empty module list is not a result.")

  # A pathway reachable from two roots would otherwise be plotted twice. These
  # subtrees are disjoint in the current release, so this is a guard, not a fix.
  multi <- grepl("|", kept$root_programme, fixed = TRUE)
  if (any(multi)) {
    warning(sum(multi), " pathway(s) sit under more than one programme; ",
            "assigning each to the first. Affected: ",
            paste(kept$stId[multi], collapse = ", "), call. = FALSE)
    kept$root_programme <- sub("\\|.*$", "", kept$root_programme)
  }

  # Largest first within a group: the legend then reads in the same order as the
  # retained-sets table in the catalogue's own stats log.
  kept <- kept %>% arrange(root_programme, desc(n_genes_panel), stId)

  keys <- module_key(kept$pathway_name)
  REACTOME_PRIMARY <- setNames(kept$stId,        keys)
  MODULE_LABELS    <- setNames(module_label(kept$pathway_name), keys)
  MODULE_NAMES     <- setNames(kept$pathway_name, keys)
  genes            <- setNames(kept$genes_panel,  keys)

  grp_order     <- intersect(names(roots), unique(kept$root_programme))
  grp_order     <- c(grp_order, setdiff(unique(kept$root_programme), grp_order))
  MODULE_GROUPS <- lapply(setNames(grp_order, grp_order),
                          function(g) keys[kept$root_programme == g])

  too_big <- lengths(MODULE_GROUPS) > length(MODULE_PALETTE)
  if (any(too_big))
    stop("Group(s) ", paste(names(MODULE_GROUPS)[too_big], collapse = ", "),
         " have more modules than MODULE_PALETTE has colours (",
         length(MODULE_PALETTE), "). Overlaid lines would reuse a colour and two ",
         "modules would be indistinguishable. Tighten CATALOGUE_RULES or split the group.")

  GROUP_LABELS <- root_labels[names(MODULE_GROUPS)]

  if (verbose) {
    message(sprintf("== %d modules from CATALOGUE_RULES ==", nrow(kept)))
    for (g in names(MODULE_GROUPS))
      message(sprintf("   %-10s %2d  [%s]: %s", g, length(MODULE_GROUPS[[g]]),
                      GROUP_LABELS[[g]], paste(MODULE_GROUPS[[g]], collapse = ", ")))
  }

  list(REACTOME_PRIMARY = REACTOME_PRIMARY,
       MODULE_GROUPS    = MODULE_GROUPS,
       MODULE_LABELS    = MODULE_LABELS,
       MODULE_NAMES     = MODULE_NAMES,
       GROUP_LABELS     = GROUP_LABELS,
       ROOT_IDS         = roots[names(MODULE_GROUPS)],
       genes            = genes,
       catalogue        = kept,
       rules            = rules)
}
