# ---------------------------------------------------------------------------
# reactome_pcd_family.R
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
#' Reactome Programmed Cell Death gene-set family -- shared resolver, plus the
#' Reactome download and name-normalisation helpers (get_reactome_file(),
#' norm_name()) used by R/reactome_catalogue_build.R.
#'
#' The family is the Reactome PCD subtree walked from the pathway hierarchy --
#' NOT a name regex. Reactome places Autophagy OUTSIDE Programmed Cell Death, and
#' "Diseases of Programmed Cell Death", the TP53-transcription terms and the
#' SARS-CoV terms sit on other branches; a regex on "apoptosis|cell death|..."
#' would wrongly pull all of those in.
#'
#' Defines only functions and constants -- sourcing runs no analysis. Requires
#' `dplyr`/`tibble` and a `gmt_terms` vector from the caller (in practice the
#' term names of the cached Enrichr Reactome_Pathways_2024 GMT).
################################################################################

PCD_ROOT       <- "R-HSA-5357801"                     # Reactome "Programmed Cell Death"
REACTOME_CACHE <- "genesets/reactome"
REACTOME_BASE  <- "https://reactome.org/download/current"

#' Fetch (and cache) one Reactome download file, recording provenance.
#'
#' reactome.org answers a bad path with 200 + an HTML page, so a non-error return is not proof of success --
#' validate the payload. The release the family was resolved against is figure
#' provenance, hence the md5 + date meta file.
get_reactome_file <- function(fname, cache_dir = REACTOME_CACHE, refresh = FALSE) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(cache_dir, fname)
  meta <- paste0(dest, ".meta.txt")
  url  <- file.path(REACTOME_BASE, fname)

  if (refresh || !file.exists(dest)) {
    message("Downloading Reactome file: ", fname, " ...")
    ok <- tryCatch({
      utils::download.file(url, dest, method = "libcurl", quiet = TRUE); TRUE
    }, error = function(e) { message("   download failed: ", e$message); FALSE })

    bad <- !ok || !file.exists(dest) || file.size(dest) < 10000
    if (!bad) {
      first <- readLines(dest, n = 1, warn = FALSE)
      bad <- length(first) == 0 || !grepl("\t", first) || grepl("^\\s*<", first)
    }
    if (bad) {
      if (file.exists(dest)) unlink(dest)
      stop("Could not fetch a usable Reactome file from ", url,
           "\n  (needs network; delete ", cache_dir, " and retry, or place the file there by hand)")
    }
    writeLines(c(
      paste("file      :", fname),
      paste("url       :", url),
      paste("downloaded:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste("md5       :", unname(tools::md5sum(dest))),
      paste("n_rows    :", length(readLines(dest, warn = FALSE)))
    ), meta)
  }
  dest
}

#' Normalise a pathway name for cross-source matching.
#'
#' Enrichr title-cases Reactome names and strips punctuation, e.g. Reactome's
#' "SMAC(DIABLO)-mediated dissociation of IAP:caspase complexes " becomes
#' "SMAC(DIABLO)-mediated Dissociation of IAP Caspase Complexes". Dropping case and
#' every non-alphanumeric character makes both collapse to the same key, and also
#' absorbs Reactome's stray trailing/double spaces.
norm_name <- function(x) gsub("[^a-z0-9]", "", tolower(as.character(x)))

#' Resolve the Reactome Programmed Cell Death family and map it onto a gene-set library.
#'
#' @param gmt_terms Character vector of term names from the Enrichr Reactome GMT.
#' @param root Reactome stable ID to walk from.
#' @param refresh TRUE to re-download the hierarchy files.
#' @return list(family = <tibble: reactome_id, reactome_name, enrichr_term>,
#'              unmatched = <tibble>, n_descendants = <int>)
get_reactome_pcd_family <- function(gmt_terms, root = PCD_ROOT, refresh = FALSE) {

  # quote="" / comment.char="": pathway names contain apostrophes and '#'.
  pth <- utils::read.delim(get_reactome_file("ReactomePathways.txt", refresh = refresh),
                           header = FALSE, quote = "", comment.char = "",
                           col.names = c("stId", "name", "species"),
                           stringsAsFactors = FALSE)
  rel <- utils::read.delim(get_reactome_file("ReactomePathwaysRelation.txt", refresh = refresh),
                           header = FALSE, quote = "", comment.char = "",
                           col.names = c("parent", "child"), stringsAsFactors = FALSE)

  pth <- pth[pth$species == "Homo sapiens", , drop = FALSE]
  rel <- rel[startsWith(rel$parent, "R-HSA-"), , drop = FALSE]

  # A release that retires or re-IDs the root must fail loudly, not silently return a
  # one-term "family".
  if (!root %in% pth$stId)
    stop("Reactome root '", root, "' not found among Homo sapiens pathways. ",
         "The release may have re-IDed it -- check ",
         file.path(REACTOME_CACHE, "ReactomePathways.txt"), " before changing PCD_ROOT.")

  kids <- split(rel$child, rel$parent)
  ids <- character(0); queue <- root
  while (length(queue)) {
    node  <- queue[1]; queue <- queue[-1]
    if (node %in% ids) next                      # the hierarchy is a DAG, not a tree
    ids   <- c(ids, node)
    queue <- c(queue, kids[[node]])
  }

  nm  <- setNames(pth$name, pth$stId)
  fam <- tibble::tibble(reactome_id = ids, reactome_name = unname(nm[ids]))
  fam <- fam[!is.na(fam$reactome_name), , drop = FALSE]

  key <- setNames(gmt_terms, norm_name(gmt_terms))
  fam$enrichr_term <- unname(key[norm_name(fam$reactome_name)])

  matched   <- fam[!is.na(fam$enrichr_term), , drop = FALSE]
  unmatched <- fam[is.na(fam$enrichr_term), , drop = FALSE]

  # A silent name-schema change on either side would shrink the family to a handful of
  # sets and quietly turn this into a different analysis.
  if (nrow(matched) < 15)
    stop("Only ", nrow(matched), " of ", nrow(fam), " Reactome PCD pathways matched the ",
         "gene-set library -- expected ~34. The name convention on one side has probably ",
         "changed; inspect the match before trusting any result.")

  list(family = matched, unmatched = unmatched, n_descendants = nrow(fam))
}

#' GMT terms carrying a cell-death keyword that the hierarchy walk did NOT pick up.
#' Everything here must be genuinely outside the PCD subtree (disease branch, TP53
#' transcription branch, signalling branches, autophagy); anything that looks like core
#' PCD means the walk is wrong. Printed into every stats log as an audit.
pcd_keyword_misses <- function(gmt_terms, family_terms) {
  grep("apopto|necro|pyropto|ferropto|caspase|cell death",
       setdiff(gmt_terms, family_terms), ignore.case = TRUE, value = TRUE)
}

#' The non-CNS / developmental term blacklist from format_res_table_enrichr(), carried
#' here for REPORTING ONLY and not applied to the family: it is built to prune a
#' genome-wide GO scan, and on a hand-specified family it can only delete pathways that
#' were deliberately chosen. In this family it would remove "Caspase-mediated Cleavage
#' of Cytoskeletal Proteins", because "skeletal" matches inside "Cytoskeletal".
PCD_TEXT_OUT <- paste(
  "cardiac|lens|Mesenchymal|myotube|glioma|melanoma|thyroid|glomerular|renal|retina|nigra|vessel|estrogen|steroid|androgen|artery|bone|skeletal|muscle|aorta|cartilage|pancreatic|myoblast|embryonic|amyotrophic|neural tube|circadian|ectoderm|stem cell|vitamin|chylomicron|coronary|osteoclast|addiction|tumor|myometrial|prolactin|glioblastoma|sensory|cancer|carcinoma|hepatitis|oocyte|cardiomyocyte|heart|cardiac|eye|kidney|ear|auditory|Allograft|lupus|graft|rett|nose",
  "endocardial|genitalia|mammary|pharyngeal|mesenchym|chondrocyte|ossification|biomineral|angiogenesis|vasoconstriction|bilateral symmetry|left/right|left / right|pattern formation|epithelial|tight junction",
  sep = "|")
