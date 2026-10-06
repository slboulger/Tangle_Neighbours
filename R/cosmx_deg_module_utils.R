#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# cosmx_deg_module_utils.R
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
# cosmx_deg_module_utils.R
#
# Cuts CosMx module gene sets out of the two CANONICAL DEG tables under one
# symmetric rule, so that a "cell-autonomous" module (PHF1+ vs PHF1-) and a
# "field" module (near vs far from a tangle) are built the same way, are the
# same size, and are therefore directly comparable when scored side by side.
#
# Pure functions only - no argparse, no setwd(), no file writing at load time.
# Sourced by:
#   R/build_cosmx_deg_modules.R   the CLI that writes the sets to disk
#
# ---------------------------------------------------------------------------
# THE TWO SOURCE TABLES, AND WHY THEY ARE NOT INTERCHANGEABLE
#
#   deg/de_cat/de_PHF1/<ct>/<ct>_PHF1TRUEVsPHF1FALSE.tsv
#     Set 1p pseudobulk limma-trend, contrast "PHF1TRUE - PHF1FALSE".
#     POSITIVE t = UP IN TANGLE-BEARING CELLS.  4464 genes tested for CBLN2
#     (a pseudobulk CPM-sum filter).
#
#   deg/de_linear_distance/<ct>/<ct>_dist_to_phf1_um_scaled.tsv
#     Set 1 cell-level dream, coefficient dist_to_phf1_um_scaled =
#     log(dist)/sd(log(dist)), fitted on PHF1-NEGATIVE cells only.
#     NEGATIVE t = expression FALLS with distance = UP NEAR TANGLES.  2468 genes
#     tested for CBLN2 (a 5% cell-level detection filter).
#
#   The two universes are DIFFERENT (4464 vs 2468) because the gene filters are
#   different, and the two contrasts are on different cells. They are aligned
#   here only in direction and in size, never pooled.
#
# SIGN TRAP, stated once so it cannot be read past. The project's own enrichR
#   output under de_linear_distance/<ct>/enrichr/DOWN/ IS the "up near tangles"
#   direction - those folders are named by the sign of logFC, not by proximity.
#   Same genes, opposite word. Confirmed against the header of
#   R/nDEG_by_celltype_linear.R.
#
# ---------------------------------------------------------------------------
# THE SELECTION RULE
#
#   1. padj < PADJ_CUT (0.1, the project DEG threshold).
#   2. Direction taken relative to the PANEL-WIDE MEDIAN t, not to zero:
#        cat :  s =  (t - median(t))     positive = up in tangle-bearing
#        dist:  s = -(t - median(t))     positive = up near tangles
#   3. Rank by s descending, take the top N, with N common to both sets.
#
# WHAT THE RECENTRING DOES AND DOES NOT DO. Measured on CBLN2: the
#   median t is +1.5592 (cat) and +1.7287 (linear), i.e. both tables carry a
#   large panel-wide POSITIVE shift. Recentring is nonetheless a NO-OP FOR SET
#   MEMBERSHIP at padj < 0.1 - the split is 1631/263 in cat and 1057/53 in
#   linear whether or not the median is subtracted, because every gene that
#   survives FDR has |t| far above the median. Subtracting a constant also
#   cannot reorder a one-sided ranking. So the recentring is retained as
#     (a) the DIRECTION rule, and
#     (b) the reported MAGNITUDE, which is what makes the two sets' effect
#         floors comparable across two tables with different medians.
#   It is NOT a re-test. diagnose_recentring() emits the counts that prove the
#   no-op, and build_deg_modules() puts them in the provenance file, so nobody
#   has to take this comment on trust.
#
# WHY THE POSITIVE MEDIAN MATTERS SCIENTIFICALLY. In the linear model - fitted
#   on PHF1-negative cells within 1000 um - the TYPICAL panel gene rises with
#   distance. The genes selected as "up near tangles" are a minority tail
#   running against that panel-wide shift.
#
# SIZE MATCHING. padj < 0.1 yields 1631 eligible cat genes but only 53 eligible
#   distance genes, so the common N is capped at 53 by the distance side. Above
#   N = 53 the distance set is RANK-EXTENDED past FDR: n_sig is carried on every
#   row so a rank-extended set can never be mistaken for a significant one.

##  ............................................................................
##  Paths                                                                   ####

## Celltype -> filesystem-safe tag. Matches the idiom in phf1_markers/ and in
## read_geneset() in R/plot_phf1_module_vs_phf1_distance_modelp.R. Identity for every
## current celltype name (they are alphanumeric plus "-").
ct_tag <- function(x) gsub("[^A-Za-z0-9_-]", "_", x)

## The DEG trees use the RAW celltype name for both the directory and the file
## stem - not the safe tag. Every consumer in this project re-implements these
## two one-liners; they are centralised here.
cat_deg_path <- function(cat_dir, celltype)
  file.path(cat_dir, celltype, paste0(celltype, "_PHF1TRUEVsPHF1FALSE.tsv"))

dist_deg_path <- function(dist_dir, celltype)
  file.path(dist_dir, celltype, paste0(celltype, "_dist_to_phf1_um_scaled.tsv"))

module_dir <- function(out_dir, celltype) file.path(out_dir, ct_tag(celltype))

module_geneset_path <- function(out_dir, celltype, set, n)
  file.path(module_dir(out_dir, celltype),
            sprintf("%s_%s_N%d_geneset.txt", ct_tag(celltype), set, n))

seg_draws_path <- function(out_dir, celltype)
  file.path(module_dir(out_dir, celltype),
            sprintf("%s_seg_draws.tsv", ct_tag(celltype)))

membership_path <- function(out_dir, celltype)
  file.path(module_dir(out_dir, celltype),
            sprintf("%s_modules_membership.tsv", ct_tag(celltype)))

provenance_path <- function(out_dir, celltype)
  file.path(module_dir(out_dir, celltype),
            sprintf("%s_modules_provenance.txt", ct_tag(celltype)))

##  ............................................................................
##  Probe names                                                             ####

## CosMx probe names can pack several genes behind one probe ("MZT2A/B"). The
## snRNA reference carries the individual symbols, so split on "/" and re-attach
## the common prefix to bare suffixes. Without the split, a packed probe silently
## fails to match the snRNA reference (a naming mismatch, not a missing gene).
expand_probe_names <- function(g) {
  parts <- strsplit(g, "/", fixed = TRUE)
  out <- unlist(lapply(parts, function(p) {
    if (length(p) == 1) return(p)
    stem <- p[1]
    c(stem, vapply(p[-1], function(s) {
      if (nchar(s) >= nchar(stem)) s
      else paste0(substr(stem, 1, nchar(stem) - nchar(s)), s)
    }, character(1)))
  }), use.names = FALSE)
  unique(out)
}

##  ............................................................................
##  Readers                                                                 ####

## Validated DEG reader. stop()s on a missing column rather than returning a
## silently truncated set - `t` in particular is easy to assume and is what the
## whole recentring rests on.
read_deg_table <- function(path, what) {
  if (!file.exists(path)) stop(what, " DEG table not found: ", path, call. = FALSE)
  d <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  need <- c("gene", "logFC", "padj", "t")
  miss <- setdiff(need, colnames(d))
  if (length(miss))
    stop(what, " DEG table missing column(s): ", paste(miss, collapse = ", "),
         "\n  file: ", path,
         "\n  has:  ", paste(colnames(d), collapse = ", "), call. = FALSE)
  d$gene <- as.character(d$gene)
  d <- d[!is.na(d$t) & nzchar(d$gene), , drop = FALSE]
  if (nrow(d) == 0) stop(what, " DEG table has no usable rows: ", path, call. = FALSE)
  if (anyDuplicated(d$gene))
    stop(what, " DEG table has duplicate gene rows: ", path, call. = FALSE)
  d
}

## Recentre on the panel-wide median t and orient so that POSITIVE s is always
## the biologically wanted direction.
##   flip = FALSE  cat      positive t  = up in tangle-bearing  -> s =  (t - med)
##   flip = TRUE   distance negative t  = up near tangles       -> s = -(t - med)
recentre_deg <- function(d, flip) {
  med <- stats::median(d$t)
  d$t_median_panel <- med
  d$s <- if (flip) -(d$t - med) else (d$t - med)
  d
}

## The proof that the recentring does not move the padj-significant split.
## Returned as a one-row data.frame and printed into the provenance file.
diagnose_recentring <- function(d, flip, padj_cut, label) {
  sig <- !is.na(d$padj) & d$padj < padj_cut
  raw <- if (flip) -d$t else d$t          # raw-sign orientation, no median
  data.frame(
    table            = label,
    n_tested         = nrow(d),
    median_t         = unname(d$t_median_panel[1]),
    mean_t           = mean(d$t),
    n_padj_sig       = sum(sig),
    n_sig_raw_pos    = sum(sig & raw > 0),
    n_sig_raw_neg    = sum(sig & raw < 0),
    n_sig_recentred_pos = sum(sig & d$s > 0),
    n_sig_recentred_neg = sum(sig & d$s < 0),
    recentring_is_noop  = sum(sig & raw > 0) == sum(sig & d$s > 0) &&
                          sum(sig & raw < 0) == sum(sig & d$s < 0),
    stringsAsFactors = FALSE
  )
}

##  ............................................................................
##  Module construction                                                     ####

## Eligible genes in the wanted direction, most-extreme recentred t first.
## padj is NOT applied here - it is applied as a top-N cap by the caller, so
## that a rank-extended set beyond FDR is possible but always labelled.
eligible_by_s <- function(d) {
  e <- d[d$s > 0, , drop = FALSE]
  e[order(-e$s), , drop = FALSE]
}

#' Build the module gene sets across a size sweep.
#'
#' @param cat_df,dist_df validated + recentred DEG tables
#' @param padj_cut       FDR threshold used for `phf1_full` and for n_sig
#' @param n_sweep        integer vector of set sizes
#' @return list(sets, meta, membership, phf1_full)
#'   sets       list keyed "<set>_N<n>" -> character vector of PROBES
#'   meta       one row per (set, N): n_genes_probe, n_sig, all_sig, s_floor, ...
#'              (the caller adds n_genes once probes are expanded to symbols)
#'   membership long data.frame, one row per probe per set per N
build_deg_modules <- function(cat_df, dist_df, padj_cut, n_sweep) {
  stopifnot(is.numeric(n_sweep), length(n_sweep) > 0, all(n_sweep > 0))
  n_sweep <- sort(unique(as.integer(n_sweep)))

  elig_cat  <- eligible_by_s(cat_df)
  elig_dist <- eligible_by_s(dist_df)

  ## The FULL padj-significant PHF1 set. This is the subtrahend for
  ## `field_strict`, and it is deliberately the full 1631-gene list rather than
  ## the size-matched one: subtracting only the top N would leave genes that ARE
  ## significantly up in tangle-bearing cells in a set called "field only".
  phf1_full <- elig_cat$gene[!is.na(elig_cat$padj) & elig_cat$padj < padj_cut]

  sets <- list(); meta <- list(); memb <- list()

  add <- function(key, n, genes, src, note) {
    sets[[sprintf("%s_N%d", key, n)]] <<- genes
    m <- src[match(genes, src$gene), , drop = FALSE]
    n_sig <- sum(!is.na(m$padj) & m$padj < padj_cut)
    meta[[length(meta) + 1L]] <<- data.frame(
      set = key, N = n, n_genes_probe = length(genes), n_sig = n_sig,
      all_sig = n_sig == length(genes),
      s_floor = if (length(genes)) min(m$s, na.rm = TRUE) else NA_real_,
      s_top   = if (length(genes)) max(m$s, na.rm = TRUE) else NA_real_,
      note = note, stringsAsFactors = FALSE)
    if (length(genes))
      memb[[length(memb) + 1L]] <<- data.frame(
        set = key, N = n, gene = genes, t = m$t, s = m$s, padj = m$padj,
        logFC = m$logFC, in_phf1_full = genes %in% phf1_full,
        stringsAsFactors = FALSE)
    invisible(NULL)
  }

  for (n in n_sweep) {
    g_phf1 <- utils::head(elig_cat$gene,  n)
    g_dist <- utils::head(elig_dist$gene, n)

    add("phf1cat",  n, g_phf1, elig_cat,
        "up in tangle-bearing (categorical PHF1+ vs PHF1-)")
    add("distance", n, g_dist, elig_dist,
        "up near tangles (linear distance, t < 0)")
    ## Derived sets. NOT size-matched by construction - the size-matching claim
    ## holds for phf1cat and distance only.
    add("field_strict",  n, setdiff(g_dist, phf1_full), elig_dist,
        "distance minus the FULL padj-significant PHF1 set")
    add("field_matched", n, setdiff(g_dist, g_phf1), elig_dist,
        "distance minus the size-matched PHF1 set")
  }

  list(sets       = sets,
       meta       = do.call(rbind, meta),
       membership = do.call(rbind, memb),
       phf1_full  = phf1_full,
       n_sweep    = n_sweep)
}

#' Size-matched SEG null draws.
#'
#' The CBLN2 SEG set is 305 genes against 53-gene signatures, and module-score
#' variance scales with set size, so the full set is not a fair null: it is
#' smoother than anything it is being compared with. Instead draw `n_draws`
#' random N-gene subsets. Scored and modelled, these give an EMPIRICAL NULL
#' DISTRIBUTION for the predictor coefficient rather than a single curve to
#' eyeball - which is the only way to say whether a signature does anything a
#' stably-expressed gene set of the same size does not also do.
#'
#' Seeded per (N, draw) so a draw is reproducible independently of how many
#' other draws or sweep levels were requested.
seg_matched_draws <- function(seg_genes, n_sweep, n_draws, seed) {
  seg_genes <- unique(as.character(seg_genes))
  n_sweep <- sort(unique(as.integer(n_sweep)))
  out <- list()
  for (n in n_sweep) {
    if (n > length(seg_genes))
      stop("SEG set has ", length(seg_genes), " genes; cannot draw N = ", n,
           ". Reduce --n_sweep or use a larger SEG set.", call. = FALSE)
    for (k in seq_len(n_draws)) {
      set.seed(seed + n * 1000L + k)
      out[[length(out) + 1L]] <- data.frame(
        N = n, draw = k, gene = sort(sample(seg_genes, n)),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

##  ............................................................................
##  Consumer-side loaders                                                   ####

## Flat newline-delimited gene list, matching the phf1_markers/ convention so
## the existing read_geneset() loaders elsewhere work on these files unchanged.
read_module_geneset <- function(out_dir, celltype, set, n) {
  f <- module_geneset_path(out_dir, celltype, set, n)
  if (!file.exists(f))
    stop("module geneset not found: ", f,
         "\n  Run R/build_cosmx_deg_modules.R first.", call. = FALSE)
  g <- trimws(readLines(f))
  unique(g[nzchar(g)])
}

## Returns a named list "seg_draw001" -> genes, for one sweep level.
read_seg_draws <- function(out_dir, celltype, n) {
  f <- seg_draws_path(out_dir, celltype)
  if (!file.exists(f))
    stop("SEG draw table not found: ", f,
         "\n  Run R/build_cosmx_deg_modules.R first.", call. = FALSE)
  d <- utils::read.delim(f, stringsAsFactors = FALSE)
  d <- d[d$N == n, , drop = FALSE]
  if (nrow(d) == 0)
    stop("SEG draw table has no rows for N = ", n, "; has N in {",
         paste(sort(unique(utils::read.delim(f)$N)), collapse = ", "), "}",
         call. = FALSE)
  sp <- split(d$gene, sprintf("seg_draw%03d", d$draw))
  lapply(sp, function(x) unique(as.character(x)))
}

## Which sweep levels are on disk, for a consumer that wants to validate --n_sweep.
available_module_sizes <- function(out_dir, celltype) {
  d <- module_dir(out_dir, celltype)
  if (!dir.exists(d)) return(integer(0))
  f <- list.files(d, pattern = "_distance_N[0-9]+_geneset\\.txt$")
  sort(as.integer(sub(".*_distance_N([0-9]+)_geneset\\.txt$", "\\1", f)))
}
