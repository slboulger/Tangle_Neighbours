# ---------------------------------------------------------------------------
# cortical_depth_utils.R
#
# Figure panels: S3B - helper
# Sourced by deg_dream_linear_distance_depth_phf1.r (the depth-adjusted distance DEG).
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
## cortical_depth_utils.R
##
## Relative cortical depth for the CosMx cohort, as a covariate for the
## depth-adjusted linear-distance DEG (deg_dream_linear_distance_depth_phf1.r).
##
## WHAT THIS IS
## ------------
## Each sample's cortex sits at a different orientation in slide space, so depth
## cannot be read off x or y directly. A manual per-sample compass direction for
## L1 (the pial surface) is stored in `l1_dir`, and each cell is projected onto
## that sample's axis. `l1_dir` / `depth_one()` / `add_cortical_depth()` below are
## identical to the copies in nn3_utils.r.
##
## WHY THE BBOX VARIANT EXISTS
## ---------------------------
## `add_cortical_depth()` normalises by min/max over *the cells it is handed*.
## The DEG scripts run per celltype, per exclusion radius, on PHF1-negative cells
## only -- handing that subset to `add_cortical_depth()` would silently give
## `depth_rel` a different meaning in every run. `cortical_depth_bbox()` freezes
## the per-sample bounding box once, over the FULL all-cell frame, and
## `add_cortical_depth_from_bbox()` applies it to any subset, so depth_rel is
## comparable across runs.
##
## ORIENTATION CHECK
## -----------------
## depth_rel is position along the manual pia -> white-matter compass axis,
## scaled 0 (pia) to 1 within each sample. `validate_depth_orientation()`
## computes the per-sample and pooled Spearman correlation between the laminar
## excitatory ordinal (L2-3 -> L3-5 -> L5 -> L6-IT -> L6-CT) and depth_rel, and
## `cat_depth_orientation()` prints it into the stats log of every run that uses
## depth.
##
## Self-test:  Rscript phf1_v2/R/cortical_depth_utils.R --selftest
## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## Compass table and projection -- identical to nn3_utils.r
## ---------------------------------------------------------------------------
l1_dir <- c(
  "BBN_24895"     = "North", "BBN00628710" = "North", "BBN00635914" = "South", "BBN00638047" = "North",
  "BBN_10231" = "South", "BBN00428931" = "South", "BBN_9889" = "North",
  "BBN_9928" = "East",  "BBN_9931" = "South"
)

depth_one <- function(x, y, dir) {
  switch(dir,
    "North" = y - min(y),      # L1 at min y
    "South" = max(y) - y,      # L1 at max y
    "East"  = max(x) - x,      # L1 at max x
    "West"  = x - min(x),
    stop(sprintf("Unknown L1 direction '%s'", dir)))
}

## Adds depth_mm + depth_rel (0 = pial, 1 = deepest) to a frame, per sample.
## NOTE: depth is normalised over ALL cells given, so pass the full all-cell
## frame, not a neuron subset -- otherwise max(depth) differs between callers.
add_cortical_depth <- function(cd, sample_col = "sample_id",
                               coord_x = "x_slide_mm", coord_y = "y_slide_mm") {
  missing <- setdiff(unique(as.character(cd[[sample_col]])), names(l1_dir))
  if (length(missing))
    stop(sprintf("No L1 orientation for sample(s): %s", paste(missing, collapse = ", ")))
  cd$depth_mm <- NA_real_; cd$depth_rel <- NA_real_
  for (smp in unique(as.character(cd[[sample_col]]))) {
    i <- which(as.character(cd[[sample_col]]) == smp)
    d <- depth_one(cd[[coord_x]][i], cd[[coord_y]][i], l1_dir[[smp]])
    cd$depth_mm[i]  <- d
    cd$depth_rel[i] <- d / max(d)
  }
  cd
}

## ---------------------------------------------------------------------------
## Frozen all-cell bounding box
## ---------------------------------------------------------------------------

## Default all-cell coordinate export. 167,176 rows covering exactly the 9
## samples in `l1_dir` (the excluded white-matter section BBN_2932 is absent).
## Written by R/export_seu_coords_for_phf1.R from the Seurat meta.data; carries
## `celltype`, so the orientation check runs off the same small CSV without
## loading the Seurat object.
##
## Paths are derived from THIS FILE's own location, never hardcoded. This file
## sits at <tree>/R/, the export at <tree>/PHF1/, so the same source resolves
## correctly from either <PROJECT_ROOT> or <PROJECT_ROOT>/phf1_v2.
.depth_tree_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  m  <- grep("^--file=", ca, value = TRUE)
  f  <- if (length(m)) sub("^--file=", "", m[1]) else ""
  # `Rscript -` (stdin) reports --file=- , which is not a path. Only trust a
  # --file= that actually exists; otherwise fall back to the working directory.
  d <- if (nzchar(f) && f != "-" && file.exists(f)) {
    dirname(normalizePath(f))
  } else {
    getwd()
  }
  # d is <tree>/R when running from the script directory; when it is already the
  # tree root (interactive session at the top of a tree), keep it.
  if (basename(d) == "R") dirname(d) else d
}
DEPTH_COORDS_DEFAULT <- file.path(.depth_tree_dir(), "PHF1", "seu_coords.csv")
DEPTH_BBOX_DEFAULT   <- file.path(.depth_tree_dir(), "cortical_depth_bbox.tsv")

## Per-sample min/max of the slide coordinates, plus the compass direction and
## the resulting depth span. One row per sample.
cortical_depth_bbox <- function(cd, sample_col = "sample_id",
                                coord_x = "x_slide_mm", coord_y = "y_slide_mm") {
  smps <- sort(unique(as.character(cd[[sample_col]])))
  missing <- setdiff(smps, names(l1_dir))
  if (length(missing))
    stop(sprintf("No L1 orientation for sample(s): %s", paste(missing, collapse = ", ")))
  x <- as.numeric(cd[[coord_x]]); y <- as.numeric(cd[[coord_y]])
  if (anyNA(x) || anyNA(y))
    stop("NA in coordinates; the bbox must come from a complete all-cell frame.")
  out <- do.call(rbind, lapply(smps, function(s) {
    i <- which(as.character(cd[[sample_col]]) == s)
    d <- depth_one(x[i], y[i], l1_dir[[s]])
    data.frame(sample_id = s, l1_dir = unname(l1_dir[[s]]),
               x_min = min(x[i]), x_max = max(x[i]),
               y_min = min(y[i]), y_max = max(y[i]),
               depth_max_mm = max(d), n_cells = length(i),
               stringsAsFactors = FALSE)
  }))
  if (any(out$depth_max_mm <= 0))
    stop("A sample has zero depth span; its cells are collinear along the compass axis.")
  rownames(out) <- NULL
  out
}

## Build (or read) the frozen bbox. Caches to `cache` so every consumer -- every
## celltype, every radius -- normalises against byte-identical numbers.
load_cortical_depth_bbox <- function(coords_csv = DEPTH_COORDS_DEFAULT,
                                     cache      = DEPTH_BBOX_DEFAULT,
                                     refresh    = FALSE) {
  if (!refresh && file.exists(cache))
    return(utils::read.delim(cache, stringsAsFactors = FALSE))
  if (!file.exists(coords_csv))
    stop("All-cell coordinate export not found: ", coords_csv)
  cd   <- utils::read.csv(coords_csv, stringsAsFactors = FALSE)
  bbox <- cortical_depth_bbox(cd)
  utils::write.table(bbox, cache, sep = "\t", quote = FALSE, row.names = FALSE)
  bbox
}

## Apply a frozen bbox to any subset of cells. depth_rel is clamped to [0,1]
## only against floating-point slop; a value genuinely outside the box means the
## subset contains cells the bbox was not built from, which is an error.
add_cortical_depth_from_bbox <- function(cd, bbox, sample_col = "sample_id",
                                         coord_x = "x_slide_mm",
                                         coord_y = "y_slide_mm",
                                         tol = 1e-9) {
  smps <- unique(as.character(cd[[sample_col]]))
  missing <- setdiff(smps, bbox$sample_id)
  if (length(missing))
    stop(sprintf("Sample(s) absent from the frozen depth bbox: %s",
                 paste(missing, collapse = ", ")))
  x <- as.numeric(cd[[coord_x]]); y <- as.numeric(cd[[coord_y]])
  cd$depth_mm <- NA_real_; cd$depth_rel <- NA_real_
  for (s in smps) {
    i <- which(as.character(cd[[sample_col]]) == s)
    b <- bbox[match(s, bbox$sample_id), ]
    d <- switch(b$l1_dir,
      "North" = y[i] - b$y_min,
      "South" = b$y_max - y[i],
      "East"  = b$x_max - x[i],
      "West"  = x[i] - b$x_min,
      stop(sprintf("Unknown L1 direction '%s'", b$l1_dir)))
    cd$depth_mm[i]  <- d
    cd$depth_rel[i] <- d / b$depth_max_mm
  }
  if (anyNA(cd$depth_rel)) stop("NA depth_rel after applying the frozen bbox.")
  rng <- range(cd$depth_rel)
  if (rng[1] < -tol || rng[2] > 1 + tol)
    stop(sprintf("depth_rel outside [0,1] (%.6f, %.6f): these cells are not the ones the bbox was built from.",
                 rng[1], rng[2]))
  cd$depth_rel <- pmin(1, pmax(0, cd$depth_rel))
  cd
}

## ---------------------------------------------------------------------------
## Orientation validation
## ---------------------------------------------------------------------------

## Laminar excitatory subtypes, shallow -> deep. Used ONLY to test the axis.
DEPTH_LAMINAR_ORDER <- c(
  "Exc-IT-L2-3-CBLN2-HOPX",
  "Exc-IT-L3-5-CHGA-IL1RAPL2",
  "Exc-ET-L5-SPON1-FGD4",
  "Exc-IT-L6-CTXN1-ERC2",
  "Exc-CT-L6-SYNPO2-SEMA3E"
)
DEPTH_RHO_WARN <- 0.2   # pooled-rho threshold for the note printed in the stats log

## Per-sample and pooled Spearman rho between the laminar ordinal and depth_rel,
## plus the mean depth_rel of each subtype. Returns a list; `warn` is TRUE when
## the pooled rho is below DEPTH_RHO_WARN.
validate_depth_orientation <- function(cd, sample_col = "sample_id",
                                       celltype_col = "celltype") {
  if (!celltype_col %in% colnames(cd) || !"depth_rel" %in% colnames(cd))
    return(list(per_sample = NULL, pooled_rho = NA_real_, warn = NA,
                note = "celltype or depth_rel absent; orientation not validated"))
  ct  <- as.character(cd[[celltype_col]])
  keep <- ct %in% DEPTH_LAMINAR_ORDER
  if (sum(keep) < 50)
    return(list(per_sample = NULL, pooled_rho = NA_real_, warn = NA,
                note = "fewer than 50 laminar Exc cells; orientation not validated"))
  ord <- match(ct[keep], DEPTH_LAMINAR_ORDER)
  dep <- cd$depth_rel[keep]
  smp <- as.character(cd[[sample_col]])[keep]

  rows <- lapply(sort(unique(smp)), function(s) {
    i <- which(smp == s)
    r <- if (length(unique(ord[i])) > 1)
           suppressWarnings(stats::cor(ord[i], dep[i], method = "spearman")) else NA_real_
    means <- vapply(seq_along(DEPTH_LAMINAR_ORDER), function(k) {
      j <- i[ord[i] == k]; if (length(j)) mean(dep[j]) else NA_real_
    }, numeric(1))
    out <- data.frame(sample_id = s, l1_dir = unname(l1_dir[[s]]),
                      n = length(i), rho = r, stringsAsFactors = FALSE)
    for (k in seq_along(DEPTH_LAMINAR_ORDER))
      out[[paste0("mean_", DEPTH_LAMINAR_ORDER[k])]] <- means[k]
    out
  })
  per_sample <- do.call(rbind, rows)
  pooled <- mean(per_sample$rho, na.rm = TRUE)
  list(per_sample = per_sample, pooled_rho = pooled,
       warn = isTRUE(pooled < DEPTH_RHO_WARN), note = NA_character_)
}

## Cohort-level validation, computed on the ALL-CELL frame and cached.
##
## This exists because the per-celltype callers CANNOT validate the axis
## themselves: a *_sce_neighbours.qs holds exactly one celltype, so the laminar
## ordinal has a single level and every rho would come back NA. The axis is a
## property of the cohort, not of the cell set being modelled, so it is measured
## once on the same all-cell frame the bbox comes from and cached beside it.
cortical_depth_validation <- function(coords_csv = DEPTH_COORDS_DEFAULT,
                                      bbox       = NULL,
                                      cache      = NULL,
                                      refresh    = FALSE) {
  if (is.null(cache))
    cache <- file.path(dirname(DEPTH_BBOX_DEFAULT), "cortical_depth_validation.tsv")
  if (!refresh && file.exists(cache)) {
    ps <- utils::read.delim(cache, stringsAsFactors = FALSE, check.names = FALSE)
    pooled <- mean(ps$rho, na.rm = TRUE)
    return(list(per_sample = ps, pooled_rho = pooled,
                warn = isTRUE(pooled < DEPTH_RHO_WARN), note = NA_character_))
  }
  if (!file.exists(coords_csv))
    stop("All-cell coordinate export not found, cannot validate the depth axis: ",
         coords_csv)
  cd <- utils::read.csv(coords_csv, stringsAsFactors = FALSE)
  if (is.null(bbox)) bbox <- cortical_depth_bbox(cd)
  cd <- add_cortical_depth_from_bbox(cd, bbox)
  v  <- validate_depth_orientation(cd)
  if (!is.null(v$per_sample))
    utils::write.table(v$per_sample, cache, sep = "\t", quote = FALSE, row.names = FALSE)
  v
}

## Print the validation into an open sink(). Reports rather than stops: depth_rel
## is used as a position covariate whatever the rho.
cat_depth_orientation <- function(v) {
  cat("\n--- Cortical-depth orientation validation ---\n")
  if (!is.na(v$note)) { cat(" ", v$note, "\n"); return(invisible(NULL)) }
  cat("Spearman rho between laminar Exc ordinal (L2-3 -> L6-CT) and depth_rel.\n\n")
  print(v$per_sample, row.names = FALSE, digits = 3)
  cat(sprintf("\nPooled mean rho: %+.3f\n", v$pooled_rho))
  if (isTRUE(v$warn)) {
    cat("\n  Note: pooled rho is below DEPTH_RHO_WARN; depth_rel is used as position along\n")
    cat("  the manual pia->WM compass axis. See the header of cortical_depth_utils.R.\n")
  }
  invisible(NULL)
}

## ---------------------------------------------------------------------------
## Self-test
## ---------------------------------------------------------------------------
if (sys.nframe() == 0L && "--selftest" %in% commandArgs(TRUE)) {
  # Tree root is derived from this file's location, so the self-test runs in
  # whichever tree it was invoked from with no hardcoded paths.
  setwd(.depth_tree_dir())

  cat("=== cortical_depth_utils.R self-test ===\n")
  cat("Tree: ", getwd(), "\n\n")
  cd <- utils::read.csv(DEPTH_COORDS_DEFAULT, stringsAsFactors = FALSE)
  cat(sprintf("All-cell frame: %d cells, %d samples\n",
              nrow(cd), length(unique(cd$sample_id))))
  stopifnot(setequal(unique(cd$sample_id), names(l1_dir)))
  cat("Sample set matches l1_dir exactly.\n\n")

  bbox <- load_cortical_depth_bbox(refresh = TRUE)
  print(bbox, row.names = FALSE, digits = 6)

  # 1. bbox path must reproduce the in-frame path exactly on the full frame
  a <- add_cortical_depth(cd)
  b <- add_cortical_depth_from_bbox(cd, bbox)
  d1 <- max(abs(a$depth_rel - b$depth_rel))
  cat(sprintf("\n[1] bbox vs in-frame depth_rel, full frame: max abs diff %.3e\n", d1))
  stopifnot(d1 < 1e-12)

  # 2. the point of the bbox: a subset must keep the same depth_rel
  set.seed(1)
  sub <- cd[sample(nrow(cd), 5000), ]
  s_b <- add_cortical_depth_from_bbox(sub, bbox)
  s_a <- add_cortical_depth(sub)                     # renormalised over the subset
  ref <- b$depth_rel[match(s_b$cell_id, b$cell_id)]
  cat(sprintf("[2] subset via bbox      vs full-frame: max abs diff %.3e\n",
              max(abs(s_b$depth_rel - ref))))
  cat(sprintf("[2] subset via add_cortical_depth() vs full-frame: max abs diff %.3e\n",
              max(abs(s_a$depth_rel - ref))))
  stopifnot(max(abs(s_b$depth_rel - ref)) < 1e-12)

  # 3. round-trip against the cached depth_rel from the nn3 analysis. That cache
  #    may sit in <PROJECT_ROOT> rather than phf1_v2, so fall back to it before skipping.
  nn3 <- "plots/nn3_neuron_spacing/source_data_nn3_cells.tsv"
  if (!file.exists(nn3)) {
    alt <- file.path(dirname(getwd()), nn3)
    if (file.exists(alt)) nn3 <- alt
  }
  if (file.exists(nn3)) {
    t3 <- utils::read.delim(nn3, stringsAsFactors = FALSE)
    if (all(c("cell_id", "depth_rel") %in% colnames(t3))) {
      m <- match(t3$cell_id, b$cell_id); ok <- !is.na(m)
      d3 <- max(abs(t3$depth_rel[ok] - b$depth_rel[m[ok]]))
      cat(sprintf("[3] vs cached nn3 depth_rel (%d shared cells): max abs diff %.3e\n",
                  sum(ok), d3))
      stopifnot(d3 < 1e-6)
    } else cat("[3] nn3 cache lacks depth_rel; skipped\n")
  } else cat("[3] nn3 cache absent; skipped\n")

  # 4. orientation validation
  cat_depth_orientation(validate_depth_orientation(b))

  cat("\n=== self-test passed ===\n")
}
