#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# nn3_utils.r
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
# nn3_utils.r
#
# Shared helper: single source of truth for the 3-NEAREST-NEURONAL-NEIGHBOUR
# SPACING metric (nn3_um) and its technical companion, the FOV edge-censoring
# distance (edge_dist_um).
#
# Deliberately mirrors the design of R/phf1_distance_utils.r: the metric lives in
# ONE place so that
#   * R/nn3_phf1_vs_neg.R        (Analysis A, cell-autonomous)
#   * R/nn3_over_phf1_distance.R (Analysis B, non-cell-autonomous)
#   * R/null_replica_nn3.R       (supplementary permutation null)
# compute it through one implementation and can never drift apart.
#
# ---------------------------------------------------------------------------
# WHY THIS METRIC
# ---------------------------------------------------------------------------
# Zwang et al. 2024 (Cell Rep., Hyman lab, with Bennett; PMC11441076) imaged
# rTg4510 cortex longitudinally and found tangle-FREE neurons died at >3x the
# rate of tangle-bearing neurons, and that before dying they became more distant
# from their neighbours -- quantified as the 3-nearest-neighbour distance among
# labelled neurons (33.9 +/- 1.9 um in dying vs 24.0 +/- 0.5 um in persisting
# neurons, p = 0.02, mixed-effects model). Neighbour spacing is therefore a
# LOCAL NEURON-LOSS readout, not merely a density descriptor.
#
# Differences from that study:
#   * 2D single section here, 3D z-stacks there. A 2D nn3 over-estimates true 3D
#     spacing because out-of-plane neighbours are invisible.
#   * Neurons are identified by transcriptomic celltype assignment, not a
#     pan-neuronal AAV label. Cells whose type could not be called
#     ("Unassigned", 40,681 cells) are NOT in the pool, so some true neurons are
#     missing from it and nn3_um is biased upward. The bias is approximately
#     uniform across PHF1 status, so it inflates the level, not the contrast.
#   * Cross-sectional, so "became more distant" cannot be observed; only the
#     standing spatial association can.
#
# ---------------------------------------------------------------------------
# CONVENTIONS
# ---------------------------------------------------------------------------
# Coordinates are x_slide_mm / y_slide_mm (slide-global centroids, mm). They are
# slide-global and NOT sample-global -- two samples can share a slide -- so every
# computation here groups by sample first, exactly as compute_dist_to_phf1_um()
# does in phf1_distance_utils.r. All returned distances are MICRONS (mm * 1000).

suppressPackageStartupMessages({
  library(cli)
})

if (!requireNamespace("RANN", quietly = TRUE)) {
  stop("Package 'RANN' is required. Install with: install.packages('RANN')")
}

## ---------------------------------------------------------------------------
## Neighbour pool definition
##
## The 10 neuron labels: neuron_order (9 confident subtypes, R/palettes.R)
## plus "Unassigned Neuron". This is IDENTICAL to the PHF1+ distance-source
## neuron set passed in run_label_phf1_neighbours.sh, so nn3_um and the
## canonical dist_to_phf1_um refer to the same underlying neuron population.
## Do not change one without changing the other.
##
## (Note the looser definition in rename_celltypes.R,
## grepl("Neuron|Exc|Inh", celltype), used only to build sce_neuron.qs. It is not
## the same set and must not be substituted here.)
## ---------------------------------------------------------------------------
NEURON_POOL_10 <- c(
  "Exc-IT-L2-3-CBLN2-HOPX",
  "Exc-IT-L3-5-CHGA-IL1RAPL2",
  "Exc-ET-L5-SPON1-FGD4",
  "Exc-IT-L6-CTXN1-ERC2",
  "Exc-CT-L6-SYNPO2-SEMA3E",
  "Inh-PVALB",
  "Inh-SST",
  "Inh-VIP",
  "Inh-LAMP5",
  "Unassigned Neuron"
)

## ---------------------------------------------------------------------------
## compute_nn3_um()
##
## Per-sample k-nearest-neuron spacing among the neighbour pool.
##
##   cd           : data.frame with cell_id, celltype, the coordinate columns and
##                  the sample column. May contain non-pool cells; they are
##                  dropped (they are neither queries nor neighbours).
##   sample_col   : name of the sample grouping column.
##   neuron_pool  : character vector of celltype labels forming the pool.
##   coord_x/y    : coordinate column names (mm; converted to um here).
##   k            : number of nearest neighbours (3 = the Zwang et al. metric).
##
## Both the QUERY set and the NEIGHBOUR set are the pool: spacing is a property
## of a neuron relative to other neurons, so it is defined for every pool cell
## and is INDEPENDENT OF PHF1 STATUS (PHF1+ neurons are in the pool, both as
## queries and as neighbours of others). That invariance is what makes the
## label-permutation null in null_replica_nn3.R nearly free -- nn3_um does not
## change when PHF1 labels are shuffled.
##
## Self-exclusion is done BY CELL INDEX, not by dropping column 1. A handful of
## cells share a centroid exactly with another pool cell (1 such cell in sample
## BBN_24895), so the distance-0 match in column 1 is not always the cell itself and
## dropping column 1 blindly would silently measure the wrong neighbours.
##
## Returns a data.frame: cell_id, nn3_um (MEAN distance to the k nearest, the
## headline metric), nn3rd_um (distance to the k-th nearest, the sensitivity
## metric), n_pool_sample, n_coincident (pool neurons at distance 0 from this
## one -- exact centroid duplicates, a segmentation artefact worth reporting).
## NA where the sample has <= k pool neurons.
## ---------------------------------------------------------------------------
compute_nn3_um <- function(cd, sample_col, neuron_pool = NEURON_POOL_10,
                           coord_x = "x_slide_mm", coord_y = "y_slide_mm",
                           k = 3L, verbose = TRUE) {

  stopifnot("cell_id" %in% colnames(cd), "celltype" %in% colnames(cd))
  stopifnot(k >= 1L)

  pool <- cd[as.character(cd$celltype) %in% neuron_pool, , drop = FALSE]
  if (nrow(pool) == 0) stop("No cells matched neuron_pool; check celltype labels.")

  out <- data.frame(
    cell_id       = pool$cell_id,
    nn3_um        = NA_real_,
    nn3rd_um      = NA_real_,
    n_pool_sample = NA_integer_,
    n_coincident  = NA_integer_,
    stringsAsFactors = FALSE
  )

  samples <- unique(pool[[sample_col]])
  for (smp in samples) {
    idx <- which(pool[[sample_col]] == smp)
    n_p <- length(idx)
    out$n_pool_sample[idx] <- n_p

    if (n_p <= k) {
      if (verbose)
        cat(sprintf("  %-10s pool neurons: %-6d  <= k = %d, nn3 left NA\n", smp, n_p, k))
      next
    }

    coords <- as.matrix(pool[idx, c(coord_x, coord_y)])
    # k + 2 gives slack so the self-match can be dropped by index even when a
    # coincident cell also sits at distance 0.
    kk <- min(k + 2L, n_p)
    nn <- RANN::nn2(data = coords, query = coords, k = kk)

    self     <- nn$nn.idx == seq_len(n_p)          # TRUE where the match is the cell itself
    has_self <- rowSums(self) == 1L
    s_col    <- max.col(self, ties.method = "first")   # column holding the self-match

    # Columns are sorted ascending, so removing one element keeps the order:
    # the j-th non-self neighbour is column j if j < s_col, else column j + 1.
    jj     <- seq_len(k)
    colsel <- outer(s_col, jj, function(s_, j_) ifelse(j_ >= s_, j_ + 1L, j_))
    ok_col <- has_self & (colsel[, k] <= kk)

    d <- matrix(NA_real_, nrow = n_p, ncol = k)
    if (any(ok_col)) {
      r  <- which(ok_col)
      cs <- colsel[r, , drop = FALSE]
      d[r, ] <- nn$nn.dists[cbind(rep(r, k), as.vector(cs))] * 1000   # mm -> um
    }

    out$nn3_um[idx]       <- rowMeans(d)
    out$nn3rd_um[idx]     <- d[, k]
    out$n_coincident[idx] <- rowSums(nn$nn.dists == 0) - 1L

    n_drop <- sum(!ok_col)
    if (n_drop > 0L)
      cli_alert_warning(paste0("Sample {smp}: {n_drop} cell(s) have more than k + 1 ",
                              "coincident pool neurons; nn3 left NA for them."))
    if (verbose)
      cat(sprintf("  %-10s pool neurons: %-6d  median nn%d = %6.1f um  (median %d-th = %6.1f um)  coincident cells: %d\n",
                  smp, n_p, k, median(out$nn3_um[idx], na.rm = TRUE), k,
                  median(out$nn3rd_um[idx], na.rm = TRUE),
                  sum(out$n_coincident[idx] > 0L, na.rm = TRUE)))
  }

  # Metric invariants (cheap, and they catch coordinate-unit mistakes).
  ok <- !is.na(out$nn3_um)
  stopifnot(all(out$nn3_um[ok] > 0))
  stopifnot(all(out$nn3_um[ok] <= out$nn3rd_um[ok] + 1e-9))   # mean of k <= k-th

  out
}

## ---------------------------------------------------------------------------
## fov_edge_dist_um()
##
## Distance from each cell to the nearest UNIMAGED region, in microns.
##
## WHY THIS IS NEEDED. CosMx FOVs here are ~0.504 mm tiles on a ~0.51 mm pitch,
## placed along the cortical ribbon -- they do NOT fill the tissue (verified:
## 40-130 FOVs per sample covering roughly two thirds of the sample bounding
## box, laid out in offset strips rather than one global grid). A neuron beside
## an unimaged FOV has some of its true nearest neighbours invisible, so its
## nn3_um is artefactually inflated. edge_dist_um quantifies that exposure and
## enters every model as a fixed effect (adjusted for rather than excluded, so no
## cells are lost).
##
## METHOD. Per sample:
##   1. Reconstruct each FOV's imaged RECTANGLE in mm. Within a FOV,
##      x_slide_mm is affine in x_FOV_px, so a robust per-sample slope (median of
##      per-FOV least-squares slopes) plus a per-FOV intercept maps pixel 0 and
##      pixel fov_px-1 onto the rectangle edges. No nominal FOV size or um/px
##      constant is assumed -- both are read off the data. Same trick the PHF1
##      mask extractor uses to go from mask pixels to slide mm.
##   2. Rasterise the union of those rectangles at `raster_um` resolution,
##      padded so the outside of the tissue counts as unimaged.
##   3. Take the raster pixels that are UNIMAGED and 4-adjacent to an imaged
##      pixel -- the boundary of the imaged region, gaps between strips included.
##   4. edge_dist_um = RANN::nn2(k = 1) from each cell to that boundary set.
##      Exact up to the raster resolution, and correct for arbitrary FOV layouts
##      including diagonal gaps, which a grid-neighbourhood rule gets wrong.
##
## Values are capped at `cap_um` because edge censoring cannot affect an nn3_um
## smaller than the distance to the nearest unimaged pixel; beyond the cap the
## covariate would only add leverage without meaning. Callers should check the
## cap sits above the upper tail of nn3_um (Analysis A prints the quantile).
##
## Returns a data.frame: cell_id, edge_dist_um (capped), edge_dist_um_raw.
## ---------------------------------------------------------------------------
fov_edge_dist_um <- function(cd, sample_col = "sample_id", fov_col = "fov",
                             coord_x = "x_slide_mm", coord_y = "y_slide_mm",
                             px_x = "x_FOV_px", px_y = "y_FOV_px",
                             fov_px = 4256, raster_um = 10, cap_um = 300,
                             verbose = TRUE) {

  stopifnot("cell_id" %in% colnames(cd))
  for (cl in c(sample_col, fov_col, coord_x, coord_y, px_x, px_y))
    if (!cl %in% colnames(cd)) stop("fov_edge_dist_um(): missing column '", cl, "'")

  raster_mm <- raster_um / 1000
  out <- data.frame(cell_id = cd$cell_id, edge_dist_um_raw = NA_real_,
                    stringsAsFactors = FALSE)

  # Robust slope of mm per pixel, estimated per sample from every FOV with
  # enough pixel spread, then median-pooled. Sign carries the axis orientation.
  robust_slope <- function(px, mm, grp) {
    sl <- vapply(split(seq_along(px), grp), function(i) {
      if (length(i) < 20L) return(NA_real_)
      vx <- var(px[i])
      if (!is.finite(vx) || vx < 1e6) return(NA_real_)   # need a real pixel spread
      unname(stats::coef(stats::lm(mm[i] ~ px[i]))[2])
    }, numeric(1))
    sl <- sl[is.finite(sl)]
    if (!length(sl)) stop("Could not estimate a pixel->mm slope for a sample.")
    stats::median(sl)
  }

  for (smp in unique(cd[[sample_col]])) {
    rows <- which(cd[[sample_col]] == smp)
    s    <- cd[rows, , drop = FALSE]
    fv   <- as.character(s[[fov_col]])

    sx <- robust_slope(s[[px_x]], s[[coord_x]], fv)
    sy <- robust_slope(s[[px_y]], s[[coord_y]], fv)

    # Per-FOV intercept: mm at pixel 0, taken as a median so stray cells with
    # bad centroids cannot shift a rectangle.
    ix <- vapply(split(seq_len(nrow(s)), fv),
                 function(i) stats::median(s[[coord_x]][i] - sx * s[[px_x]][i]), numeric(1))
    iy <- vapply(split(seq_len(nrow(s)), fv),
                 function(i) stats::median(s[[coord_y]][i] - sy * s[[px_y]][i]), numeric(1))

    # Rectangle = pixel 0 .. fov_px - 1 mapped through the affine.
    ex <- cbind(ix, ix + sx * (fov_px - 1)); ey <- cbind(iy, iy + sy * (fov_px - 1))
    rect <- data.frame(x0 = pmin(ex[, 1], ex[, 2]), x1 = pmax(ex[, 1], ex[, 2]),
                       y0 = pmin(ey[, 1], ey[, 2]), y1 = pmax(ey[, 1], ey[, 2]))

    # Raster over the FOV union, padded by 2 pixels of guaranteed "unimaged".
    pad <- 2L
    ox  <- min(rect$x0) - pad * raster_mm
    oy  <- min(rect$y0) - pad * raster_mm
    nx  <- as.integer(ceiling((max(rect$x1) - ox) / raster_mm)) + pad
    ny  <- as.integer(ceiling((max(rect$y1) - oy) / raster_mm)) + pad
    occ <- matrix(FALSE, nrow = nx, ncol = ny)

    for (r in seq_len(nrow(rect))) {
      i0 <- max(1L, as.integer(floor((rect$x0[r] - ox) / raster_mm)) + 1L)
      i1 <- min(nx, as.integer(ceiling((rect$x1[r] - ox) / raster_mm)))
      j0 <- max(1L, as.integer(floor((rect$y0[r] - oy) / raster_mm)) + 1L)
      j1 <- min(ny, as.integer(ceiling((rect$y1[r] - oy) / raster_mm)))
      if (i1 >= i0 && j1 >= j0) occ[i0:i1, j0:j1] <- TRUE
    }

    # Unimaged pixels 4-adjacent to an imaged one = boundary of the imaged region.
    shift <- function(m, di, dj) {
      z <- matrix(FALSE, nrow(m), ncol(m))
      i_src <- max(1L, 1L - di):min(nrow(m), nrow(m) - di)
      j_src <- max(1L, 1L - dj):min(ncol(m), ncol(m) - dj)
      z[i_src + di, j_src + dj] <- m[i_src, j_src]
      z
    }
    near_occ <- shift(occ, 1, 0) | shift(occ, -1, 0) | shift(occ, 0, 1) | shift(occ, 0, -1)
    bnd <- which(!occ & near_occ, arr.ind = TRUE)
    if (!nrow(bnd)) {
      out$edge_dist_um_raw[rows] <- cap_um   # no boundary at all: fully interior
      next
    }
    # Pixel centres in mm.
    bxy <- cbind(ox + (bnd[, 1] - 0.5) * raster_mm, oy + (bnd[, 2] - 0.5) * raster_mm)

    q <- RANN::nn2(data = bxy, query = as.matrix(s[, c(coord_x, coord_y)]), k = 1)
    out$edge_dist_um_raw[rows] <- q$nn.dists[, 1] * 1000

    if (verbose)
      cat(sprintf("  %-10s FOVs: %-4d  mm/px: %+.6f, %+.6f  imaged px: %-8d median edge = %6.1f um\n",
                  smp, nrow(rect), sx, sy, sum(occ),
                  stats::median(out$edge_dist_um_raw[rows])))
  }

  out$edge_dist_um <- pmin(out$edge_dist_um_raw, cap_um)
  out
}

## ---------------------------------------------------------------------------
## tissue_edge_dist_um()
##
## Distance from each cell to the edge of the TISSUE, in microns, using only
## base R + RANN. No sf, no GDAL, no concaveman.
##
## WHY THIS EXISTS. The obvious way to get a tissue boundary is a concave hull
## (concaveman + sf::st_buffer), which needs sf/GDAL; this avoids that dependency. A
## polygon hull also cannot represent HOLES, and this cohort has them: CosMx FOVs
## are 0.51 mm tiles laid along the cortical ribbon and do not tile the tissue, so
## there are genuine unimaged gaps in the interior. Neurons beside those gaps have
## artefactually inflated neighbour distances exactly like neurons at the outer
## edge, and a hull misses every one of them.
##
## METHOD -- morphological closing on an occupancy raster, then a distance
## transform, which is the same technique fov_edge_dist_um() below already uses:
##   1. Rasterise the sample bounding box at `raster_um`, padded so the outside
##      is guaranteed non-tissue.
##   2. A pixel is TISSUE if any cell centroid lies within `close_um` of it. This
##      is a dilation by close_um and is what closes the gaps between individual
##      cells so the interior is solid rather than perforated.
##   3. Boundary = non-tissue pixels 4-adjacent to tissue -- the outer edge AND
##      the rim of every interior hole.
##   4. Distance from each cell to the nearest boundary pixel, minus close_um to
##      undo the dilation, so a cell sitting on the outermost centroid comes out
##      at ~0 rather than at +close_um. Exact to the raster resolution.
##
## CHOOSING close_um. It must exceed the typical spacing between ADJACENT CELLS
## (not neurons) or the interior perforates. Measured on this cohort, all-cell
## nearest-neighbour distance is 10.3 um median, 20.8 um at the 90th percentile
## and 49.2 um at the 99.9th, so the 50 um default closes essentially every
## interior gap while remaining far smaller than a 510 um FOV gap.
##
## Returns a data.frame: cell_id, tissue_edge_dist_um. Positive = inside the
## tissue; values near 0 sit on the boundary. Compare to a buffer to exclude.
## ---------------------------------------------------------------------------
tissue_edge_dist_um <- function(cd, sample_col = "sample_id",
                                coord_x = "x_slide_mm", coord_y = "y_slide_mm",
                                raster_um = 10, close_um = 50, verbose = TRUE) {

  stopifnot("cell_id" %in% colnames(cd))
  for (cl in c(sample_col, coord_x, coord_y))
    if (!cl %in% colnames(cd)) stop("tissue_edge_dist_um(): missing column '", cl, "'")

  r  <- raster_um / 1000
  cl <- close_um  / 1000
  out <- data.frame(cell_id = cd$cell_id, tissue_edge_dist_um = NA_real_,
                    stringsAsFactors = FALSE)

  shift <- function(m, di, dj) {
    z <- matrix(FALSE, nrow(m), ncol(m))
    i <- max(1L, 1L - di):min(nrow(m), nrow(m) - di)
    j <- max(1L, 1L - dj):min(ncol(m), ncol(m) - dj)
    z[i + di, j + dj] <- m[i, j]
    z
  }

  for (smp in unique(cd[[sample_col]])) {
    rows <- which(cd[[sample_col]] == smp)
    axy  <- as.matrix(cd[rows, c(coord_x, coord_y)])

    ox <- min(axy[, 1]) - 3 * cl; oy <- min(axy[, 2]) - 3 * cl
    nx <- as.integer(ceiling((max(axy[, 1]) + 3 * cl - ox) / r))
    ny <- as.integer(ceiling((max(axy[, 2]) + 3 * cl - oy) / r))
    gx <- ox + (seq_len(nx) - 0.5) * r
    gy <- oy + (seq_len(ny) - 0.5) * r
    grid <- cbind(rep(gx, times = ny), rep(gy, each = nx))

    occ <- matrix(RANN::nn2(data = axy, query = grid, k = 1)$nn.dists[, 1] <= cl, nx, ny)
    bnd <- which(!occ & (shift(occ, 1, 0) | shift(occ, -1, 0) |
                         shift(occ, 0, 1) | shift(occ, 0, -1)), arr.ind = TRUE)
    if (!nrow(bnd)) {
      out$tissue_edge_dist_um[rows] <- Inf   # no boundary at all: everything interior
      next
    }
    bxy <- cbind(ox + (bnd[, 1] - 0.5) * r, oy + (bnd[, 2] - 0.5) * r)
    d   <- RANN::nn2(data = bxy, query = axy, k = 1)$nn.dists[, 1] * 1000 - close_um
    out$tissue_edge_dist_um[rows] <- d

    if (verbose)
      cat(sprintf("  %-10s cells: %-6d  tissue px: %-8d boundary px: %-7d median edge = %6.1f um\n",
                  smp, length(rows), sum(occ), nrow(bnd), stats::median(d)))
  }
  out
}

## ---------------------------------------------------------------------------
## Local neuron count -- CONSTRUCT-VALIDITY DIAGNOSTIC ONLY, never an outcome.
##
## Neurons within `radius_um`, self excluded. nn3_um should be strongly
## NEGATIVELY correlated with it; if it is not, the metric is wrong.
## The count saturates at k_max and is uncorrected for that cap, which is why it
## is a diagnostic and not a model term.
## ---------------------------------------------------------------------------
local_neuron_count <- function(cd, sample_col, neuron_pool = NEURON_POOL_10,
                               coord_x = "x_slide_mm", coord_y = "y_slide_mm",
                               radius_um = 50, k_max = 60) {
  pool <- cd[as.character(cd$celltype) %in% neuron_pool, , drop = FALSE]
  out  <- data.frame(cell_id = pool$cell_id, n_local = NA_real_,
                     stringsAsFactors = FALSE)
  for (smp in unique(pool[[sample_col]])) {
    idx <- which(pool[[sample_col]] == smp)
    if (length(idx) < 2L) next
    coords <- as.matrix(pool[idx, c(coord_x, coord_y)])
    q <- RANN::nn2(data = coords, query = coords, k = min(k_max, length(idx)))
    out$n_local[idx] <- rowSums(q$nn.dists * 1000 <= radius_um) - 1
  }
  out
}

## ---------------------------------------------------------------------------
## Cortical depth
##
## L1 (pial) compass direction per sample -> normalised cortical-depth axis.
##
## Used only for the depth-ADJUSTED model variants, which are reported but are
## not the headline: neuron spacing varies strongly with cortical layer and
## PHF1+ neurons are layer-biased (L2-3 and L3-5 carry 255 of 390 PHF1+ neurons).
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
