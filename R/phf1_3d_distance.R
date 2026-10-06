#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# phf1_3d_distance.R
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
# phf1_3d_distance.R
#
# STAGE 1 of the out-of-plane tangle correction. Produces:
#
#   1. SIMULATED 3D TANGLE COORDINATES, one set per donor, in the real FOV
#      mosaic -> results/tangle_3d/simulated_coords.tsv.gz (projection validation)
#   2. The per-donor 2D -> 3D INFLATION FACTOR, from lambda2 measured by COUNTS
#      and imaged area, not by inverting a pooled nearest-neighbour median.
#   3. UNSEEN-PLANE CORRECTED DISTANCE TABLES keyed by cell_id, one wide TSV per
#      configuration -> dist3d/ (plus the identity table and, with --arm_b, the
#      in-plane augmented-anchor tables). These are what the Stage 2 DEG consumes.
#
# ============================================================================
# WHY THIS IS A SEPARATE STAGE FROM THE DEG
# ============================================================================
# spatstat is not installed in the HPC environment where dream runs, so everything
# spatial happens HERE, locally, and Stage 2 reads nothing but plain TSVs.
#
# ============================================================================
# WHAT THE SIMULATION IS AND IS NOT FOR
# ============================================================================
# p_detect is ANALYTIC -- (D_i + t)/(D_s + t), see t3_geometry(). The simulation
# does NOT estimate it. Its three jobs are:
#   (a) validation: can a projected 3D process reproduce the observed anchor
#       pattern at all? If not, the correction rests on a model that cannot make
#       the data, and the run aborts.
#   (b) the clustered scenario, which the analytic Poisson formula cannot give.
#   (c) projection-validation coordinates (simulated_coords.tsv.gz).
#
# The distances that reach the DEG are computed for REAL cells (section 5);
# simulated tangles serve only as additional anchors. They have no cell_id and
# cannot carry expression.
#
# ============================================================================
# ARM B (--arm_b) IS BIASED TOWARD NULL BY CONSTRUCTION
# ============================================================================
# Hidden anchors are a position-weighted random draw, not an identification of
# which cells are mislabelled, so augmentation attenuates the distance slope
# through anchor-set dilution WHETHER OR NOT the field effect is real:
#
#     survival under B  -> robustness
#     weakening under B -> UNINFORMATIVE
#
# The size of that attenuation is measured, not guessed -- see
# results/tangle_3d/armB_dilution_constant.tsv from R/test_tangle_3d_utils.R.
# This script refuses to run without it.
#
# Usage (local, needs spatstat):
#   Rscript R/phf1_3d_distance.R                       # headline + corners
#   Rscript R/phf1_3d_distance.R --R 3 --skip_sim      # quick, distances only

suppressPackageStartupMessages({
  library(argparse); library(dplyr)
})

.this_dir <- local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else "R"
})
source(file.path(.this_dir, "phf1_distance_utils.r"))
source(file.path(.this_dir, "nn3_utils.r"))
source(file.path(.this_dir, "kernel_field_utils.R"))
source(file.path(.this_dir, "tangle_field_utils.R"))
source(file.path(.this_dir, "tangle_3d_utils.R"))
source(file.path(.this_dir, "phf1_hidden_anchors.R"))

parser <- ArgumentParser()
parser$add_argument("--coords",   default = "PHF1/seu_coords.csv")
parser$add_argument("--celltype", default = "Exc-IT-L2-3-CBLN2-HOPX")
parser$add_argument("--R",        type = "integer", default = 10L,
                    help = "Replicate augmented draws at the headline cell [default 10]")
parser$add_argument("--R_corner", type = "integer", default = 3L,
                    help = "Replicates at the envelope corners [default 3]")
parser$add_argument("--sigma_um", type = "double", default = HIDDEN_SIGMA_UM,
                    help = "Kernel sigma for intensity weighting (um)")
parser$add_argument("--sigma_xy", type = "double", default = 150)
parser$add_argument("--sigma_z",  type = "double", default = 150)
parser$add_argument("--mu",       type = "double", default = 100,
                    help = "Offspring per cluster parent for the Thomas scenario")
parser$add_argument("--max_dist_um", type = "double", default = T3_MAX_DIST_UM)
parser$add_argument("--seed",     type = "integer", default = 42L)
parser$add_argument("--z_max_um", type = "double", default = 500,
                    help = "Depth above/below the section for unseen tangles (um)")
parser$add_argument("--arm_b", action = "store_true", default = FALSE,
                    help = "Also emit the in-plane hidden-anchor (stress-test) tables")
parser$add_argument("--skip_sim", action = "store_true", default = FALSE,
                    help = "Skip the 3D simulation; produce distance tables only")
parser$add_argument("--sce_dir",  default = "celltype_sce_neighbours",
                    help = "Only used to source the identity table from the SCE's own stored column")
parser$add_argument("--distdir",  default = "dist3d")
parser$add_argument("--resdir",   default = "results/tangle_3d")
parser$add_argument("--outdir",   default = "plots/phf1_3d_distance")
args <- parser$parse_args()

CT <- args$celltype
for (d in c(args$distdir, args$resdir, args$outdir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

sink_path <- file.path(args$outdir, "stats_phf1_3d_distance.txt")
con <- file(sink_path, open = "wt"); sink(con, split = TRUE)

cat("=== Stage 1: out-of-plane tangle correction ===\n")
cat("Celltype for the DEG arm:", CT, "\n")
cat("Seed:", args$seed, "| max_dist_um:", args$max_dist_um, "\n\n")

## ---------------------------------------------------------------------------
## 0. The Arm B dilution constant must exist before anything else
## ---------------------------------------------------------------------------
DIL <- arm_b_dilution_constant(file.path(args$resdir, "armB_dilution_constant.tsv"))
cat(sprintf("Arm B dilution constant: %.3f [%.3f, %.3f]\n",
            DIL[["dilution_slope_aug_over_obs"]], DIL[["dilution_ci_lo"]],
            DIL[["dilution_ci_hi"]]))
cat("  Every Arm B slope is reported against this. A weakening of this size is\n")
cat("  expected from the procedure itself (anchor-set dilution).\n\n")

## ---------------------------------------------------------------------------
## 1. Coordinates, anchor ledger, and the distance identity check
## ---------------------------------------------------------------------------
cd <- read.csv(args$coords, stringsAsFactors = FALSE)
cd$PHF1     <- as.character(cd$PHF1)
cd$celltype <- as.character(cd$celltype)
stopifnot(all(c("cell_id", "sample_id", "fov", "x_slide_mm", "y_slide_mm",
                "x_FOV_px", "y_FOV_px", "PHF1", "celltype") %in% colnames(cd)))

anchors <- phf1_source_idx(cd, NEURON_POOL_10)
n_phf1  <- sum(cd$PHF1 == "TRUE")
cat("=== Anchor ledger ===\n")
cat("cells:", nrow(cd), "| PHF1+:", n_phf1, "| anchors (NEURON_POOL_10):",
    length(anchors), "\n")
per_donor <- sort(table(cd$sample_id[anchors]))
print(per_donor)
# Assert against the published anchor counts rather than trusting the file.
if (n_phf1 != 601L || length(anchors) != 398L)
  stop("Anchor ledger disagrees with the published counts (expected 601 PHF1+ / 398 ",
       "anchors, got ", n_phf1, " / ", length(anchors), "). Do not proceed.")
cat("matches the published counts (601 PHF1+ / 398 anchors)\n\n")

d_can <- compute_dist_to_phf1_um(cd, "sample_id", NEURON_POOL_10, verbose = FALSE)
cd$dist_to_phf1_um <- d_can[match(cd$cell_id, names(d_can))]

# The identity check: p_detect = 1 draws H = 0 hidden anchors, so the override
# path must reproduce the canonical k = 1 path EXACTLY. This also proves the
# k = 2 self-match branch is a no-op when no hidden anchor is added.
dr_id <- draw_hidden_anchors(cd, NEURON_POOL_10, p_detect = 1)
d_id  <- compute_dist_to_phf1_um(cd, "sample_id", NEURON_POOL_10,
                                 source_idx_override = dr_id, verbose = FALSE)
cmp   <- is.finite(d_can) & is.finite(d_id)
dmax  <- max(abs(d_can[cmp] - d_id[cmp]))
cat(sprintf("IDENTITY CHECK (p_detect = 1): max abs diff = %.3g um over %d cells\n",
            dmax, sum(cmp)))
if (dmax > 1e-6)
  stop("The override path does not reproduce the canonical distance. Do not proceed.")
cat("\n")

## ---------------------------------------------------------------------------
## 2. Per-donor imaged area and lambda2 FROM COUNTS
## ---------------------------------------------------------------------------
cat("=== Per-donor anchor intensity (from counts and imaged area) ===\n")
samples <- sort(unique(cd$sample_id))
wins <- list(); areas <- numeric(0)
for (s in samples) {
  cds <- cd[cd$sample_id == s, , drop = FALSE]
  rect <- fov_rectangles(cds)
  w    <- fov_owin(rect, verbose = FALSE)
  wins[[s]]  <- w
  areas[[s]] <- spatstat.geom::area(w)
}
n_anch <- vapply(samples, function(s)
  sum(cd$PHF1 == "TRUE" & cd$celltype %in% NEURON_POOL_10 & cd$sample_id == s),
  numeric(1))
lam <- t3_lambda2_from_counts(n_anch, areas[samples], samples)
print(lam, row.names = FALSE, digits = 4)
LAM2 <- lam$lambda2_per_mm2[lam$sample_id == "POOLED"]

d_obs_all <- cd$dist_to_phf1_um[is.finite(cd$dist_to_phf1_um) &
                                cd$dist_to_phf1_um <= args$max_dist_um]
lam2_inv <- 1e6 * log(2) / (pi * median(d_obs_all)^2)
cat(sprintf("\npooled lambda2 from COUNTS          : %.3f /mm2\n", LAM2))
cat(sprintf("pooled lambda2 from median inversion: %.3f /mm2  (observed median %.1f um)\n",
            lam2_inv, median(d_obs_all)))
cat("The two disagree because cells beside unimaged FOV gaps have UPWARD-CENSORED\n")
cat("nearest-anchor distances, which inflates the pooled median and deflates the\n")
cat("inverted intensity. The count-based figure is the sounder one and is used below.\n\n")

## ---------------------------------------------------------------------------
## 3. Geometry ledger and the inflation table
## ---------------------------------------------------------------------------
grid <- t3_geometry_grid()
cells <- grid[grid$corner != "interior", , drop = FALSE]
cat("=== Geometry cells to be carried through (derived from the triple) ===\n")
print(cells[, c("t_um", "d_incl_um", "d_soma_um", "p_detect", "h_eff", "h_cell", "corner")],
      row.names = FALSE, digits = 4)

infl <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  x <- t3_median_d3d(LAM2, cells$h_eff[i], d2_observed_um = median(d_obs_all))
  cbind(corner = cells$corner[i], p_detect = cells$p_detect[i], x)
}))
cat("\n=== 2D -> 3D inflation (Poisson, lambda2 from counts) ===\n")
print(infl[, c("corner", "p_detect", "h_eff_um", "median_d2d_ref_um",
               "median_d3d_um", "inflation")], row.names = FALSE, digits = 4)
write.table(infl, file.path(args$resdir, "inflation_table.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Per-donor inflation, since lambda2 spans a 23x range across donors.
infl_d <- do.call(rbind, lapply(samples, function(s) {
  l <- lam$lambda2_per_mm2[lam$sample_id == s]
  ds <- cd$dist_to_phf1_um[cd$sample_id == s & is.finite(cd$dist_to_phf1_um) &
                           cd$dist_to_phf1_um <= args$max_dist_um]
  x <- t3_median_d3d(l, T3_D_INCL_UM + T3_T_UM, d2_observed_um = median(ds))
  cbind(sample_id = s, n_anchor = lam$n_anchor[lam$sample_id == s], x)
}))
cat("\n=== Per-donor inflation at the headline h_eff ===\n")
print(infl_d[, c("sample_id", "n_anchor", "lambda2_per_mm2", "median_d2d_ref_um",
                 "median_d3d_um", "inflation")], row.names = FALSE, digits = 4)
write.table(infl_d, file.path(args$resdir, "inflation_by_donor.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## ---------------------------------------------------------------------------
## 4. Simulated 3D coordinates
## ---------------------------------------------------------------------------
if (!args$skip_sim) {
  cat("\n=== Simulating tangles in continuous 3D, per donor ===\n")
  geom_h <- t3_geometry()          # headline triple
  sim_rows <- list(); val_rows <- list()

  for (s in samples) {
    w      <- wins[[s]]
    n_a    <- as.integer(n_anch[[s]])
    if (n_a < 1L) { cat(sprintf("  %-9s no anchors -- skipped\n", s)); next }
    padz   <- t3_pad_z(lam$lambda2_per_mm2[lam$sample_id == s], geom_h$h_eff,
                       args$sigma_z)
    n_tg   <- t3_n_tangle_for(n_a, padz, geom_h$h_eff)
    n_neur <- max(n_tg * 4L, 20000L)

    for (mdl in c("poisson", "thomas")) {
      sim <- t3_simulate_3d(w, n_neuron = n_neur, n_tangle = n_tg, model = mdl,
                            sigma_xy = args$sigma_xy, sigma_z = args$sigma_z,
                            mu = args$mu, neuron_model = "uniform",
                            pad_z_um = padz, seed = args$seed + which(samples == s))
      # Query from the simulated neurons that are in the object, matching the
      # real analysis (query cells are cells, not arbitrary points).
      qq <- sim$neurons[abs(sim$neurons$z_um) < geom_h$h_cell / 2,
                        c("x_mm", "y_mm"), drop = FALSE]
      pr <- t3_project(sim, geom_h, query_xy_mm = qq, max_dist_um = args$max_dist_um)

      tg <- pr$tangles
      sim_rows[[paste(s, mdl)]] <- data.frame(
        sample_id = s, model = mdl,
        x_um = tg$x_mm * 1000, y_um = tg$y_mm * 1000, z_um = tg$z_um,
        class = tg$class, z_source = "simulated", z_measured = FALSE,
        stringsAsFactors = FALSE)

      d_obs_s <- cd$dist_to_phf1_um[cd$sample_id == s &
                                    is.finite(cd$dist_to_phf1_um) &
                                    cd$dist_to_phf1_um <= args$max_dist_um]
      val_rows[[paste(s, mdl)]] <- data.frame(
        sample_id = s, model = mdl, n_anchor_obs = n_a,
        n_detected = pr$n_detected,
        n_missed_recoverable = pr$n_missed_recoverable,
        n_missed_deep = pr$n_missed_deep,
        p_detect_realised = pr$p_detect_realised,
        median_d2d_sim = if (length(pr$d2_um)) median(pr$d2_um) else NA_real_,
        median_d2d_obs = median(d_obs_s),
        median_d3d_sim = if (length(pr$d3_um)) median(pr$d3_um) else NA_real_,
        stringsAsFactors = FALSE)
    }
    cat(sprintf("  %-9s anchors %3d | pad_z %5.0f um | n_tangle %6d\n",
                s, n_a, padz, n_tg))
  }

  sim_df <- bind_rows(sim_rows); val_df <- bind_rows(val_rows)

  # Observed anchors, carried into the same table so the figure has one source.
  # z is NOT MEASURED for these -- there is a single confocal plane. They are
  # placed at the slab mid-plane because that is the only defensible position,
  # and the flag says so on every row.
  obs_df <- cd[cd$PHF1 == "TRUE" & cd$celltype %in% NEURON_POOL_10, , drop = FALSE]
  obs_df <- data.frame(
    sample_id = obs_df$sample_id, model = "observed",
    x_um = obs_df$x_slide_mm * 1000, y_um = obs_df$y_slide_mm * 1000, z_um = 0,
    class = "observed_anchor", z_source = "observed_plane_assumed",
    z_measured = FALSE, stringsAsFactors = FALSE)

  out_sim <- bind_rows(sim_df, obs_df)
  gz <- gzfile(file.path(args$resdir, "simulated_coords.tsv.gz"), "wt")
  write.table(out_sim, gz, sep = "\t", row.names = FALSE, quote = FALSE)
  close(gz)

  cat("\n=== Projection validation, per donor ===\n")
  print(val_df, row.names = FALSE, digits = 4)
  write.table(val_df, file.path(args$resdir, "projection_validation.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  vp <- val_df[val_df$model == "poisson", ]
  rel <- abs(vp$median_d2d_sim - vp$median_d2d_obs) / pmax(vp$median_d2d_obs, 1)
  cat(sprintf("\nPoisson: simulated vs observed median d2D, median rel. diff = %.1f%%\n",
              100 * median(rel, na.rm = TRUE)))
  cat(sprintf("realised p_detect %.3f (analytic %.3f)\n",
              mean(vp$p_detect_realised, na.rm = TRUE), geom_h$p_detect))
  cat(sprintf("implied 3D inflation from the simulation: %.2fx\n",
              median(vp$median_d2d_sim / vp$median_d3d_sim, na.rm = TRUE)))
  cat("\nNOTE: z is NOT MEASURED for observed anchors -- one confocal plane. They\n")
  cat("carry z_measured = FALSE and z_source = 'observed_plane_assumed'.\n")
}

## ---------------------------------------------------------------------------
## 5. UNSEEN-PLANE distance tables -- what the DEG consumes
##
## WHAT THIS DOES NOT ASSUME. Every tangle in the imaged plane is taken as
## detected. Nothing is relabelled and no misclassification is posited. The
## in-plane hidden-anchor arm, which does posit unlabelled tangle-bearing
## neurons, is kept behind --arm_b as a stress test only.
##
## What a single section cannot see is tissue it never cut. Tangles are
## simulated ONLY at |z| > h_eff/2, at the volumetric intensity implied by
## the observed in-plane density, and each cell's corrected distance is the 3D
## minimum over the observed in-plane anchors and those.
##
## The correction shortens distances mostly in the TAIL: near cells keep their
## distance, far cells move closer. The Spearman rho(d2D, d3D) for each
## configuration is printed below and written to MANIFEST.tsv.
## ---------------------------------------------------------------------------
geom_h_eff <- t3_geometry()$h_eff
cat("\n=== Unseen-plane corrected distances ===\n")
cat(sprintf("Unseen tangles simulated over +/-%g um, outside the %g um sampled slice.\n",
            args$z_max_um, geom_h_eff))

ct_ids <- cd$cell_id[cd$celltype == CT]
cat("Cells in", CT, ":", length(ct_ids), "\n")
man <- list()

unseen_table <- function(h_eff, clustered, R, seed0) {
  M <- matrix(NA_real_, nrow = length(ct_ids), ncol = R,
              dimnames = list(ct_ids, NULL))
  d2ref <- rep(NA_real_, length(ct_ids)); nun <- integer(R)
  for (i in seq_len(R)) {
    acc <- list()
    for (s in samples) {
      cs  <- cd[cd$sample_id == s, , drop = FALSE]
      w   <- wins[[s]]
      anc <- cs[cs$PHF1 == "TRUE" & cs$celltype %in% NEURON_POOL_10,
                c("x_slide_mm", "y_slide_mm"), drop = FALSE]
      q   <- cs[cs$cell_id %in% ct_ids & cs$PHF1 != "TRUE", , drop = FALSE]
      if (!nrow(q) || !nrow(anc)) next
      lam2 <- nrow(anc) / spatstat.geom::area(w)
      u <- t3_unseen_plane_anchors(w, lam2, h_eff, z_max_um = args$z_max_um,
                                   clustered = clustered,
                                   anchor_xy_mm = as.matrix(anc),
                                   sigma_um = args$sigma_um,
                                   seed = seed0 + i * 1000L + which(samples == s))
      dd <- t3_dist3d_unseen(q[, c("x_slide_mm", "y_slide_mm")], anc, u)
      nun[i] <- nun[i] + attr(dd, "n_unseen")
      acc[[s]] <- data.frame(cell_id = q$cell_id, d3 = as.numeric(dd),
                             d2 = as.numeric(attr(dd, "d2d")))
    }
    a <- bind_rows(acc)
    M[, i] <- a$d3[match(ct_ids, a$cell_id)]
    if (i == 1L) d2ref <- a$d2[match(ct_ids, a$cell_id)]
  }
  list(M = M, d2 = d2ref, n_unseen = round(mean(nun)))
}

cfgs <- list(
  list(h = NA, cl = TRUE,  R = args$R,        tag = "headline", nm = "unseen_clustered"),
  list(h = NA, cl = FALSE, R = args$R,        tag = "bound",    nm = "unseen_uniform"),
  list(h = 16, cl = TRUE,  R = args$R_corner, tag = "heff16",   nm = "unseen_clustered_h16"),
  list(h = 20, cl = TRUE,  R = args$R_corner, tag = "heff20",   nm = "unseen_clustered_h20"))

for (cf in cfgs) {
  h <- if (is.na(cf$h)) geom_h_eff else cf$h
  cat(sprintf("\n-- %s (h_eff %g um, %s, R = %d) --\n", cf$nm, h,
              if (cf$cl) "clustered" else "uniform", cf$R))
  r <- unseen_table(h, cf$cl, cf$R, args$seed)
  keep <- is.finite(r$d2) & r$d2 <= args$max_dist_um
  if (any(r$M[keep, ] <= 0, na.rm = TRUE))
    stop("Non-positive corrected distance in ", cf$nm, " -- log() requires > 0.")
  rho <- suppressWarnings(cor(r$d2[keep], r$M[keep, 1], method = "spearman"))
  cat(sprintf("   unseen tangles simulated: %d (observed anchors 398)\n", r$n_unseen))
  cat(sprintf("   median d2D %.1f -> d3D %.1f um (%.1f%% shorter) | Spearman rho = %.3f\n",
              median(r$d2[keep]), median(r$M[keep, ], na.rm = TRUE),
              100 * (1 - median(r$M[keep, ], na.rm = TRUE) / median(r$d2[keep])), rho))
  fmtq <- function(v) ifelse(is.na(v), "NA", sprintf("%.17g", v))
  colnames(r$M) <- sprintf("rep_%02d", seq_len(ncol(r$M)))
  out <- data.frame(cell_id = ct_ids, dist_observed_um = fmtq(r$d2),
                    lapply(as.data.frame(r$M), fmtq), check.names = FALSE)
  f <- file.path(args$distdir, paste0(cf$nm, ".tsv"))
  write.table(out, f, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("   written:", f, "\n")
  man[[cf$nm]] <- data.frame(config = cf$nm, tag = cf$tag, p_detect = NA_real_,
    weight = if (cf$cl) "clustered" else "uniform", R = cf$R, h_eff = h, file = f,
    n_cells = nrow(out), median_obs_um = median(r$d2[keep]),
    median_aug_um = median(r$M[keep, ], na.rm = TRUE),
    spearman_rho = round(rho, 4), n_unseen = r$n_unseen, stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## 5a. Coordinates for the per-sample 3D figure
##
## One headline replicate, exported so the figure shows EXACTLY what the method
## does: the observed in-plane tangles (taken as complete) and the simulated
## tangles in the planes the section never cut. There are only TWO classes.
## (results/tangle_3d/simulated_coords.tsv.gz, from section 4, is the projection
## validation output and carries different classes; it is not the figure source.)
## ---------------------------------------------------------------------------
cat("\n=== Coordinates for the 3D figure (headline settings, one replicate) ===\n")
coord_rows <- list()
for (s in samples) {
  cs  <- cd[cd$sample_id == s, , drop = FALSE]
  w   <- wins[[s]]
  anc <- cs[cs$PHF1 == "TRUE" & cs$celltype %in% NEURON_POOL_10,
            c("x_slide_mm", "y_slide_mm"), drop = FALSE]
  if (!nrow(anc)) next
  lam2 <- nrow(anc) / spatstat.geom::area(w)
  u <- t3_unseen_plane_anchors(w, lam2, geom_h_eff, z_max_um = args$z_max_um,
                               clustered = TRUE, anchor_xy_mm = as.matrix(anc),
                               sigma_um = args$sigma_um, seed = args$seed)
  coord_rows[[paste0(s, "_obs")]] <- data.frame(
    sample_id = s, x_um = anc$x_slide_mm * 1000, y_um = anc$y_slide_mm * 1000,
    z_um = 0, class = "observed_in_section", z_measured = FALSE,
    z_source = "observed_plane_assumed", stringsAsFactors = FALSE)
  if (nrow(u))
    coord_rows[[paste0(s, "_uns")]] <- data.frame(
      sample_id = s, x_um = u$x_mm * 1000, y_um = u$y_mm * 1000, z_um = u$z_um,
      class = "simulated_unseen", z_measured = FALSE, z_source = "simulated",
      stringsAsFactors = FALSE)
}
uc <- bind_rows(coord_rows)
gz <- gzfile(file.path(args$resdir, "unseen_coords.tsv.gz"), "wt")
write.table(uc, gz, sep = "\t", row.names = FALSE, quote = FALSE); close(gz)
print(table(uc$sample_id, uc$class))
cat("written:", file.path(args$resdir, "unseen_coords.tsv.gz"), "\n")

## ---------------------------------------------------------------------------
## 5b. In-plane hidden-anchor tables (--arm_b, stress test)
## ---------------------------------------------------------------------------
if (args$arm_b) {

cat("\n=== Arm B (stress test): in-plane augmented-anchor distances ===\n")
configs <- list()
hl <- cells[cells$corner == "headline", ]
configs[[1]] <- list(p = hl$p_detect, w = "intensity", R = args$R,       tag = "headline")
configs[[2]] <- list(p = hl$p_detect, w = "uniform",   R = max(5L, args$R %/% 2L), tag = "headline")
for (i in which(cells$corner != "headline"))
  configs[[length(configs) + 1]] <- list(p = cells$p_detect[i], w = "intensity",
                                         R = args$R_corner, tag = cells$corner[i])

man <- list()
for (cfg in configs) {
  ptag <- gsub("\\.", "p", sprintf("%.3f", cfg$p))
  nm   <- sprintf("aug_p%s_%s", ptag, cfg$w)
  cat(sprintf("\n-- %s (%s, R = %d) --\n", nm, cfg$tag, cfg$R))

  M <- hidden_dist_matrix(cd, ct_ids, NEURON_POOL_10, p_detect = cfg$p,
                          R = cfg$R, weight = cfg$w, sigma_um = args$sigma_um,
                          seed = args$seed, verbose = FALSE)
  colnames(M) <- sprintf("rep_%02d", seq_len(ncol(M)))

  ref <- cd$dist_to_phf1_um[match(ct_ids, cd$cell_id)]
  keep <- is.finite(ref) & ref <= args$max_dist_um
  cat(sprintf("   median d2D %.1f -> %.1f um (%.1f%% shorter) over %d modelled cells\n",
              median(ref[keep]), median(M[keep, ], na.rm = TRUE),
              100 * (1 - median(M[keep, ], na.rm = TRUE) / median(ref[keep])),
              sum(keep)))
  st <- attr(M, "draw_stats")
  cat(sprintf("   H per replicate: %d | uniform fallbacks: %d\n",
              sum(st$H[st$rep == 1]), sum(st$weight == "uniform_fallback")))
  if (any(M[keep, ] <= 0, na.rm = TRUE))
    stop("Non-positive augmented distance in ", nm,
         " -- the log transform requires dist > 0. Fix upstream, do not offset.")

  # Full precision. write.table's default 15 significant digits loses ~5e-11 um
  # on the round trip, which is physically nothing but breaks the identity check
  # from being exact -- and an identity check that is only approximately exact
  # invites the question of how approximate is acceptable. %.17g round-trips a
  # double exactly.
  fmt <- function(v) ifelse(is.na(v), "NA", sprintf("%.17g", v))
  out <- data.frame(cell_id = ct_ids, dist_observed_um = fmt(ref),
                    lapply(as.data.frame(M), fmt), check.names = FALSE)
  f <- file.path(args$distdir, paste0(nm, ".tsv"))
  write.table(out, f, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("   written:", f, "\n")

  man[[nm]] <- data.frame(
    config = nm, tag = cfg$tag, p_detect = cfg$p, weight = cfg$w, R = cfg$R,
    h_eff = cells$h_eff[which.min(abs(cells$p_detect - cfg$p))],
    file = f, n_cells = nrow(out),
    median_obs_um = median(ref[keep]),
    median_aug_um = median(M[keep, ], na.rm = TRUE),
    spearman_rho = NA_real_, n_unseen = NA_integer_,
    stringsAsFactors = FALSE)
}
}   # end --arm_b

# Identity config: p_detect = 1, H = 0. The DEG run on this MUST reproduce the
# canonical table exactly, which is the proof that the fork changed nothing but
# the distance.
# The identity table is sourced from the SCE's OWN STORED COLUMN, not recomputed
# from seu_coords.csv. Recomputing agrees only to ~5e-11 um, because
# seu_coords.csv is a 15-significant-digit CSV export of coordinates that the
# stored distances were computed from at full precision. 5e-11 um is 5e-17 m and
# means nothing physically -- but the identity arm exists to prove the DEG fork
# changed NOTHING but the distance, and a check that is only approximately exact
# invites the question of how approximate is acceptable. Taking the stored column
# makes it exact. The distance MACHINERY is separately proven faithful by the
# override-path check above (max abs diff 0 um), so nothing is left untested.
sce_id <- qs::qread(file.path(args$sce_dir, sprintf("%s_sce_neighbours.qs", CT)))
sce_cid <- if ("cell_id" %in% colnames(SummarizedExperiment::colData(sce_id)))
             sce_id$cell_id else colnames(sce_id)
d_ct <- unname(sce_id$dist_to_phf1_um[match(ct_ids, sce_cid)])
rm(sce_id); gc(verbose = FALSE)
.fmt17 <- function(v) ifelse(is.na(v), "NA", sprintf("%.17g", v))
write.table(data.frame(cell_id = ct_ids, dist_observed_um = .fmt17(d_ct),
                       rep_01 = .fmt17(d_ct)),
            file.path(args$distdir, "aug_p1p000_identity.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
man[["identity"]] <- data.frame(
  config = "aug_p1p000_identity", tag = "identity", p_detect = 1, weight = "none",
  R = 1L, h_eff = NA_real_, file = file.path(args$distdir, "aug_p1p000_identity.tsv"),
  n_cells = length(ct_ids), median_obs_um = median(d_ct[is.finite(d_ct)]),
  median_aug_um = median(d_ct[is.finite(d_ct)]),
  spearman_rho = 1, n_unseen = 0L, stringsAsFactors = FALSE)

manifest <- bind_rows(man)
write.table(manifest, file.path(args$distdir, "MANIFEST.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n=== Distance table manifest ===\n")
print(manifest[, c("config", "tag", "weight", "R", "h_eff", "median_obs_um",
                   "median_aug_um", "spearman_rho", "n_unseen")],
      row.names = FALSE, digits = 4)
cat("\nTotal DEG fits implied:", sum(manifest$R), "\n")

cat("\n"); print(sessionInfo())
sink(); close(con)
cat("\nStage 1 complete. Distance tables in", args$distdir, "\n")
cat("Stats log:", sink_path, "\n")
