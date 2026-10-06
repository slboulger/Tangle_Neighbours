#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# test_tangle_3d_utils.R
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
# test_tangle_3d_utils.R
#
# Calibration gate for the out-of-plane tangle correction. Pure simulation --
# no SCE, no Seurat, no cluster. Runs locally in a couple of minutes and exits 1
# on any failure.
#
#   Rscript R/test_tangle_3d_utils.R [--reps 40] [--seed 1] [--quick]
#
# ============================================================================
# WHAT EACH SECTION IS FOR, AND WHICH ONES ARE GATES
# ============================================================================
#
#  1. GEOMETRY. p_detect = (D_i + t)/(D_s + t) against Monte Carlo, including
#     the eccentric offset. Asserts the delta-invariance that the header of
#     tangle_3d_utils.R claims: p_detect does not move for any |delta| within
#     (D_s - D_i)/2. If this fails the closed form is wrong and every downstream
#     number is wrong with it.
#
#  2. POISSON INVERSION. The 2D and 3D median nearest-neighbour closed forms
#     against direct simulation. Cheap, and it is the arithmetic that produced
#     the headline "in-plane distance overstates true separation by ~2.8x".
#
#  3. PROJECTION. Simulate, project, and check that the realised detection
#     fraction equals p_detect and that the realised in-plane anchor intensity
#     equals the target. This is what makes t3_n_tangle_for() trustworthy; get
#     it wrong and the simulated cohort is a hundredfold too sparse.
#
#  4. THE ARM B DILUTION CONSTANT -- THE GATE THAT MATTERS.
#     Hidden anchors are a position-weighted random draw, NOT an identification
#     of which cells are mislabelled, so the expected consequence of
#     augmentation is ATTENUATION OF THE DISTANCE SLOPE THROUGH ANCHOR-SET
#     DILUTION, independently of the underlying effect.
#
#     This section measures how much. It plants a known field, deletes a known
#     fraction of anchors, adds the same number back by weighted draw, and
#     reports the slope ratio. THAT RATIO IS THE CONSTANT EVERY ARM B NUMBER IS
#     REPORTED AGAINST; it separates the attenuation produced by the procedure
#     itself from any change in the field.
#
#     Note what is deliberately NOT asserted here: that augmentation RECOVERS
#     the true slope. It does not, and it is not supposed to.
#
#  5. GUARDS. Non-positive distance, inclusion larger than soma, grid centre on
#     the grid, n_tangle inversion round-trip.

suppressPackageStartupMessages({
  library(argparse)
})

.this_dir <- local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else "R"
})
source(file.path(.this_dir, "tangle_3d_utils.R"))

parser <- ArgumentParser()
parser$add_argument("--reps",  type = "integer", default = 40L)
parser$add_argument("--seed",  type = "integer", default = 1L)
parser$add_argument("--quick", action = "store_true", default = FALSE,
                    help = "Smaller n for a fast smoke run; NOT a valid gate.")
args <- parser$parse_args()
if (args$quick) args$reps <- max(8L, args$reps %/% 4L)

FAILURES <- character(0)
check <- function(ok, msg) {
  ok <- isTRUE(ok)
  cat(if (ok) "  [PASS] " else "  [FAIL] ", msg, "\n", sep = "")
  if (!ok) FAILURES <<- c(FAILURES, msg)
  invisible(ok)
}
hdr <- function(x) cat("\n", strrep("=", 72), "\n", x, "\n", strrep("=", 72), "\n", sep = "")

set.seed(args$seed)


## ===========================================================================
hdr("1. GEOMETRY -- p_detect closed form, and delta-invariance")
## ===========================================================================

g <- t3_geometry(t_um = 7, d_incl_um = 11, d_soma_um = 17)
check(abs(g$p_detect - 0.75) < 1e-12,
      sprintf("p_detect(7, 11, 17) = %.6f, expected exactly 0.75", g$p_detect))
check(g$h_eff == 18 && g$h_cell == 24,
      sprintf("h_eff = %g (18), h_cell = %g (24)", g$h_eff, g$h_cell))

# Monte Carlo, at delta = 0 and at the extreme admissible offset.
for (dl in c(0, 1.5, 3)) {
  mc <- t3_p_detect_mc(7, 11, 17, delta_um = dl, n = 1e6L, seed = args$seed + dl)
  check(abs(mc - 0.75) < 0.005,
        sprintf("MC p_detect at delta = %.1f um: %.4f (analytic 0.750)", dl, mc))
}

# Beyond the admissible range p_detect genuinely DOES fall. Asserting this
# stops a future reader "fixing" the invariance into a universal claim.
mc_wide <- t3_p_detect_mc(7, 11, 17, delta_um = 8, n = 1e6L, seed = args$seed)
check(mc_wide < 0.74,
      sprintf("beyond |delta| <= (Ds-Di)/2 the invariance ends: %.4f at delta = 8 um",
              mc_wide))

grid <- t3_geometry_grid()
check(sum(grid$corner == "headline") == 1L,
      "grid contains exactly one 'headline' cell (centre is ON the grid)")
check(abs(min(grid$p_detect) - 0.615) < 0.002 && abs(max(grid$p_detect) - 0.833) < 0.002,
      sprintf("grid p_detect spans %.3f - %.3f (expected 0.615 - 0.833)",
              min(grid$p_detect), max(grid$p_detect)))
check(min(grid$h_eff) == 16 && max(grid$h_eff) == 20,
      sprintf("grid h_eff spans %g - %g um (expected 16 - 20)",
              min(grid$h_eff), max(grid$h_eff)))
check(all(grid$d_incl_um <= grid$d_soma_um),
      "no grid cell has the inclusion larger than the soma")


## ===========================================================================
hdr("2. POISSON INVERSION -- 2D and 3D median NN against direct simulation")
## ===========================================================================

# 2D
L2 <- 4000; n2 <- 20000
p2 <- cbind(runif(n2, 0, L2), runif(n2, 0, L2))
core2 <- p2[, 1] > 800 & p2[, 1] < 3200 & p2[, 2] > 800 & p2[, 2] < 3200
nn2d <- RANN::nn2(p2, p2[core2, , drop = FALSE], k = 2)$nn.dists[, 2]
lam2_um <- n2 / L2^2
med2_cf <- sqrt(log(2) / (pi * lam2_um))
check(abs(median(nn2d) / med2_cf - 1) < 0.05,
      sprintf("2D median NN: sim %.2f vs closed form %.2f (ratio %.3f)",
              median(nn2d), med2_cf, median(nn2d) / med2_cf))

# 3D
L3 <- 1000; n3 <- 60000
p3 <- cbind(runif(n3, 0, L3), runif(n3, 0, L3), runif(n3, 0, L3))
core3 <- apply(p3, 1, function(r) all(r > 200 & r < 800))
nn3d <- RANN::nn2(p3, p3[core3, , drop = FALSE], k = 2)$nn.dists[, 2]
lam3_um <- n3 / L3^3
med3_cf <- (3 * log(2) / (4 * pi * lam3_um))^(1 / 3)
check(abs(median(nn3d) / med3_cf - 1) < 0.05,
      sprintf("3D median NN: sim %.2f vs closed form %.2f (ratio %.3f)",
              median(nn3d), med3_cf, median(nn3d) / med3_cf))

# The headline inflation, from the count-based intensity.
inf <- t3_median_d3d(2.84, 18, d2_observed_um = 304.6)
check(inf$inflation > 2.5 && inf$inflation < 3.5,
      sprintf("headline inflation at lambda2 = 2.84/mm2, h_eff = 18: %.2fx",
              inf$inflation))


## ===========================================================================
hdr("3. PROJECTION -- realised detection fraction and anchor intensity")
## ===========================================================================

if (!requireNamespace("spatstat.geom", quietly = TRUE)) {
  check(FALSE, "spatstat.geom available (Stage 1 requires it)")
} else {
  w      <- spatstat.geom::owin(c(0, 4), c(0, 4))          # 16 mm^2
  area_w <- spatstat.geom::area(w)
  geom   <- t3_geometry()
  lam2_t <- 2.84                                            # target, /mm^2
  n_det_t <- round(lam2_t * area_w)
  padz   <- t3_pad_z(lam2_t, geom$h_eff, sigma_z = 150)
  n_tg   <- t3_n_tangle_for(n_det_t, padz, geom$h_eff)

  cat(sprintf("  target %d detected anchors over %.1f mm2; pad_z = %.0f um -> n_tangle = %d\n",
              n_det_t, area_w, padz, n_tg))

  det_frac <- numeric(0); det_n <- numeric(0)
  for (r in seq_len(min(args$reps, 15L))) {
    sim <- t3_simulate_3d(w, n_neuron = max(n_tg * 3L, 5000L), n_tangle = n_tg,
                          model = "poisson", sigma_xy = 150, sigma_z = 150,
                          pad_z_um = padz, seed = args$seed + r)
    pr  <- t3_project(sim, geom)
    det_frac <- c(det_frac, pr$p_detect_realised)
    det_n    <- c(det_n, pr$n_detected)
  }
  check(abs(mean(det_frac) - geom$p_detect) < 0.03,
        sprintf("realised detection fraction %.3f vs p_detect %.3f (over %d reps)",
                mean(det_frac), geom$p_detect, length(det_frac)))
  check(abs(mean(det_n) / n_det_t - 1) < 0.12,
        sprintf("realised detected count %.0f vs target %d (t3_n_tangle_for round-trip)",
                mean(det_n), n_det_t))

  # ---- clustering: assert the mechanism, MEASURE the consequence ------------
  #
  # The assertion is on the ANCHOR-TO-ANCHOR nearest-neighbour distance, which
  # is unambiguous: a clustered anchor pattern has a shorter one. If it does
  # not, the thinning weight is doing nothing and the "scenario sweep" is a
  # no-op dressed as a sensitivity.
  #
  # The median QUERY distance is only reported, never asserted, because its
  # sign is not a property of the code. It depends on whether the query cells
  # share the anchors' clustering, and both couplings are physically arguable:
  # CBLN2 neurons sit in pre-alpha islands (shared), but tangles may also
  # cluster beyond the neuron process (not shared). The table below is printed
  # so the direction is on record for whichever coupling the calibrated model
  # ends up using.
  cl <- NULL
  for (nm in c("uniform", "clustered")) for (tm in c("poisson", "thomas")) {
    ann <- medq <- numeric(0)
    for (r in seq_len(min(args$reps, 8L))) {
      s <- t3_simulate_3d(w, n_neuron = max(n_tg * 4L, 20000L), n_tangle = n_tg,
                          model = tm, sigma_xy = 150, sigma_z = 150, mu = 100,
                          neuron_model = nm, neuron_mu = 200,
                          pad_z_um = padz, seed = 100 + r)
      qq <- s$neurons[abs(s$neurons$z_um) < geom$h_cell / 2, c("x_mm", "y_mm")]
      pr <- t3_project(s, geom, query_xy_mm = qq)
      a  <- pr$tangles[pr$tangles$in_slab_incl, c("x_mm", "y_mm")]
      if (nrow(a) > 5)
        ann <- c(ann, median(RANN::nn2(a, a, k = 2)$nn.dists[, 2] * 1000))
      if (length(pr$d2_um) > 10) medq <- c(medq, median(pr$d2_um))
    }
    cl <- rbind(cl, data.frame(neuron_model = nm, tangle_model = tm,
                               anchor_nn_um = mean(ann), median_query_um = mean(medq)))
  }
  cat("\n  Clustering coupling (query = simulated neurons, as in the real analysis):\n")
  print(cl, row.names = FALSE, digits = 3)

  for (nm in c("uniform", "clustered")) {
    a_p <- cl$anchor_nn_um[cl$neuron_model == nm & cl$tangle_model == "poisson"]
    a_t <- cl$anchor_nn_um[cl$neuron_model == nm & cl$tangle_model == "thomas"]
    check(a_t < 0.9 * a_p,
          sprintf("[%s neurons] Thomas genuinely clusters anchors: anchor NN %.0f -> %.0f um",
                  nm, a_p, a_t))
  }

  q_uni <- cl$median_query_um[cl$neuron_model == "uniform"  & cl$tangle_model == "poisson"]
  q_clu <- cl$median_query_um[cl$neuron_model == "clustered" & cl$tangle_model == "poisson"]
  check(q_clu < q_uni,
        sprintf("shared (neuron) clustering SHORTENS the query distance: %.0f -> %.0f um",
                q_uni, q_clu))
  cat(sprintf("\n  NOTE: tangle clustering BEYOND the neuron process LENGTHENS it\n"))
  cat(sprintf("        (uniform neurons %.0f -> %.0f um; clustered neurons %.0f -> %.0f um).\n",
              q_uni, cl$median_query_um[cl$neuron_model == "uniform"   & cl$tangle_model == "thomas"],
              q_clu, cl$median_query_um[cl$neuron_model == "clustered" & cl$tangle_model == "thomas"]))
  cat("        So \"clustering shortens d3D, therefore Poisson is conservative\" is\n")
  cat("        conditional, not general. The observed d2D distribution is what\n")
  cat("        pins down which scenarios are admissible -- see t3_validate().\n")
}


## ===========================================================================
hdr("4. ARM B DILUTION CONSTANT -- the gate")
## ===========================================================================
#
# Plant a known log-distance field, delete a known fraction of anchors, add the
# same number back by intensity-weighted draw, and measure what happens to the
# fitted slope. Three slopes per replicate:
#
#   slope_true  : distance from the FULL anchor set          (ground truth)
#   slope_obs   : distance from the THINNED anchor set       (what we measure)
#   slope_aug   : distance from thinned + weighted add-back  (what Arm B does)
#
# The reported constant is slope_aug / slope_obs. If it is below 1, augmentation
# attenuates the slope on its own and every Arm B result must be read against
# that number.

sim_field <- function(n_cell = 6000, n_anchor = 120, L = 4, beta = -0.30,
                      p_keep = 0.75, sigma_um = 200, seed = 1) {
  set.seed(seed)
  cx <- runif(n_cell, 0, L); cy <- runif(n_cell, 0, L)
  ax <- runif(n_anchor, 0, L); ay <- runif(n_anchor, 0, L)

  d_full <- RANN::nn2(cbind(ax, ay), cbind(cx, cy), k = 1)$nn.dists[, 1] * 1000
  d_full[d_full <= 0] <- 0.1
  # Outcome generated from the TRUE distance.
  y <- beta * (log(d_full) / sd(log(d_full))) + rnorm(n_cell, 0, 1)

  keep <- sort(sample.int(n_anchor, round(n_anchor * p_keep)))
  d_obs <- RANN::nn2(cbind(ax, ay)[keep, , drop = FALSE], cbind(cx, cy),
                     k = 1)$nn.dists[, 1] * 1000
  d_obs[d_obs <= 0] <- 0.1

  # Add back H cells drawn from the CELL pool, weighted by kernel-smoothed
  # observed-anchor intensity -- exactly what draw_hidden_anchors() will do.
  H <- n_anchor - length(keep)
  kx <- cbind(ax, ay)[keep, , drop = FALSE]
  dk <- RANN::nn2(kx, cbind(cx, cy), k = min(10L, nrow(kx)))$nn.dists * 1000
  wt <- rowSums(exp(-0.5 * (dk / sigma_um)^2))
  wt[!is.finite(wt) | wt < 0] <- 0
  if (all(wt == 0)) wt <- rep(1, n_cell)
  hid <- sample.int(n_cell, H, prob = wt)

  aug <- rbind(kx, cbind(cx, cy)[hid, , drop = FALSE])
  # k = 2 self-match handling, as compute_dist_to_phf1_um() does on the
  # override path: a hidden anchor drawn from the query pool would otherwise
  # return a spurious zero for itself.
  nnA <- RANN::nn2(aug, cbind(cx, cy), k = 2)
  d_aug <- nnA$nn.dists[, 1] * 1000
  self  <- nnA$nn.dists[, 1] < 1e-12
  d_aug[self] <- nnA$nn.dists[self, 2] * 1000
  d_aug[d_aug <= 0] <- 0.1

  sl <- function(d) unname(coef(lm(y ~ I(log(d) / sd(log(d)))))[2])
  c(true = sl(d_full), obs = sl(d_obs), aug = sl(d_aug),
    med_obs = median(d_obs), med_aug = median(d_aug))
}

res <- t(vapply(seq_len(args$reps), function(r)
  sim_field(seed = args$seed * 1000 + r), numeric(5)))

ratio_aug_obs  <- res[, "aug"]  / res[, "obs"]
ratio_obs_true <- res[, "obs"]  / res[, "true"]
DILUTION <- mean(ratio_aug_obs)
dil_ci <- DILUTION + c(-1.96, 1.96) * sd(ratio_aug_obs) / sqrt(nrow(res))

cat(sprintf("\n  slope_true      %+.4f\n", mean(res[, "true"])))
cat(sprintf("  slope_obs       %+.4f   (thinned anchors; obs/true = %.3f)\n",
            mean(res[, "obs"]), mean(ratio_obs_true)))
cat(sprintf("  slope_aug       %+.4f   (augmented)\n", mean(res[, "aug"])))
cat(sprintf("  median d2D      %.0f um -> %.0f um on augmentation\n",
            mean(res[, "med_obs"]), mean(res[, "med_aug"])))
cat(sprintf("\n  *** ARM B DILUTION CONSTANT  slope_aug/slope_obs = %.3f [%.3f, %.3f] ***\n",
            DILUTION, dil_ci[1], dil_ci[2]))
cat("  Arm B slopes are compared against this number: a weakening of this size\n")
cat("  is produced by the PROCEDURE itself.\n\n")

check(is.finite(DILUTION) && DILUTION > 0,
      sprintf("dilution constant is finite and positive (%.3f)", DILUTION))
check(mean(res[, "med_aug"]) < mean(res[, "med_obs"]),
      sprintf("augmentation shortens the median distance (%.0f -> %.0f um)",
              mean(res[, "med_obs"]), mean(res[, "med_aug"])))
# Augmentation is not expected to restore the true slope: Arm B is a robustness
# check, not an estimator.
check(abs(mean(res[, "aug"]) - mean(res[, "true"])) >
      0.2 * abs(mean(res[, "true"]) - mean(res[, "obs"])) ||
      DILUTION < 1.05,
      "augmentation does not recover the true slope (it is a robustness check, not an estimator)")

# Write the constant where the driver can read it rather than re-deriving it.
const_path <- file.path(dirname(.this_dir), "results", "tangle_3d")
dir.create(const_path, recursive = TRUE, showWarnings = FALSE)
write.table(
  data.frame(quantity = c("dilution_slope_aug_over_obs", "dilution_ci_lo",
                          "dilution_ci_hi", "reps", "p_keep", "seed"),
             value = c(DILUTION, dil_ci[1], dil_ci[2], args$reps, 0.75, args$seed)),
  file.path(const_path, "armB_dilution_constant.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
cat("  written:", file.path(const_path, "armB_dilution_constant.tsv"), "\n")


## ===========================================================================
hdr("5. GUARDS")
## ===========================================================================

check(inherits(tryCatch(t3_geometry(7, 20, 17), error = function(e) e), "error"),
      "t3_geometry() stops when the inclusion is larger than the soma")
check(inherits(tryCatch(t3_geometry(0, 11, 17), error = function(e) e), "error"),
      "t3_geometry() stops on a non-positive thickness")
check(inherits(tryCatch(t3_hidden_count(100, 0), error = function(e) e), "error"),
      "t3_hidden_count() stops at p_detect = 0")
check(t3_hidden_count(398, 1) == 0L,
      "t3_hidden_count() returns 0 at p_detect = 1 (the identity arm)")
check(t3_hidden_count(398, 0.75) == 133L,
      sprintf("t3_hidden_count(398, 0.75) = %d (expected 133)",
              t3_hidden_count(398, 0.75)))

lam <- t3_lambda2_from_counts(c(4, 120), c(13.9, 17.8), c("a", "b"))
check(nrow(lam) == 3 && lam$sample_id[3] == "POOLED",
      "t3_lambda2_from_counts() appends a POOLED row")
check(abs(lam$lambda2_per_mm2[3] - 124 / 31.7) < 1e-9,
      "POOLED intensity is total anchors / total area, not a mean of ratios")


## ===========================================================================
hdr("RESULT")
## ===========================================================================
if (length(FAILURES)) {
  cat("\nFAILED", length(FAILURES), "check(s):\n")
  for (f in FAILURES) cat("  - ", f, "\n", sep = "")
  cat("\nDo NOT run this on real data until it passes.",
      "Do not loosen a tolerance to make it pass.\n")
  quit(save = "no", status = 1)
}
cat("\nAll checks passed.\n")
cat(sprintf("Arm B dilution constant: %.3f [%.3f, %.3f]\n",
            DILUTION, dil_ci[1], dil_ci[2]))
cat("\n"); print(sessionInfo())
