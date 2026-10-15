#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# imc_tangle_depth_density_shift.R
#
# Figure panels: 4C (default, excitatory), 4D (`inh` argument)
# 4C/4D are the PHF1 intensity panels (plot_density_shift_vs_phf1_intensity.pdf). The annotation,
# plaque and distance panels are written alongside them.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# imc_tangle_depth_density_shift.R
#
# IMC (EC cohort), EXCITATORY neurons (default) or INHIBITORY neurons (`inh`). Each ROI is cut into
# 10 relative cortical-depth deciles (0 = pia; scaled_Y over ALL cells of the ROI). For every decile:
#
#   x = % difference in columnar density of excitatory neurons, AD vs control -- the x axis of
#       plot_cluster_loss_vs_tangle_distance_raw.pdf (imc_cluster_loss_vs_tangle_distance.R, computed
#       as in the published IMC analysis code, 12b-Cluster_Percentage_Change.Rmd:39-103), with the
#       decile in place of the cluster:
#         columnar density = exc neurons in the decile / ROI width (per mm), per ROI;
#         donor value      = mean over the donor's 3 ROIs;
#         x                = 100 * (mean AD donor value / mean control donor value - 1),
#         AD = Braak_3_4 + Braak_6 (Braak >= 2, the Rmd cutoff), control = Braak_0_1.
#   y = the decile's PHF1 load in AD donors, pooled across the cohort (not a mean of donors):
#         annotation -- % of exc neurons that are PHF1_Otsu+;
#         intensity  -- mean within-donor z-scored PHF1 intensity of ALL neurons (exc + inh), which
#                       does NOT use the PHF1+/- call. Assay "exprs" = asinh(counts/1), z-scored
#                       within each donor over that donor's neurons (as imc_cluster_loss_vs_phf1_intensity.R),
#                       so the decile ranking is the same in the exc and inh runs;
#         plaque     -- % of exc neurons inside the 40 um-dilated 4G8+ plaque mask;
#         distance   -- geometric-mean distance from PHF1- exc neurons to the nearest PHF1+ neuron.
#   Points are coloured by the decile's rank on that panel's y (rank 1 = most PHF1 = darkest), so the
#   panels rank deciles independently.
#
# Uncertainty: donor-level STRATIFIED bootstrap (control and AD donors resampled separately), both
# axes recomputed from the same resampled donors -- as in imc_cluster_loss_vs_tangle_distance.R.
# Percentile 95% CIs. Association across deciles -- that script's test, copied verbatim
# (spear / boot_p / perm_p, B = 9999, B_PERM = 199999):
#   Spearman rho, percentile 95% CI from the stratified bootstrap, CI-inversion bootstrap p, and a
#   decile-label permutation p (Monte Carlo; k = 10 > 8 so not enumerated).
#
# ROI x decile cells with 0 exc neurons are counted as 0 (for a depth band a zero is a real
# observation); the drop-zero version is in the stats log as a sensitivity. The per-donor x decile
# table is written to source_data_density_shift_donor.tsv.
#
# Outputs: plots/imc_tangle_depth_density_shift/ (exc) or plots/imc_tangle_depth_density_shift_inh/,
# triple per panel.
#   Rscript phf1_v2/R/imc_tangle_depth_density_shift.R          # excitatory (default)
#   Rscript phf1_v2/R/imc_tangle_depth_density_shift.R inh      # inhibitory -> ..._inh/
# The inhibitory run is identical with "inhibitory" for "excitatory" throughout: its own columnar
# density, PHF1- inhibitory targets for the distance panel (anchors stay PHF1+ cells of any neuron
# cluster); the intensity y is the same all-neuron measure in both runs.

suppressPackageStartupMessages({
  library(SpatialExperiment); library(dplyr); library(tibble); library(tidyr); library(ggplot2)
})

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")
source("R/imc_tangle_depth_deciles_plot_utils.R")   # plot_density_shift()

CLASS <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(CLASS)) CLASS <- "exc"
stopifnot(CLASS %in% c("exc", "inh"))
NOUN     <- c(exc = "excitatory", inh = "inhibitory")[[CLASS]]
NOUN_CAP <- c(exc = "Excitatory", inh = "Inhibitory")[[CLASS]]
out_dir <- if (CLASS == "exc") "plots/imc_tangle_depth_density_shift" else "plots/imc_tangle_depth_density_shift_inh"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE, row.names = FALSE)

SEED   <- 20261013
B      <- 9999
B_PERM <- 199999   # decile-label permutations, as imc_cluster_loss_vs_tangle_distance.R (k = 10; Monte Carlo)
N_DEC  <- 10
EXC    <- paste(NOUN_CAP, "Neurons")   # the analysed class (exc or inh)
DEC_LEVELS <- sprintf("D%02d", seq_len(N_DEC))
DEC_LABELS <- setNames(sprintf("%d-%d%%", (seq_len(N_DEC) - 1) * 10, seq_len(N_DEC) * 10), DEC_LEVELS)

# -------------------------------------------------------------------
# Load + cohort (donors with all 3 ROIs)
# -------------------------------------------------------------------
localwd <- "<IMC_ROOT>/"
spe <- readRDS(paste0(localwd, "spe.rds"))
cd <- as.data.frame(colData(spe)); cd$cell_id <- colnames(spe)
sc <- spatialCoords(spe); cd$x <- sc[, "Pos_X"]; cd$y <- sc[, "Pos_Y"]
cd$phf1_raw <- as.numeric(assay(spe, "exprs")["PHF1", ])
rm(spe); invisible(gc())

cv_removal <- cd %>% distinct(patient_id, sample_id) %>% count(patient_id) %>%
  filter(n < 3) %>% pull(patient_id)
cells <- cd %>%
  filter(!patient_id %in% cv_removal) %>%
  mutate(patient_id = as.character(patient_id), sample_id = as.character(sample_id),
         BraakGroup = as.character(BraakGroup), layer = as.character(layer)) %>%
  group_by(sample_id) %>%
  mutate(scaled_Y = 100 * (y - min(y)) / (max(y) - min(y))) %>%
  ungroup() %>%
  mutate(dec_i  = pmin(pmax(ceiling(round(scaled_Y, 9) / 10), 1L), N_DEC),
         decile = factor(DEC_LEVELS[dec_i], levels = DEC_LEVELS),
         ad     = BraakGroup != "Braak_0_1")
stopifnot(sum(cells$decile == "D01") == sum(cells$layer %in% "1"),
          all((cells$decile == "D01") == (cells$layer %in% "1")))

roi <- cells %>% distinct(sample_id, patient_id, BraakGroup, Braak, ad, width_px)
stopifnot(!any(duplicated(roi$sample_id)), all(table(roi$patient_id) == 3),
          all(roi$ad == (as.numeric(as.character(roi$Braak)) >= 2)))

exc <- cells %>% filter(broad_celltypes == EXC)
stopifnot(all(grepl(paste0("^", NOUN_CAP, " neuron"), exc$celltype_clusters)), all(is.finite(exc$phf1_raw)))
# within-donor z over the donor's excitatory neurons (house convention, analysed population = exc)
exc <- exc %>% group_by(patient_id) %>%
  mutate(phf1_z = as.numeric(scale(phf1_raw))) %>% ungroup() %>%
  mutate(phf1_pos = PHF1_Otsu == "PHF1_pos",
         plq = Matched_4G8_40)   # nucleus inside the 40 um-dilated 4G8+ amyloid mask (logical)
stopifnot(is.logical(exc$plq), !anyNA(exc$plq))
zchk <- exc %>% group_by(patient_id) %>% summarise(m = mean(phf1_z), s = sd(phf1_z), .groups = "drop")
stopifnot(max(abs(zchk$m)) < 1e-9, max(abs(zchk$s - 1)) < 1e-9)
# the Otsu call is made on this channel, so PHF1+ exc neurons must be brighter in every donor with any
dz <- exc %>% group_by(patient_id) %>%
  summarise(n_pos = sum(phf1_pos), d = mean(phf1_z[phf1_pos]) - mean(phf1_z[!phf1_pos]), .groups = "drop")
stopifnot(all(dz$d[dz$n_pos > 0] > 0))

# Intensity panel y: PHF1 intensity across ALL NEURONS (exc + inh), z-scored within
# donor over the donor's neurons -- the imc_cluster_loss_vs_phf1_intensity.R population -- so the
# decile ranking is the SAME in the exc and inh runs; only the x (density) is class-specific.
neu <- cells %>% filter(broad_celltypes %in% c("Excitatory Neurons", "Inhibitory Neurons")) %>%
  group_by(patient_id) %>% mutate(phf1_zn = as.numeric(scale(phf1_raw))) %>% ungroup()
stopifnot(all(grepl("neuron", neu$celltype_clusters, ignore.case = TRUE)), all(is.finite(neu$phf1_zn)))
zn_chk <- neu %>% group_by(patient_id) %>% summarise(m = mean(phf1_zn), s = sd(phf1_zn), .groups = "drop")
stopifnot(max(abs(zn_chk$m)) < 1e-9, max(abs(zn_chk$s - 1)) < 1e-9)

# -------------------------------------------------------------------
# Donor x decile matrices
# -------------------------------------------------------------------
roi_dec <- exc %>%
  group_by(sample_id, decile) %>%
  summarise(n = n(), n_pos = sum(phf1_pos), z_sum = sum(phf1_z), n_plq = sum(plq), .groups = "drop") %>%
  right_join(tidyr::expand_grid(sample_id = roi$sample_id,
                                decile = factor(DEC_LEVELS, levels = DEC_LEVELS)),
             by = c("sample_id", "decile")) %>%
  mutate(across(c(n, n_pos, z_sum, n_plq), ~ coalesce(.x, 0))) %>%
  left_join(roi, by = "sample_id") %>%
  mutate(col_dens = n / width_px * 1e3)                  # exc neurons per mm of ROI width
n_zero_roi_dec <- sum(roi_dec$n == 0)

don <- roi_dec %>% group_by(patient_id, ad, decile) %>%
  # col_dens_nz first: summarise() overwrites col_dens with its mean for every later expression
  summarise(col_dens_nz = if (any(n > 0)) mean(col_dens[n > 0]) else NA_real_,
            col_dens = mean(col_dens),
            n = sum(n), n_pos = sum(n_pos), z_sum = sum(z_sum), n_plq = sum(n_plq), .groups = "drop")
wide <- function(col, keep) {
  m <- don %>% filter(ad == keep) %>% select(patient_id, decile, all_of(col)) %>%
    pivot_wider(names_from = decile, values_from = all_of(col)) %>% arrange(patient_id)
  as.matrix(m[, DEC_LEVELS])
}
Dc <- wide("col_dens", FALSE); Da <- wide("col_dens", TRUE)
Dc_nz <- wide("col_dens_nz", FALSE); Da_nz <- wide("col_dens_nz", TRUE)
N  <- wide("n", TRUE); P <- wide("n_pos", TRUE); Z <- wide("z_sum", TRUE); Q <- wide("n_plq", TRUE)
ad_order <- don %>% filter(ad) %>% distinct(patient_id) %>% arrange(patient_id) %>% pull(patient_id)
neu_w <- function(col) {
  m <- neu %>% filter(ad) %>% group_by(patient_id, decile) %>%
    summarise(n = n(), z_sum = sum(phf1_zn), .groups = "drop") %>%
    select(patient_id, decile, all_of(col)) %>%
    complete(patient_id = ad_order, decile = factor(DEC_LEVELS, levels = DEC_LEVELS), fill = setNames(list(0), col)) %>%
    pivot_wider(names_from = decile, values_from = all_of(col)) %>% arrange(patient_id)
  stopifnot(identical(m$patient_id, ad_order))   # same donor rows as N, so bootstrap indices line up
  as.matrix(m[, DEC_LEVELS])
}
NN <- neu_w("n"); ZN <- neu_w("z_sum")
stopifnot(nrow(Dc) == 18, nrow(Da) == 44)

pct_change <- function(Dc, Da) 100 * (colMeans(Da, na.rm = TRUE) / colMeans(Dc, na.rm = TRUE) - 1)
y_ann <- function(i) 100 * colSums(P[i, , drop = FALSE]) / colSums(N[i, , drop = FALSE])
y_int <- function(i) colSums(ZN[i, , drop = FALSE]) / colSums(NN[i, , drop = FALSE])   # all neurons
y_plq <- function(i) 100 * colSums(Q[i, , drop = FALSE]) / colSums(N[i, , drop = FALSE])

# distance to the nearest PHF1+ neuron -- the y of plot_cluster_loss_vs_tangle_distance_raw.pdf,
# constructed exactly as imc_cluster_loss_vs_tangle_distance.R (decile in place of cluster):
# anchors = PHF1_Otsu+ cells of ANY neuron cluster (exc + inh), targets = PHF1- exc neurons, nearest
# anchor within the ROI (NA if the ROI has none), donor x decile mean of natural-log distance where
# >= MIN_CELLS targets, y = exp(mean over AD donors). Uncapped, as imc_cluster_loss_vs_tangle_distance.R.
MIN_CELLS <- 5
anc <- cells %>% filter(ad, !celltype_clusters %in% c("Artefact cluster", "Unassigned cluster"),
                        PHF1_Otsu == "PHF1_pos", grepl("neuron", celltype_clusters, ignore.case = TRUE))
tgt <- exc %>% filter(ad, !phf1_pos)
tgt$d <- NA_real_
for (s_ in unique(tgt$sample_id)) {
  i <- which(tgt$sample_id == s_); j <- which(anc$sample_id == s_)
  if (length(j)) tgt$d[i] <- RANN::nn2(as.matrix(anc[j, c("x", "y")]), as.matrix(tgt[i, c("x", "y")]),
                                         k = 1)$nn.dists[, 1]
}
# GATE: the cohort of imc_cluster_loss_vs_tangle_distance.R (its gate 2)
stopifnot(nrow(anc) == 2467, CLASS != "exc" || nrow(tgt) == 34563,   # 34563 = the exc targets of imc_cluster_loss_vs_tangle_distance.R
          all(tgt$d > 0, na.rm = TRUE))                                # log(), never log1p
L <- tgt %>% filter(!is.na(d)) %>% group_by(patient_id, decile) %>%
  summarise(n_t = n(), ld = mean(log(d)), .groups = "drop") %>%
  mutate(ld = ifelse(n_t >= MIN_CELLS, ld, NA_real_)) %>%
  select(patient_id, decile, ld) %>%
  complete(patient_id = sort(unique(don$patient_id[don$ad])), decile = DEC_LEVELS) %>%
  pivot_wider(names_from = decile, values_from = ld) %>% arrange(patient_id)
stopifnot(identical(L$patient_id, sort(unique(don$patient_id[don$ad]))))   # same row order as N
L <- as.matrix(L[, DEC_LEVELS])
y_dst <- function(i) exp(colMeans(L[i, , drop = FALSE], na.rm = TRUE))

obs_x <- pct_change(Dc, Da)
obs_x_nz <- pct_change(Dc_nz, Da_nz)
obs <- list(annotation = y_ann(seq_len(nrow(N))), intensity = y_int(seq_len(nrow(N))),
            plaque = y_plq(seq_len(nrow(N))), distance = y_dst(seq_len(nrow(N))))

set.seed(SEED)
# adding a panel draws no extra random numbers, so earlier panels reproduce exactly
bx <- matrix(NA_real_, B, N_DEC); by <- list(annotation = bx, intensity = bx, plaque = bx, distance = bx)
rho_b <- matrix(NA_real_, B, length(obs), dimnames = list(NULL, names(obs)))
for (b in seq_len(B)) {
  ic <- sample.int(nrow(Dc), replace = TRUE); ia <- sample.int(nrow(Da), replace = TRUE)
  bx[b, ] <- pct_change(Dc[ic, , drop = FALSE], Da[ia, , drop = FALSE])
  by$annotation[b, ] <- y_ann(ia); by$intensity[b, ] <- y_int(ia); by$plaque[b, ] <- y_plq(ia); by$distance[b, ] <- y_dst(ia)
  for (k in names(obs)) rho_b[b, k] <- cor(bx[b, ], by[[k]][b, ], method = "spearman")
}
q <- function(m, p) apply(m, 2, quantile, p, na.rm = TRUE)

# verbatim from imc_cluster_loss_vs_tangle_distance.R (spear, boot_p, perm_p)
spear <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))
boot_p <- function(r) { r <- r[!is.na(r)]; min(1, 2 * min((sum(r <= 0) + 1), (sum(r >= 0) + 1)) / (length(r) + 1)) }
perm_p <- function(x, y, n_perm) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; r0 <- spear(x, y); k <- length(x)
  if (k <= 8) {                                    # exact enumeration
    perms <- function(v) if (length(v) <= 1) list(v) else
      do.call(c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(p) c(v[i], p))))
    r <- sapply(perms(seq_len(k)), function(p) spear(x, y[p]))
    return(c(p = mean(abs(r) >= abs(r0) - 1e-12), n = length(r)))
  }
  rx <- rank(x); ry <- rank(y)
  r <- vapply(seq_len(n_perm), function(i) cor(rx, sample(ry)), 0)
  c(p = (sum(abs(r) >= abs(r0) - 1e-12) + 1) / (n_perm + 1), n = n_perm)
}

# per-donor x decile table (cf. the _donor.tsv of imc_cluster_loss_vs_tangle_distance.R)
wr(don %>% mutate(group = ifelse(ad, "AD", "Control"),
                  pct_phf1_pos = ifelse(n > 0, 100 * n_pos / n, NA_real_),
                  phf1_z_mean = ifelse(n > 0, z_sum / n, NA_real_),
                  pct_plq = ifelse(n > 0, 100 * n_plq / n, NA_real_),
                  depth = DEC_LABELS[as.character(decile)]) %>%
     left_join(roi %>% distinct(patient_id, BraakGroup), by = "patient_id") %>%
     select(patient_id, group, BraakGroup, decile, depth, columnar_density_per_mm = col_dens,
            n_exc = n, n_phf1_pos = n_pos, pct_phf1_pos, phf1_z_mean,
            n_plaque_40um = n_plq, pct_plaque_40um = pct_plq),
   "source_data_density_shift_donor.tsv")

# -------------------------------------------------------------------
# Panels
# -------------------------------------------------------------------
PANELS <- list(
  annotation = list(ylab = paste0("PHF1+ ", NOUN, " neurons in AD (%)"), legend = "Cortical depth decile,\nranked by % PHF1+",
                    hline = NULL, y_name = "pct_phf1_pos_ad",
                    what = "% of excitatory neurons that are PHF1_Otsu+ (pooled AD donors)",
                    file = "phf1_annotation", short = "PHF1"),
  intensity  = list(ylab = "PHF1 intensity (within-donor z)", legend = "Cortical depth decile,\nranked by PHF1 intensity",
                    hline = 0, y_name = "phf1_z_ad",
                    what = "mean within-donor z-scored PHF1 intensity of ALL neurons (exc + inh; pooled AD donors)",
                    file = "phf1_intensity", short = "PHF1"),
  plaque     = list(# two lines: one line overruns the 6 cm panel height
                    ylab = bquote(atop(.(paste(NOUN_CAP, "neurons within")), "40 " * mu * "m of a plaque (%)")),
                    legend = "Cortical depth decile,\nranked by % near plaque",
                    hline = NULL, y_name = "pct_plaque_40um_ad",
                    what = paste("% of excitatory neurons whose nucleus lies inside the 40 um-dilated 4G8+",
                                 "amyloid plaque mask (Matched_4G8_40; pooled AD donors)"),
                    file = "plaque", short = "plaque proximity"),
  distance   = list(ylab = expression("Distance to nearest PHF1+ neuron (" * mu * "m)"),
                    legend = "Cortical depth decile,\nranked by distance\nto nearest PHF1+ neuron",
                    hline = NULL, y_name = "dist_phf1_um_ad", log_y = TRUE, rank_desc = FALSE,
                    tags = c("closest", "farthest"),
                    what = paste("geometric-mean distance from PHF1- exc neurons to the nearest PHF1+ neuron",
                                 "(exc or inh), exp(mean over AD donors of donor-mean log distance), as the",
                                 "y of plot_cluster_loss_vs_tangle_distance_raw.pdf"),
                    file = "distance", short = "tangle proximity"))

ledger <- function() {
  cat("=== Cohort ===\n")
  cat(sprintf("Excitatory neurons, 3-ROI rule: %d control donors (Braak_0_1, %d ROIs), %d AD donors\n",
              nrow(Dc), 3 * nrow(Dc), nrow(Da)))
  cat(sprintf("(Braak_3_4 + Braak_6, %d ROIs). %d exc neurons in all, %d in AD, %d PHF1+ in AD.\n",
              3 * nrow(Da), nrow(exc), sum(N), sum(P)))
  cat(sprintf("ROI x decile cells with 0 exc neurons: %d of %d (counted as 0).\n\n",
              n_zero_roi_dec, nrow(roi_dec)))
}

# the log text is written for excitatory neurons; for the inhibitory run swap the class words only
relabel_log <- function(f) {
  x <- readLines(f)
  for (r in list(c("excitatory", "inhibitory"), c("Excitatory", "Inhibitory"), c("exc neurons", "inh neurons"),
                 c("exc columnar", "inh columnar"), c("AD exc", "AD inh")))
    x <- gsub(r[1], r[2], x, fixed = TRUE)
  writeLines(x, f)
}
res <- list(); tests <- list()
for (k in names(PANELS)) {
  pn <- PANELS[[k]]
  d <- tibble(decile = DEC_LEVELS, depth = DEC_LABELS[DEC_LEVELS],
              pct_change = obs_x, pct_lo = q(bx, 0.025), pct_hi = q(bx, 0.975),
              y = obs[[k]], y_lo = q(by[[k]], 0.025), y_hi = q(by[[k]], 0.975),
              rank = rank(if (isFALSE(pn$rank_desc)) obs[[k]] else -obs[[k]], ties.method = "first"),
              ctrl_col_dens_per_mm = colMeans(Dc), ad_col_dens_per_mm = colMeans(Da),
              n_ad_exc = colSums(N), n_ad_phf1_pos = colSums(P),
              pct_change_dropzero = obs_x_nz) %>%
    arrange(rank)
  res[[k]] <- d
  rho <- spear(d$pct_change, d$y)
  rho_ci <- quantile(rho_b[, k], c(0.025, 0.975), na.rm = TRUE)
  set.seed(SEED + 1)
  pp <- perm_p(d$pct_change, d$y, B_PERM)
  tests[[k]] <- tibble(panel = k, n_deciles = nrow(d), rho = rho, ci_lo = rho_ci[[1]], ci_hi = rho_ci[[2]],
                       p_boot = boot_p(rho_b[, k]), p_perm = unname(pp["p"]), n_perm = unname(pp["n"]),
                       frac_boot_neg = mean(rho_b[, k] < 0, na.rm = TRUE),
                       n_boot_na = sum(is.na(rho_b[, k])))

  ggsave(file.path(out_dir, paste0("plot_density_shift_vs_", pn$file, ".pdf")),
         plot_density_shift(d, pn$ylab, pn$legend, pn$hline, log_y = isTRUE(pn$log_y),
                            xlab = paste0("Difference in ", NOUN, " neuron\ncolumnar density, AD vs control (%)"),
                            tags = if (is.null(pn$tags)) c("highest", "lowest") else pn$tags),
         width = 9.5, height = 6, units = "cm", device = "pdf")   # = the size of the plot_cluster_loss_vs_tangle_distance_raw.pdf panel
  wr(d %>% rename(!!pn$y_name := y, !!paste0(pn$y_name, "_lo") := y_lo,
                  !!paste0(pn$y_name, "_hi") := y_hi),
     paste0("source_data_density_shift_vs_", pn$file, ".tsv"))

  sink(file.path(out_dir, paste0("stats_density_shift_vs_", pn$file, ".txt")))
  cat("Depth-decile excitatory neuron loss vs", pn$short, "--", k, "\n\n")
  ledger()
  cat("=== Definitions ===\n")
  cat("Decile: ceiling(scaled_Y/10), scaled_Y = 100*(y-min y)/(max y-min y) per ROI over ALL cells,\n")
  cat("0 = pia (D01 == stored layer '1', asserted).\n")
  cat("x: columnar density = exc neurons in decile / ROI width * 1000 (per mm), mean over the donor's\n")
  cat("   3 ROIs; x = 100*(mean AD donor / mean control donor - 1). Same quantity as the x axis of\n")
  cat("   plot_cluster_loss_vs_tangle_distance_raw.pdf, per decile instead of per cluster.\n")
  cat("y: ", pn$what, ".\n", sep = "")
  if (k == "intensity")
    cat("   PHF1 = assay exprs (asinh(counts/1)), z-scored within each donor over all of that donor's\n",
        "   neurons, exc + inh (asserted mean 0 / SD 1 per donor). Does not use the PHF1+/- annotation.\n",
        "   Same y and decile ranking in the exc and inh runs; only the x (density) is class-specific.\n",
        sep = "")
  cat("Colour/legend: each point is a cortical depth decile, ranked by its y (1 = highest, darkest).\n")
  cat("CIs: donor-level stratified bootstrap, B = ", B, ", control and AD donors resampled separately,\n",
      "both axes from the same resampled donors; percentile 95%.\n\n", sep = "")
  cat("=== Per decile (effect sizes with 95% CI), ordered by rank on y ===\n")
  print(as.data.frame(d %>% transmute(rank, decile, depth,
                                      x_pct = round(pct_change, 1),
                                      x_CI = sprintf("[%.1f, %.1f]", pct_lo, pct_hi),
                                      y = round(y, 3), y_CI = sprintf("[%.3f, %.3f]", y_lo, y_hi),
                                      ctrl_per_mm = round(ctrl_col_dens_per_mm, 1),
                                      ad_per_mm = round(ad_col_dens_per_mm, 1),
                                      x_dropzero = round(pct_change_dropzero, 1))), row.names = FALSE)
  cat(sprintf("\nWhole-ROI exc columnar density: control %.1f, AD %.1f per mm (%+.1f%%).\n",
              sum(colMeans(Dc)), sum(colMeans(Da)), 100 * (sum(colMeans(Da)) / sum(colMeans(Dc)) - 1)))
  tt <- tests[[k]]
  cat("\n=== Test (as plot_cluster_loss_vs_tangle_distance_raw.pdf) ===\n")
  cat(sprintf("Spearman rho(x, y) across 10 deciles = %.3f, 95%% CI [%.3f, %.3f]; p(boot, CI inversion) = %s;\n",
              tt$rho, tt$ci_lo, tt$ci_hi, signif(tt$p_boot, 3)))
  cat(sprintf("p(perm, decile labels, %d draws) = %s. %.1f%% of bootstrap rho < 0.\n",
              tt$n_perm, signif(tt$p_perm, 3), 100 * tt$frac_boot_neg))
  if (k == "distance") cat("rho > 0 = deciles whose neurons sit CLOSER to tangles lost more excitatory neurons.\n")
  else cat("rho < 0 = deciles with more ", pn$short, " lost more excitatory neurons.\n", sep = "")
  cat("p(boot) = donor sampling: would a new set of donors give the same sign? (both axes resampled\n")
  cat("together, control and AD separately, so the shared AD donors and the control denominator are\n")
  cat("propagated). p(perm) permutes the decile labels across the 10 bands.\n\n")
  cat("=== Notes ===\n")
  cat(" - Relative depth is a fraction of each ROI's own extent (0 = pia).\n")
  cat(" - Columnar density is per mm of ROI width.\n")
  cat(" - y pools AD exc neurons across donors (cohort-level).\n")
  if (k == "plaque")
    cat(" - Matched_4G8_40 is a binary flag from a 40 um dilation of the 4G8 mask.\n", sep = "")
  if (k == "intensity")
    cat(" - Mean z is over all neurons, PHF1+ and PHF1-.\n", sep = "")
  cat("\n"); print(sessionInfo())
  sink()
  if (CLASS == "inh") relabel_log(file.path(out_dir, paste0("stats_density_shift_vs_", pn$file, ".txt")))
}

cmp <- res$annotation %>% select(decile, depth, rank_annotation = rank, pct_phf1_pos_ad = y) %>%
  left_join(res$intensity %>% select(decile, rank_intensity = rank, phf1_z_ad = y), by = "decile") %>%
  left_join(res$plaque %>% select(decile, rank_plaque = rank, pct_plaque_40um_ad = y), by = "decile") %>%
  left_join(res$distance %>% select(decile, rank_distance = rank, dist_phf1_um_ad = y), by = "decile") %>%
  arrange(decile)
wr(cmp, "stats_density_shift_rank_comparison.tsv")
wr(bind_rows(tests), "stats_density_shift_spearman.tsv")
print(as.data.frame(bind_rows(tests)))
cat(sprintf("Rank agreement annotation vs intensity: Spearman %.3f; annotation vs plaque: %.3f\n",
            cor(cmp$rank_annotation, cmp$rank_intensity, method = "spearman"),
            cor(cmp$rank_annotation, cmp$rank_plaque, method = "spearman")))
cat("Done:", out_dir, "\n")
