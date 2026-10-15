#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# pseudotime_by_phf1_distance.R
#
# Figure panels: 3C
# 3C is plot_pseudotime_by_phf1_distance_rollmean.pdf. The raw-micron hinge model is a
# documented exception to the log-distance transform (see docs/MODELS.md).
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# pseudotime_by_phf1_distance.R
#
# Where do tangle-free neurons sit on the tangle-associated Slingshot trajectory (2B) as a
# function of distance to the nearest PHF1+ neuron, and how does that compare with the
# tangle-bearing neurons themselves?
#
# Rolling-mean panel (3C): pooled 75 um sliding-window mean of pseudotime in tangle-free cells
# (donor-aligned to the 500-1000 um baseline), 95% CI = 1.96 x SEM; dotted line = mean pseudotime
# of PHF1+ neurons; dashed grey line + band = breakpoint of the single-hinge model on raw distance
# with its donor-bootstrap 95% CI. Only the single-hinge breakpoint is drawn.
#
# Binned panel: x = distance bin, far -> near, then PHF1+; y = pseudotime shift relative to the
# same donor's far (500-1000 um) neurons, in SD of pseudotime. One faint point per donor, mean and
# 95% t CI across donors. PHF1+ is shown for donors with >= MIN_PHF1_PER_DONOR tangle-bearing cells
# (the trajectory's own threshold).
#
# Stats:
#   (1) Donor-paired: <50 um vs far, and PHF1+ vs far -- mean shift, t CI, dz, donors in the same
#       direction; plus the fraction of the PHF1+ shift reached at <50 um.
#   (2) Set 1 LMM on the tangle-free cells (canonical distance-DEG cell set: PHF1-, <= 1000 um):
#         pt_scaled ~ dist_to_phf1_um_scaled + nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)
#       dist_to_phf1_um_scaled = log(d) / sd(log(d)) (canonical, never log1p); pt_scaled = pt / SD(pt)
#       over all trajectory cells, exactly as deg_dream_linear_pseudotime_phf1.r.
#   (3) Categorical LMM over the plotted groups (all trajectory cells <= 1000 um + PHF1+), same
#       covariates; emmeans treatment-vs-far contrasts (Dunnett), estimate + 95% CI + Cohen's d.
#   (4) Hinge breakpoints on raw um (single hinge and plateau-ramp-plateau), donor bootstrap CIs.
#   (5) Decay length constant (R/decay_length_utils.R), stats log only.
#   Sensitivity: the donor-paired table without the 1000 um cap (far bin = >= 500 um).
#
# Inputs: celltype_sce_neighbours/<ct>_sce_neighbours.qs,
#   plots/<tag>_slingshot_trajectory/source_data_<tag>_traj_percell.tsv (cbln2_slingshot_trajectory.R)
# Output: plots/pseudotime_by_phf1_distance/<ct>/  (plot / source_data / stats triple)
# Usage:  Rscript R/pseudotime_by_phf1_distance.R [celltype]

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
  library(SingleCellExperiment); library(lme4); library(lmerTest); library(emmeans)
})

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")

args <- commandArgs(trailingOnly = TRUE)
CT   <- if (length(args) >= 1) args[1] else "Exc-IT-L2-3-CBLN2-HOPX"
TAG  <- c("Exc-IT-L2-3-CBLN2-HOPX" = "cbln2", "Exc-IT-L3-5-CHGA-IL1RAPL2" = "chga")[[CT]]

MAX_DIST_UM        <- 1000                                 # canonical distance-DEG cap
BREAKS             <- c(0, 50, 100, 200, 350, 500, MAX_DIST_UM)
FAR_LABEL          <- "500-1000"
MIN_PHF1_PER_DONOR <- 5                                    # trajectory tangle threshold
CONF_VARS          <- c("Sex", "Age", "PMI")

sce_path <- file.path("celltype_sce_neighbours", paste0(CT, "_sce_neighbours.qs"))
pt_path  <- file.path("plots", paste0(TAG, "_slingshot_trajectory"), paste0("source_data_", TAG, "_traj_percell.tsv"))
for (f in c(sce_path, pt_path)) if (!file.exists(f)) stop("missing input: ", f)
out_dir <- file.path("plots/pseudotime_by_phf1_distance", CT)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------ cells ---------------------------------------------------------
sce <- qs::qread(sce_path)
cd  <- as.data.frame(colData(sce))
cd$cell_id     <- colnames(sce)
cd$nUMI_log    <- log2(colSums(counts(sce)) + 1)
cd$percent_neg <- as.numeric(cd$percent.neg)
cd$Age         <- as.numeric(scale(as.numeric(cd$Age)))
cd$PMI         <- as.numeric(scale(as.numeric(cd$PMI)))
cd$sample_id   <- droplevels(as.factor(cd$sample_id))
cd$PHF1        <- as.logical(as.character(cd$PHF1))
rm(sce); invisible(gc())

pt <- read.delim(pt_path, check.names = FALSE)[, c("cell_id", "pt")]
x  <- inner_join(cd[, c("cell_id", "sample_id", "PHF1", "dist_to_phf1_um", "nUMI_log", "percent_neg", CONF_VARS)],
                 pt, by = "cell_id")
if (nrow(x) != nrow(cd) || nrow(x) != nrow(pt))
  stop("cell join mismatch: SCE ", nrow(cd), ", pseudotime ", nrow(pt), ", joined ", nrow(x))
if (any(!x$PHF1 & !(x$dist_to_phf1_um > 0), na.rm = TRUE) || any(!x$PHF1 & is.na(x$dist_to_phf1_um)))
  stop("non-positive or missing distance on a PHF1-negative cell")

pt_sd <- sd(x$pt)                                          # over ALL trajectory cells, as the DEG
x$pt_scaled <- x$pt / pt_sd

grp_levels <- c(FAR_LABEL, "350-500", "200-350", "100-200", "50-100", "<50", "PHF1+")
bin_of <- function(d, far_label) {
  b <- cut(d, BREAKS, right = FALSE, labels = c("<50", "50-100", "100-200", "200-350", "350-500", far_label))
  as.character(b)
}
x$grp <- ifelse(x$PHF1, "PHF1+", bin_of(x$dist_to_phf1_um, FAR_LABEL))
x$grp <- factor(x$grp, levels = grp_levels)                # NA = PHF1- beyond the cap

# ------------------------------ donor-level shifts --------------------------------------------
donor_shift <- function(dat, far) {
  dat %>% filter(!is.na(grp)) %>%
    group_by(sample_id) %>% mutate(far_mean = mean(pt_scaled[grp == far])) %>%
    group_by(sample_id, grp) %>%
    summarise(n_cells = n(), pt_mean = mean(pt_scaled), shift = mean(pt_scaled - far_mean), .groups = "drop")
}
dd <- donor_shift(x, FAR_LABEL)
n_phf1 <- x %>% filter(PHF1) %>% dplyr::count(sample_id, name = "n_phf1")
dd <- dd %>% left_join(n_phf1, by = "sample_id") %>% mutate(n_phf1 = coalesce(n_phf1, 0L),
  shown = grp != "PHF1+" | n_phf1 >= MIN_PHF1_PER_DONOR)

tci <- function(v) { n <- length(v); m <- mean(v); h <- if (n > 1) qt(0.975, n - 1) * sd(v) / sqrt(n) else NA
  c(mean = m, lo = m - h, hi = m + h, n = n) }
summ <- dd %>% filter(shown) %>% group_by(grp) %>%
  summarise(mean = tci(shift)["mean"], lo = tci(shift)["lo"], hi = tci(shift)["hi"],
            n_donors = n(), n_cells = sum(n_cells), .groups = "drop")

paired <- function(dat, a, far, min_phf1 = 0) {
  w <- dat %>% filter(grp %in% c(a, far)) %>%
    filter(!(grp == "PHF1+" & n_phf1 < min_phf1)) %>%
    select(sample_id, grp, pt_mean) %>% pivot_wider(names_from = grp, values_from = pt_mean) %>%
    filter(!is.na(.data[[a]]), !is.na(.data[[far]]))
  d <- w[[a]] - w[[far]]; tt <- t.test(d)
  tibble(contrast = paste(a, "vs", far), n_donors = length(d), mean_shift = mean(d),
         ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2], dz = mean(d) / sd(d),
         n_same_direction = sum(d > 0), p_paired_t = tt$p.value)
}
pair_tbl <- bind_rows(paired(dd, "<50", FAR_LABEL),
                      paired(dd, "PHF1+", FAR_LABEL, MIN_PHF1_PER_DONOR))
frac_tbl <- dd %>% filter(grp %in% c("<50", "PHF1+"), n_phf1 >= MIN_PHF1_PER_DONOR) %>%
  select(sample_id, grp, shift) %>% pivot_wider(names_from = grp, values_from = shift) %>%
  mutate(frac_of_phf1 = `<50` / `PHF1+`)
frac_of_means <- with(frac_tbl, mean(`<50`) / mean(`PHF1+`))

# Sensitivity: no 1000 um cap (far = >= 500 um)
BREAKS_UNCAP <- c(0, 50, 100, 200, 350, 500, Inf)
xu <- x; xu$grp <- ifelse(xu$PHF1, "PHF1+",
  as.character(cut(xu$dist_to_phf1_um, BREAKS_UNCAP, right = FALSE,
                   labels = c("<50", "50-100", "100-200", "200-350", "350-500", ">=500"))))
xu$grp <- factor(xu$grp, levels = c(">=500", grp_levels[-1]))
ddu <- donor_shift(xu, ">=500") %>% left_join(n_phf1, by = "sample_id") %>% mutate(n_phf1 = coalesce(n_phf1, 0L))
pair_uncap <- bind_rows(paired(ddu, "<50", ">=500"), paired(ddu, "PHF1+", ">=500", MIN_PHF1_PER_DONOR))

# ------------------------------ (2) Set 1 LMM, continuous log distance ------------------------
xn <- x %>% filter(!PHF1, dist_to_phf1_um <= MAX_DIST_UM)
dist_sd <- sd(log(xn$dist_to_phf1_um))
xn$dist_to_phf1_um_scaled <- log(xn$dist_to_phf1_um) / dist_sd
form_c <- as.formula(paste("pt_scaled ~ dist_to_phf1_um_scaled + nUMI_log + percent_neg +",
                           paste(CONF_VARS, collapse = " + "), "+ (1|sample_id)"))
m_c  <- lmer(form_c, data = xn, REML = TRUE)
co   <- coef(summary(m_c))["dist_to_phf1_um_scaled", ]
ci_c <- confint(m_c, parm = "dist_to_phf1_um_scaled", method = "Wald")
vc   <- as.data.frame(VarCorr(m_c)); sd_tot <- sqrt(sum(vc$vcov))
lmm_tbl <- tibble(term = "dist_to_phf1_um_scaled", estimate = co[["Estimate"]], se = co[["Std. Error"]],
                  ci_lo = ci_c[1], ci_hi = ci_c[2], df = co[["df"]], p = co[["Pr(>|t|)"]],
                  cohens_d = co[["Estimate"]] / sd_tot, n_cells = nrow(xn), n_donors = nlevels(droplevels(xn$sample_id)),
                  dist_sd = dist_sd, pt_sd = pt_sd, donor_re_var = vc$vcov[vc$grp == "sample_id"])

# ------------------------------ (3) categorical LMM over the plotted groups -------------------
xc <- x %>% filter(!is.na(grp))
form_g <- as.formula(paste("pt_scaled ~ grp + nUMI_log + percent_neg +",
                           paste(CONF_VARS, collapse = " + "), "+ (1|sample_id)"))
m_g  <- lmer(form_g, data = xc, REML = TRUE)
vc_g <- as.data.frame(VarCorr(m_g)); sd_tot_g <- sqrt(sum(vc_g$vcov))
emm  <- emmeans(m_g, ~ grp, lmer.df = "asymptotic")        # large n: asymptotic df (no KR/Satterthwaite)
ctr  <- contrast(emm, method = "trt.vs.ctrl", ref = 1, adjust = "dunnettx")
emm_tbl <- as.data.frame(confint(ctr)) %>%
  left_join(as.data.frame(ctr)[, c("contrast", "z.ratio", "p.value")], by = "contrast") %>%
  mutate(cohens_d = estimate / sd_tot_g)

# ------------------------------ plot ----------------------------------------------------------
dd_plot <- dd %>% filter(shown) %>% mutate(grp = factor(grp, levels = grp_levels))
summ    <- summ %>% mutate(grp = factor(grp, levels = grp_levels), is_phf1 = grp == "PHF1+")
p <- ggplot(mapping = aes(grp)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_line(data = dd_plot %>% filter(grp != "PHF1+"), aes(y = shift, group = sample_id),
            colour = "grey80", linewidth = 0.25) +
  geom_point(data = dd_plot, aes(y = shift), colour = "grey65", size = 0.6) +
  geom_line(data = filter(summ, !is_phf1), aes(y = mean, group = 1), colour = "#D55E00", linewidth = 0.5) +
  geom_errorbar(data = summ, aes(ymin = lo, ymax = hi), width = 0, colour = "#D55E00", linewidth = 0.4) +
  geom_point(data = summ, aes(y = mean, shape = is_phf1), colour = "#D55E00", fill = "#D55E00", size = 1.4) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 23), guide = "none") +
  labs(x = expression("Distance to PHF1+ neuron (" * mu * "m)"),
       y = "Pseudotime shift (SD)") +
  theme_classic(base_size = 8) + fig_theme +
  theme(plot.margin = margin(4, 4, 4, 4))
ggsave(file.path(out_dir, "plot_pseudotime_by_phf1_distance.pdf"), p, width = 2.2, height = 2, device = "pdf")

write.table(bind_rows(dd_plot %>% mutate(row_type = "donor") %>%
                        select(row_type, grp, sample_id, n_cells, n_phf1, shift),
                      summ %>% mutate(row_type = "mean_ci") %>%
                        select(row_type, grp, n_donors, n_cells, mean, lo, hi)),
            file.path(out_dir, "source_data_pseudotime_by_phf1_distance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ------------------------------ rolling-mean variant (module-score convention) ----------------
# No bins: a pooled WINDOW_UM sliding-window mean + 95% CI (1.96 x SEM) over the tangle-free cells, as
# the integrated module-score panels (plot_phf1_module_vs_phf1_distance_modelp.R, 75 um default), with
# the PHF1+ neurons as a dotted horizontal reference (mean, 1.96 x SEM, all PHF1+ cells) as there.
# y is in PSEUDOTIME UNITS, donor-aligned: each cell's pt minus its donor's 500-1000 um mean, plus
# the mean of those donor far means -- so the far baseline reads as an actual pseudotime value and
# donors whose cells sit mostly far away cannot bend the pooled curve. The ribbon is descriptive;
# inference is the LMM / donor-paired dz.
ROLLMEAN_MIN_IN_WINDOW <- 20
WINDOW_UM <- 75
far_mean  <- x %>% filter(grp == FAR_LABEL) %>% group_by(sample_id) %>%
  summarise(far_mean_pt = mean(pt), .groups = "drop")
baseline_pt <- mean(far_mean$far_mean_pt)
x <- x %>% left_join(far_mean, by = "sample_id") %>% mutate(pt_aligned = pt - far_mean_pt + baseline_pt)
xr <- x %>% filter(!PHF1, dist_to_phf1_um <= MAX_DIST_UM)

grid <- seq(0, MAX_DIST_UM, length.out = 200)
roll_pooled <- bind_rows(lapply(grid, function(g) {
  k <- which(abs(xr$dist_to_phf1_um - g) <= WINDOW_UM / 2); n <- length(k)
  if (n < ROLLMEAN_MIN_IN_WINDOW) return(tibble(dist_to_phf1_um = g, roll_mean = NA_real_, sem = NA_real_, n_window = n))
  tibble(dist_to_phf1_um = g, roll_mean = mean(xr$pt_aligned[k]), sem = sd(xr$pt_aligned[k]) / sqrt(n), n_window = n)
}))
pv <- x$pt_aligned[x$PHF1]
phf1_ref <- tibble(ref_mean = mean(pv), ci_lo = mean(pv) - 1.96 * sd(pv) / sqrt(length(pv)),
                   ci_hi = mean(pv) + 1.96 * sd(pv) / sqrt(length(pv)), n_phf1_pos = length(pv))

# ------------------------------ hinge breakpoints (raw um; documented exception) ---------------
# Where does pseudotime leave the far baseline, and where does the increase level off?
# Plateau-ramp-plateau on RAW MICRONS (docs/MODELS.md, documented raw-um exception):
#   pt = P_donor + beta * [ (psi2 - d)+ - (psi1 - d)+ ] + nUMI_log + percent_neg      psi1 < psi2
# i.e. flat (baseline) beyond psi2, a linear ramp from psi2 in to psi1, flat again inside psi1.
# Also the single hinge  pt = P_donor + beta * (psi - d)+ + ...  (no near plateau). Raw um because the
# estimand is a distance where a flat baseline ends, and that elbow is linear in um.
# Donor = fixed intercept (Sex/Age/PMI aliased with it; see docs/MODELS.md). Exact grid
# search: the residualised hinge basis R is computed once, then every (psi1, psi2) pair's RSS follows
# from crossprod(R) and R'y. Donor bootstrap (integer weights = duplicated donors) for the CIs.
HINGE_GRID   <- seq(10, 900, by = 10)
HINGE_MINSEP <- 30
HINGE_B      <- 999
Zh   <- model.matrix(~ nUMI_log + percent_neg + sample_id, data = transform(xn, sample_id = droplevels(sample_id)))
Hraw <- vapply(HINGE_GRID, function(p) pmax(p - xn$dist_to_phf1_um, 0), numeric(nrow(xn)))
hinge_fit <- function(w = rep(1, nrow(xn))) {
  sw <- sqrt(w); q <- qr(sw * Zh); yt <- qr.resid(q, sw * xn$pt); R <- qr.resid(q, sw * Hraw)
  yy <- sum(yt^2); G <- crossprod(R); r <- drop(crossprod(R, yt))
  rss_h <- yy - r^2 / diag(G)
  num <- outer(-r, r, "+"); den <- outer(diag(G), diag(G), "+") - 2 * G
  ok  <- outer(HINGE_GRID, HINGE_GRID, function(a, b) b - a >= HINGE_MINSEP)
  rss_pr <- matrix(Inf, length(HINGE_GRID), length(HINGE_GRID)); rss_pr[ok] <- yy - num[ok]^2 / den[ok]
  k <- which(rss_pr == min(rss_pr), arr.ind = TRUE)[1, ]; jh <- which.min(rss_h)
  list(psi_h = HINGE_GRID[jh], rss_h = rss_h[jh], beta_h = r[jh] / G[jh, jh],
       psi1 = HINGE_GRID[k[1]], psi2 = HINGE_GRID[k[2]], rss_pr = rss_pr[k[1], k[2]],
       beta_pr = num[k[1], k[2]] / den[k[1], k[2]], rss0 = yy)
}
hf <- hinge_fit()
# Exactness check: projection RSS == lm() at the selected breakpoints
lm_rss <- sum(resid(lm(pt ~ I(pmax(hf$psi2 - dist_to_phf1_um, 0) - pmax(hf$psi1 - dist_to_phf1_um, 0)) +
                         nUMI_log + percent_neg + sample_id, data = xn))^2)
if (abs(lm_rss - hf$rss_pr) > 1e-6 * lm_rss) stop("hinge projection RSS != lm() RSS")
nh <- nrow(xn)
ic <- function(rss, k, pen) nh * log(rss / nh) + k * pen
rss_lin <- sum(qr.resid(qr(cbind(Zh, xn$dist_to_phf1_um)), xn$pt)^2)
rss_log <- sum(qr.resid(qr(cbind(Zh, log(xn$dist_to_phf1_um))), xn$pt)^2)
# extra parameters beyond the covariate/donor block: slope 1; hinge slope + psi 2; plateau model 3
hinge_cmp <- tibble(model = c("no distance term", "linear in um", "log-linear (canonical axis)",
                              "single hinge (flat beyond psi)", "plateau-ramp-plateau"),
                    rss = c(hf$rss0, rss_lin, rss_log, hf$rss_h, hf$rss_pr), k_extra = c(0, 1, 1, 2, 3)) %>%
  mutate(BIC = ic(rss, k_extra, log(nh)), AIC = ic(rss, k_extra, 2),
         dBIC = BIC - min(BIC), dAIC = AIC - min(AIC))
set.seed(42)
Dh <- levels(droplevels(xn$sample_id))
hb <- t(vapply(seq_len(HINGE_B), function(b) {
  m <- table(factor(sample(Dh, replace = TRUE), levels = Dh)); w <- as.numeric(m[as.character(xn$sample_id)])
  f <- hinge_fit(w); c(psi1 = f$psi1, psi2 = f$psi2, psi_h = f$psi_h, ramp = f$beta_pr, plateau = f$beta_pr * (f$psi2 - f$psi1))
}, numeric(5)))
qci <- function(v) unname(quantile(v, c(0.025, 0.975), type = 1))
hinge_tbl <- tibble(
  parameter = c("psi2_leaves_baseline_um", "psi1_increase_levels_off_um", "psi_single_hinge_um",
                "ramp_pt_per_100um", "near_plateau_excess_pt"),
  estimate  = c(hf$psi2, hf$psi1, hf$psi_h, 100 * hf$beta_pr, hf$beta_pr * (hf$psi2 - hf$psi1)),
  ci_lo = c(qci(hb[, "psi2"])[1], qci(hb[, "psi1"])[1], qci(hb[, "psi_h"])[1], 100 * qci(hb[, "ramp"])[1], qci(hb[, "plateau"])[1]),
  ci_hi = c(qci(hb[, "psi2"])[2], qci(hb[, "psi1"])[2], qci(hb[, "psi_h"])[2], 100 * qci(hb[, "ramp"])[2], qci(hb[, "plateau"])[2]),
  n_boot = HINGE_B)
write.table(hinge_tbl, file.path(out_dir, "stats_pseudotime_by_phf1_distance_hinge.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(hinge_cmp, file.path(out_dir, "stats_pseudotime_by_phf1_distance_hinge_modelcomparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.data.frame(hb), file.path(out_dir, "source_data_pseudotime_by_phf1_distance_hinge_bootstrap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# Plotted: ONLY the single-hinge breakpoint. The plateau-ramp-plateau model is reported in the
# stats log.
bp_lines <- tibble(which = "psi_single_hinge", psi = hf$psi_h,
                   lo = qci(hb[, "psi_h"])[1], hi = qci(hb[, "psi_h"])[2])

p_r <- ggplot() +
  # breakpoints + donor-bootstrap 95% CIs, drawn first so they sit behind everything else
  geom_rect(data = bp_lines, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf), fill = "grey50", alpha = 0.12, inherit.aes = FALSE) +
  geom_vline(data = bp_lines, aes(xintercept = psi), colour = "grey35", linetype = "dashed", linewidth = 0.35) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = phf1_ref$ci_lo, ymax = phf1_ref$ci_hi, fill = "#D55E00", alpha = 0.10) +
  geom_hline(yintercept = phf1_ref$ref_mean, colour = "#D55E00", linetype = "dotted", linewidth = 0.4) +
  annotate("text", x = MAX_DIST_UM, y = phf1_ref$ref_mean, label = "PHF1+", hjust = 1, vjust = -0.6, size = 2, colour = "#D55E00") +
  geom_hline(yintercept = baseline_pt, linewidth = 0.25, colour = "grey60") +
  geom_ribbon(data = filter(roll_pooled, is.finite(roll_mean)),
              aes(dist_to_phf1_um, ymin = roll_mean - 1.96 * sem, ymax = roll_mean + 1.96 * sem), fill = "#D55E00", alpha = 0.25) +
  geom_line(data = filter(roll_pooled, is.finite(roll_mean)), aes(dist_to_phf1_um, roll_mean), colour = "#D55E00", linewidth = 0.6) +
  coord_cartesian(xlim = c(0, MAX_DIST_UM)) +
  labs(x = expression("Distance to PHF1+ neuron (" * mu * "m)"), y = "Pseudotime") +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), plot.margin = margin(4, 6, 4, 4))
ggsave(file.path(out_dir, "plot_pseudotime_by_phf1_distance_rollmean.pdf"), p_r, width = 2.2, height = 2, device = "pdf")

write.table(bind_rows(roll_pooled %>% mutate(row_type = "pooled"),
                      phf1_ref %>% transmute(row_type = "phf1_ref", roll_mean = ref_mean,
                                             sem = (ci_hi - ref_mean) / 1.96, n_window = n_phf1_pos),
                      tibble(row_type = "far_baseline", roll_mean = baseline_pt),
                      bp_lines %>% transmute(row_type = paste0("breakpoint_", which), dist_to_phf1_um = psi,
                                             ci_lo_um = lo, ci_hi_um = hi)) %>%
              mutate(window_um = WINDOW_UM, y_units = "pseudotime (donor-aligned to the 500-1000 um baseline)"),
            file.path(out_dir, "source_data_pseudotime_by_phf1_distance_rollmean.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ------------------------------ decay (length constant) ---------------------------------------
# Length constant of the approach to baseline, stats log only. Decay engine R/decay_length_utils.R
# (raw-um exception, docs/MODELS.md): per-donor
# pt = P_j + A_j exp(-d / lambda_j) + nUMI_log + percent_neg, DL + Knapp-Hartung meta-analysis of
# log lambda, d95 = 3 lambda, with the engine's reporting gates.
# Outcome = raw pt (P_j absorbs donor level, so no alignment needed; lambda does not depend on the
# outcome scale).
source("R/decay_length_utils.R")
dec <- decay_analyse(as.data.frame(xn), "pt", cap_um = MAX_DIST_UM)
write.table(decay_summary_row(dec, celltype = CT, outcome = "pseudotime"),
            file.path(out_dir, "stats_pseudotime_by_phf1_distance_decay_lambda.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------ stats ---------------------------------------------------------
eff <- bind_rows(
  pair_tbl %>% transmute(analysis = "donor_paired", contrast, estimate = mean_shift, ci_lo, ci_hi,
                         effect_size = dz, effect_size_type = "dz", n_donors, n_same_direction, p = p_paired_t),
  lmm_tbl %>% transmute(analysis = "set1_lmm_log_distance", contrast = term, estimate, ci_lo, ci_hi,
                        effect_size = cohens_d, effect_size_type = "cohens_d", n_donors, p),
  emm_tbl %>% transmute(analysis = "categorical_lmm_emmeans", contrast, estimate, ci_lo = asymp.LCL, ci_hi = asymp.UCL,
                        effect_size = cohens_d, effect_size_type = "cohens_d", p = p.value))
write.table(eff, file.path(out_dir, "stats_pseudotime_by_phf1_distance_effectsize.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

r4 <- function(d) as.data.frame(mutate(d, across(where(~ is.numeric(.x) && any(.x != round(.x), na.rm = TRUE)), ~ signif(.x, 4))))   # counts stay exact
sink(file.path(out_dir, "stats_pseudotime_by_phf1_distance.txt"))
cat("Pseudotime position by distance to the nearest PHF1+ neuron --", CT, "\n")
cat("Date:", format(Sys.time()), "\n\n")
cat("Inputs:", sce_path, "|", pt_path, "\n")
cat("Cells joined:", nrow(x), "(PHF1+", sum(x$PHF1), "); PHF1- beyond", MAX_DIST_UM, "um (not plotted, not modelled):",
    sum(!x$PHF1 & x$dist_to_phf1_um > MAX_DIST_UM), "\n")
cat(sprintf("pt_scaled = pt / %.4f (SD over all trajectory cells, as deg_dream_linear_pseudotime_phf1.r)\n", pt_sd))
cat("y = per-donor mean pt_scaled in the group minus that donor's", FAR_LABEL, "um mean.\n")
cat("PHF1+ plotted for donors with >=", MIN_PHF1_PER_DONOR, "PHF1+ cells:\n"); print(as.data.frame(n_phf1))
cat("\n=== Per-group summary across donors (plotted) ===\n"); print(r4(summ))
cat("\n=== (1) Donor-paired (headline) ===\n"); print(r4(pair_tbl))
cat("\nFraction of the PHF1+ shift reached at <50 um (donors with >=", MIN_PHF1_PER_DONOR, "PHF1+):\n")
print(r4(frac_tbl))
cat(sprintf("ratio of means = %.3f; median per-donor fraction = %.3f\n", frac_of_means, median(frac_tbl$frac_of_phf1)))
cat("\nSensitivity, no 1000 um cap (far = >= 500 um):\n"); print(r4(pair_uncap))
cat("\n=== Rolling-mean panel (plot_pseudotime_by_phf1_distance_rollmean) ===\n")
cat(sprintf("Window %g um, grid 200 points over 0-%g um, windows < %d cells blank; CI = 1.96 x SEM (cell-level,\n",
            WINDOW_UM, MAX_DIST_UM, ROLLMEAN_MIN_IN_WINDOW))
cat("descriptive). y = pt - donor's", FAR_LABEL, "um mean + mean of donor far means\n")
cat(sprintf("(far baseline = %.3f pseudotime units; 1 SD of pseudotime = %.3f units).\n", baseline_pt, pt_sd))
cat("PHF1+ reference (all PHF1+ cells, same alignment):\n"); print(r4(phf1_ref))
cat("Pooled rolling mean at selected distances:\n")
print(r4(roll_pooled %>% dplyr::slice(vapply(c(0, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 750, 1000),
                                             function(d) which.min(abs(grid - d)), integer(1)))))
cat("\n=== Hinge breakpoints (raw um; docs/MODELS.md distance-transform exception) ===\n")
cat("Plateau-ramp-plateau: pt = P_donor + beta*[(psi2 - d)+ - (psi1 - d)+] + nUMI_log + percent_neg (donor fixed);\n")
cat("single hinge: pt = P_donor + beta*(psi - d)+ + ... . Exact grid search", min(HINGE_GRID), "-", max(HINGE_GRID),
    "um step", diff(HINGE_GRID)[1], "um, psi2 - psi1 >=", HINGE_MINSEP, "um; donor bootstrap B =", HINGE_B, "(seed 42).\n")
cat(sprintf("Projection RSS checked against lm(): %.6g vs %.6g\n", hf$rss_pr, lm_rss))
cat("Model comparison (identical rows; k_extra = parameters beyond donor + covariates):\n"); print(r4(hinge_cmp))
cat("\nEstimates (95% CI = donor-bootstrap percentile):\n"); print(r4(hinge_tbl))
cat("PLOTTED: the single-hinge psi only (dashed line + grey CI band). psi1 and psi2 of the plateau model are\n")
cat("reported for comparison.\n")
cat("\n=== Decay (length constant), canonical engine R/decay_length_utils.R, outcome = raw pseudotime ===\n")
decay_write_stats_block(dec, label = "pseudotime")
cat("\n=== (2) Set 1 LMM, tangle-free cells <=", MAX_DIST_UM, "um ===\n")
cat("Formula:", deparse(form_c, width.cutoff = 500), "\n")
cat(sprintf("dist_to_phf1_um_scaled = log(d) / %.4f\n", dist_sd))
cat("Sign: estimate < 0 = later in pseudotime NEARER tangles.\n")
print(summary(m_c)); cat("\n"); print(r4(lmm_tbl))
cat("Cohen's d = estimate / sqrt(sum of VarCorr variances). The donor-paired dz above is the donor-level\n")
cat("effect size.\n")
cat("\n=== (3) Categorical LMM over the plotted groups, emmeans vs", FAR_LABEL, "um (Dunnett, asymptotic df) ===\n")
cat("Formula:", deparse(form_g, width.cutoff = 500), "\n")
print(anova(m_g)); cat("\n"); print(r4(emm_tbl))
cat("\n"); print(sessionInfo())
sink()
cat("done ->", out_dir, "\n")
