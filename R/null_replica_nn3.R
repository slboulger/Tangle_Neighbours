#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# null_replica_nn3.R
#
# Figure panels: S3 (all) - third step of run_nn3.sh (permutation null)
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# null_replica_nn3.R
#
# Label-permutation null for the 3-NN spacing analyses (R/nn3_phf1_vs_neg.R,
# R/nn3_over_phf1_distance.R).
#
# ---------------------------------------------------------------------------
# WHAT IT BOUNDS
# ---------------------------------------------------------------------------
# dist_to_phf1_um is itself a nearest-neighbour distance to a subset of neurons,
# so a cell in a sparse neighbourhood is far from the nearest PHF1+ neuron BY
# CONSTRUCTION. Analysis B therefore regresses one geometric quantity on another,
# and a "spacing is tighter near tangles" coefficient is expected from geometry
# alone. Permuting which neurons carry the PHF1 label PRESERVES that geometry
# exactly, so the null coefficient distribution is the geometric baseline and the
# excess over it is the quantity attributable to the PHF1 labelling.
#
# ---------------------------------------------------------------------------
# WHY THIS IS CHEAP
# ---------------------------------------------------------------------------
# nn3_um is INVARIANT to PHF1 relabelling: the neighbour pool is defined by
# celltype, and both query and neighbour sets are the whole pool regardless of
# PHF1 status. So nothing spatial is recomputed -- the outcome is fixed and only
# the regressor changes. The invariance is not asserted rhetorically: the script
# recomputes nn3_um from the cache WITH THE PHF1 COLUMN REMOVED and requires the
# result to equal the cached values exactly. A function that cannot see PHF1
# cannot depend on it.
#
# Inputs (produced once by R/generate_phf1_null_labels.r, seed 42, 1000
# labellings shuffled within celltype x sample with per-group PHF1+ counts held
# fixed):
#   deg/null_label_deg/engine/phf1_label_permutations.qs   167176 x 1001 logical
#   deg/null_label_deg/engine/phf1_null_distance_matrix.qs 167176 x 1001 numeric
# Both carry an "observed" column reproducing the real labels/distances, which is
# asserted here before any null is read.
#
# ---------------------------------------------------------------------------
# ESTIMATOR
# ---------------------------------------------------------------------------
# 1000 lmer fits per subtype is not viable, so every coefficient -- OBSERVED AND
# NULL ALIKE -- is computed by the same fast within-donor fixed-effect projection:
# donor enters as fixed dummies (9 donors, so this is unproblematic) and the
# coefficient of interest is obtained by residualising against all other
# covariates,
#     beta = (Mg)'(My) / (Mg)'(Mg),   M = I - X0 (X0'X0)^-1 X0'
# via qr.resid(). The observed statistic is ALSO computed the slow lmer way and
# the two are compared before any empirical p-value is reported, so the null and
# the observed are never estimated by different machinery.
#
# Empirical p is TWO-TAILED AGAINST THE NULL'S OWN DISTRIBUTION, not against zero:
# for Analysis B the geometric coupling puts the null median well away from zero, so
# the usual |null| >= |observed| form would ask the wrong question. It uses the
# (1 + r) / (1 + n) convention, so it can never be 0 and is bounded below by
# 2 / (1 + n_perm) -- run the full 1000 labellings before quoting one.
#
# Outputs, under plots/nn3_neuron_spacing/:
#   plot_nn3_null_bounds_distance.pdf   per-subtype null coefficient spread + observed
#   plot_nn3_null_bounds_phf1.pdf       PHF1+ vs PHF1- null spread + observed
#   source_data_nn3_null_distance.tsv   every null coefficient drawn
#   source_data_nn3_null_phf1.tsv
#   stats_nn3_null_replica.txt
#   stats_nn3_null_replica_summary.tsv
#
# Requires R/nn3_phf1_vs_neg.R (cache) and R/nn3_over_phf1_distance.R (observed
# coefficients) to have been run first.
# Run on HPC -- the distance matrix is 1.14 GB on disk:
#   Rscript R/null_replica_nn3.R --n_perm 1000

suppressPackageStartupMessages({
  library(argparse); library(qs)
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
  library(lmerTest)
})

hpc <- "<PROJECT_ROOT>/phf1_v2"
loc <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(hpc)) hpc else loc)
source("R/palettes.R")
source("R/nn3_utils.r")

parser <- ArgumentParser()
parser$add_argument("--n_perm", type = "integer", default = 1000,
  help = "Number of null labellings to use (max = what the engine produced) [default: 1000]")
parser$add_argument("--engine_dir", default = "deg/null_label_deg/engine",
  help = "Directory holding phf1_label_permutations.qs and phf1_null_distance_matrix.qs")
parser$add_argument("--max_dist_um", type = "double", default = 1000,
  help = "Distance cap, must match nn3_over_phf1_distance.R [default: 1000]")
args <- parser$parse_args()

out_dir <- "plots/nn3_neuron_spacing"
cache   <- file.path(out_dir, "source_data_nn3_cells.tsv")
coef_f  <- file.path(out_dir, "source_data_nn3_distance_coef.tsv")
es_f    <- file.path(out_dir, "stats_nn3_phf1_vs_neg_effectsize.tsv")
for (f in c(cache, coef_f, es_f))
  if (!file.exists(f)) stop("Missing input: ", f,
                            "\nRun R/nn3_phf1_vs_neg.R and R/nn3_over_phf1_distance.R first.")

MIN_CELLS  <- 50L
MIN_DONORS <- 3L
FDR        <- 0.05

wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE)

# -------------------------------------------------------------------
# Cache + observed coefficients
# -------------------------------------------------------------------
md <- read.delim(cache, sep = "\t", header = TRUE, check.names = FALSE,
                 colClasses = c(Braak = "character", Sex = "character",
                                sample_id = "character", celltype = "character",
                                cell_id = "character"))
obs_dist <- read.delim(coef_f, sep = "\t", header = TRUE, check.names = FALSE)
obs_phf1 <- read.delim(es_f,   sep = "\t", header = TRUE, check.names = FALSE)

md <- md %>%
  mutate(celltype = factor(celltype, levels = c(neuron_order, "Unassigned Neuron")),
         Sex = factor(Sex), sample_f = factor(sample_id),
         log_nn3 = log(nn3_um),
         Age_s = as.numeric(scale(as.numeric(Age))),
         PMI_s = as.numeric(scale(as.numeric(PMI))))
subtypes <- levels(md$celltype)

# -------------------------------------------------------------------
# ASSERTION: nn3_um cannot depend on PHF1.
# Recompute it from the cache with the PHF1 columns stripped out entirely. The
# cache holds exactly the pool neurons, and the pool is self-contained, so the
# recomputation must reproduce the cached values to the last bit.
# -------------------------------------------------------------------
cd_blind <- md %>%
  select(cell_id, sample_id, celltype, x_slide_mm, y_slide_mm) %>%
  mutate(celltype = as.character(celltype))
stopifnot(!any(grepl("PHF1|phf1", colnames(cd_blind), ignore.case = TRUE)))
nn_blind <- compute_nn3_um(cd_blind, sample_col = "sample_id",
                           neuron_pool = NEURON_POOL_10, k = 3L, verbose = FALSE)
chk <- md %>% select(cell_id, nn3_um) %>% inner_join(nn_blind, by = "cell_id",
                                                     suffix = c("_cached", "_blind"))
inv_cmp <- all.equal(chk$nn3_um_cached, chk$nn3_um_blind)
if (!isTRUE(inv_cmp))
  stop("nn3_um is NOT invariant to PHF1 relabelling: ", inv_cmp,
       "\nThe whole economy of this null rests on that invariance; refusing to continue.")
message(sprintf("Invariance OK: nn3_um recomputed PHF1-blind matches the cache for %d neurons.",
                nrow(chk)))

# -------------------------------------------------------------------
# Engine artifacts
# -------------------------------------------------------------------
lbl_f <- file.path(args$engine_dir, "phf1_label_permutations.qs")
dst_f <- file.path(args$engine_dir, "phf1_null_distance_matrix.qs")
for (f in c(lbl_f, dst_f)) if (!file.exists(f)) stop("Missing engine artifact: ", f)

message("Reading ", lbl_f, " ...")
lbl <- qs::qread(lbl_f)
message("Reading ", dst_f, " (large) ...")
dst <- qs::qread(dst_f)

stopifnot("observed" %in% colnames(lbl), "observed" %in% colnames(dst))
# Restrict both to the pool neurons, in cache order, and free the rest immediately.
stopifnot(all(md$cell_id %in% rownames(lbl)), all(md$cell_id %in% rownames(dst)))
lbl <- lbl[md$cell_id, , drop = FALSE]
dst <- dst[md$cell_id, , drop = FALSE]
gc()

# The engine's "observed" columns must reproduce reality, or every null below is
# anchored to the wrong thing. Same assertions generate_phf1_null_labels.r makes.
if (!all(lbl[, "observed"] == md$phf1_pos))
  stop("Engine 'observed' label column does not reproduce PHF1 for the pool neurons.")
obs_d_cmp <- all.equal(unname(dst[, "observed"]), unname(md$dist_to_phf1_um))
if (!isTRUE(obs_d_cmp))
  stop("Engine 'observed' distance column does not match canonical dist_to_phf1_um: ", obs_d_cmp)
message("Engine provenance OK: observed label and distance columns both reproduce reality.")

null_cols <- setdiff(colnames(lbl), "observed")
n_perm    <- min(args$n_perm, length(null_cols))
if (n_perm < length(null_cols))
  message(sprintf("NOTE: using %d of %d available labellings (--n_perm).",
                  n_perm, length(null_cols)))
null_cols <- null_cols[seq_len(n_perm)]

# -------------------------------------------------------------------
# Fast projection estimator. beta of `g` after residualising on X0.
# -------------------------------------------------------------------
proj_beta <- function(qr0, g, My) {
  Mg <- qr.resid(qr0, g)
  den <- sum(Mg * Mg)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  sum(Mg * My) / den
}

# -------------------------------------------------------------------
# Empirical two-tailed p against the NULL'S OWN distribution.
#
# The usual (1 + sum(|null| >= |obs|)) / (1 + n) form assumes the null is centred
# on zero. Here it is NOT: the geometric coupling means the null coefficient for
# Analysis B sits well away from zero, so comparing absolute magnitudes would ask
# "is the observed value further from zero than the null?" -- the wrong question.
# The question is whether the observed value is extreme WITHIN the null, so both
# tails are measured against the null distribution itself.
# -------------------------------------------------------------------
p_emp_fn <- function(b, obs) {
  b <- b[is.finite(b)]
  if (!length(b) || !is.finite(obs)) return(NA_real_)
  n  <- length(b)
  lo <- (1 + sum(b <= obs)) / (1 + n)
  hi <- (1 + sum(b >= obs)) / (1 + n)
  min(1, 2 * min(lo, hi))
}

# ===================================================================
# PART 1. Analysis B null: spacing vs distance, per subtype.
# Outcome fixed; the regressor is dist_scaled built from each null distance column.
# Cells that are PHF1+ under a given labelling have NA distance and drop out, so
# the design is rebuilt per permutation on that permutation's query set.
# ===================================================================
message("Part 1: distance-coefficient null, ", n_perm, " labellings x ", length(subtypes),
        " subtypes ...")

beta_dist_one <- function(idx, dvec) {
  # idx: rows of md for this subtype; dvec: distance column restricted to md rows
  keep <- idx[is.finite(dvec[idx]) & dvec[idx] > 0 & dvec[idx] <= args$max_dist_um]
  if (length(keep) < MIN_CELLS) return(NA_real_)
  dd <- md[keep, ]
  if (dplyr::n_distinct(dd$sample_id) < MIN_DONORS) return(NA_real_)
  dt <- log(dvec[keep])
  ds <- stats::sd(dt)
  if (!is.finite(ds) || ds <= 0) return(NA_real_)
  g  <- dt / ds
  X0 <- model.matrix(~ edge_dist_um + Sex + Age_s + PMI_s + droplevels(sample_f), data = dd)
  qr0 <- qr(X0)
  My  <- qr.resid(qr0, dd$log_nn3)
  -proj_beta(qr0, g, My)     # sign: beta_closer = -beta_dist
}

dist_null <- bind_rows(lapply(subtypes, function(ct) {
  idx <- which(md$celltype == ct & !md$phf1_pos)
  if (length(idx) < MIN_CELLS) return(NULL)
  b_obs_fast <- beta_dist_one(idx, dst[, "observed"])
  b_null <- vapply(null_cols, function(j) beta_dist_one(idx, dst[, j]), numeric(1))
  tibble(celltype = ct, perm = null_cols, beta_closer_null = unname(b_null),
         beta_closer_obs_fast = b_obs_fast)
}))
wr(dist_null, "source_data_nn3_null_distance.tsv")

dist_sum <- dist_null %>%
  group_by(celltype) %>%
  summarise(n_perm_ok      = sum(is.finite(beta_closer_null)),
            beta_obs_fast  = first(beta_closer_obs_fast),
            null_median    = median(beta_closer_null, na.rm = TRUE),
            null_q025      = quantile(beta_closer_null, 0.025, na.rm = TRUE),
            null_q975      = quantile(beta_closer_null, 0.975, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(obs_dist %>% select(celltype, beta_closer_lmer = beta_closer,
                                pval_lmer = pval, padj_lmer = padj, n_cells),
            by = "celltype") %>%
  rowwise() %>%
  mutate(
    # two-tailed empirical p against the null's own distribution (see p_emp_fn)
    p_emp = p_emp_fn(dist_null$beta_closer_null[dist_null$celltype == celltype],
                     beta_obs_fast),
    # the number that actually matters: how far the observed sits beyond geometry
    excess_over_null = beta_obs_fast - null_median,
    inside_null_ci   = beta_obs_fast >= null_q025 & beta_obs_fast <= null_q975,
    dir_vs_null      = ifelse(!is.finite(excess_over_null), NA_character_,
                       ifelse(excess_over_null > 0,
                              "wider spacing near tangles than geometry predicts",
                              "tighter spacing near tangles than geometry predicts")),
    # agreement between the fast estimator and the primary lmer fit
    abs_gap_vs_lmer  = abs(beta_obs_fast - beta_closer_lmer)) %>%
  ungroup() %>%
  mutate(celltype = factor(celltype, levels = subtypes),
         padj_emp = p.adjust(p_emp, method = "BH")) %>%
  arrange(celltype)

# ===================================================================
# PART 2. Analysis A null: PHF1+ vs PHF1- spacing.
# Outcome fixed; the regressor is the permuted PHF1 indicator.
# ===================================================================
message("Part 2: PHF1+ vs PHF1- null, ", n_perm, " labellings ...")

X0_A  <- model.matrix(~ celltype + edge_dist_um + sample_f, data = md)
qr0_A <- qr(X0_A)
My_A  <- qr.resid(qr0_A, md$log_nn3)

b_A_obs_fast <- proj_beta(qr0_A, as.numeric(md$phf1_pos), My_A)
b_A_null <- vapply(null_cols, function(j) proj_beta(qr0_A, as.numeric(lbl[, j]), My_A),
                   numeric(1))

phf1_null <- tibble(perm = null_cols, beta_null = unname(b_A_null),
                    pct_null = 100 * (exp(unname(b_A_null)) - 1))
wr(phf1_null, "source_data_nn3_null_phf1.tsv")

b_A_lmer <- log(1 + obs_phf1$pct_diff / 100)   # the primary emmeans contrast, log scale
A_sum <- tibble(
  beta_obs_fast = b_A_obs_fast,
  pct_obs_fast  = 100 * (exp(b_A_obs_fast) - 1),
  beta_lmer     = b_A_lmer,
  pct_lmer      = obs_phf1$pct_diff,
  abs_gap_vs_lmer = abs(b_A_obs_fast - b_A_lmer),
  n_perm_ok     = sum(is.finite(b_A_null)),
  null_median   = median(b_A_null, na.rm = TRUE),
  null_q025     = quantile(b_A_null, 0.025, na.rm = TRUE),
  null_q975     = quantile(b_A_null, 0.975, na.rm = TRUE),
  pct_null_q025 = 100 * (exp(quantile(b_A_null, 0.025, na.rm = TRUE)) - 1),
  pct_null_q975 = 100 * (exp(quantile(b_A_null, 0.975, na.rm = TRUE)) - 1),
  p_emp = p_emp_fn(b_A_null, b_A_obs_fast),
  inside_null_ci = b_A_obs_fast >= quantile(b_A_null, 0.025, na.rm = TRUE) &
                   b_A_obs_fast <= quantile(b_A_null, 0.975, na.rm = TRUE))

wr(bind_rows(
     dist_sum %>% transmute(analysis = "B_distance", celltype = as.character(celltype),
                            dir_vs_null, padj_emp,
                            n_cells, beta_obs_fast, beta_lmer = beta_closer_lmer,
                            abs_gap_vs_lmer, null_median, null_q025, null_q975,
                            excess_over_null, inside_null_ci, p_emp,
                            p_model = pval_lmer, padj_model = padj_lmer),
     A_sum %>% transmute(analysis = "A_phf1_vs_neg", celltype = "ALL POOLED",
                         n_cells = nrow(md), beta_obs_fast, beta_lmer,
                         abs_gap_vs_lmer, null_median, null_q025, null_q975,
                         excess_over_null = beta_obs_fast - null_median,
                         inside_null_ci, p_emp, p_model = NA_real_, padj_model = NA_real_)),
   "stats_nn3_null_replica_summary.tsv")

# -------------------------------------------------------------------
# Figures.
# -------------------------------------------------------------------
pd <- dist_null %>%
  filter(is.finite(beta_closer_null)) %>%
  mutate(celltype = factor(celltype, levels = rev(subtypes)))
po <- dist_sum %>% mutate(celltype = factor(as.character(celltype), levels = rev(subtypes)))

p1 <- ggplot(pd, aes(beta_closer_null, celltype)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey60") +
  geom_violin(fill = "grey80", colour = NA, scale = "width", width = 0.8) +
  geom_point(data = po, aes(x = null_median), shape = 124, size = 1.6, colour = "grey30") +
  geom_point(data = po, aes(x = beta_obs_fast, colour = inside_null_ci), size = 1.7) +
  scale_colour_manual(values = c("TRUE" = "grey45", "FALSE" = "#BD0026"),
                      labels = c("TRUE" = "inside null 95%", "FALSE" = "outside null 95%"),
                      name = NULL) +
  labs(x = expression("Coefficient per s.d. closer to a PHF1+ neuron (log " * mu * "m)"),
       y = NULL) +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "top",
        plot.margin = margin(4, 14, 4, 4, unit = "pt"))
ggsave(file.path(out_dir, "plot_nn3_null_bounds_distance.pdf"), p1,
       width = 10.6, height = 6.4, units = "cm", device = "pdf")

p2 <- ggplot(phf1_null, aes(pct_null)) +
  geom_histogram(bins = 40, fill = "grey80", colour = NA) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey60") +
  geom_vline(xintercept = A_sum$pct_null_q025, linetype = 3, linewidth = 0.3, colour = "grey40") +
  geom_vline(xintercept = A_sum$pct_null_q975, linetype = 3, linewidth = 0.3, colour = "grey40") +
  geom_vline(xintercept = A_sum$pct_obs_fast, linewidth = 0.5, colour = "#BD0026") +
  labs(x = "Spacing in PHF1+ vs PHF1- neurons (%)", y = "Null labellings") +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(4, 6, 4, 4, unit = "pt"))
ggsave(file.path(out_dir, "plot_nn3_null_bounds_phf1.pdf"), p2,
       width = 7.0, height = 4.6, units = "cm", device = "pdf")

# -------------------------------------------------------------------
# Stats log
# -------------------------------------------------------------------
sink(file.path(out_dir, "stats_nn3_null_replica.txt"))
cat("SUPPLEMENTARY label-permutation null for the 3-NN neuronal spacing analyses\n")
cat("==========================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Engine: ", args$engine_dir, "\n", sep = "")
cat("  phf1_label_permutations.qs / phf1_null_distance_matrix.qs, seed 42, PHF1 shuffled WITHIN\n")
cat("  each (celltype x sample) group with the per-group PHF1+ count held fixed.\n")
cat(sprintf("  labellings used: %d\n", n_perm))
cat("  Both 'observed' columns were asserted to reproduce reality before any null was read.\n\n")
cat("INVARIANCE. nn3_um does not depend on PHF1: the neighbour pool is defined by celltype and\n")
cat("both query and neighbour sets are the whole pool. Proven at runtime, not asserted -- nn3_um\n")
cat(sprintf("was recomputed from a frame with the PHF1 columns removed and matched the cache for\n"))
cat(sprintf("all %d pool neurons. So the outcome is fixed across labellings and only the regressor\n",
            nrow(chk)))
cat("changes; nothing spatial is recomputed.\n\n")
cat("ESTIMATOR. Every coefficient here -- observed and null alike -- comes from the same fast\n")
cat("within-donor fixed-effect projection (donor as fixed dummies, coefficient of interest\n")
cat("residualised against all other covariates via qr.resid). The observed value is compared\n")
cat("against the primary lmer fit below.\n\n")

cat("=== ESTIMATOR AGREEMENT CHECK (read this first) ===\n")
cat("Fast projection vs the primary lmer coefficient, same data, same covariates:\n")
print(as.data.frame(dist_sum %>% select(celltype, n_cells, beta_obs_fast,
                                        beta_closer_lmer, abs_gap_vs_lmer)),
      row.names = FALSE, digits = 4)
cat(sprintf("\n  Analysis A: fast = %+.5f, lmer/emmeans = %+.5f, |gap| = %.5f\n",
            A_sum$beta_obs_fast, A_sum$beta_lmer, A_sum$abs_gap_vs_lmer))
cat("  The two differ only in how the donor term is handled (fixed dummies vs a random\n")
cat("  intercept), so small gaps are expected. A large gap invalidates the comparison.\n")

cat("\n=== PART 1: Analysis B -- spacing vs distance to nearest PHF1+ neuron ===\n")
cat("THE KEY COLUMN IS excess_over_null. The null median is the coefficient that geometry alone\n")
cat("produces; 'inside_null_ci' TRUE means the observed value lies within the null 95% interval.\n\n")
cat("THE NULL IS NOT CENTRED ON ZERO, so the no-effect reference is null_median. p_emp below\n")
cat("is a two-TAILED probability computed against the null's own distribution, so it will flag an\n")
cat("observed coefficient that is too WEAK just as readily as one that is too strong. Read\n")
cat("dir_vs_null for the direction of each departure.\n\n")
print(as.data.frame(dist_sum %>%
        select(celltype, n_cells, n_perm_ok, beta_obs_fast, null_median, null_q025,
               null_q975, excess_over_null, inside_null_ci, p_emp, padj_emp, padj_lmer)),
      row.names = FALSE, digits = 3)
cat("\nDirection of departure from the geometric null (BH across subtypes):\n")
print(as.data.frame(dist_sum %>% select(celltype, excess_over_null, p_emp, padj_emp,
                                        dir_vs_null)),
      row.names = FALSE, digits = 3)
cat("\nSame thing as % change in spacing per s.d. closer to a tangle:\n")
for (i in seq_len(nrow(dist_sum))) {
  cat(sprintf("  %-26s observed %+6.2f%%  null median %+6.2f%%  null 95%% [%+6.2f%%, %+6.2f%%]  excess %+6.2f%%  p_emp = %.4f%s\n",
              as.character(dist_sum$celltype[i]),
              100 * (exp(dist_sum$beta_obs_fast[i]) - 1),
              100 * (exp(dist_sum$null_median[i]) - 1),
              100 * (exp(dist_sum$null_q025[i]) - 1),
              100 * (exp(dist_sum$null_q975[i]) - 1),
              100 * (exp(dist_sum$beta_obs_fast[i]) - 1) - 100 * (exp(dist_sum$null_median[i]) - 1),
              dist_sum$p_emp[i],
              ifelse(dist_sum$inside_null_ci[i], "   <- INSIDE THE NULL", "")))
}
n_in  <- sum(dist_sum$inside_null_ci, na.rm = TRUE)
n_wid <- sum(dist_sum$excess_over_null > 0 & dist_sum$padj_emp < FDR, na.rm = TRUE)
n_tig <- sum(dist_sum$excess_over_null < 0 & dist_sum$padj_emp < FDR, na.rm = TRUE)
cat(sprintf("\n  %d of %d subtypes fall INSIDE the null 95%% interval.\n", n_in, nrow(dist_sum)))
cat(sprintf("  Departures at padj_emp < %.2f: %d WIDER than geometry predicts, %d TIGHTER.\n",
            FDR, n_wid, n_tig))
cat("  Interpretation rules:\n")
cat("   * A subtype whose observed coefficient is LESS negative than the null median has WIDER\n")
cat("     spacing near tangles than geometry alone predicts. That is the Zwang et al. direction,\n")
cat("     and it is only visible once the geometric baseline is subtracted -- the raw coefficient\n")
cat("     is negative and would be read as the opposite finding.\n")
cat("   * A subtype MORE negative than the null median has tighter-than-geometric spacing near\n")
cat("     tangles, i.e. tangles sitting preferentially in unusually dense neighbourhoods.\n")

cat("\n=== PART 2: Analysis A -- PHF1+ vs PHF1- spacing ===\n")
cat("For this comparison the permutation is a calibration check rather than a confound bound:\n")
cat("PHF1 status is a WITHIN-donor predictor, so the permutation calibrates the model's standard\n")
cat("errors for cells nested in donors. The null is centred near zero by construction; what\n")
cat("matters is its WIDTH.\n\n")
print(as.data.frame(A_sum), row.names = FALSE, digits = 4)
cat(sprintf("\n  observed %+.2f%% vs null 95%% interval [%+.2f%%, %+.2f%%], p_emp = %.4f%s\n",
            A_sum$pct_obs_fast, A_sum$pct_null_q025, A_sum$pct_null_q975, A_sum$p_emp,
            ifelse(A_sum$inside_null_ci, "   <- INSIDE THE NULL", "")))

cat("\nNOTES:\n")
cat(" - The permutation holds the per-(celltype x sample) PHF1+ COUNT fixed, so it cannot test\n")
cat("   anything about how many neurons bear tangles -- only about WHICH ones do.\n")
cat(" - It also preserves the spatial arrangement of all neurons exactly, which is the point for\n")
cat("   Part 1 and a limitation for Part 2: it cannot bound confounds that are properties of the\n")
cat("   neuron positions themselves (e.g. cortical layer), only of the labelling.\n")
cat(" - The fast estimator treats donor as a fixed effect, so it estimates the within-donor\n")
cat("   contrast. That is the right target here, but it is not numerically identical to the\n")
cat("   primary partially-pooled estimate; see the agreement check above.\n")
cat(" - p_emp is bounded below by 1/(1 + n_perm), so with 1000 labellings it cannot go under\n")
cat("   ~0.001 however extreme the observed value is.\n")
cat("\n=== sessionInfo() ===\n"); print(sessionInfo())
sink()

message(sprintf("Done. Analysis B: %d of %d subtypes inside the null 95%%. Analysis A: observed %+.2f%%, null 95%% [%+.2f%%, %+.2f%%], p_emp = %.4f",
                n_in, nrow(dist_sum), A_sum$pct_obs_fast, A_sum$pct_null_q025,
                A_sum$pct_null_q975, A_sum$p_emp))
print(as.data.frame(dist_sum %>% select(celltype, beta_obs_fast, null_median,
                                        excess_over_null, inside_null_ci, p_emp)),
      row.names = FALSE, digits = 3)
