# ---------------------------------------------------------------------------
# cbln2_trajectory_utils.R
#
# Figure panels: 2B - shared utility (clustering, slingshot and AUC helpers)
# Sourced by the trajectory scripts; no panel of its own.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# cbln2_trajectory_utils.R
#
# Analysis helpers shared by
#   R/cbln2_trajectory_optimise.R        (search over trajectory settings)
#   R/cbln2_slingshot_trajectory.R       (fit and description of the selected trajectory)
#   R/cbln2_trajectory_selection_null.R  (selection-adjusted permutation test)
# so the search, the final fit and the null use identical clustering, slingshot calls
# and AUCs.

SEED_KMEANS   <- 100      # set.seed(100) before kmeans
APPROX_POINTS <- 300      # slingshot approx_points
`%||%` <- function(a, b) if (is.null(a)) b else a

# rank AUC = P(pseudotime of a PHF1+ cell > a PHF1- cell); identical to
# plot_phf1_intensity_validation*.R's auc_fn().
auc_fn <- function(score, pos) {
  ok <- !is.na(score) & !is.na(pos)
  score <- score[ok]; pos <- as.logical(pos[ok])
  n1 <- sum(pos); n0 <- sum(!pos)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
spearman <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "pairwise.complete.obs"))

kmeans_clusters <- function(rd, k) {
  set.seed(SEED_KMEANS)                       # fixed seed and controls for every clustering
  cl <- as.character(stats::kmeans(rd, centers = k, nstart = 10, iter.max = 50)$cluster)
  names(cl) <- rownames(rd)
  cl
}

## Per-cluster PHF1 table: raw fraction, enrichment over the cohort, and the
## donor-standardised O/E (expected = sum over the cluster's cells of their donor's
## PHF1+ rate), so a cluster that is PHF1-poor only because it is built from
## low-burden donors is visible as such.
cluster_table <- function(cl, pos, donor) {
  donor_rate <- tapply(pos, donor, mean)
  expct <- donor_rate[donor]
  tibble(cluster = cl, pos = pos, expct = as.numeric(expct), donor = donor) %>%
    group_by(cluster) %>%
    summarise(n = n(), n_phf1 = sum(pos), frac_phf1 = mean(pos),
              expected_phf1 = sum(expct), n_donors_with_phf1 = n_distinct(donor[pos]),
              .groups = "drop") %>%
    mutate(enrich = frac_phf1 / mean(pos), oe_donor_std = n_phf1 / expected_phf1)
}
pick_cluster <- function(tab, rule, min_n, which = c("min", "max")) {
  which <- match.arg(which)
  t <- tab[tab$n >= min_n, ]
  s <- if (rule == "fraction") cbind(t$frac_phf1, t$oe_donor_std) else cbind(t$oe_donor_std, t$frac_phf1)
  if (which == "max") s <- -s
  t$cluster[order(s[, 1], s[, 2], -t$n)[1]]
}

## slingshot with a fallback chain: if the default call fails, retry with dist.method = 'mnn',
## then on the first 10 dimensions with 'mnn'. The method used is stored as attr 'fit_note'.
run_slingshot <- function(rd, cl, start, label) {
  t0 <- Sys.time()
  call_sling <- function(rd_use, extra = list())
    do.call(slingshot::slingshot, c(list(data = rd_use, clusterLabels = cl, start.clus = start,
                                         approx_points = APPROX_POINTS), extra))
  attempts <- list(
    list(note = "slingshot cluster distances", fn = function() call_sling(rd)),
    list(note = "dist.method = 'mnn'", fn = function() call_sling(rd, list(dist.method = "mnn"))),
    list(note = "first 10 dimensions, dist.method = 'mnn'",
         fn = function() call_sling(rd[, seq_len(min(10, ncol(rd))), drop = FALSE],
                                    list(dist.method = "mnn"))))
  for (a in attempts) {
    sds <- try(a$fn(), silent = TRUE)
    if (!inherits(sds, "try-error")) {
      message(sprintf("  slingshot [%s] start %s: %d lineage(s) in %.1f s [%s]", label, start,
                      length(slingshot::slingLineages(sds)),
                      as.numeric(difftime(Sys.time(), t0, units = "secs")), a$note))
      attr(sds, "fit_note") <- a$note
      return(sds)
    }
    message(sprintf("  [%s] %s failed: %s", label, a$note, sub("^Error[^:]*: ", "", as.character(sds))))
  }
  NULL
}
avg_pt <- function(sds) {
  if (is.null(sds)) return(NULL)
  if ("slingAvgPseudotime" %in% getNamespaceExports("slingshot")) return(slingshot::slingAvgPseudotime(sds))
  pt <- slingshot::slingPseudotime(sds, na = TRUE); w <- slingshot::slingCurveWeights(sds)
  pt[is.na(pt)] <- 0; out <- rowSums(pt * w) / rowSums(w); out[rowSums(w) == 0] <- NA_real_; out
}

## Per-donor AUC table + donor-level summary rows (>= 1 and >= MIN_POS_DONOR PHF1+).
donor_auc <- function(pt, pos, donor) {
  tibble(pt = pt, pos = pos, donor = donor) %>%
    group_by(donor) %>%
    summarise(n = n(), n_phf1 = sum(pos), auc = auc_fn(pt, pos),
              median_pt_phf1 = if (any(pos)) median(pt[pos]) else NA_real_,
              median_pt_neg = median(pt[!pos]), .groups = "drop")
}
summ_auc <- function(tab, min_pos, label) {
  a <- tab$auc[tab$n_phf1 >= min_pos & !is.na(tab$auc)]
  if (length(a) < 3) return(tibble(donor_set = label, n_donors = length(a)))
  tt <- t.test(a, mu = 0.5)
  tibble(donor_set = label, n_donors = length(a), mean_auc = mean(a),
         ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2], t_p_vs_0.5 = tt$p.value,
         n_donors_auc_gt_0.5 = sum(a > 0.5), dz = (mean(a) - 0.5) / sd(a))
}

