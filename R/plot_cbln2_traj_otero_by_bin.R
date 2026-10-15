#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_cbln2_traj_otero_by_bin.R
#
# Figure panels: 2D
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# plot_cbln2_traj_otero_by_bin.R
#
# Otero-Garcia Exc1/2 (L2/3) AT8-UP module score along the CBLN2 tangle-associated trajectory,
# in the style of plot_cbln2_traj_phf1_by_bin_mean_5.08cm.pdf (same x-axis length, 6 cm tall here),
# with the score in tangle-bearing and tangle-free neurons alongside.
#
# SCORE: genesets/fig2_contrast/otero_l23_up.txt (written by prepare_fig2_contrast_genesets.R; the
# 319 panel genes of otero_signatures/otero_at8_signatures.rds $otero_L23_up, matched by plain
# intersect), Seurat AddModuleScore on the SCT assay of seu_PHF1.rds subset to CBLN2 cells,
# ctrl = 100, nbin = 24 (falls back to 15 if 24 bins fail), seed 42, the set in its own call.
#
# STAGES
#   1. score (needs seu_PHF1.rds, 8.4 GB): writes the per-cell score cache
#      plots/cbln2_slingshot_trajectory/source_data_cbln2_traj_otero_score_percell.tsv
#   2. plot: reads the cache + the trajectory's per-cell pseudotime and bins.
#   The script does stage 1 only if the cache is missing (or --rescore).
#
# PLOT
#   left   tangle-status panel: per donor, the mean score over all that donor's tangle-free and
#          tangle-bearing cells (NOT stratified by pseudotime), as a box of the donor means with
#          each donor's two means joined; Wilcoxon signed-rank p (exact, two-sided) on the
#          donor-paired means. Donors = those with >= min_pos_donor PHF1+ cells, as the % tangle plot.
#   right  x = the trajectory's 10 equal-width pseudotime bins, same as the % tangle plot: DONOR
#          mean of per-donor bin means, 95% t CI, black line, grey ribbon; a bar over bins 1..10
#          with the Friedman p across bins (bins within donor).
#   Both panels share one y axis.
#
# RUN:  Rscript R/plot_cbln2_traj_otero_by_bin.R
# OUTPUT: plots/cbln2_slingshot_trajectory/{plot,source_data,stats}_cbln2_traj_otero_by_bin.*
#   + plot_/source_data_cbln2_traj_otero_by_bin_donor.*: per-donor curves coloured by Braak.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(argparse)
})
setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")
source("R/cbln2_embedding_utils.R")
source("R/fixed_panel_utils.R")      # BIN_PANEL_W_CM: same x-axis length as the % tangle bin plot

parser <- ArgumentParser(description = "Otero Exc1/2 UP score along the CBLN2 trajectory")
parser$add_argument("--seu", default = "seu_PHF1.rds")
parser$add_argument("--rescore", action = "store_true", default = FALSE)
parser$add_argument("--min_pos_donor", type = "integer", default = 5L,
                    help = "donors for the front line (>= this many PHF1+), as the % tangle plot")
args <- parser$parse_args()

CT          <- CBLN2_CT
TAG         <- ct_tag(CT)
DIR         <- file.path("plots", paste0(TAG, "_slingshot_trajectory"))
NAME        <- paste0(TAG, "_traj_otero_by_bin")
SCORE_CACHE <- file.path(DIR, paste0("source_data_", TAG, "_traj_otero_score_percell.tsv"))
GENESET     <- "genesets/fig2_contrast/otero_l23_up.txt"
CTRL <- 100; NBIN <- 24; SEED <- 42       # AddModuleScore settings
H_CM <- 6                    # height (cm)
BOX_PANEL_W_CM <- 1.70       # width of the tangle-status side panel (cm); canvas width follows
COLS <- c("All neurons" = "black", "Tangle-free" = "#0072B2", "Tangle-bearing" = "#D55E00")

## 1. Score ----------------------------------------------------------------------------------
if (args$rescore || !file.exists(SCORE_CACHE)) {
  suppressPackageStartupMessages(library(Seurat))
  genes <- unique(readLines(GENESET)); genes <- genes[nzchar(genes)]
  message("Loading ", args$seu, " ...")
  seu <- readRDS(args$seu)
  DefaultAssay(seu) <- "SCT"
  present <- intersect(genes, rownames(seu[["SCT"]]))
  stopifnot(length(present) == length(genes))      # the file IS the on-panel set (319)
  seu_g <- seu[, which(seu$celltype == CT)]
  rm(seu); invisible(gc())
  DefaultAssay(seu_g) <- "SCT"
  set.seed(SEED)
  seu_g <- tryCatch(
    AddModuleScore(seu_g, features = list(otero_l23_up = present), name = "OTERO_",
                   assay = "SCT", ctrl = CTRL, nbin = NBIN, seed = SEED),
    error = function(e) { message("nbin = ", NBIN, " failed (", conditionMessage(e), "); retrying at 15")
      set.seed(SEED)
      AddModuleScore(seu_g, features = list(otero_l23_up = present), name = "OTERO_",
                     assay = "SCT", ctrl = CTRL, nbin = 15, seed = SEED) })
  write.table(tibble(cell_id = colnames(seu_g), otero_l23_up = seu_g$OTERO_1,
                     PHF1 = as.logical(as.character(seu_g$PHF1)), sample_id = as.character(seu_g$sample_id)),
              SCORE_CACHE, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Scored ", ncol(seu_g), " cells on ", length(present), " genes -> ", SCORE_CACHE)
  rm(seu_g); invisible(gc())
}

## 2. Join to pseudotime ---------------------------------------------------------------------
sc <- read.delim(SCORE_CACHE, colClasses = c(cell_id = "character", sample_id = "character"))
pc <- read.delim(ct_percell_file(CT), colClasses = c(cell_id = "character", sample_id = "character"))
df <- pc %>% select(cell_id, sample_id, Braak, PHF1, pt, bin) %>%
  inner_join(sc %>% select(cell_id, otero_l23_up, PHF1_sc = PHF1), by = "cell_id")
stopifnot(nrow(df) == nrow(pc), all(as.logical(df$PHF1) == as.logical(df$PHF1_sc)))
df$PHF1 <- as.logical(df$PHF1)
n_bins <- max(df$bin)
n_pos_d <- tapply(df$PHF1, df$sample_id, sum)
head_donors <- names(n_pos_d)[n_pos_d >= args$min_pos_donor]

tci <- function(x) { n <- length(x); m <- mean(x); h <- if (n > 1) qt(0.975, n - 1) * sd(x) / sqrt(n) else NA_real_
  c(mean = m, lo = m - h, hi = m + h, n = n) }
front <- df %>% filter(sample_id %in% head_donors) %>% group_by(sample_id, bin) %>%
  summarise(m = mean(otero_l23_up), .groups = "drop") %>% group_by(bin) %>%
  summarise(mean = mean(m), lo = tci(m)[["lo"]], hi = tci(m)[["hi"]], n = n(), .groups = "drop") %>%
  mutate(group = "All neurons", level = "donor mean (95% t CI)")
## Reference levels across the whole distribution (independent of pseudotime), donor as unit.
donor_grp <- df %>% filter(sample_id %in% head_donors) %>%
  mutate(group = ifelse(PHF1, "Tangle-bearing", "Tangle-free")) %>%
  group_by(sample_id, group) %>% summarise(m = mean(otero_l23_up), n_cells = n(), .groups = "drop")
ref <- donor_grp %>% group_by(group) %>%
  summarise(mean = mean(m), lo = tci(m)[["lo"]], hi = tci(m)[["hi"]], n = n(), .groups = "drop") %>%
  mutate(level = "donor mean over all cells (95% t CI)")

## 3. Plot: curve (fixed 4.20 cm x-axis, = the % tangle plot) + tangle-status panel, same y axis ------
## Side panel = the donor-level comparison, annotated with the Wilcoxon signed-rank p (scores are
## bounded, so nonparametric): a box of the donor means per group, each donor's two means joined.
## Donor means (not cells) so it shares the curve's y range.
paired_w <- donor_grp %>% select(sample_id, group, m) %>% pivot_wider(names_from = group, values_from = m)
pdiff <- paired_w[["Tangle-bearing"]] - paired_w[["Tangle-free"]]
# Exact two-sided Wilcoxon signed-rank test on the donor-paired means.
ptest <- wilcox.test(paired_w[["Tangle-bearing"]], paired_w[["Tangle-free"]], paired = TRUE, exact = TRUE)
p_txt <- sprintf("p = %s", formatC(ptest$p.value, format = "g", digits = 2))
## Tests across pseudotime bins: unit = donor x bin mean (the points the curve averages).
## Friedman (bins within donor) is annotated; Kruskal-Wallis (treats the donors' bins as
## independent) is kept in the stats log.
donor_bin <- df %>% filter(sample_id %in% head_donors) %>% group_by(sample_id, Braak, bin) %>%
  summarise(mean = mean(otero_l23_up), n_cells = n(), n_phf1 = sum(PHF1), .groups = "drop") %>%
  mutate(Braak = factor(as.character(Braak), levels = braak_levels))
stopifnot(all(table(donor_bin$bin) == length(head_donors)))       # every donor in every bin
kw <- kruskal.test(mean ~ factor(bin), data = donor_bin)
kw_eps2 <- unname(kw$statistic) / (nrow(donor_bin) - 1)            # epsilon^2 = H / (n - 1)
fr <- friedman.test(mean ~ bin | sample_id, data = donor_bin)
# Annotated test = Friedman, label just "p = " (test named in the caption); KW in the log.
# plotmath label: small p as "1.6 x 10^-7", otherwise plain; drawn with parse = TRUE.
p_plotmath <- function(p) {
  if (p >= 1e-3) return(sprintf('p == "%s"', formatC(p, format = "g", digits = 2)))
  m <- strsplit(formatC(p, format = "e", digits = 1), "e")[[1]]
  sprintf('p == "%s" %%*%% 10^%d', m[1], as.integer(m[2]))
}
kw_txt <- p_plotmath(fr$p.value)
y_top <- max(front$hi, donor_grp$m, na.rm = TRUE)
y_lim <- c(min(0, front$lo, na.rm = TRUE), y_top * 1.14)
yb <- y_top * 1.05                                                  # bracket height
dm <- donor_grp %>% mutate(group = factor(group, levels = c("Tangle-free", "Tangle-bearing")))

p_curve <- ggplot(front, aes(bin, mean)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey75", alpha = 0.6) +
  geom_line(linewidth = 0.6) + geom_point(size = 1) +
  annotate("segment", x = 1, xend = n_bins, y = yb, yend = yb, linewidth = 0.3) +
  # no end ticks: a plain bar spans all bins, so it does not read as a first-vs-last contrast
  annotate("text", x = (1 + n_bins) / 2, y = yb, label = kw_txt, parse = TRUE, vjust = -0.4, size = 2.1) +
  scale_x_continuous(breaks = seq_len(n_bins)) +
  scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0.02, 0))) +
  labs(x = "Pseudotime bin", y = NULL) +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        plot.margin = margin(4, 4, 4, 2))
p_box <- ggplot(dm, aes(group, m)) +
  geom_boxplot(aes(fill = group), width = 0.6, outlier.shape = NA, alpha = 0.45, linewidth = 0.3, colour = "grey20") +
  geom_line(aes(group = sample_id), colour = "grey40", linewidth = 0.3) +
  geom_point(aes(fill = group), shape = 21, colour = "black", size = 1.2, stroke = 0.25) +
  annotate("segment", x = 1, xend = 2, y = yb, yend = yb, linewidth = 0.3) +
  annotate("segment", x = c(1, 2), xend = c(1, 2), y = yb, yend = yb - 0.025 * diff(y_lim), linewidth = 0.3) +
  annotate("text", x = 1.5, y = yb, label = p_txt, vjust = -0.4, size = 2.1) +
  scale_fill_manual(values = COLS[c("Tangle-free", "Tangle-bearing")], guide = "none") +
  scale_x_discrete(labels = c("Tangle-free" = "Tangle-\nfree", "Tangle-bearing" = "Tangle-\nbearing")) +
  scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0.02, 0))) +
  labs(x = NULL, y = "Otero-Garcia Exc1/2 UP score") +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 6, lineheight = 0.85),
        plot.margin = margin(4, 2, 4, 4))
## Box panel left, curve right; shared y axis drawn once, on the box panel.
pw <- patchwork::wrap_plots(p_box, p_curve, nrow = 1) +
  patchwork::plot_layout(widths = grid::unit(c(BOX_PANEL_W_CM, BIN_PANEL_W_CM), "cm"))
W_USED <- save_patchwork_exact_width(pw, file.path(DIR, paste0("plot_", NAME, ".pdf")), H_CM)
message(sprintf("Curve x-axis %.2f cm (= %% tangle plot), box panel %.2f cm, canvas %.2f x %.2f cm",
                BIN_PANEL_W_CM, BOX_PANEL_W_CM, W_USED, H_CM))
plot_df <- bind_rows(front %>% mutate(bin = as.numeric(bin), panel = "curve"),
                     dm %>% transmute(group = as.character(group), level = "donor mean over all cells",
                                      sample_id, mean = m, panel = "tangle_status"))
write.table(plot_df %>% select(panel, group, level, any_of("sample_id"), bin, mean, lo, hi, n),
            file.path(DIR, paste0("source_data_", NAME, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

## 3a. Donor-by-bin: every donor's own curve (thin, coloured by Braak) + the donor mean (black) --------
## Same donors as the curve above (donor_bin, built for the bin tests), so the black line is exactly
## their mean; same 4.20 cm x-axis, 6 cm tall.
p_donor <- ggplot() +
  geom_line(data = donor_bin, aes(bin, mean, group = sample_id, colour = Braak), linewidth = 0.4, alpha = 0.9) +
  geom_point(data = donor_bin, aes(bin, mean, colour = Braak), size = 0.6, alpha = 0.9) +
  geom_line(data = front, aes(bin, mean), colour = "black", linewidth = 0.7) +
  geom_point(data = front, aes(bin, mean), colour = "black", size = 1) +
  scale_colour_manual(values = braak_line_palette, name = "Braak", drop = FALSE) +
  scale_x_continuous(breaks = seq_len(n_bins)) +
  labs(x = "Pseudotime bin", y = "Otero-Garcia Exc1/2 UP score") +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        legend.key.height = grid::unit(7, "pt"), legend.key.width = grid::unit(10, "pt"),
        legend.box.spacing = grid::unit(3, "pt"), plot.margin = margin(4, 4, 4, 4))
W_DONOR <- save_patchwork_exact_width(
  patchwork::wrap_plots(p_donor) + patchwork::plot_layout(widths = grid::unit(BIN_PANEL_W_CM, "cm")),
  file.path(DIR, paste0("plot_", NAME, "_donor.pdf")), H_CM)
message(sprintf("Donor-by-bin plot: x-axis %.2f cm, canvas %.2f x %.2f cm", BIN_PANEL_W_CM, W_DONOR, H_CM))
write.table(bind_rows(donor_bin %>% mutate(Braak = as.character(Braak), line = "donor"),
                      front %>% transmute(bin, mean, line = "donor mean (black)")),
            file.path(DIR, paste0("source_data_", NAME, "_donor.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

cells <- df %>% filter(sample_id %in% head_donors) %>%
  mutate(group = factor(ifelse(PHF1, "Tangle-bearing", "Tangle-free"), levels = c("Tangle-free", "Tangle-bearing")))

## 4. Stats ---------------------------------------------------------------------------------------
first_last <- df %>% filter(sample_id %in% head_donors, bin %in% c(1, n_bins)) %>%
  group_by(sample_id, bin) %>% summarise(m = mean(otero_l23_up), .groups = "drop") %>%
  pivot_wider(names_from = bin, values_from = m)
dd <- first_last[[as.character(n_bins)]] - first_last[["1"]]
tt <- t.test(dd)
paired <- donor_grp %>% select(sample_id, group, m) %>% pivot_wider(names_from = group, values_from = m)
pd <- paired[["Tangle-bearing"]] - paired[["Tangle-free"]]
pt_t <- t.test(pd)
sink(file.path(DIR, paste0("stats_", NAME, ".txt")))
cat("Otero-Garcia Exc1/2 UP score along the", CT, "trajectory\n", strrep("=", 72), "\n", sep = "")
cat("Run:", format(Sys.time()), "\nScore cache:", SCORE_CACHE, "\nPseudotime:", ct_percell_file(CT), "\n")
cat(sprintf("Score: AddModuleScore, SCT, CBLN2 cells, %s (%d genes), ctrl %d, nbin %d, seed %d\n",
            GENESET, length(readLines(GENESET)), CTRL, NBIN, SEED))
cat("Front line: donor mean of per-donor bin means over", length(head_donors), "donors (>=",
    args$min_pos_donor, "PHF1+):", paste(head_donors, collapse = ", "), "\n")
cat("Reference lines: per donor, the mean score over ALL that donor's tangle-free / tangle-bearing cells\n")
cat("(not stratified by pseudotime), then the mean over the same donors with a 95% t CI.\n\n")
cat("Front line (by pseudotime bin):\n")
print(as.data.frame(front %>% select(bin, mean, lo, hi, n)), digits = 3, row.names = FALSE)
cat("\nReference levels:\n")
print(as.data.frame(ref %>% select(group, mean, lo, hi, n)), digits = 3, row.names = FALSE)
cat("\nPer donor:\n"); print(as.data.frame(donor_grp), digits = 3, row.names = FALSE)
cat("\nDonor-by-bin (plot_", NAME, "_donor.pdf): per-donor bin means; smallest donor x bin cell count ",
    min(donor_bin$n_cells), "\n", sep = "")
print(as.data.frame(donor_bin %>% select(sample_id, Braak, bin, mean) %>%
                      pivot_wider(names_from = bin, values_from = mean)), digits = 3, row.names = FALSE)
cat("\nCell-level distribution (descriptive -- the test is donor-paired):\n")
print(as.data.frame(cells %>% group_by(group) %>%
        summarise(n_cells = n(), median = median(otero_l23_up), q25 = quantile(otero_l23_up, .25),
                  q75 = quantile(otero_l23_up, .75), sd = sd(otero_l23_up), .groups = "drop")), digits = 3, row.names = FALSE)
w_t <- suppressWarnings(wilcox.test(pd, mu = 0))
cat(sprintf("Wilcoxon signed-rank on the donor differences (THE p ANNOTATED ON THE BOX PANEL): V = %g, p = %.3g\n",
            unname(w_t$statistic), w_t$p.value))
cat(sprintf("\nTangle-bearing minus tangle-free, donor-paired: %+.4f [%.4f, %.4f], dz = %.2f, %d/%d donors higher, paired t p = %.3g\n",
            mean(pd), pt_t$conf.int[1], pt_t$conf.int[2], mean(pd) / sd(pd), sum(pd > 0), length(pd), pt_t$p.value))
cat(sprintf("\nKruskal-Wallis across the %d pseudotime bins (log only; unit = donor x bin mean, n = %d):\n  H = %.2f, df = %d, p = %.3g, epsilon^2 = %.2f\n",
            n_bins, nrow(donor_bin), unname(kw$statistic), unname(kw$parameter), kw$p.value, kw_eps2))
cat(sprintf("  Friedman, bins within donor (THE p ANNOTATED ON THE BIN CURVE; KW ignores that each donor contributes every bin):\n  chi2 = %.2f, df = %d, p = %.3g, Kendall's W = %.2f\n",
            unname(fr$statistic), unname(fr$parameter), fr$p.value, unname(fr$statistic) / (length(head_donors) * (n_bins - 1))))
cat(sprintf("\nEffect size, last (bin %d) minus first (bin 1), donor-paired: %+.4f [%.4f, %.4f], dz = %.2f, %d/%d donors higher\n",
            n_bins, mean(dd), tt$conf.int[1], tt$conf.int[2], mean(dd) / sd(dd), sum(dd > 0), length(dd)))
cat("\n\n")
print(sessionInfo())
sink()
message("Written: ", file.path(DIR, paste0("plot_", NAME, ".pdf")))
