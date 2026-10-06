# ---------------------------------------------------------------------------
# nDEG_by_celltype_linear.R
#
# Figure panels: 3A
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
library(tidyverse)
library(stringr)
library(ggplot2)

# Root directory with linear distance DE results
base_dir <- "<PROJECT_ROOT>/phf1_v2/deg/de_linear_distance"
# base_dir <- "<PROJECT_ROOT>/phf1_v2/deg/de_linear_distance_raw"

# Shared Nature-compact theme + celltype ordering (match the distance plots)
source(file.path(dirname(dirname(base_dir)), "R", "palettes.R"))  # fig_theme, neuron_order, glia_order

# Significance threshold
PADJ_THRESHOLD <- 0.1

# logFC threshold (per SD of log(um)).
# Unlike categorical DEG there is no natural fold-change cutoff here —
# logFC is a continuous distance coefficient, not a group contrast.
# Set to 0 to count all significant genes regardless of effect size,
# or raise to a small value (e.g. 0.05) to require a minimum slope.
LOGFC_THRESHOLD <- 0

celltypes <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

# Main linear DEG results file: *_dist_to_phf1_um_scaled.tsv
deg_files <- unlist(lapply(celltypes, function(ct) {
  list.files(ct, pattern = "_dist_to_phf1_um_scaled\\.tsv$", full.names = TRUE, recursive = FALSE)
}))

if (length(deg_files) == 0) stop("No *_dist_to_phf1_um_scaled.tsv files found under ", base_dir)

# Function to count DEGs per file
# Direction: negative logFC = lower expression further from PHF1+ = upregulated near PHF1+
#            positive logFC = higher expression further from PHF1+ = downregulated near PHF1+
count_degs <- function(file) {
  df <- read_tsv(file, show_col_types = FALSE)

  df <- df %>% dplyr::filter(!is.na(padj), padj < PADJ_THRESHOLD, abs(logFC) > LOGFC_THRESHOLD)

  down   <- sum(df$logFC > LOGFC_THRESHOLD)   # higher expression with increasing distance
  up <- sum(df$logFC < -LOGFC_THRESHOLD)  # lower expression with increasing distance (near-PHF1+ upregulated)

  tibble(
    CellType = basename(dirname(file)),
    File     = basename(file),
    Up       = up,
    Down     = down
  )
}

# Combine all counts
deg_counts <- bind_rows(lapply(deg_files, count_degs))

# Reshape for plotting — Down counts negated to plot on the left of the axis.
# Neurons only: ordering on neuron_order sets glia/Other/Unassigned to NA, dropped.
deg_counts_long <- deg_counts %>%
  tidyr::pivot_longer(cols = c("Up", "Down"),
                      names_to = "Direction",
                      values_to = "Count") %>%
  mutate(Count = ifelse(Direction == "Down", -Count, Count),
         Direction = factor(Direction, levels = c("Up", "Down")),
         CellType  = factor(CellType, levels = neuron_order)) %>%
  dplyr::filter(!is.na(CellType))

# Round `x` up to the nearest 1/2/5 x 10^n; never below 1, since counts are
# integers and a fractional tick step would be meaningless.
nice_step <- function(x) {
  if (!is.finite(x) || x <= 0) return(1)
  mag  <- 10^floor(log10(x))
  mult <- c(1, 2, 5, 10)
  max(1, mag * mult[which(mult >= x / mag)[1]])
}

# x axis whose outermost labelled tick on each side lies strictly beyond the
# largest bar on that side, so every bar ends inside the numbered range and the
# reader can read a count off the axis rather than extrapolating past the last
# label. Step is shared by both sides and targets ~5 intervals overall.
count_axis <- function(counts) {
  up   <- max(c(counts, 0))
  down <- -min(c(counts, 0))
  step <- nice_step((up + down) / 5)
  hi   <-  step * (floor(up   / step) + 1)
  lo   <- -step * (floor(down / step) + 1)
  scale_x_continuous(limits = c(lo, hi), breaks = seq(lo, hi, by = step),
                     labels = abs, expand = expansion(mult = 0.02))
}

# Bar plot (celltypes on the y axis, horizontal bars; matches the nPathway plots).
# neuron_order has its first entry at the top (scale_y_discrete limits = rev).
# Override fig_theme's rotated x text + wide left margin -- both were for
# celltypes on the x axis; here the x axis is just the numeric count.
ggplot(deg_counts_long, aes(x = Count, y = CellType, fill = Direction)) +
  geom_col(position = "stack") +
  geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.4) +
  scale_fill_manual(values = c("Up" = "#DC0000FF", "Down" = "#3C5488FF"),
                    labels = c("Up" = "Up near PHF1+", "Down" = "Down near PHF1+"),
                    name = NULL) +
  scale_y_discrete(limits = rev(neuron_order), expand = expansion(add = 0.6)) +
  count_axis(deg_counts_long$Count) +
  theme_classic(base_size = 7) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        plot.margin = margin(t = 4, r = 6, b = 4, l = 2, unit = "pt")) +
  labs(x = "Number of DEGs", y = NULL)

ggsave(filename = file.path(base_dir, "nDEGs_by_celltype_linear.pdf"),
       width = 9.9, height = 5.3, units = "cm", device = "pdf")
