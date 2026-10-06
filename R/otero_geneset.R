# ---------------------------------------------------------------------------
# otero_geneset.R
#
# Produces: Table S9
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
#
# WHAT THIS DOES
#   Builds the published Otero-Garcia et al. neurofibrillary-tangle gene sets (layer 2-3
#   AT8-positive neurons, UP and DOWN) from the paper's supplementary workbook, resolving
#   symbols against the CosMx 6k panel so the sets can be scored on this data.
#
#   Used as the external comparator signature in Figure 4B and Figure S4, and as the
#   reference set for the overlap test reported with Table S9.
#
# INPUTS   OteroGarcia2026_S6.xlsx (supplementary table of the source publication)
# OUTPUTS  otero_signatures/otero_at8_signatures.tsv and .rds
# ---------------------------------------------------------------------------
library(tidyverse)
library(readxl)

setwd("<PROJECT_ROOT>/phf1_v2")

# ---- config ----
xlsx      <- "OteroGarcia2026_S6.xlsx"
primary_tabs <- c("Ex1", "Ex2")                 # CBLN2-HOPX-matched L2/3
panex_tabs   <- paste0("Ex", 1:13)              # pan-excitatory robustness
padj_cut  <- 0.05
lfc_cut   <- 0.2                                 # matches Otero-Garcia's threshold
recur_min <- 2                                   # pan-Ex: significant in >= this many clusters
strip_confounds <- TRUE                          # remove RP/MT/lncRNA/MAPT from final sets
outdir    <- "otero_signatures"
dir.create(outdir, showWarnings = FALSE)

lfc_col  <- "AD-AT8+ vs AD-AT8- avg_log2FC"
padj_col <- "AD-AT8+ vs AD-AT8- p_val_adj"

# ---- repair Excel-corrupted gene symbols -------------------------------
# Excel autoconverts some gene symbols to dates, read back as either date
# strings ("2021-09-06") or numeric serials (44446). Map the known ones.
# Verify against the paper's non-Excel list if you can; unrecognised
# corruptions are dropped rather than mis-assigned.
date_fixes <- c(
  "44446"      = "SEPT9",   "2021-09-06" = "SEPT9",  "6-Sep" = "SEPT9",
  "44261"      = "MARCH6",  "2021-03-06" = "MARCH6", "6-Mar" = "MARCH6"
)
repair_symbol <- function(x) {
  x <- as.character(x)
  hit <- x %in% names(date_fixes)
  x[hit] <- date_fixes[x[hit]]
  x
}

# ---- read one tab -> tidy sig genes for the AT8+ vs AT8- contrast ----
read_tab <- function(tab) {
  read_excel(xlsx, sheet = tab, col_types = "text") %>%   # force text: no silent coercion
    transmute(
      cluster = tab,
      gene    = repair_symbol(Gene),
      lfc     = as.numeric(.data[[lfc_col]]),
      padj    = as.numeric(.data[[padj_col]])
    ) %>%
    filter(!is.na(lfc), !is.na(padj)) %>%              # drops the empty AT8 rows
    filter(padj < padj_cut, abs(lfc) > lfc_cut) %>%
    mutate(dir = if_else(lfc > 0, "up", "down"))
}

# =========================================================
# 1. PRIMARY: Ex1 U Ex2, union with sign-consistency
# =========================================================
prim <- map_dfr(primary_tabs, read_tab)

# resolve sign per gene: keep only genes with consistent direction where significant
prim_signed <- prim %>%
  group_by(gene) %>%
  summarise(
    n_sig    = n(),
    n_up     = sum(dir == "up"),
    n_down   = sum(dir == "down"),
    mean_lfc = mean(lfc),
    .groups  = "drop"
  ) %>%
  filter(!(n_up > 0 & n_down > 0)) %>%                 # drop opposite-sign conflicts
  mutate(
    direction   = if_else(n_up > 0, "up", "down"),
    shared_core = n_sig == length(primary_tabs)        # sig in BOTH clusters
  )

primary_up   <- prim_signed %>% filter(direction == "up")   %>% pull(gene)
primary_down <- prim_signed %>% filter(direction == "down") %>% pull(gene)
core_up      <- prim_signed %>% filter(direction == "up",   shared_core) %>% pull(gene)
core_down    <- prim_signed %>% filter(direction == "down", shared_core) %>% pull(gene)

# =========================================================
# 2. PAN-EX ROBUSTNESS: recurrent across >= recur_min Ex clusters, same sign
# =========================================================
panex <- map_dfr(panex_tabs, read_tab)

panex_signed <- panex %>%
  group_by(gene, dir) %>%
  summarise(n_clusters = n(), .groups = "drop") %>%
  pivot_wider(names_from = dir, values_from = n_clusters,
              values_fill = 0, names_prefix = "n_") %>%
  # ensure both columns exist even if one direction never occurs
  { if (!"n_up"   %in% names(.)) mutate(., n_up = 0)   else . } %>%
  { if (!"n_down" %in% names(.)) mutate(., n_down = 0) else . } %>%
  filter(xor(n_up >= recur_min, n_down >= recur_min)) %>%   # recurrent & unambiguous
  mutate(direction = if_else(n_up >= recur_min, "up", "down"))

panex_up   <- panex_signed %>% filter(direction == "up")   %>% pull(gene)
panex_down <- panex_signed %>% filter(direction == "down") %>% pull(gene)

# =========================================================
# 3. OPTIONAL confound strip (activity/depth-driven + circularity)
# =========================================================
strip_set <- function(g) {
  g[!(str_detect(g, "^RP[LS]") |                       # ribosomal proteins
        str_detect(g, "^MT-")    |                       # mito-encoded
        g %in% c("MALAT1","XIST","MEG3","NEAT1","KCNQ1OT1") |  # abundant lncRNAs
        g == "MAPT")]                                    # circularity pre-empt vs tau anchor
}
if (strip_confounds) {
  primary_up   <- strip_set(primary_up)
  primary_down <- strip_set(primary_down)
  core_up      <- strip_set(core_up)
  core_down    <- strip_set(core_down)
  panex_up     <- strip_set(panex_up)
  panex_down   <- strip_set(panex_down)
}

# =========================================================
# report + write
# =========================================================
cat("PRIMARY (Ex1 U Ex2, sign-consistent):\n")
cat(sprintf("  up   = %d  (shared core %d)\n", length(primary_up),   length(core_up)))
cat(sprintf("  down = %d  (shared core %d)\n", length(primary_down), length(core_down)))
n_conflict <- prim %>% group_by(gene) %>%
  summarise(u = any(dir=="up"), d = any(dir=="down"), .groups="drop") %>%
  filter(u & d) %>% nrow()
cat(sprintf("  dropped opposite-sign conflicts: %d\n", n_conflict))
cat(sprintf("PAN-EX (recurrent >= %d clusters): up = %d, down = %d\n",
            recur_min, length(panex_up), length(panex_down)))
cat(sprintf("confound strip (RP/MT/lncRNA/MAPT): %s\n", strip_confounds))

sig_list <- list(
  otero_L23_up        = primary_up,
  otero_L23_down      = primary_down,
  otero_L23_core_up   = core_up,
  otero_L23_core_down = core_down,
  otero_panEx_up      = panex_up,
  otero_panEx_down    = panex_down
)
saveRDS(sig_list, file.path(outdir, "otero_at8_signatures.rds"))

# also flat tsv (gene, set) for portability
imap_dfr(sig_list, ~ tibble(gene = .x, set = .y)) %>%
  write_tsv(file.path(outdir, "otero_at8_signatures.tsv"))

cat("Written:", file.path(outdir, "otero_at8_signatures.rds"), "and .tsv\n")