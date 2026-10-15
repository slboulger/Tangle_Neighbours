#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# prepare_fig2_contrast_genesets.R
#
# Figure panels: 2D - upstream step
# Writes genesets/fig2_contrast/otero_l23_up.txt, the Otero-Garcia L2/3 AT8-up set
# (319 panel genes) that plot_cbln2_traj_otero_by_bin.R scores for 2D.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# prepare_fig2_contrast_genesets.R
#
# Builds several pre-specified gene sets as flat lists of CosMx panel probes, plus a
# membership table and a provenance file. No cells are read: inputs are the Otero-Garcia
# supplementary workbooks, the cached Reactome GMT and one small per-celltype SCE for the
# panel rownames.
#
# ONLY `otero_l23_up` IS USED IN THE PAPER. It is the same 319-gene Otero-Garcia set scored
# in 3B, and is the score drawn in 2D. The other sets below are written alongside it but are
# not drawn in any figure panel.
#
# SETS (key -> source)
#   otero_l23_up        Otero-Garcia 2022 Table S6, L2/3 AT8-up (Ex1 u Ex2). Copied
#                       verbatim from otero_signatures/otero_at8_signatures.rds
#                       [["otero_L23_up"]] (built by R/otero_geneset.R) so it is the
#                       exact set scored in 3B. That file carries "SEPT9" where the
#                       source cell is SEPT7; both are off the panel, so the scored set
#                       is unaffected.
#   otero_s2_core_up    Table S2: genes in the DEG block (padj < 0.05, |logFC| > 0.1,
#   otero_s2_core_down  >= 20% cells) of ALL FIVE susceptible subtypes Ex1/2/3/7/10
#                       (the five named in Table S5), with the same sign in all five.
#   otero_s5_presyn     Table S5, Sheet5: SynGO SYNGO:presynprocess "process in the
#                       presynapse", the 63-gene headline set. The SynGO input was
#                       genes dysregulated in either direction, so this set is not
#                       uniformly up: its per-gene direction (from Sheet1 logFCs across
#                       the five clusters) is written to the membership table and the
#                       up/down/mixed split is printed and asserted.
#   otero_s5_presyn_up  The consistently-up subset of the above (logFC > 0 in every
#                       cluster where it is DE).
#   reactome_apoptosis  Reactome "Apoptosis" (R-HSA-109581) from the cached Enrichr
#                       Reactome_Pathways_2024 GMT (168 genes).
#   vatpase_foiani      Foiani et al. 2026 V-ATPase list (11 subunits).
#
# SYMBOLS
#   Excel date-mangled cells (serial numbers or "7-Sep" strings) are converted to the
#   current HGNC symbol (SEPTIN*/MARCHF*/DELEC1). A symbol not on the panel is then
#   resolved through org.Hs.eg.db ALIAS -> SYMBOL (accepted only if UNIQUE) and kept
#   if the resolved symbol is on the panel (e.g. TCEB2 -> ELOB). Panel matching is at
#   the PROBE level via expand_probe_names() (CosMx packs genes, e.g. "FKBP1A/C").
#   EXCEPTION otero_l23_up: matched by plain intersect, exactly as scored in 3B
#   (319 probes). The full rule would add FKBP1A/C, H2AZ1 (H2AFZ) and H3-3A (H3F3A);
#   those are logged as not_used.
#
# OUTPUT  genesets/fig2_contrast/
#   <key>.txt                     PANEL PROBES to score, newline-delimited
#   fig2_contrast_membership.tsv  set, source symbol, cleaned symbol, fix, probe,
#                                 on_panel, source direction
#   fig2_contrast_summary.tsv     per set: n source, n on panel (probes)
#   fig2_contrast_provenance.txt
#
# RUN  Rscript R/prepare_fig2_contrast_genesets.R

suppressPackageStartupMessages({
  library(readxl); library(qs); library(SingleCellExperiment)
})

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/cosmx_deg_module_utils.R")  # expand_probe_names()

S2_XLSX   <- "../OteroGarcia_S2.xlsx"
S5_XLSX   <- "../OteroGarcia_S5.xlsx"
OTERO_RDS <- "otero_signatures/otero_at8_signatures.rds"
GMT       <- "genesets/enrichr_gmt/Reactome_Pathways_2024.gmt"
PANEL_SCE <- "celltype_sce_neighbours/Exc-ET-L5-SPON1-FGD4_sce_neighbours.qs"
OUT_DIR   <- "genesets/fig2_contrast"
S2_CORE_SUBTYPES <- c("Ex1", "Ex2", "Ex3", "Ex7", "Ex10")
VATPASE <- c("ATP6V1G1", "ATP6V1G2", "ATP6V1B2", "ATP6V1C1", "ATP6V1E1", "ATP6V1A",
             "ATP6V0D1", "ATP6V1F", "ATP6V0A1", "ATP6V1H", "ATP6V1D")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
for (f in c(S2_XLSX, S5_XLSX, OTERO_RDS, GMT, PANEL_SCE))
  if (!file.exists(f)) stop("missing input: ", f, call. = FALSE)

## ---- panel (probe level) -------------------------------------------------------
panel <- rownames(qs::qread(PANEL_SCE))
stopifnot("panel should be the 6175-probe CosMx 6k panel" = length(panel) == 6175)
probe_members <- setNames(lapply(panel, expand_probe_names), panel)
sym2probe <- function(s) {
  hit <- names(probe_members)[vapply(probe_members, function(m) s %in% m, logical(1))]
  if (length(hit) == 0) NA_character_ else hit[1]
}
panel_syms <- unique(unlist(probe_members))

## ---- symbol cleaning -----------------------------------------------------------
DATE_PREFIX <- c(Sep = "SEPTIN", Mar = "MARCHF", Dec = "DELEC")
fix_date <- function(g) {
  out <- g; fix <- rep("none", length(g))
  ser <- grepl("^[0-9]{5}$", g)                         # Excel serial read as text
  if (any(ser)) {
    d <- as.Date(as.numeric(g[ser]), origin = "1899-12-30")
    mon <- format(d, "%b"); day <- as.integer(format(d, "%d"))
    ok <- mon %in% names(DATE_PREFIX)
    out[ser][ok] <- paste0(DATE_PREFIX[mon[ok]], day[ok]); fix[ser][ok] <- "excel_date"
  }
  txt <- grepl("^[0-9]{1,2}-(Sep|Mar|Dec)$", g)          # "7-Sep" stored as text
  if (any(txt)) {
    day <- sub("-.*", "", g[txt]); mon <- sub(".*-", "", g[txt])
    out[txt] <- paste0(DATE_PREFIX[mon], as.integer(day)); fix[txt] <- "excel_date"
  }
  bad <- grepl("^[0-9]+$", out)
  if (any(bad)) stop("unresolved numeric gene name(s): ", paste(out[bad], collapse = ", "))
  list(sym = out, fix = fix)
}

suppressPackageStartupMessages(library(org.Hs.eg.db))
official <- unique(AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL"))
alias_to_symbol <- function(s) {
  m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = s, keytype = "ALIAS",
                                              columns = "SYMBOL"))
  m <- m[!is.na(m$SYMBOL), , drop = FALSE]
  u <- tapply(m$SYMBOL, m$ALIAS, function(x) if (length(unique(x)) == 1) x[1] else NA)
  out <- setNames(rep(NA_character_, length(s)), s)
  out[names(u)] <- u
  out
}

# One row per source symbol: date fix -> panel match -> (if off panel) unique alias rescue.
# exact = TRUE reproduces a published set's plain `intersect(genes, panel)` (no date fix,
# no alias rescue, no packed-probe expansion); what the full rule WOULD add is recorded
# in `fix` as "not_used:<reason>" with on_panel = FALSE.
clean_set <- function(key, genes, direction = NA_character_, exact = FALSE) {
  genes <- unique(as.character(genes[!is.na(genes) & nzchar(genes)]))
  if (length(direction) == 1) direction <- rep(direction, length(genes))
  if (exact) {
    full <- clean_set(key, genes, direction, exact = FALSE)
    plain <- full$symbol_source %in% panel
    full$probe[!plain] <- NA_character_
    full$probe[plain]  <- full$symbol_source[plain]
    would <- !plain & full$on_panel
    full$fix[would] <- paste0("not_used:", ifelse(full$fix[would] == "none", "packed_probe",
                                                  full$fix[would]))
    full$on_panel <- plain
    return(full)
  }
  df <- fix_date(genes)
  sym <- df$sym; fix <- df$fix
  off <- !(sym %in% panel_syms)
  if (any(off)) {
    al <- alias_to_symbol(sym[off])
    resc <- !is.na(al) & al != sym[off] & al %in% panel_syms
    idx <- which(off)[resc]
    sym[idx] <- al[resc]
    fix[idx] <- ifelse(fix[idx] == "none", "alias", paste0(fix[idx], "+alias"))
  }
  probe <- vapply(sym, sym2probe, character(1))
  data.frame(set = key, symbol_source = genes, symbol = sym, fix = fix,
             probe = unname(probe), on_panel = !is.na(probe),
             source_direction = direction, stringsAsFactors = FALSE)
}

## ---- S6: Otero L2/3 AT8-up (verbatim, = the 3B set) ----------------------------
otero <- readRDS(OTERO_RDS)
stopifnot("otero_L23_up missing" = "otero_L23_up" %in% names(otero))
l23_up <- otero[["otero_L23_up"]]
stopifnot("otero_L23_up should be 623 genes (Fig. 3B)" = length(l23_up) == 623)

## ---- S2: five-subtype shared NFT core ------------------------------------------
read_s2_deg <- function(ct) {
  sh <- sprintf("DGE - %s", ct)
  x <- suppressMessages(read_excel(S2_XLSX, sheet = sh, skip = 5, col_names = TRUE,
                                   col_types = "text"))
  x <- as.data.frame(x)[, 1:6]
  names(x) <- c("gene", "p_val", "avg_logFC", "pct1", "pct2", "p_val_adj")
  x <- x[!is.na(x$gene), , drop = FALSE]
  x$avg_logFC <- as.numeric(x$avg_logFC); x$p_val_adj <- as.numeric(x$p_val_adj)
  stopifnot("S2 DEG block: header misread" = all(is.finite(x$avg_logFC)),
            "S2 DEG block should hold padj < 0.05 genes only" = all(x$p_val_adj < 0.05))
  x$gene <- fix_date(x$gene)$sym          # so a date-mangled gene matches across sheets
  x$subtype <- ct
  x
}
s2 <- do.call(rbind, lapply(S2_CORE_SUBTYPES, read_s2_deg))
cat("S2 DEG rows per subtype:\n"); print(table(s2$subtype))
sgn <- tapply(sign(s2$avg_logFC), s2$gene, function(v) if (length(v) == 5 && length(unique(v)) == 1) v[1] else 0)
core_up   <- sort(names(sgn)[sgn ==  1])
core_down <- sort(names(sgn)[sgn == -1])
cat(sprintf("S2 five-subtype core: up %d, down %d\n", length(core_up), length(core_down)))
stopifnot("S2 core up should be 102"  = length(core_up)   == 102,
          "S2 core down should be 22" = length(core_down) == 22)

## ---- S5: SynGO process in the presynapse ---------------------------------------
s5 <- as.data.frame(suppressMessages(read_excel(S5_XLSX, sheet = "Sheet5", skip = 1)))
i5 <- which(s5[[1]] == "SYNGO:presynprocess")
stopifnot("SYNGO:presynprocess not found once in S5 Sheet5" = length(i5) == 1)
presyn <- strsplit(s5[i5, 10], ";", fixed = TRUE)[[1]]
stopifnot("S5 presynprocess should be 63 genes" = length(presyn) == 63)
# per-gene direction from Sheet1 (logFC in the five clusters; NA where not DE)
s5r <- as.data.frame(suppressMessages(read_excel(S5_XLSX, sheet = "Sheet1", skip = 3,
                                                 col_names = FALSE, col_types = "text")))
lfc <- sapply(c(2, 4, 6, 8, 10), function(j) as.numeric(s5r[[j]]))
rownames(lfc) <- s5r[[1]]
s5_dir <- vapply(presyn, function(g) {
  if (!g %in% rownames(lfc)) return(NA_character_)
  v <- lfc[g, ]; v <- v[!is.na(v)]
  if (all(v > 0)) "up" else if (all(v < 0)) "down" else "mixed"
}, character(1))
cat("S5 presynprocess direction (Otero logFC, 5 clusters):\n"); print(table(s5_dir, useNA = "ifany"))
stopifnot("S5: every presynprocess gene should be in Sheet1" = !anyNA(s5_dir))
presyn_up <- presyn[s5_dir == "up"]

## ---- Reactome Apoptosis --------------------------------------------------------
gmt <- strsplit(readLines(GMT, warn = FALSE), "\t", fixed = TRUE)
names(gmt) <- vapply(gmt, `[`, "", 1)
stopifnot("Reactome 'Apoptosis' not found once" = sum(names(gmt) == "Apoptosis") == 1)
apop <- gmt[["Apoptosis"]][-(1:2)]
apop <- apop[nzchar(apop)]
stopifnot("Reactome Apoptosis should be 168 genes" = length(apop) == 168)

## ---- assemble ------------------------------------------------------------------
mem <- rbind(
  clean_set("otero_l23_up",       l23_up,    "up", exact = TRUE),   # = the 3B set, 319
  clean_set("otero_s2_core_up",   core_up,   "up"),
  clean_set("otero_s2_core_down", core_down, "down"),
  clean_set("otero_s5_presyn",    presyn,    unname(s5_dir)),
  clean_set("otero_s5_presyn_up", presyn_up, "up"),
  clean_set("reactome_apoptosis", apop),
  clean_set("vatpase_foiani",     VATPASE)
)
summ <- do.call(rbind, lapply(split(mem, factor(mem$set, unique(mem$set))), function(d)
  data.frame(set = d$set[1], n_source = nrow(d),
             n_on_panel_probes = length(unique(d$probe[d$on_panel])),
             n_date_fixed = sum(grepl("excel_date", d$fix) & d$on_panel),
             n_alias_rescued = sum(grepl("alias", d$fix) & d$on_panel),
             off_panel = paste(d$symbol[!d$on_panel], collapse = ","),
             stringsAsFactors = FALSE)))
print(summ[, 1:5], row.names = FALSE)

# Gates: expected set sizes.
g <- setNames(summ$n_on_panel_probes, summ$set)
cat("\nV-ATPase off panel:", summ$off_panel[summ$set == "vatpase_foiani"], "\n")
stopifnot("V-ATPase: 9 of 11 expected on panel" = g[["vatpase_foiani"]] == 9,
          "otero_l23_up must reproduce Fig. 3B's 319 panel genes" = g[["otero_l23_up"]] == 319,
          "S2 core down must reach MIN_GENES_PRESENT = 10 (alias rescue failed?)" =
            g[["otero_s2_core_down"]] >= 10)

# <key>.txt = the PANEL PROBES to score (all symbol decisions are made here, so the
# consumer only has to check they are rownames of the SCT assay).
for (k in unique(mem$set))
  writeLines(sort(unique(mem$probe[mem$set == k & mem$on_panel])),
             file.path(OUT_DIR, paste0(k, ".txt")))
write.table(mem,  file.path(OUT_DIR, "fig2_contrast_membership.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(summ, file.path(OUT_DIR, "fig2_contrast_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

rescued <- mem[mem$fix != "none", c("set", "symbol_source", "symbol", "fix", "on_panel")]
prov <- c(
  "Fig. 2 contrast gene sets - provenance",
  sprintf("Generated : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "Script    : R/prepare_fig2_contrast_genesets.R",
  sprintf("Panel     : %s rownames (%d probes, %d packed)", PANEL_SCE, length(panel),
          sum(lengths(probe_members) > 1)),
  sprintf("org.Hs.eg.db %s (alias rescue, unique mappings only)",
          as.character(utils::packageVersion("org.Hs.eg.db"))),
  "",
  "Sources:",
  sprintf("  otero_l23_up        %s [[otero_L23_up]] (verbatim; = Fig. 3B)", OTERO_RDS),
  sprintf("  otero_s2_core_*     %s DGE - %s, DEG block, same sign in all five",
          S2_XLSX, paste(S2_CORE_SUBTYPES, collapse = "/")),
  sprintf("  otero_s5_presyn     %s Sheet5 SYNGO:presynprocess (63)", S5_XLSX),
  sprintf("                      direction from Sheet1: %s",
          paste(sprintf("%s %d", names(table(s5_dir)), table(s5_dir)), collapse = ", ")),
  "                      -> the 63-gene set is NOT uniformly up (SynGO input was",
  "                         dysregulated genes); plot it without an expected-direction",
  "                         arrow; otero_s5_presyn_up is the consistently-up sensitivity.",
  sprintf("  reactome_apoptosis  %s 'Apoptosis'", GMT),
  "  vatpase_foiani      Foiani et al. 2026, 11 subunits",
  "", "Per-set counts (probe level):",
  capture.output(print(summ, row.names = FALSE)),
  "", "Symbols changed by date fix / alias rescue:",
  if (nrow(rescued)) capture.output(print(rescued, row.names = FALSE)) else "  none",
  "", "Consumed by R/plot_cbln2_traj_otero_by_bin.R (otero_l23_up)")
writeLines(prov, file.path(OUT_DIR, "fig2_contrast_provenance.txt"))
cat("\nWrote", OUT_DIR, "\n")
