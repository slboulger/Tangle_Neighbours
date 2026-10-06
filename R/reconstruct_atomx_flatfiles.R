# ---------------------------------------------------------------------------
# reconstruct_atomx_flatfiles.R
#
# Data release - rebuilds the per-slide AtoMx exprMat and fov_positions flat files
# deposited on GEO from the vendor tx_file, verified against the AtoMx Seurat export.
# No figure panel.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# Reconstruct the two AtoMx flat files that were never exported for GEO:
#   <slide>_exprMat_file.csv.gz        cell x target counts
#   <slide>_fov_positions_file.csv.gz  per-FOV global position
#
# Both are rebuilt from the vendor files we DO hold (tx_file, metadata_file, polygons, the AtoMx
# Seurat export, the instrument RunSummary) and are verified against them before anything is
# written. Nothing is approximated: the script stop()s on any disagreement.
#
# exprMat  = per-cell transcript counts from tx_file (cell_ID 0 = unassigned, excluded as in
#            AtoMx SIP >= v1.3.2). Must equal the AtoMx Seurat counts (RNA + negprobes +
#            falsecode) exactly, and the metadata nCount_/nFeature_ columns exactly.
#            Column order follows AtoMx v1.3.2: fov, cell_ID, genes (alphabetical),
#            Negative1..N, SystemControl1..N.
# fov_pos  = AtoMx v1.3.2 columns FOV, x_global_px, y_global_px, x_global_mm, y_global_mm.
#            In these exports x_global = x_local + x_off and y_global = y_off - y_local (y is
#            flipped; x_off/y_off are the exact integer offsets of the cell/polygon frame; transcripts
#            carry an extra constant sub-pixel shift, reported in the log).
#            The file reports the FOV's MINIMUM global corner, i.e.
#              x_global_px = x_off,  y_global_px = y_off - FOV_PX
#            so a FOV spans [x_global_px, x_global_px + FOV_PX] x [y_global_px, y_global_px + FOV_PX]
#            and  x_global = x_global_px + x_local,  y_global = y_global_px + FOV_PX - y_local.
#            mm uses the px -> mm map AtoMx itself applied to cells (x_slide_mm / y_slide_mm in
#            the Seurat export), fitted here and required to be exact.
#
# Run (local mount or HPC interactive session):
#   COSMX_BASE=<project dir> Rscript R/reconstruct_atomx_flatfiles.R

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(Seurat)
})
BASE <- Sys.getenv("COSMX_BASE",
  "<PROJECT_ROOT>")
setwd(BASE)

FOV_PX <- 4256L
SLIDES <- data.table(
  folder = c("1b", "2a"),
  run    = c("20260114_151706_S1", "20260114_151706_S2"),
  seurat = c("seuratObject_IGFQ002102_matthews_13.1.2026_CosMx_RNA_1b.RDS",
             "seuratObject_IGFQ002102_matthews_13.1.2026_CosMx_RNA_2a.RDS")
)
FLAT <- "flatFiles/flatFiles"
OVERWRITE <- identical(Sys.getenv("OVERWRITE"), "1")

log_file <- "flatFiles/reconstruct_atomx_flatfiles.log"
sink(log_file, split = TRUE)
cat("reconstruct_atomx_flatfiles.R  ", format(Sys.time()), "\nBASE:", BASE, "\n\n")

check <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(ok)) "PASS" else "FAIL", msg))
  if (!isTRUE(ok)) { sink(); stop(msg, call. = FALSE) }
}
gz_read <- function(f, ...) fread(cmd = paste("gzip -dc", shQuote(f)), ...)
num_suffix <- function(x, prefix) x[order(as.integer(sub(prefix, "", x)))]

new_files <- character()

for (i in seq_len(nrow(SLIDES))) {
  sl <- SLIDES[i]
  tag <- paste0("IGFQ002102_matthews_1312026_CosMx_RNA_", sl$folder)
  dir <- file.path(FLAT, tag)
  f_tx   <- file.path(dir, paste0(tag, "_tx_file.csv.gz"))
  f_meta <- file.path(dir, paste0(tag, "_metadata_file.csv.gz"))
  f_poly <- file.path(dir, paste0(tag, "-polygons.csv.gz"))
  f_expr <- file.path(dir, paste0(tag, "_exprMat_file.csv.gz"))
  f_fov  <- file.path(dir, paste0(tag, "_fov_positions_file.csv.gz"))
  run_dir <- file.path("DecodedFiles/DecodedFiles", tag, sl$run, "RunSummary")
  f_stage <- list.files(run_dir, pattern = "_RNA_FOV_Locations\\.csv$", full.names = TRUE)

  cat("==================================================================\n")
  cat("Slide folder", sl$folder, "\n")
  if (!OVERWRITE && any(file.exists(c(f_expr, f_fov))))
    stop("Output exists; set OVERWRITE=1 to replace: ", f_expr)
  check(length(f_stage) == 1, paste("one RunSummary FOV_Locations file:", f_stage))

  # ---- inputs -------------------------------------------------------------------------------
  meta <- gz_read(f_meta, select = c("fov", "cell_ID", "cell", "slide_ID",
    "CenterX_local_px", "CenterY_local_px", "CenterX_global_px", "CenterY_global_px",
    "nCount_RNA", "nFeature_RNA", "nCount_negprobes", "nFeature_negprobes",
    "nCount_falsecode", "nFeature_falsecode"))
  setorder(meta, fov, cell_ID)
  slide_id <- unique(meta$slide_ID)
  cat("  metadata cells:", nrow(meta), " AtoMx slide_ID:", slide_id, "\n")
  check(length(slide_id) == 1, "single slide_ID in metadata")
  check(!anyDuplicated(meta$cell), "metadata cell keys unique")

  t0 <- Sys.time()
  tx <- gz_read(f_tx, select = c("fov", "cell_ID", "x_local_px", "y_local_px",
                                 "x_global_px", "y_global_px", "target"))
  cat(sprintf("  tx_file transcripts: %s (read %.1f min)\n", format(nrow(tx), big.mark = ","),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  seu <- readRDS(sl$seurat)
  genes <- rownames(seu[["RNA"]])
  negs  <- num_suffix(rownames(seu[["negprobes"]]), "Negative")
  sysc  <- num_suffix(rownames(seu[["falsecode"]]), "SystemControl")
  features <- c(sort(genes, method = "radix"), negs, sysc)
  cat("  features: genes", length(genes), "| Negative", length(negs),
      "| SystemControl", length(sysc), "\n")
  check(!anyDuplicated(features), "feature names unique across the three assays")
  check(all(unique(tx$target) %in% features), "every tx_file target is a known feature")

  # ---- FOV offsets: metadata, polygons and transcripts must agree -------------------------------
  # Cells and polygons use ONE integer offset per FOV (exact, zero spread). Transcript local and
  # global px are rounded per point from sub-pixel positions, so their (global - local) straddles
  # the integer offset plus a constant sub-pixel shift; its per-FOV mean must round to the
  # integer offset. The integer offset is what fov_positions reports.
  fov_offsets <- function(d, xg, xl, yg, yl) {
    d[, .(x_off = mean(get(xg) - get(xl)), y_off = mean(get(yg) + get(yl)),
          x_spread = diff(range(get(xg) - get(xl))), y_spread = diff(range(get(yg) + get(yl))),
          max_local = max(get(xl), get(yl)), n = .N), by = fov]
  }
  off_meta <- fov_offsets(meta, "CenterX_global_px", "CenterX_local_px", "CenterY_global_px", "CenterY_local_px")
  check(max(off_meta$x_spread, off_meta$y_spread) == 0, "metadata: exactly one integer offset per FOV")

  poly <- gz_read(f_poly, select = c("fov", "x_local_px", "y_local_px", "x_global_px", "y_global_px"))
  off_poly <- fov_offsets(poly, "x_global_px", "x_local_px", "y_global_px", "y_local_px")
  rm(poly)
  check(max(off_poly$x_spread, off_poly$y_spread) == 0, "polygons: exactly one integer offset per FOV")

  off_tx <- fov_offsets(tx, "x_global_px", "x_local_px", "y_global_px", "y_local_px")
  check(max(off_tx$x_spread, off_tx$y_spread) <= 1, "tx: (global - local) spans <= 1 px within every FOV")
  check(max(off_tx$max_local) <= FOV_PX, sprintf("tx local px <= %d", FOV_PX))

  fov <- merge(off_meta[, .(fov, x_off, y_off)],
               off_poly[, .(fov, x_off_p = x_off, y_off_p = y_off)], by = "fov", all = TRUE)
  fov <- merge(fov, off_tx[, .(fov, x_off_t = x_off, y_off_t = y_off, n_tx = n)], by = "fov", all = TRUE)
  check(all(!is.na(fov$x_off)) && all(!is.na(fov$x_off_p)) && all(!is.na(fov$x_off_t)),
        "metadata, polygons and tx cover the same FOVs")
  check(with(fov, all(x_off == x_off_p & y_off == y_off_p)), "polygon offsets identical to metadata offsets")
  check(with(fov, all(round(x_off_t) == x_off & round(y_off_t) == y_off)),
        "per-FOV transcript offset rounds to the metadata offset")
  fov[, `:=`(sx = x_off_t - x_off, sy = y_off_t - y_off)]
  cat(sprintf(paste0("  transcript sub-pixel shift vs cell frame: x %.3f [%.3f, %.3f] px, ",
                     "y %.3f [%.3f, %.3f] px (mean [min, max] over FOVs; min tx per FOV %d)\n"),
              mean(fov$sx), min(fov$sx), max(fov$sx), mean(fov$sy), min(fov$sy), max(fov$sy), min(fov$n_tx)))
  steps <- c(diff(sort(unique(fov$x_off))), diff(sort(unique(fov$y_off))))
  cat("  distinct FOV grid steps (px):", paste(sort(unique(steps)), collapse = ", "), "\n")

  # Instrument stage positions: the FOV list, and an independent linear check on the offsets.
  stage <- fread(f_stage)
  cat("  FOVs: tx", nrow(off_tx), "| metadata", nrow(off_meta), "| RunSummary", nrow(stage), "\n")
  check(setequal(stage$FOV, fov$fov), "FOV set identical to the instrument RunSummary")
  fov <- merge(fov, stage[, .(fov = FOV, X_mm, Y_mm)], by = "fov")
  fx <- lm(x_off ~ X_mm + Y_mm, fov); fy <- lm(y_off ~ X_mm + Y_mm, fov)
  cat(sprintf("  offsets vs stage mm: x_off R2 %.8f max|resid| %.2f px; y_off R2 %.8f max|resid| %.2f px\n",
              summary(fx)$r.squared, max(abs(resid(fx))), summary(fy)$r.squared, max(abs(resid(fy)))))
  print(round(rbind(x_off = coef(fx), y_off = coef(fy)), 3))
  check(max(abs(resid(fx)), abs(resid(fy))) < 5, "offsets linear in instrument stage position (< 5 px)")

  # ---- px -> mm map AtoMx used for cells --------------------------------------------------------
  sm <- data.table(cell = colnames(seu), x_mm = seu$x_slide_mm, y_mm = seu$y_slide_mm)
  sm <- merge(sm, meta[, .(cell, CenterX_global_px, CenterY_global_px)], by = "cell")
  check(nrow(sm) == nrow(meta), "Seurat x/y_slide_mm available for every metadata cell")
  mx <- lm(x_mm ~ CenterX_global_px, sm); my <- lm(y_mm ~ CenterY_global_px, sm)
  cat(sprintf("  px->mm: x = %.10g + %.10g*px (max|resid| %.2e mm); y = %.10g + %.10g*px (max|resid| %.2e mm)\n",
              coef(mx)[1], coef(mx)[2], max(abs(resid(mx))), coef(my)[1], coef(my)[2], max(abs(resid(my)))))
  check(max(abs(resid(mx)) / coef(mx)[2], abs(resid(my)) / abs(coef(my)[2])) < 1,
        "AtoMx px->mm map linear to within the 1 px rounding of CenterX/Y_global_px")

  fov_out <- fov[, .(FOV = fov, x_global_px = as.integer(x_off), y_global_px = as.integer(y_off - FOV_PX))]
  fov_out[, x_global_mm := round(predict(mx, data.frame(CenterX_global_px = x_global_px)), 6)]
  fov_out[, y_global_mm := round(predict(my, data.frame(CenterY_global_px = y_global_px)), 6)]
  setorder(fov_out, FOV)

  # ---- counts ----------------------------------------------------------------------------------
  n_unassigned <- tx[cell_ID == 0, .N]
  agg <- tx[cell_ID != 0, .N, by = .(fov, cell_ID, target)]
  rm(tx); gc(verbose = FALSE)
  meta[, row := .I]
  agg <- merge(agg, meta[, .(fov, cell_ID, row)], by = c("fov", "cell_ID"), all.x = TRUE)
  check(!anyNA(agg$row), "every cell-assigned transcript belongs to a metadata cell")
  agg[, col := match(target, features)]
  M <- sparseMatrix(i = agg$col, j = agg$row, x = as.numeric(agg$N),
                    dims = c(length(features), nrow(meta)), dimnames = list(features, meta$cell))
  rm(agg)
  cat(sprintf("  transcripts: assigned %s | unassigned (cell_ID 0, not in exprMat) %s\n",
              format(sum(M), big.mark = ","), format(n_unassigned, big.mark = ",")))

  # Cross-check 1: AtoMx Seurat export, exactly.
  S <- rbind(LayerData(seu, assay = "RNA", layer = "counts"),
             LayerData(seu, assay = "negprobes", layer = "counts"),
             LayerData(seu, assay = "falsecode", layer = "counts"))
  check(setequal(colnames(S), meta$cell), "Seurat cells == metadata cells")
  S <- S[features, meta$cell]
  D <- drop0(M - S)
  check(length(D@x) == 0, sprintf("tx-derived counts identical to AtoMx Seurat counts (%d differing entries)", length(D@x)))
  rm(S, D, seu); gc(verbose = FALSE)

  # Cross-check 2: metadata totals.
  is_g <- features %in% genes; is_n <- features %in% negs; is_s <- features %in% sysc
  check(all(colSums(M[is_g, ]) == meta$nCount_RNA),       "colSums(genes) == nCount_RNA")
  check(all(colSums(M[is_g, ] > 0) == meta$nFeature_RNA), "genes detected == nFeature_RNA")
  check(all(colSums(M[is_n, ]) == meta$nCount_negprobes), "colSums(Negative) == nCount_negprobes")
  check(all(colSums(M[is_s, ]) == meta$nCount_falsecode), "colSums(SystemControl) == nCount_falsecode")

  # ---- write exprMat (dense, chunked through a local temp file) ---------------------------------
  tmp <- tempfile(fileext = ".csv")
  chunks <- split(seq_len(ncol(M)), ceiling(seq_len(ncol(M)) / 10000))
  for (k in seq_along(chunks)) {
    j <- chunks[[k]]
    dt <- data.table(fov = meta$fov[j], cell_ID = meta$cell_ID[j])
    dense <- as.matrix(t(M[, j, drop = FALSE])); storage.mode(dense) <- "integer"
    dt <- cbind(dt, as.data.table(dense))
    fwrite(dt, tmp, append = k > 1, col.names = k == 1)
  }
  stopifnot(system(paste("gzip -c", shQuote(tmp), ">", shQuote(f_expr))) == 0)
  unlink(tmp)
  fwrite(fov_out, f_fov, compress = "gzip")

  # ---- verify what landed on disk ---------------------------------------------------------------
  back <- gz_read(f_expr)
  check(identical(names(back), c("fov", "cell_ID", features)), "exprMat header = fov, cell_ID, features")
  check(nrow(back) == nrow(meta) && all(back$fov == meta$fov & back$cell_ID == meta$cell_ID),
        "exprMat rows = metadata cells, same order")
  bm <- as.matrix(back[, -(1:2)])
  check(all(colSums(bm) == rowSums(M)) && all(rowSums(bm) == colSums(M)),
        "written exprMat row/column sums equal the verified matrix")
  rm(back, bm); gc(verbose = FALSE)
  fback <- gz_read(f_fov)
  check(nrow(fback) == nrow(stage) && isTRUE(all.equal(fback$x_global_px, fov_out$x_global_px)),
        "fov_positions written and read back")
  cat("  wrote", f_expr, sprintf("(%.0f MB)", file.size(f_expr) / 2^20), "\n")
  cat("  wrote", f_fov, "\n")
  new_files <- c(new_files, f_expr, f_fov)
  rm(M, meta, fov, fov_out); gc(verbose = FALSE)
}

# ---- checksums for the reconstructed files (vendor manifest left untouched) ----------------------
md5 <- data.table(md5sum = unname(tools::md5sum(new_files)),
                  file = sub(paste0("^", FLAT, "/"), "", new_files))
fwrite(md5, "flatFiles/md5sum/md5sum_reconstructed_flatFiles.csv")
cat("\nmd5 -> flatFiles/md5sum/md5sum_reconstructed_flatFiles.csv\n"); print(md5)
cat("\n"); print(sessionInfo())
sink()
