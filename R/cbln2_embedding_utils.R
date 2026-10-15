# ---------------------------------------------------------------------------
# cbln2_embedding_utils.R
#
# Figure panels: 2B-2D, 2G, 3C - shared utility (embedding definition and cache)
# Sourced by the trajectory scripts; no panel of its own.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# cbln2_embedding_utils.R
#
# The single definition of the CBLN2-only embedding, shared by
#   R/plot_cbln2_umap_phf1_p95.R       (builds the cache; UMAP coloured by PHF1 p95, PHF1+ ringed)
#   R/cbln2_trajectory_optimise.R      (trajectory settings search)
#   R/cbln2_slingshot_trajectory.R     (slingshot pseudotime on the same UMAP)
#   R/cbln2_trajectory_selection_null.R
# so the trajectory is drawn on exactly the same UMAP.
#
# Recipe: gene filter (>= 0.5% and >= 3 cells), LogNormalize + vst HVGs, SCTransform
# glmGamPoi regressing nCount_RNA + percent.neg + percent.false, PCA (30), PCs by the
# project's min(co1, co2) elbow rule unless dims are given, harmony on sample_id, UMAP
# n.epochs 200, seed 42. Keep the operation order as it is: a rebuild must reproduce the
# same RNG stream.
#
# Cache: subcluster_phf1/cbln2_umap_<norm>_dims-<dims>.qs, a list with
#   coords   cell_id, UMAP_1, UMAP_2
#   pca      cells x 30 PC scores
#   harmony  cells x length(use_dims)
#   stdev    PC standard deviations, for the elbow plot
#   prep     provenance (norm, genes, HVGs, vars regressed, elbow numbers, dims, ...)
# A coords-only cache still serves the UMAP plot. A caller that needs pca/harmony/
# stdev (need_full = TRUE) rebuilds it; the previous cache is kept as *.prev.qs and
# the UMAP concordance between the two builds is reported.

CBLN2_CT            <- "Exc-IT-L2-3-CBLN2-HOPX"
CBLN2_SCE_PATH      <- file.path("celltype_sce", paste0(CBLN2_CT, "_sce.qs"))

## Other celltypes go through the same code via --celltype. ct_tag() names every
## per-celltype file and folder; CBLN2 is "cbln2". Unknown celltypes get a sanitised
## lower-case name.
CT_TAGS <- c("Exc-IT-L2-3-CBLN2-HOPX"    = "cbln2",
             "Exc-IT-L3-5-CHGA-IL1RAPL2" = "chga")
ct_tag <- function(ct) {
  if (ct %in% names(CT_TAGS)) return(unname(CT_TAGS[ct]))
  gsub("[^a-z0-9]+", "_", tolower(ct))
}
ct_sce_path <- function(ct) file.path("celltype_sce", paste0(ct, "_sce.qs"))
## Per-cell pseudotime written by R/cbln2_slingshot_trajectory.R for a celltype.
ct_percell_file <- function(ct) file.path("plots", paste0(ct_tag(ct), "_slingshot_trajectory"),
                                          paste0("source_data_", ct_tag(ct), "_traj_percell.tsv"))
EMB_CACHE_DIR       <- "subcluster_phf1"
MIN_CELLS_GENE_FRAC <- 0.005   # R/SCTransform.R:7
MIN_CELLS_GENE_ABS  <- 3
NFEATURES           <- 2000
NPCS                <- 30
UMAP_EPOCHS         <- 200     # as in the whole-cohort clustering
SCT_METHOD          <- "glmGamPoi"
EMB_SEED            <- 42

## PC selection -- the project's elbow rule, as used for neuron subclustering. Also
## returns the spectrum so the elbow can be plotted.
pick_pcs <- function(stdv) {
  percent_stdv <- (stdv / sum(stdv)) * 100
  cumulative   <- cumsum(percent_stdv)
  # 1st PC that cumulatively explains >90% of variance while itself explaining <5%
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1]
  # PC with the last drop in variance >0.1 (the knee)
  co2 <- sort(which((percent_stdv[1:(length(percent_stdv) - 1)] -
                       percent_stdv[2:length(percent_stdv)]) > 0.1),
              decreasing = TRUE)[1] + 1
  min_pc <- min(co1, co2, na.rm = TRUE)
  if (!is.finite(min_pc) || min_pc < 2) min_pc <- min(10, length(stdv))
  list(min_pc = min_pc, co1 = co1, co2 = co2,
       percent_stdv = percent_stdv, cumulative = cumulative)
}

parse_dims <- function(x) {
  if (identical(x, "elbow")) return(NULL)
  if (grepl("^[0-9]+$", x)) return(seq_len(as.integer(x)))
  if (grepl("^[0-9]+:[0-9]+$", x)) {
    ab <- as.integer(strsplit(x, ":")[[1]]); return(ab[1]:ab[2])
  }
  stop("--dims must be 'elbow', N or A:B; got: ", x)
}

emb_cache_file <- function(norm, dims, ct = CBLN2_CT) {
  file.path(EMB_CACHE_DIR, sprintf("%s_umap_%s_dims-%s.qs", ct_tag(ct), norm, gsub(":", "-", dims)))
}

## Build the embedding from a counts matrix + colData frame (rownames = cell ids).
build_cbln2_embedding <- function(counts, meta, dims_req, norm) {
  stopifnot(identical(rownames(meta), colnames(counts)))
  min_cells <- max(MIN_CELLS_GENE_ABS, ceiling(MIN_CELLS_GENE_FRAC * ncol(counts)))
  keep      <- Matrix::rowSums(counts > 0) >= min_cells

  meta <- meta[, setdiff(colnames(meta), c("nCount_RNA", "nFeature_RNA",
                                           "ident", "orig.ident")), drop = FALSE]
  obj  <- Seurat::CreateSeuratObject(counts = counts[keep, , drop = FALSE], meta.data = meta,
                                     min.cells = 0, min.features = 0)
  obj  <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4)
  obj  <- Seurat::FindVariableFeatures(obj, selection.method = "vst",
                                       nfeatures = min(NFEATURES, nrow(obj)))

  vtr <- intersect(c("nCount_RNA", "percent.neg", "percent.false"), colnames(obj@meta.data))
  bad <- vtr[vapply(vtr, function(v) anyNA(obj@meta.data[[v]]), logical(1))]
  if (length(bad)) {
    warning("Dropping vars.to.regress with NAs: ", paste(bad, collapse = ", "), call. = FALSE)
    vtr <- setdiff(vtr, bad)
  }

  if (identical(norm, "sct")) {
    obj <- Seurat::SCTransform(obj, assay = "RNA", new.assay.name = "SCT", method = SCT_METHOD,
                               vars.to.regress = if (length(vtr)) vtr else NULL,
                               return.only.var.genes = TRUE, seed.use = EMB_SEED)
    Seurat::DefaultAssay(obj) <- "SCT"
  } else {
    obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(obj),
                             vars.to.regress = if (length(vtr)) vtr else NULL)
  }

  npcs <- min(NPCS, ncol(obj) - 1, length(Seurat::VariableFeatures(obj)) - 1)
  obj  <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(obj), npcs = npcs,
                         seed.use = EMB_SEED)
  stdv <- obj[["pca"]]@stdev
  pcs  <- pick_pcs(stdv)
  use_dims <- if (is.null(dims_req)) 1:pcs$min_pc else dims_req
  use_dims <- use_dims[use_dims <= ncol(Seurat::Embeddings(obj, "pca"))]
  stopifnot(length(use_dims) >= 2)

  obj <- harmony::RunHarmony(obj, group.by.vars = "sample_id", reduction.use = "pca",
                             dims.use = use_dims, reduction.save = "harmony")
  obj <- Seurat::RunUMAP(obj, reduction = "harmony", dims = use_dims, n.epochs = UMAP_EPOCHS,
                         seed.use = EMB_SEED, reduction.name = "umap")

  emb <- as.data.frame(Seurat::Embeddings(obj, "umap"))
  names(emb) <- c("UMAP_1", "UMAP_2")
  list(coords  = tibble::rownames_to_column(emb, "cell_id"),
       pca     = Seurat::Embeddings(obj, "pca"),
       harmony = Seurat::Embeddings(obj, "harmony"),
       stdev   = stdv,
       prep = list(norm = norm, n_genes = sum(keep), min_cells_gene = min_cells,
                   n_hvg = length(Seurat::VariableFeatures(obj)), vars_regressed = vtr,
                   pcs_elbow = pcs[c("min_pc", "co1", "co2")], use_dims = use_dims,
                   harmony = "sample_id", umap_epochs = UMAP_EPOCHS, seed = EMB_SEED,
                   seurat_version = as.character(packageVersion("Seurat"))))
}

## Reuse the cache if it has what the caller needs, else build it (SCT needs
## glmGamPoi). `get_inputs` is a function returning list(counts, meta), so the counts
## matrix is only pulled when a build actually happens.
load_or_build_cbln2_embedding <- function(get_inputs, norm = "sct", dims = "elbow",
                                          need_full = FALSE, force_rebuild = FALSE, ct = CBLN2_CT) {
  dir.create(EMB_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
  cache_file <- emb_cache_file(norm, dims, ct)
  old <- if (file.exists(cache_file)) qs::qread(cache_file) else NULL
  has_full <- !is.null(old) && all(c("pca", "harmony", "stdev") %in% names(old))

  if (!force_rebuild && !is.null(old) && (!need_full || has_full)) {
    message("Reusing cached embedding: ", cache_file)
    old$prep$cache <- cache_file
    return(old)
  }
  if (!is.null(old) && need_full && !has_full)
    message("Cached embedding ", cache_file, " has coordinates only; rebuilding to add ",
            "PCA / harmony / PC spectrum (the previous cache is kept as *.prev.qs).")
  if (identical(norm, "sct") && !requireNamespace("glmGamPoi", quietly = TRUE))
    stop("glmGamPoi is not installed in this R, so the SCT embedding cannot be built here. ",
         "Run once in an R with glmGamPoi installed to write ", cache_file,
         "; after that every caller reuses it. The script will not substitute a ",
         "normalisation method for you.", call. = FALSE)

  inp <- get_inputs()
  set.seed(EMB_SEED)   # same RNG state as the top-level set.seed in plot_cbln2_umap_phf1_p95.R
  new <- build_cbln2_embedding(inp$counts, inp$meta, parse_dims(dims), norm)

  if (!is.null(old)) {
    prev_file <- sub("\\.qs$", ".prev.qs", cache_file)
    qs::qsave(old, prev_file)
    m <- match(new$coords$cell_id, old$coords$cell_id)
    new$prep$rebuild_vs_previous <- list(
      previous_file = prev_file,
      max_abs_diff_umap = max(abs(as.matrix(new$coords[, c("UMAP_1", "UMAP_2")]) -
                                  as.matrix(old$coords[m, c("UMAP_1", "UMAP_2")]))),
      cor_umap1 = cor(new$coords$UMAP_1, old$coords$UMAP_1[m]),
      cor_umap2 = cor(new$coords$UMAP_2, old$coords$UMAP_2[m]),
      same_use_dims = identical(new$prep$use_dims, old$prep$use_dims))
    message(sprintf("Rebuilt embedding vs previous cache: max |dUMAP| = %.3g, r = %.4f / %.4f%s",
                    new$prep$rebuild_vs_previous$max_abs_diff_umap,
                    new$prep$rebuild_vs_previous$cor_umap1,
                    new$prep$rebuild_vs_previous$cor_umap2,
                    if (new$prep$rebuild_vs_previous$max_abs_diff_umap > 1e-6)
                      " -- the UMAP CHANGED; re-run plot_cbln2_umap_phf1_p95.R so the figures agree" else ""))
  }
  qs::qsave(new, cache_file)
  message("Cached embedding: ", cache_file)
  new$prep$cache <- cache_file
  new
}
