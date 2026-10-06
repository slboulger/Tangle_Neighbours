#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# deg_pb_category_phf1.r
#
# Produces: Table S4
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# Perform differential expression analysis 
# Nurun Fancy <n.fancy@imperial.ac.uk>
# Use run_pb_deg_category.sh

##  ............................................................................
##  Load packages                                                           ####
library(tidyverse)
library(argparse)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(caret)
library(Matrix)

source("<PROJECT_ROOT>/phf1_v2/R/generate_pseudobulk_limma_trend_deg_gazestani.r")

##  ............................................................................
##  Parse command-line arguments                                            ####

# create parser object
parser <- ArgumentParser()

# specify options
required <- parser$add_argument_group("Required", "required arguments")
optional <- parser$add_argument_group("Optional", "required arguments")

required$add_argument(
  "--sce",
  help = "path to sce.qs object",
  metavar = "<PROJECT_ROOT>/phf1_v2/celltype_sce/celltype_sce/Exc-IT-L2-3_sce.qs",
  required = TRUE
)

required$add_argument(
  "--dependent_var",
  help = "dependent variable",
  metavar = "PHF1",
  required = TRUE
)

required$add_argument(
  "--ref_class",
  help = "reference class within dependent variable",
  metavar = "NULL",
  required = TRUE
)

required$add_argument(
  "--confounding_vars",
  help = "confounding variables",
  metavar = "Age,Sex,PostMortemInterval",
  required = TRUE
)

required$add_argument(
  "--output_dir",
  help = "output dir path",
  metavar = "<PROJECT_ROOT>/phf1_v2/deg/de_cat",
  required = TRUE
)


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
args <- parser$parse_args()

# set working directory

celltype <- gsub("_sce.qs", "", basename(args$sce))

args$confounding_vars <- strsplit(args$confounding_vars, ",")[[1]]


outdir <- sprintf("%s/de_%s/%s", args$output_dir, args$dependent_var, celltype)
cli::cli_text(c(
  c("Creating ", cli::col_green(c(outdir, " \n")))
))
dir.create(outdir, recursive = TRUE)


#### preparing sce ####

cli::cli_text(c(
  c("Reading {.strong {celltype}} from ", cli::col_green(c(args$sce, " \n")))
))

sce <- qs::qread(args$sce)
rownames(sce) <- rowData(sce)$gene
sce$sample_id <- as.factor(sce$sample_id)
sce$sample_id <- droplevels(sce$sample_id)



cli::cli_text(c(
  c("Generating pseudobulk"))
)


sce$celltype_sample_phf1 <- paste(sce$celltype, sce$sample_id, colData(sce)[[args$dependent_var]], sep = "_")

mean_pseudocell_size=100000
min_pseudocell_size=3

pseudocell_data_all = 
  .sconline.PseudobulkGeneration(
    argList = NULL,
    n_clusters = NULL,
    parsing.col.names = c("celltype_sample_phf1"), # must make this column in sce
    use.sconline.cluster4parsing = FALSE,
    cluster_obj = NULL,
    pseudocell.size = mean_pseudocell_size, # NULL or mean_pseudocell_size
    inputExpData = sce,
    min_size_limit = min_pseudocell_size, # NULL or min_pseudocell_size
    inputPhenoData = as.data.frame(colData(sce)),
    inputEmbedding = NULL,
    tol_level = 0.9,
    use.sconline.embeddings = FALSE,
    nPCs = 20,
    ncores = 1,
    rand_pseudobulk_mod = TRUE,
    organism = "Human"
  )


#### Check 
Y <- counts(pseudocell_data_all)          # genes x pseudocells
ph <- colData(pseudocell_data_all)

# Indices for each group (of the dependent variable, not hardcoded PHF1)
ix_true  <- which(ph[[args$dependent_var]] == "TRUE")
ix_false <- which(ph[[args$dependent_var]] == "FALSE")

N_true  <- length(ix_true)
N_false <- length(ix_false)

cat("N pseudocells: TRUE =", N_true, " FALSE =", N_false, "\n")

# Prevalence thresholds: at least 50% of pseudocells/samples in that group (but at least 2)
T2_true  <- max(2, ceiling(0.5 * N_true))
T2_false <- max(2, ceiling(0.5 * N_false))

# Expression thresholds (total CPM within group)
# Start with something modest; you can tune this
T1_true  <- 100 * N_true    # total CPM in PHF1 TRUE
T1_false <- 100 * N_false    # total CPM in PHF1 FALSE

cat("True group thresholds:  T1_true =", T1_true,  " T2_true  =", T2_true,  "\n")
cat("False group thresholds: T1_false =", T1_false," T2_false =", T2_false, "\n")



# 1) Prevalence: # of pseudocells with >0 counts per gene in each group
tmpCount2_true  <- apply(Y[, ix_true,  drop = FALSE], 1, function(x) sum(x > 0))
tmpCount2_false <- apply(Y[, ix_false, drop = FALSE], 1, function(x) sum(x > 0))

# 2) Total CPM per gene within each group
cpm_true  <- edgeR::cpm(as.matrix(Y[, ix_true,  drop = FALSE]))
cpm_false <- edgeR::cpm(as.matrix(Y[, ix_false, drop = FALSE]))

tmpCount_true  <- rowSums(cpm_true)
tmpCount_false <- rowSums(cpm_false)


keep_true  <- (tmpCount_true  > T1_true)  & (tmpCount2_true  >= T2_true)
keep_false <- (tmpCount_false > T1_false) & (tmpCount2_false >= T2_false)

keep <- keep_true | keep_false

bkg_genes <- rownames(pseudocell_data_all)[keep]

cat("Total genes:", nrow(pseudocell_data_all), "\n")
cat("Genes kept:", length(bkg_genes), "\n")
cat("Proportion kept:", length(bkg_genes) / nrow(pseudocell_data_all), "\n")

####


# tmpCount2=apply(counts(pseudocell_data_all),1,function(x) sum(x>0))
# tmpCount=rowSums(edgeR::cpm(as.matrix(counts(pseudocell_data_all))))
# keep=tmpCount>max(0.1*ncol(pseudocell_data_all),min(20,ncol(pseudocell_data_all)/3))
# keep=keep & tmpCount2>max(0.05*ncol(pseudocell_data_all),min(20,ncol(pseudocell_data_all)/2))
# bkg_genes <- rownames(pseudocell_data_all)[keep]



pseudocell_data <- pseudocell_data_all
  
  pseudocell_data$pseudocell_size_scale=as.numeric(scale(pseudocell_data$pseudocell_size))
  pseudocell_data$nUMI=colSums(counts(pseudocell_data))
  pseudocell_data$nUMI_scaled=log2(pseudocell_data$nUMI)
  pseudocell_data$nGene=apply(counts(pseudocell_data),2,function(x) sum(x>0))
  pseudocell_data$nGene_scaled=log2(pseudocell_data$nGene)
  pseudocell_data$QC_MT.pct <- log2(pseudocell_data$mean_pc_mito+1)
  nGeneUMI=prcomp(as.matrix(colData(pseudocell_data)[,c("nGene","nUMI")]),scale. = T)
  pseudocell_data$nGeneUMI=nGeneUMI$x[,1]
  pseudocell_data$Age=as.numeric(scale(pseudocell_data$Age))
  pseudocell_data$PMI=as.numeric(scale(pseudocell_data$PMI))
  
  
  
  default_vars <- c("nUMI_scaled") # Set 1p: pseudocell_size_scale and QC_MT.pct (mean percent.neg) are not used
  
  
  ## --- Sanity checks before limma ---
  
  Y   <- counts(pseudocell_data)
  cd  <- colData(pseudocell_data)
  
  cat("Before limma: dim(counts) =", paste(dim(Y), collapse = " x "), "\n")
  cat("Before limma: dim(colData) =", paste(dim(cd), collapse = " x "), "\n")
  
  if (ncol(Y) != nrow(cd)) {
    stop("Mismatch: ncol(counts(pseudocell_data)) != nrow(colData(pseudocell_data))")
  }

  model_vars   <- c(args$dependent_var, default_vars, args$confounding_vars)
  
  cat("Model covariates:", paste(model_vars, collapse = ", "), "\n")
  cat("Available colData columns:\n")
  print(colnames(cd))
  
  missing <- setdiff(model_vars, colnames(cd))
  if (length(missing) > 0) {
    stop("Missing model covariates in colData(pseudocell_data): ",
         paste(missing, collapse = ", "))
  }
  
  # Show first few rows of the phenotype used for the design
  cat("Head of model covariates in colData(pseudocell_data):\n")
  print(head(as.data.frame(cd[, model_vars, drop = FALSE])))
  ## --------------------------------
  
  
  # quantile_norm <- TRUE
  
    fit = .sconline.fitLimmaFn(
      inputExpData = pseudocell_data,
      covariates = c(args$dependent_var, 
                     default_vars, 
                     args$confounding_vars),
      randomEffect = "sample_id", 
      normalization = "CPM", # "TMM" or "CPM", etc.
      DEmethod = "Trend",
      prior.count = 1, 
      bkg_genes = bkg_genes
    )
    
    colnames(fit$fit)

    
    # Contrast = (positive) - (negative) for the chosen dependent variable.
    # Coefficient names are "<dependent_var>TRUE"/"<dependent_var>FALSE".
    .coef_true  <- paste0(args$dependent_var, "TRUE")
    .coef_false <- paste0(args$dependent_var, "FALSE")
    if (!all(c(.coef_true, .coef_false) %in% colnames(fit$model))) {
      cat(sprintf("Skipping %s: '%s' does not have both TRUE and FALSE pseudocells (design coefs: %s).\n",
                  celltype, args$dependent_var, paste(colnames(fit$model), collapse = ", ")))
      quit(save = "no", status = 0)
    }
    contrast1 = makeContrasts(contrasts = paste0(.coef_true, " - ", .coef_false),
                              levels = fit$model)

    
    
    fit2=contrasts.fit(fit$fit,contrasts=contrast1)
    fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)
    
    colnames(fit2)
    
    
    for(i in colnames(fit2)){
    jk_res=topTable(fit2,
                    number=dim(fit2)[1],
                    adjust.method = "BH", 
                    coef=i, 
                    confint = TRUE
    )
    
    jk_res <- jk_res %>%
      mutate(gene = rownames(.),
             contrast = i,
             model = paste0(c("~", args$dependent_var, 
                              paste0(default_vars, collapse = "+"), 
                              paste0(args$confounding_vars, collapse = "+")),
                            collapse = "+")
      ) %>%
      dplyr::rename(pval = P.Value,
                    padj = adj.P.Val,
      ) %>%
      dplyr::select(c(gene, logFC, pval, padj, AveExpr, CI.L, CI.R, t, B, contrast, model))
    
    
    
    file_name <- sprintf("%s/%s_%s.tsv", outdir, celltype, gsub(" ", "", i) %>% gsub("-", "Vs", .))
    write.table(
      x = jk_res, 
      file = file_name,
      col.names = TRUE, row.names = FALSE, sep = "\t")
    
}
    
    
    cat("duplicateCorrelation consensus r =", fit$dc$consensus.correlation, "\n")
    cat("Blocked analysis applied:", fit$blocked_analysis, "\n")
