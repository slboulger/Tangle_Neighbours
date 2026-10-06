# ---------------------------------------------------------------------------
# generate_pseudobulk_limma_trend_deg_gazestani.r
#
# Shared utility - sourced by the scripts above, no panel of its own
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
.sconline.PseudobulkGeneration=function(argList=NULL,n_clusters=NULL, parsing.col.names = c("anno_batch"), use.sconline.cluster4parsing=T,cluster_obj=NULL,
                                        pseudocell.size=40,inputExpData=NULL,min_size_limit=20,inputPhenoData=NULL,inputEmbedding=NULL,tol_level=0.9,use.sconline.embeddings=F,nPCs=NULL,ncores=5,calculate.outlier.score=F,rand_pseudobulk_mod=T,organism,cols_to_sum=NULL){
  
  
  #parsing.col.names: The columns in the pheonData that will be used to parse the expression data and generate the pseudocell/pseudobulk data
  #use.sconline.cluster4parsing: use the sconline cluster as factor for the parsing
  #inputPhenoData: in case we want to run the function outside sconline space
  #min_size_limit: minimum acceptable size (ie, #cells) for each pseudobulk
  #pseudocell.size: average pseudocell size.
  #if pseudocell.size=null, function turns to a tranditional pseudobulk method, ie, all cells at each parsing level are combined together
  #nPCs: the dimension of the embedding space for the construction of pseudobulk data
  #inputEmbedding: the embedding space to be used for the generation of the pseudobulk. only needed when pseudocell.size is not null
  #use.sconline.embeddings: use the embedding space used for sconline to generate the pseudobulk samples. Only needed when the pseudocell.size is not set to null
  
  #Function adds/modifies three annotation columns: pseuodcell_size, QC_Gene_total_count, QC_Gene_unique_count
  #QC_Gene_total_count: equivalant to nUMI for the pseudobulk samples
  #QC_Gene_unique_count: equivalant to nGene for the pseudobulk samples
  #use scale(tst$QC_Gene_total_count) and scale(tst$pseudocell_size) as additional covariates for the DE analysis
  
  
  if(is.null(inputPhenoData)&is.null(argList)){
    stop("Both argList and inputPhenoData cannot be null!")
  }
  
  if(!is.null(pseudocell.size)){
    if(pseudocell.size<2){
      stop("pseudocell.size cannot be less than 2!")
    }
  }
  
  if(is.null(inputExpData)&!rand_pseudobulk_mod){
    warning("It's advised to provide the inputEmbedding")
  }
  
  if(is.null(nPCs)&!is.null(pseudocell.size)&!rand_pseudobulk_mod){
    if(!is.null(argList)){
      nPCs=argList$nPCs
    } else if(!is.null(inputEmbedding)){
      warning(paste0("Setting the nPCs to ", nrow(inputEmbedding)," based on the inputEmbedding"))
    } else {
      stop("nPCs argument (number of PCs for the generation of embedding) need to be provided")
    }
    
    if(nPCs>pseudocell.size){
      warning(paste0("nPCs larger than pseuodcell.size is not advised. setting nPCs to ",pseudocell.size-5))
      nPCs=pseudocell.size-5
    }
  }
  
  
  
  if(is.null(argList)){
    use.sconline.cluster4parsing=F
    use.sconline.embeddings=F
  }
  
  
  if(is.null(inputPhenoData)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else {
    pd=inputPhenoData
  }
  
  
  if(use.sconline.cluster4parsing){
    parsing.col.names=c(parsing.col.names,"cluster_anno_res")
    prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
    prop_mat=Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
    
    if(!is.null(n_clusters)){
      cat(paste0("Analysis based on ",n_clusters," clusters"))
      if(is.null(cluster_obj)){
        stop("Cluster object needs to be provided!")
      }
      diff_clust=cluster_obj$cluster_object
      d_conMat=cutree(diff_clust,k=n_clusters)
      prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
      prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
      prop_merged=Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
      #prop anno
      
      
    } else {
      cat(paste0("Analysis at the pseudocell level"))
      prop_merged=prop_mat
    }
    
    colMax_vals_m=qlcMatrix::colMax(prop_merged)
    colMax_vals_m=prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
    prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
    prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
    prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
    prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
    prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
    pd$cluster_anno_res=paste0("C",as.character(prop_m_hardCluster$i))
  }
  
  if(is.null(inputExpData)){
    if(!file.exists(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))){
      stop("Expression data is missing!")
    }
    inputExpData=qread(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))
    inputExpData=.extraExport2ExpressionSetFn(counts=inputExpData@assays$RNA@counts,pd=as.data.frame(inputExpData@meta.data),fd=as.data.frame(inputExpData@assays$RNA@meta.features))
  } else if(tolower(class(inputExpData))=="seurat"){
    inputExpData=.sconline.convertSeuratToExpressionSet(object=inputExpData)
  } else if(class(inputExpData)!="SingleCellExperiment") {
    stop("Unrecognized inputExpression data!")
  }
  
  inputExpData=inputExpData[,colnames(inputExpData) %in% row.names(pd)]
  if(ncol(inputExpData)==0){
    stop("Expression data doesn't match with the phenoData!")
  }
  pd=pd[match(colnames(inputExpData), row.names(pd)),]
  
  for(icol in parsing.col.names){
    if(sum(colnames(colData(inputExpData))==icol)==0){
      if(sum(colnames(pd)==icol)==0){
        stop(paste0(icol," column was not identified!"))
      } else {
        if(class(pd[,icol])[1]==class(factor())){
          colData(inputExpData)[,icol]=as.character(pd[,icol])
        } else {
          colData(inputExpData)[,icol]=pd[,icol]
        }
        
      }
    }
  }
  
  if(length(unique(parsing.col.names))==1){
    inputExpData$lib_anno=colData(inputExpData)[,unique(parsing.col.names)]
  } else {
    inputExpData$lib_anno=apply(as.data.frame(colData(inputExpData)[,unique(parsing.col.names)]),1,function(x)paste(x,collapse="_"))
  }
  
  
  if(is.null(pseudocell.size)){
    #inputData=inputExpData;colName="lib_anno";mode="sum";cols_to_sum=cols_to_sum
    sl_data=.extra_sconline.PseudobulkFn(inputData=inputExpData,colName="lib_anno",mode="sum",min_size_limit=min_size_limit,cols_to_sum=cols_to_sum,ncores=ncores)
  } else {
    inputExpList=.mySplitObject_v2(inputExpData,colName="lib_anno",min_dataset_size=min_size_limit,ncores=ncores)
    #inputExpList2=.mySplitObject_v2(inputExpData,colName="library_anno2",min_dataset_size=min_size_limit,ncores=ncores)
    if(sum(unlist(lapply(inputExpList,class))=="SingleCellExperiment")!=length(inputExpList)){
      cat("Error: consider increasing RAM! re-trying with lower number of cores")
      inputExpList=.mySplitObject_v2(inputExpData,colName="lib_anno",min_dataset_size=min_size_limit,ncores=1)
      if(sum(unlist(lapply(inputExpList,class))=="SingleCellExperiment")!=length(inputExpList)){
        stop("Persistent error: increase RAM!")
      }
    }
    if(use.sconline.embeddings&!rand_pseudobulk_mod){
      load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
      inputEmbedding=harmony_embeddings[,1:nPCs,drop=F]
    }
    
    data_size=unlist(lapply(inputExpList,ncol))
    inputExpList=inputExpList[data_size>=min_size_limit]
    data_size=unlist(lapply(inputExpList,ncol))
    
    if(F){
      for(icheck in inputExpList){
        #inputExp=icheck;include.outlier.score=calculate.outlier.score
        tst=.extra_sconline.FixedSizeFn(inputExp=icheck,inputEmbedding=inputEmbedding,pseudocell_size=pseudocell.size,nPCs=nPCs,include.outlier.score=calculate.outlier.score)
      }
    }
    
    
    sl_data=suppressWarnings(parallel::mclapply(inputExpList,.extra_sconline.FixedSizeFn,inputEmbedding=inputEmbedding,pseudocell_size=pseudocell.size,nPCs=nPCs,include.outlier.score=calculate.outlier.score,rand_pseudobulk_mod=rand_pseudobulk_mod,mc.cores = ncores,cols_to_sum=cols_to_sum))
    if(sum(unlist(lapply(sl_data,class))=="SingleCellExperiment")!=length(inputExpList)){
      cat("Error: consider increasing RAM! re-trying with lower number of cores")
      sl_data=suppressWarnings(parallel::mclapply(inputExpList,.extra_sconline.FixedSizeFn,inputEmbedding=inputEmbedding,pseudocell_size=pseudocell.size,nPCs=nPCs,include.outlier.score=calculate.outlier.score,cols_to_sum=cols_to_sum,mc.cores = 1))
      if(sum(unlist(lapply(sl_data,class))=="SingleCellExperiment")!=length(inputExpList)){
        stop("Persistent error: increase RAM!")
      }
    }
    sl_data_size=lapply(sl_data,function(x) max(x$pseudocell_size))
    sl_data=sl_data[sl_data_size>=min_size_limit]
    sl_data=lapply(sl_data,function(x) x[,which(x$pseudocell_size>=min_size_limit),drop=F])
    
    sl_data=.mycBindFn(sl_data)
    sl_data$QC_Gene_total_count=apply(counts(sl_data),2,sum)
    sl_data$QC_Gene_unique_count=apply(counts(sl_data),2,function(x) sum(x>0))
  }
  
  sl_data$QC_MT.pct=.extraMitoPctFn(inputData = sl_data,organism = organism)
  
  
  return(sl_data)
}

.extra_sconline.PseudobulkFn=function(inputData,colName,mode="sum",min_size_limit=20,ncores=1,cols_to_sum=NULL){
  
  #colName: the column in the meta/pheno data that specifies the pseudobulk level.
  #colName: usually a column that is a combination of subjectId/libraryId + cellType
  #min_size_limit: the minimum acceptable size of the pseudobulk data.
  #min_size_limit: pseudobulks with cells less than this threshold are excluded from the analysis 
  
  #pdlist=split(as.data.frame(colData(inputData)),colData(inputData)[,colName])
  #res=matrix(0,nrow=nrow(inputData),ncol=length(pdlist))
  #colnames(res)=names(pdlist)
  
  require(Matrix)
  
  design_mat=as.matrix(.myOneHotFn(colData(inputData)[,colName]))
  if(!is.null(min_size_limit)){
    design_mat=design_mat[,colSums(design_mat)>=min_size_limit,drop=F]
  }
  
  design_mat=t(design_mat)
  design_mat=as(design_mat,"dgCMatrix")
  
  if(mode=="sum"){
    agg_mat=design_mat %*% t(counts(inputData))
    agg_mat=t(agg_mat)
    
  } else if (mode=="mean"){
    design_mat=Matrix::diag(x=1/rowSums(design_mat)) %*% design_mat
    agg_mat=design_mat %*% t(counts(inputData))
    agg_mat=t(agg_mat)
  }
  
  if(is.null(cols_to_sum)&sum(colnames(colData(inputData)) %in% cols_to_sum)==0){
    res_pd=as.data.frame(colData(inputData))
    res_pd=res_pd[!duplicated(res_pd[,colName]),]
    row.names(res_pd)=res_pd[,colName]
    res_pd=res_pd[match(colnames(agg_mat),row.names(res_pd)),]
    
  } else {
    if(length(setdiff(cols_to_sum,colnames(colData(inputData))))>0){
      warning("some of provided cols_to_sum cols were not identified in the dataset")
    }
    res_pd=as.data.frame(colData(inputData)[,!colnames(colData(inputData)) %in% cols_to_sum])
    res_pd=res_pd[!duplicated(res_pd[,colName]),]
    row.names(res_pd)=res_pd[,colName]
    res_pd=res_pd[match(colnames(agg_mat),row.names(res_pd)),]
    
    sums_res=design_mat %*% as.matrix(as.data.frame(colData(inputData)[,cols_to_sum]))
    if(any(row.names(sums_res)!=row.names(res_pd),na.rm = F)){
      print("Error in summing the cols!")
    }
    res_pd=cbind(res_pd,sums_res)
  }
  
  
  res=SingleCellExperiment(assays = list(counts = agg_mat),colData = res_pd,rowData=as.data.frame(rowData(inputData)))
  
  res$QC_Gene_total_count=apply(counts(res),2,sum)
  res$QC_Gene_unique_count=apply(counts(res),2,function(x) sum(x>0))
  lib_sizes=as.data.frame(table(colData(inputData)[,colName]))
  lib_sizes=lib_sizes[match(colData(res)[,colName],lib_sizes[,1]),]
  res$pseudocell_size=lib_sizes[,2]
  
  return(res)
  
}


.mySplitObject_v2=function(object,colName,min_dataset_size,ncores=5){
  if(class(object)=="Seurat"){
    res=Seurat::SplitObject(object,split.by = colName)
  } else{
    pd=colData(object)[,colName]
    res=as.data.frame(summary(counts(object)))
    res$anno=pd[res[,2]]
    res=split(res,pd[res[,2]])
    res_pd=split(as.data.frame(colData(object)),pd)
    res=parallel::mclapply(1:length(res),function(x){
      
      y=Matrix::sparseMatrix(i = res[[x]][,1],
                             j = as.numeric(as.factor(as.character(res[[x]][,2]))),
                             x = res[[x]][,3],dims = c(nrow(object),nrow(res_pd[[x]])))
      
      #res_pd[[x]]$pseudocell_size=nrow(res_pd[[x]])
      if(ncol(y)==1){#.5*min_dataset_size){
        #y=as(matrix(rowSums(y),ncol=1),"dgCMatrix")
        
        y=SingleCellExperiment(assays = list(counts = y),colData = res_pd[[x]][1,],rowData=as.data.frame(rowData(object)))
        
      } else {
        row.names(y)=row.names(object)
        colnames(y)=row.names(res_pd[[x]])
        y=SingleCellExperiment(assays = list(counts = y),colData = res_pd[[x]],rowData=as.data.frame(rowData(object)))
        
      }
      return(y)
      
    },mc.cores=ncores)
    
    if(F){
      for(i in unique(pd)){
        res=c(res,list(object[,which(pd==i)]))
        names(res)[length(res)]=i
      }
    }
    
  }
  
  size_dist=unlist(lapply(res,ncol))
  res=res[which(size_dist>=min_dataset_size)]
  
  return(res)
}


.extra_sconline.FixedSizeFn=function(inputExp,inputEmbedding=NULL,pseudocell_size=40,n.neighbors=20,n.trees=50,nPCs=30,k.param = 20,include.outlier.score=F,rand_pseudobulk_mod=T,cols_to_sum=NULL){
  #inputExp=datalist[[i]];inputEmbedding=NULL;pseudocell_size=50;nPCs=30
  #n.neighbors=20;n.trees=50;k.param = 20
  
  if(ncol(inputExp)==1){
    inputExp$pseudocell_size=1
    return(inputExp)
  }
  
  pseudocell_count=round(ncol(inputExp)/pseudocell_size)
  if(pseudocell_count<=1){
    expData=matrix(as.numeric(rowSums(counts(inputExp))),ncol=1)
    row.names(expData)=row.names(inputExp)
    colnames(expData)=colnames(inputExp)[1]
    
    pd=as.data.frame(colData(inputExp))
    pd=pd[1,]
    pd$pseudocell_outlier_score=NA
    pd$pseudocell_size=ncol(inputExp)
    
    if(!is.null(cols_to_sum)){
      for(icoltosum in cols_to_sum){
        pd[1,icoltosum]=sum(as.data.frame(colData(inputExp))[,icoltosum])
      }
    }
    
    pd$mean_pc_mito = mean(inputExp$pc_mito)
    
    
    fd=as.data.frame(rowData(inputExp))
    
    res=SingleCellExperiment(assays = list(counts = expData),colData = pd,rowData=fd)
    return(res)
  }
  
  if(rand_pseudobulk_mod){
    set.seed(42)
    prop_m_hardCluster=sample(colnames(inputExp))
    prop_m_hardCluster=split(prop_m_hardCluster,floor(ecdf(seq_along(prop_m_hardCluster))(seq_along(prop_m_hardCluster))*pseudocell_count-0.001))
    prop_m_hardCluster=lapply(prop_m_hardCluster,function(x){
      x=data.frame(pseudocell=x[1],cells=x,stringsAsFactors = F)
      return(x)
    })
    prop_m_hardCluster=do.call("rbind",prop_m_hardCluster)
    
    prop_m_hardCluster$name=prop_m_hardCluster[,1]
    prop_m_hardCluster$name=gsub(" ",".",as.character(prop_m_hardCluster$name))
    
    prop_m_hardCluster$name=gsub("[[:punct:]]+", ".", prop_m_hardCluster$name)
    
    m2=.myOneHotFn(inputVector=prop_m_hardCluster$name)
    row.names(m2)=prop_m_hardCluster[,2]
    m2_colname=prop_m_hardCluster[!duplicated(prop_m_hardCluster[,1]),]
    m2_colname=m2_colname[match(colnames(m2),m2_colname$name),]
    colnames(m2)=m2_colname[,1]
    prop_m_hardCluster=t(as.matrix(m2))
    rm(m2)
  } else {
    if(is.null(inputEmbedding)){
      tmpData=suppressWarnings(.extraExport2SeuratFn(inputExp))
      tmpData = NormalizeData(tmpData,verbose =F)
      tmpData = FindVariableFeatures(tmpData, selection.method = "vst", nfeatures = 2000,verbose =F)
      tmpData <- ScaleData(tmpData,verbose =F)
      tmpData <- RunPCA(tmpData,npcs =nPCs,verbose =F)
      inputEmbedding=tmpData@reductions$pca@cell.embeddings[,1:nPCs]
    } else {
      inputEmbedding=inputEmbedding[colnames(inputExp),1:nPCs]
    }
    
    harmony_embeddings=inputEmbedding
    
    
    pseudocell_names=.extra_sconline.FixedSizekmeansFn(harmony_embeddings=harmony_embeddings,nPCs = nPCs,pseudocell_count = pseudocell_count)#,kmeansMethod=kmeans_method)
    
    
    #pca_centroid=res_clust$centers
    
    pca_centroid=harmony_embeddings[pseudocell_names,,drop=F]
    row.names(pca_centroid)=paste0("ps_",1:nrow(pca_centroid))
    sl_pseudo=data.frame(cluster=paste0("ps_",1:nrow(pca_centroid)),pseudocell=pseudocell_names,stringsAsFactors = F)
    
    
    idx=Seurat:::AnnoyBuildIndex(data = harmony_embeddings, metric = "euclidean", n.trees = n.trees)
    nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = harmony_embeddings,k=k.param,include.distance = T,search.k = -1)
    
    affinities=Matrix::Diagonal(x=1/(nn.ranked.1$nn.dists[,2]+0.000001)) %*% nn.ranked.1$nn.dists
    affinities@x=-1*affinities@x^2
    affinities@x=exp(affinities@x)
    affinities[,1]=affinities[,2]
    j <- as.numeric(t(nn.ranked.1$nn.idx))
    i <- ((1:length(j)) - 1)%/%n.neighbors + 1
    x=as.numeric(t(affinities))
    adj = sparseMatrix(i = i, j = j, x = x, dims = c(nrow(inputEmbedding),nrow(inputEmbedding)))
    rownames(adj) <- row.names(inputEmbedding)
    colnames(adj)=c(row.names(inputEmbedding))
    adj= Matrix::Diagonal(x=1/rowSums(adj)) %*% adj
    
    
    adj_t=t(adj)
    adj_t=Matrix::Diagonal(x=1/rowSums(adj_t)) %*% adj_t
    adj=adj[sl_pseudo$pseudocell,]
    #row.names(adj)=sl_pseudo$cluster
    adj=adj %*% adj_t
    
    colMax_vals_m=qlcMatrix::colMax(adj)
    colMax_vals_m=adj %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
    itol=0.5
    tst_mat=Matrix::drop0(colMax_vals_m,tol=itol)
    tst_mat@x=rep(1,length(tst_mat@x))
    while(all(rowSums(tst_mat)>=pseudocell_size)&itol<0.95){
      itol=itol+0.05
      tst_mat=Matrix::drop0(colMax_vals_m,tol=itol)
      tst_mat@x=rep(1,length(tst_mat@x))
    }
    prop_m_hardCluster=Matrix::drop0(colMax_vals_m,tol=itol)
    prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
    
    o=table(prop_m_hardCluster$i)
    o=o[order(as.numeric(o),decreasing = F)]
    
    prop_m_hardCluster$groups=factor(as.character(prop_m_hardCluster$i),levels=names(o))
    prop_m_hardCluster=prop_m_hardCluster[order(prop_m_hardCluster$groups,decreasing = F),]
    prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
    
    prop_m_hardCluster=sparseMatrix(i = prop_m_hardCluster$i, j = prop_m_hardCluster$j, x = rep(1,nrow(prop_m_hardCluster)), dims = c(nrow(adj), ncol(adj)))
    row.names(prop_m_hardCluster)=row.names(adj)
    colnames(prop_m_hardCluster)=row.names(inputEmbedding)
    prop_m_hardCluster=prop_m_hardCluster[which(rowSums(prop_m_hardCluster)>0),]
    
  }
  expData=t(counts(inputExp))
  expData=t(prop_m_hardCluster %*% expData[colnames(prop_m_hardCluster),])
  
  #inputData
  #sl_pseudo=sl_pseudo[match(row.names(prop_m_hardCluster),sl_pseudo$pseudocell),]
  
  if(is.null(cols_to_sum)&sum(colnames(colData(inputExp)) %in% cols_to_sum)==0){
    pd=as.data.frame(colData(inputExp))
    pd=pd[match(row.names(prop_m_hardCluster),colnames(inputExp)),]
    pd$pseudocell_size=rowSums(prop_m_hardCluster)
    
    
  } else {
    if(length(setdiff(cols_to_sum,colnames(colData(inputExp))))>0){
      warning("some of provided cols_to_sum cols were not identified in the dataset")
    }
    
    pd=as.data.frame(colData(inputExp))[,!colnames(colData(inputExp)) %in% cols_to_sum]
    pd=pd[match(row.names(prop_m_hardCluster),colnames(inputExp)),]
    pd$pseudocell_size=rowSums(prop_m_hardCluster)
    
    
    sums_res=prop_m_hardCluster %*% as.matrix(as.data.frame(colData(inputExp)[,cols_to_sum]))
    if(any(row.names(sums_res)!=row.names(pd),na.rm = F)){
      print("Error in summing the cols!")
    }
    pd=cbind(pd,sums_res)
  }
  
  
  mean_pc_mito_df = .getMeanPcMitoPerPseudocell(prop_m_hardCluster, 
                                                as.data.frame(colData(inputExp)),
                                                mito_col = "percent.neg")
  pd$mean_pc_mito = mean_pc_mito_df$mean_pc_mito[match(rownames(pd), mean_pc_mito_df$pseudocell_barcode)]
  
  
  fd=as.data.frame(rowData(inputExp))
  
  pd$pseudocell_outlier_score=NA
  
  if(include.outlier.score&ncol(expData)>3){
    bkg_genes=counts(inputExp)
    bkg_genes=rowSums(bkg_genes>0)/max(ncol(bkg_genes),10)
    if(sum(bkg_genes>0.1)>100){
      bkg_genes=row.names(expData)[bkg_genes>0.1]
      
      logCPM=edgeR::cpm(expData[bkg_genes,],normalized.lib.sizes = F, log = TRUE, prior.count = 1)
      
      tocheck=require(WGCNA,quietly = T)
      if(!tocheck){
        stop("WGCNA package is missing!")
      }
      
      normadj <- (0.5+0.5*bicor(logCPM, use='pairwise.complete.obs'))^2
      netsummary <- fundamentalNetworkConcepts(normadj)
      outlier_score=scale(netsummary$Connectivity)
      pd$pseudocell_outlier_score=as.numeric(outlier_score)
    }
  }
  
  res=SingleCellExperiment(assays = list(counts = expData),colData = pd,rowData=fd)
  
  return(res)
}

.getMeanPcMitoPerPseudocell <- function(prop_m_hardCluster, coldata_df, mito_col = "percent.neg") {
  if (!is.matrix(prop_m_hardCluster)) stop("prop_m_hardCluster must be a matrix")
  if (!all(colnames(prop_m_hardCluster) %in% rownames(coldata_df))) {
    stop("Not all cell barcodes in prop_m_hardCluster are present in coldata_df rownames")
  }
  if (!mito_col %in% colnames(coldata_df)) {
    stop(paste("Column", mito_col, "not found in coldata_df"))
  }
  
  
  coldata_df <- coldata_df[colnames(prop_m_hardCluster), ]
  
  mito_vec <- coldata_df[colnames(prop_m_hardCluster), mito_col]
  mito_sum <- prop_m_hardCluster %*% mito_vec
  cell_counts <- rowSums(prop_m_hardCluster)
  mean_pc_mito <- mito_sum / cell_counts
  
  data.frame(
    pseudocell_barcode = row.names(prop_m_hardCluster),
    mean_pc_mito = as.vector(mean_pc_mito),
    stringsAsFactors = FALSE
  )
}

.myOneHotFn=function (inputVector) {
  if(T){
    #using caret package
    require(caret)
    if(sum(is.na(inputVector))>0){
      inputVector[is.na(inputVector)]="NA"
    }
    formula="~."
    data=data.frame(data=inputVector,stringsAsFactors = F)
    sep = "."
    levelsOnly = FALSE
    fullRank = FALSE
    
    formula <- as.formula(formula)
    if (!is.data.frame(data)) 
      data <- as.data.frame(data, stringsAsFactors = FALSE)
    vars <- all.vars(formula)
    if (any(vars == ".")) {
      vars <- vars[vars != "."]
      vars <- unique(c(vars, colnames(data)))
    }
    isFac <- unlist(lapply(data[, vars, drop = FALSE], is.factor))
    if (sum(isFac) > 0) {
      facVars <- vars[isFac]
      lvls <- lapply(data[, facVars, drop = FALSE], levels)
      if (levelsOnly) {
        tabs <- table(unlist(lvls))
        if (any(tabs > 1)) {
          stop(paste("You requested `levelsOnly = TRUE` but", 
                     "the following levels are not unique", "across predictors:", 
                     paste(names(tabs)[tabs > 1], collapse = ", ")))
        }
      }
    } else {
      facVars <- NULL
      lvls <- NULL
    }
    trms <- attr(model.frame(formula, data), "terms")
    out <- list(call = match.call(), form = formula, vars = vars, 
                facVars = facVars, lvls = lvls, sep = sep, terms = trms, 
                levelsOnly = levelsOnly, fullRank = fullRank)
    class(out) <- "dummyVars"
    trsf <- data.frame(predict(out, newdata = data))
    
    tmp=gsub("^data","",colnames(trsf))
    tmp=gsub("^\\.","",tmp)
    colnames(trsf)=tmp
  } else {
    if(length(names(inputVector))==0){
      trsf=reshape2::dcast(data=data.frame(var=1:length(inputVector),outcome=as.character(inputVector)),var ~ outcome, length)
      row.names(trsf)=trsf[,1]
      trsf=trsf[,-1]
    } else {
      inputVector=inputVector[order(inputVector)]
      inputVector=names(inputVector)
      trsf=reshape2::dcast(data=data.frame(var=1:length(inputVector),outcome=as.character(inputVector)),var ~ outcome, length)
      trsf=as.matrix(trsf)
      row.names(trsf)=trsf[,1]
      trsf=trsf[,-1]
    }
    
  }
  
  
  return(trsf)
}

.mycBindFn=function(inputList,batchNames=NULL,verbose=F){
  #cbinds multiple singleCellExpression datasets with differring number of rows.
  #inputList: the list of datasets to be merged
  #batchNames: the batch name to be assigned to each dataset. length(batchNames)==length(inputList)
  res_m=""
  if(!is.null(batchNames)){
    if(length(inputList)==1){
      res_m=inputList[[1]]
    } else if(length(inputList)==2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = batchNames[1:2])
    } else if(length(inputList)>2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = batchNames[1:2])
      for(i in 3:length(inputList)){
        res_m=.extracBindDetailFn(x1=res_m,x2=inputList[[i]],batchNames=c("",batchNames[i]))
      }
    }
  } else {
    if(length(inputList)==1){
      res_m=inputList[[1]]
    } else if(length(inputList)==2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = c("",""))
    } else if(length(inputList)>2){
      #x1=inputList[[1]];x2=inputList[[2]];batchNames = c("","")
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = c("",""))
      for(i in 3:length(inputList)){
        
        #x1=res_m;x2=inputList[[i]];batchNames=c("","")
        res_m=.extracBindDetailFn(x1=res_m,x2=inputList[[i]],batchNames=c("",""))
        if(verbose){
          print(paste("dataset:",i,"; nrow:",nrow(res_m),"; ncol",ncol(res_m)))
        }
      }
    }
  }
  
  print("batch information is in the anno_batch variable")
  return(res_m)
  
}

.extracBindDetailFn=function(x1,x2,batchNames){
  
  tmp1=setdiff(row.names(x1),row.names(x2))
  tmpAnno1=rowData(x1)[which(row.names(x1) %in% tmp1),]
  tmpMat1=sparseMatrix(i=NULL,j=NULL,dims = c(length(tmp1),ncol(x2)))
  row.names(tmpMat1)=tmp1
  
  x2_c=rbind(counts(x2),tmpMat1)
  x2_c_pd=as.data.frame(colData(x2))
  x2_c_fd=as.data.frame(plyr::rbind.fill(as.data.frame(rowData(x2)),as.data.frame(tmpAnno1)))
  x2_c=SingleCellExperiment(assays = list(counts = x2_c),colData = x2_c_pd,rowData=x2_c_fd)
  
  tmp2=setdiff(row.names(x2),row.names(x1))
  tmpAnno2=rowData(x2)[which(row.names(x2) %in% tmp2),]
  tmpMat2=sparseMatrix(i=NULL,j=NULL,dims = c(length(tmp2),ncol(x1)))
  row.names(tmpMat2)=tmp2
  
  x1_c=rbind(counts(x1),tmpMat2)
  x1_c_pd=as.data.frame(colData(x1))
  x1_c_fd=plyr::rbind.fill(as.data.frame(rowData(x1)),as.data.frame(tmpAnno2))
  x1_c=SingleCellExperiment(assays = list(counts = x1_c),colData = x1_c_pd,rowData=x1_c_fd)
  
  x2_c=x2_c[match(row.names(x1_c),row.names(x2_c)),]
  
  if(all(row.names(x2_c)==row.names(x1_c))){
    
    if(batchNames[1]!=""){
      row.names(colData(x1_c))=paste0(batchNames[1],"_",colnames(x1_c))
      colnames(x1_c)=paste0(batchNames[1],"_",colnames(x1_c))
      x1_c$anno_batch=batchNames[1]
    }
    
    if(batchNames[2]!=""){
      colnames(x2_c)=paste0(batchNames[2],"_",colnames(x2_c))
      row.names(colData(x2_c))=paste0(batchNames[2],"_",colnames(x2_c))
      x2_c$anno_batch=batchNames[2]
    }
    
    #tst1=counts(x1_c)
    #tst2=counts(x2_c)
    #tst1=summary(tst1)
    #tst2=summary(tst2)
    #tst2$j=tst2$j+max(tst1$j)
    #tst=rbind(tst1,tst2)
    #tst=sparseMatrix(i = tst,j = tstt,x = tsttt)
    
    x_m_exp=cbind(counts(x1_c),counts(x2_c))
    pd_m_exp=plyr::rbind.fill(as.data.frame(colData(x1_c)),as.data.frame(colData(x2_c)))
    fd=as.data.frame(rowData(x1_c))
    x_m=SingleCellExperiment(assays=list(counts=x_m_exp),colData=pd_m_exp,rowData=fd)
  } else {
    stop("Error in the merging!")
  }
  return(x_m)
}

.extraMitoPctFn=function(inputData,organism,inputMTgenes=NULL,inputGeneName=NULL,redownload_files=T,recalculate_nUMI=T){
  
  rwNames=tolower(row.names(inputData))
  x=.extraMitoGenes(organism = organism,redownload_files=redownload_files)
  if(sum(grepl("\\.",rwNames))>0){
    rwNames=strsplit(rwNames,"\\.")
    rwNames=unlist(lapply(rwNames,function(x)x[1]))
  }
  if(is.null(inputGeneName)){
    tmpCols=c("ensembl_gene_id","entrezgene_id",'gene_name')
    
    slCounts=0
    slCol=""
    
    for(i in tmpCols){
      tmp=sum(rwNames %in% tolower(x[,i]))
      if(tmp>slCounts){
        slCounts=tmp
        slCol=i
      }
    }
    if(slCol!=""){
      inputMTgenes=tolower(x[,slCol])
    } else {
      inputMTgenes=""
    }
    
  } else {
    if(inputGeneName=="ensembl_gene_id"){
      inputMTgenes=tolower(x$ensembl_gene_id)
    } else {
      if(inputGeneName=="entrezgene_id"){
        x=.extraMitochondrialGenes()
        inputMTgenes=tolower(x$entrezgene_id)
      }
    }
  }
  
  print(paste("Number of MT genes in the dataset:",length(which(rwNames %in% inputMTgenes)),"/",sum(x$gene_biotype=="protein_coding")))
  
  inputMTgenes=inputMTgenes[inputMTgenes!=""]
  if(length(inputMTgenes)>0){
    tmpColSums=c()
    if(sum(colnames(colData(inputData))=='QC_Gene_total_count')==0|recalculate_nUMI){
      if(ncol(inputData)>10000){
        for(i in seq(1,ncol(inputData),10000)){
          tmpColSums=c(tmpColSums,apply(counts(inputData)[,i:min(i+10000-1,ncol(inputData))],2,sum))
        }
      } else {
        tmpColSums=apply(counts(inputData),2,sum)
      }
    } else {
      tmpColSums=inputData$QC_Gene_total_count
    }
    
    
    
    res=colSums(x = counts(inputData)[which(rwNames %in% inputMTgenes), , drop = FALSE])/tmpColSums
    res=res*100
  } else {
    res=rep(NA,ncol(inputData))
  }
  
  return(res)
}

.extraMitoGenes=function(organism,redownload_files=T){
  #organism: Human, Mouse
  
  if(tolower(organism)=="human"){
    gns=.extraHumanGeneAnnoAdderFn()
  } else if (tolower(organism)=="mouse"){
    gns=.extraMouseGeneAnnoAdderFn()
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
  } else {
    stop("Wrong organism name!")
  }
  
  gns=gns[gns$seqnames=="MT",]
  gns2=unlist(gns$entrezid)
  gns$entrezgene_id=gns2
  gns$ensembl_gene_id=gns$gene_id
  return(gns)
}


.extraHumanGeneAnnoAdderFn=function(inputGeneNames=NULL,server=T,redownload_files=T){
  #require(EnsDb.Hsapiens.v75)
  require(EnsDb.Hsapiens.v86)
  
  if(!dir.exists("~/serverFiles")){
    dir.create("~/serverFiles",recursive = T)
  }
  
  gns <- as.data.frame(genes(EnsDb.Hsapiens.v86))
  gns$gene_short_name=gns$gene_name
  gns$symbol=toupper(gns$symbol)
  gns$ensembl_gene_id=row.names(gns)
  
  if(!is.null(inputGeneNames)){
    rwNames=toupper(inputGeneNames)
    psCols=c("gene_short_name","ensembl_gene_id")
    slCounts=0
    slCol=""
    if(sum(grepl("\\.",rwNames)&grepl("^ENS",rwNames))>0){
      rwNames=strsplit(rwNames,"\\.")
      rwNames=unlist(lapply(rwNames,function(x)x[1]))
    }
    if(server){
      #library(googleCloudStorageR)
      if(!file.exists("~/serverFiles/human_map_to_ensembl.rda")){
        system(paste0("gsutil -m cp gs://macosko_data/vgazesta/serverFiles/orthologsFeb3/human_map_to_ensembl.rda ~/serverFiles/human_map_to_ensembl.rda"))
        #gcs_get_object("vgazesta/serverFiles/orthologsFeb3/human_map_to_ensembl.rda", saveToDisk = "~/serverFiles/human_map_to_ensembl.rda",overwrite=T)
      }
      
      load("~/serverFiles/human_map_to_ensembl.rda")
    } else {
      load("~/Desktop/human_map_to_ensembl.rda")
    }
    
    map_to_ensmbl$source=toupper(map_to_ensmbl$source)
    
    if(!file.exists("~/serverFiles/human_mapping_hg19.rda")){
      system(paste0("gsutil -m cp gs://macosko_data/vgazesta/serverFiles/orthologsFeb3/human_mapping_hg19.rda ~/serverFiles/human_mapping_hg19.rda"))
      #gcs_get_object("vgazesta/serverFiles/orthologsFeb3/human_mapping_hg19.rda", saveToDisk = "~/serverFiles/human_mapping_hg19.rda",overwrite=T)
    }
    
    load("~/serverFiles/human_mapping_hg19.rda")
    human_hg19$source=toupper(human_hg19$source)
    
    if(sum(toupper(rwNames) %in% human_hg19$source) > sum(toupper(rwNames) %in% map_to_ensmbl$source)){
      map_to_ensmbl=merge(human_hg19,data.frame(source=toupper(rwNames),stringsAsFactors = F),by="source",all.y=T)
    } else {
      map_to_ensmbl=merge(map_to_ensmbl,data.frame(source=toupper(rwNames),stringsAsFactors = F),by="source",all.y=T)
    }
    
    gns=merge(gns,map_to_ensmbl,by.x="ensembl_gene_id",by.y="target",all.y=T)
    gns=gns[match(rwNames,gns$source),]
    row.names(gns)=inputGeneNames
    gns$gene_id=inputGeneNames
    gns=gns[,-which(colnames(gns) %in% c("source","target"))]
  }
  
  return(gns)
}

.sconline.fitLimmaFn=function(inputExpData,covariates,randomEffect, DEmethod = "Trend", normalization = "CPM",
                              quantile.norm = F, bkg_genes = NULL,
                              VST_fitType = "parametric", prior.count = 1,
                              include.malat1.as.covariate = F,
                              dc.object = NULL, dream_ncores = 4) {
  
  normalization = match.arg(normalization, c("CPM", "TMM", "VST", "rmTop50", "none"))
  DEmethod      = match.arg(DEmethod,      c("Trend", "Voom", "VoomSampleWeights", "Dream"))
  
  # checks on inputExpData omitted (keep as is)
  
  if (sum(colnames(colData(inputExpData)) %in% covariates) < length(covariates)) {
    stop(paste0("Covariates ",
                paste(setdiff(covariates, colnames(colData(inputExpData))), collapse = ", "),
                " were not identified in the inputExpData!"))
  }
  
  # Turn colData into a data.frame for safe model.matrix usage
  pd <- as.data.frame(colData(inputExpData))
  
  if (DEmethod != "Dream") {
    # Build formula string
    model_formula <- "~0"
    active_covariates <- c()
    
    for (icov in covariates) {
      if (length(unique(pd[[icov]])) > 1) {
        if (is.factor(pd[[icov]])) {
          pd[[icov]] <- as.character(pd[[icov]])
        }
        model_formula <- paste0(model_formula, "+", icov)
        active_covariates <- c(active_covariates, icov)
      } else {
        warning(paste0("Excluding ", icov, " covariate as it has only one level!"))
      }
    }
    
    # optional: drop samples with NA in covariates actually used
    if (length(active_covariates) > 0) {
      cc <- stats::complete.cases(pd[, active_covariates, drop = FALSE])
    } else {
      cc <- rep(TRUE, nrow(pd))
    }
    
    if (!all(cc)) {
      warning("Dropping ", sum(!cc),
              " samples with NA in model covariates: ",
              paste(active_covariates, collapse = ", "))
    }
    
    pd_use  <- pd[cc, , drop = FALSE]
    sl_use  <- inputExpData[, cc]
    
    model_matrix <- model.matrix(as.formula(model_formula), data = pd_use)
    
  } else {
    # Dream branch – can leave as is for now
    model_formula <- "~"
    for (icov in covariates) {
      if (length(unique(pd[[icov]])) > 1) {
        if (is.factor(pd[[icov]])) {
          pd[[icov]] <- as.character(pd[[icov]])
        }
        if (model_formula == "~") {
          model_formula <- paste0(model_formula, icov)
        } else {
          model_formula <- paste0(model_formula, "+", icov)
        }
      } else {
        warning(paste0("Excluding ", icov, " covariate as it has only one level!"))
      }
    }
    model_matrix <- as.formula(paste0(model_formula, " + (1|", randomEffect, ")"))
    pd_use <- pd
    sl_use <- inputExpData
  }
  
  # Optional rmTop50 block – keep as is, but apply to sl_use if you use it
  if (normalization == "rmTop50") {
    if (sum(colnames(rowData(sl_use)) == "QC_top50_expressed") > 0) {
      sl_use <- sl_use[rowData(sl_use)$QC_top50_expressed == "No", ]
    } else {
      warning("Column QC_top50_expressed was not identified in the gene attribute dataframe. Skipping removal of top 50 expressed genes!")
    }
  }
  
  # Now call methods using sl_use and model_matrix
  switch(DEmethod,
         Trend =.extra_sconline.Fit_LimmaTrendFn(
           sl_data       = sl_use,
           model         = model_matrix,
           random_effect = randomEffect,
           quantile.norm = quantile.norm,
           TMMnorm       = (normalization == "TMM"),
           VSTnorm       = (normalization == "VST"),
           prior.count   = prior.count,
           bkg_genes     = bkg_genes,
           no_normalization = (normalization == "none"),
           dc.object     = dc.object,
           include.malat1.as.covariate = include.malat1.as.covariate
         ),
         Dream =.extra_sconline.Fit_LimmaDreamFn(
           sl_data       = sl_use,
           pd            = pd_use,
           model         = model_matrix,
           random_effect = randomEffect,
           quantile.norm = quantile.norm,
           TMMnorm       = (normalization == "TMM"),
           VSTnorm       = (normalization == "VST"),
           prior.count   = prior.count,
           bkg_genes     = bkg_genes,
           no_normalization = (normalization == "none"),
           dc.object     = dc.object,
           include.malat1.as.covariate = include.malat1.as.covariate,
           ncores        = dream_ncores
         ),
         Voom =.extra_sconline.Fit_LimmaVoomFn(
           sl_data       = sl_use,
           model         = model_matrix,
           random_effect = randomEffect,
           quantile.norm = quantile.norm,
           sample.weights = FALSE,
           TMMnorm       = (normalization == "TMM"),
           bkg_genes     = bkg_genes,
           dc.object     = dc.object
         ),
         VoomSampleWeights =.extra_sconline.Fit_LimmaVoomFn(
           sl_data       = sl_use,
           model         = model_matrix,
           random_effect = randomEffect,
           quantile.norm = quantile.norm,
           sample.weights = TRUE,
           TMMnorm       = (normalization == "TMM"),
           bkg_genes     = bkg_genes,
           dc.object     = dc.object
         )
  )
}


.extra_sconline.Fit_LimmaTrendFn=function(sl_data,model,random_effect=NULL,TMMnorm = FALSE, VSTnorm = FALSE,
                                          prior.count = 1, quantile.norm = FALSE,
                                          bkg_genes = NULL, no_normalization = FALSE,
                                          dc.object = NULL,
                                          include.malat1.as.covariate = FALSE) {
  require(edgeR)
  require(limma)
  library(Matrix)
  
  ## --- Build logCPM ------------------------------------------------------
  if (VSTnorm) {
    logCPM <-.myRNAseqNormVSTfn(inputCountData = sl_data, fitType = "parametric")
    logCPM <- logCPM$originalData
    logCPM <- counts(logCPM)
    if (!is.null(bkg_genes)) {
      keep <- row.names(sl_data) %in% bkg_genes
    } else {
      keep <- rowSums(logCPM > 3) > (0.05 * ncol(logCPM))
    }
    logCPM <- logCPM[keep, , drop = FALSE]
    
  } else if (no_normalization) {
    logCPM <- counts(sl_data)
    if (quantile.norm) {
      logCPM <- limma::normalizeQuantiles(logCPM)
    }
    
  } else {
    if (!is.null(bkg_genes)) {
      keep <- row.names(sl_data) %in% bkg_genes
    } else {
      tmpCount2 <- apply(counts(sl_data), 1, function(x) sum(x > 0))
      tmpCount  <- rowSums(edgeR::cpm(as.matrix(counts(sl_data))))
      keep      <- tmpCount > max(0.01 * ncol(sl_data), min(15, ncol(sl_data) / 3))
      keep      <- keep & tmpCount2 > max(0.01 * ncol(sl_data), min(10, ncol(sl_data) / 3))
    }
    
    cat("Number of expressed genes:", sum(keep), "\n")
    tmpexp <- counts(sl_data)[keep, , drop = FALSE]
    dge    <- DGEList(tmpexp)
    
    if (TMMnorm) {
      dge <- calcNormFactors(dge)
    }
    
    logCPM <- new("EList")
    logCPM$E <- edgeR::cpm(dge, normalized.lib.sizes = TMMnorm,
                           log = TRUE, prior.count = prior.count)
    if (quantile.norm) {
      logCPM$E <- limma::normalizeQuantiles(logCPM$E)
    }
  }
  ## ----------------------------------------------------------------------
  
  ## --- Dimension check: model vs logCPM ---------------------------------
  if (inherits(logCPM, "EList")) {
    n_arrays <- ncol(logCPM$E)
  } else {
    n_arrays <- ncol(logCPM)
  }
  
  if (nrow(model) != n_arrays) {
    stop("In.extra_sconline.Fit_LimmaTrendFn: nrow(model) (",
         nrow(model), ") != ncol(logCPM) (", n_arrays, ")")
  }
  ## ----------------------------------------------------------------------
  
  ## --- Optional MALAT1 covariate ----------------------------------------
  if (include.malat1.as.covariate) {
    if (inherits(logCPM, "EList")) {
      rn <- rownames(logCPM$E)
    } else {
      rn <- rownames(logCPM)
    }
    malat_idx <- grepl("malat1", tolower(rn))
    malat_gene <- rn[malat_idx]
    if (length(malat_gene) > 0) {
      if (inherits(logCPM, "EList")) {
        model <- cbind(model,
                       malat = as.numeric(logCPM$E[malat_gene[1], ]))
      } else {
        model <- cbind(model,
                       malat = as.numeric(logCPM[malat_gene[1], ]))
      }
    }
  }
  ## ----------------------------------------------------------------------
  
  ## --- duplicateCorrelation & lmFit -------------------------------------
  dc <- NULL
  if (!is.null(random_effect)) {
    if (is.null(dc.object)) {
      dc <-.extra_sconline.duplicateCorrelation(
        logCPM,
        design = model,
        block  = colData(sl_data)[, random_effect]
      )
    } else {
      dc <- dc.object
    }
  }
  
  blocked_analysis <- FALSE
  if (!is.null(dc)) {
    if (!is.nan(dc$consensus.correlation)) {
      if (abs(dc$consensus.correlation) < 0.9) {
        fit <- lmFit(
          logCPM,
          model,
          block = as.character(colData(sl_data)[, random_effect]),
          correlation = dc$consensus.correlation
        )
        blocked_analysis <- TRUE
      } else {
        fit <- lmFit(logCPM, model)
      }
    } else {
      fit <- lmFit(logCPM, model)
    }
  } else {
    fit <- lmFit(logCPM, model)
  }
  ## ----------------------------------------------------------------------
  
  list(fit = fit, dc = dc, model = model,
       normData = logCPM, blocked_analysis = blocked_analysis)
}

.extra_sconline.Fit_LimmaDreamFn=function(sl_data,pd,model,random_effect=NULL,TMMnorm=F,VSTnorm=F,prior.count=1,quantile.norm=F,bkg_genes=NULL,no_normalization=F,dc.object=NULL,include.malat1.as.covariate=F,ncores=4){
  require(edgeR)
  require(limma)
  require(variancePartition)
  
  param = SnowParam(ncores, "SOCK", progressbar=TRUE)
  
  # estimate weights using linear mixed model of dream
  
  sl_data=sl_data[rowSds(as.matrix(counts(sl_data))) > 0, ]
  
  dge <- DGEList(counts=counts(sl_data))
  if(is.null(bkg_genes)){
    tmpCount=apply(counts(sl_data),1,function(x) sum(x>0))
    keep=tmpCount>10
    
  } else {
    keep=row.names(sl_data) %in% bkg_genes
  }
  
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  
  if(TMMnorm){
    dge <- calcNormFactors(dge)
  }
  
  
  
  # estimate weights using linear mixed model of dream
  #model=as.formula("~status + anno_sex +nUMI_scaled  + pseudocell_size_scale +(1 | anno_batch)")
  
  if(quantile.norm){
    vobjDream = voomWithDreamWeights( dge, model,as.data.frame(pd), BPPARAM=param, normalize.method="quantile" )
  } else {
    vobjDream = voomWithDreamWeights( dge, model,as.data.frame(pd), BPPARAM=param )
  }
  
  
  
  
  fit = dream( vobjDream, model, as.data.frame(pd) )
  
  
  return(list(fit=fit,dc=NULL,model=model,normData=NULL,blocked_analysis=NULL))
}

.extra_sconline.Fit_LimmaVoomFn=function(sl_data,model,random_effect=NULL,quantile.norm=F,sample.weights=F,TMMnorm=F,bkg_genes=NULL,dc.object=NULL){
  require(edgeR)
  require(limma)
  
  
  dge <- DGEList(counts=counts(sl_data))
  if(is.null(bkg_genes)){
    tmpCount=apply(counts(sl_data),1,function(x) sum(x>0))
    keep=tmpCount>10
    
  } else {
    keep=row.names(sl_data) %in% bkg_genes
  }
  
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  
  if(TMMnorm){
    dge <- calcNormFactors(dge)
  }
  
  if(sample.weights){
    if(quantile.norm){
      logCPM <- voomWithQualityWeights(dge, model, normalize.method="quantile",plot = T,save.plot=T)
    } else {
      logCPM <- voomWithQualityWeights(dge, model,plot = T,save.plot=T)
    }
  } else {
    if(quantile.norm){
      logCPM <- voom(dge, model, normalize.method="quantile",plot = T,save.plot=T)
    } else {
      logCPM <- voom(dge, model,plot = T,save.plot=T)
    }
  }
  
  dc=NULL
  if(!is.null(random_effect)){
    if(is.null(dc.object)){
      dc <- .extra_sconline.duplicateCorrelation(logCPM,design=model, block=colData(sl_data)[,random_effect])
    } else {
      dc=dc.object
    }
  }
  
  
  blocked_analysis=F
  if(!is.null(dc)){
    if(!is.nan(dc$consensus.correlation)){
      if(abs(dc$consensus.correlation)<0.9){
        fit <- lmFit(logCPM, model,block = colData(sl_data)[,random_effect], correlation=dc$consensus.correlation)
        blocked_analysis=T
      } else {
        fit <- lmFit(logCPM, model)
      }
    } else {
      fit <- lmFit(logCPM, model)
    }
  } else {
    fit <- lmFit(logCPM, model)
  }
  
  
  #logCPM=preprocessCore::normalize.quantiles(logCPM)
  
  
  return(list(fit=fit,dc=dc,model=model,blocked_analysis=blocked_analysis,normData=logCPM))
}



.extra_sconline.duplicateCorrelation=function (object, design = NULL, ndups = 2, spacing = 1, block = NULL, 
                                               trim = 0.15, weights = NULL) {
  require(limma)
  y <- limma:::getEAWP(object)
  M <- y$exprs
  ngenes <- nrow(M)
  narrays <- ncol(M)

  if (is.null(design)) 
    design <- y$design
  if (is.null(design)) 
    design <- matrix(1, ncol(y$exprs), 1)
  else {
    design <- as.matrix(design)
    if (mode(design) != "numeric") 
      stop("design must be a numeric matrix")
  }
  if (nrow(design) != narrays) 
    stop("Number of rows of design matrix does not match number of arrays")
  ne <- limma:::nonEstimable(design)
  if (!is.null(ne)) 
    cat("Coefficients not estimable:", paste(ne, collapse = " "), 
        "\n")
  nbeta <- ncol(design)
  if (missing(ndups) && !is.null(y$printer$ndups)) 
    ndups <- y$printer$ndups
  if (missing(spacing) && !is.null(y$printer$spacing)) 
    spacing <- y$printer$spacing
  if (missing(weights) && !is.null(y$weights)) 
    weights <- y$weights
  if (!is.null(weights)) {
    weights <- asMatrixWeights(weights, dim(M))
    weights[weights <= 0] <- NA
    M[!is.finite(weights)] <- NA
  }
  if (is.null(block)) {
    if (ndups < 2) {
      warning("No duplicates: correlation between duplicates not estimable")
      return(list(cor = NA, cor.genes = rep(NA, nrow(M))))
    }
    if (is.character(spacing)) {
      if (spacing == "columns") 
        spacing <- 1
      if (spacing == "rows") 
        spacing <- object$printer$nspot.c
      if (spacing == "topbottom") 
        spacing <- nrow(M)/2
    }
    Array <- rep(1:narrays, rep(ndups, narrays))
  }
  else {
    ndups <- 1
    nspacing <- 1
    Array <- block
  }
  if (is.null(block)) {
    M <- limma:::unwrapdups(M, ndups = ndups, spacing = spacing)
    ngenes <- nrow(M)
    if (!is.null(weights)) 
      weights <- limma:::unwrapdups(weights, ndups = ndups, spacing = spacing)
    design <- design %x% rep(1, ndups)
  }
  if (!requireNamespace("statmod", quietly = TRUE)) 
    stop("statmod package required but is not installed")
  rho <- rep(NA, ngenes)
  nafun <- function(e) NA
  for (i in 1:ngenes) {
    y <- drop(M[i, ])
    o <- is.finite(y)
    A <- factor(Array[o])
    nobs <- sum(o)
    nblocks <- length(levels(A))
    if (nobs > (nbeta + 2) && nblocks > 1 && nblocks < nobs - 1) {
      y <- y[o]
      X <- design[o, , drop = FALSE]
      Z <- model.matrix(~0 + A)
      if (!is.null(weights)) {
        w <- drop(weights[i, ])[o]
        s <- tryCatch(statmod::mixedModel2Fit(y, X, Z, 
                                              w, only.varcomp = TRUE, maxit = 20)$varcomp, error = nafun)
      } else{
        s <- tryCatch(statmod::mixedModel2Fit(y, X, 
                                              Z, only.varcomp = TRUE, maxit = 20)$varcomp, 
                      error = nafun)
      } 
      if (!is.na(s[1])) 
        rho[i] <- s[2]/sum(s)
    }
  }
  arho <- atanh(pmax(-1, rho))
  mrho <- tanh(mean(arho, trim = trim, na.rm = TRUE))
  list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho,value.list=arho)
}

