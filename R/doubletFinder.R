doubletFinder <- function(object, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE, idents=NULL,
                             batch.size=Inf, get.neighbor.doublets=T) {
                              
  require(Seurat); require(fields); require(KernSmooth)

  ## Generate new list of doublet classificatons from existing pANN vector to save time
  if (reuse.pANN) {
    pANN.old <- object@meta.data[ , reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing=T)[1:nExp]] <- "Doublet"
    object@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    return(object)
  }

  # if (reuse.pANN == FALSE) {
  #   ## Make merged real-artifical data
  #   real.cells <- rownames(object@meta.data)
  #   data <- seu@assays$RNA$counts[, real.cells]
  #   n_real.cells <- length(real.cells)
  #   n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
  #   print(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
  #   real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  #   real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  #   doublets <- (data[, real.cells1] + data[, real.cells2])/2
  #   colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
  #   data_wdoublets <- cbind(data, doublets)

    # data <- object@assays$RNA@layers$counts ## M
    # colnames(data) <- colnames(object) ## M
    # rownames(data) <- rownames(object) ## M
    # data <- data[, real.cells] ## M

  ## Make merged real-artifical data
  real.cells <- rownames(object@meta.data)
  # data <- object@assays$RNA@layers$data[, real.cells]
    data <- object@assays$RNA@layers$counts ## M
    colnames(data) <- colnames(object) ## M
    rownames(data) <- rownames(object) ## M
    data <- data[, real.cells] ## M
  n_real.cells <- length(real.cells)
  n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
  print(paste("Patch: Creating",n_doublets,"artificial doublets...",sep=" "))
  if(!is.null(idents) & length(idents) != n_real.cells) {
    print(paste("Number of passed idents not as number of cells in object. Ignoring passed idents"))
    idents <- NULL
  }

  if(!is.null(idents)) {

    # Keep track of the types of the simulated doublets
    if(length(idents) != n_real.cells) {
      print("Number of passed idents not as number of cells in object. Ignoring passed idents")
      idents <- NULL
    }
    if(length(unique(idents)) == 1) {
      print("All cells share the same passed idents. Ignoring passed idents.")
      idents <- NULL
    }

    # TODO: revisit setting ident names if needed

    # if(!is.null(idents)) idents <- setNames(idents, real.cells)

    stopifnot(typeof(idents)=="character")
    stopifnot(length(idents)==length(Cells(object)))
    stopifnot(!any(is.na(idents)))
    idents <- factor(idents)
    names(idents) <- Cells(object)
    doublet_types <- idents[real.cells]
  }

  real.cells1 <- c(); real.cells2 <- c();
  while(length(real.cells1) < n_doublets) {
    rc1 <- sample(real.cells, n_doublets - length(real.cells1), replace = T)
    rc2 <- sample(real.cells, n_doublets - length(real.cells1), replace = T)

    if(is.null(idents)) {
      real.cells1 <- rc1; real.cells2 <- rc2;
    } else {
      real.cells1 <- c(real.cells1, rc1[idents[rc1] != idents[rc2]])
      real.cells2 <- c(real.cells2, rc2[idents[rc1] != idents[rc2]])
    }
  }
  rm(rc1, rc2)
  gc()
  Sys.sleep(10)

  doublets <- (data[, real.cells1] + data[, real.cells2])/2
  colnames(doublets) <- paste("X", 1:n_doublets, sep = "")

  ## Store important pre-processing information
  orig.commands <- object@commands

  ## Pre-process Seurat object
  print("Creating Seurat object...")
  seu_wdoublets <- CreateSeuratObject(counts = cbind(data, doublets))
  rm(data, doublets)

  if (!sct) {
    print("Normalizing Seurat object...")
    seu_wdoublets1 <- NormalizeData(seu_wdoublets,
                                   normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                   scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                   margin = orig.commands$NormalizeData.RNA@params$margin)
    rm(seu_wdoublets); 
    gc() # Free up memory
    Sys.sleep(10)

    print("Finding variable genes...")
    seu_wdoublets2 <- FindVariableFeatures(seu_wdoublets1,
                                          selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                                          loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                                          clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                                          mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                                          dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                                          num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                                          binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                                          nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                                          mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                                          dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
    rm(seu_wdoublets1); 
    gc() # Free up memory
    Sys.sleep(10)

    print("Scaling data...")
    seu_wdoublets3 <- ScaleData(seu_wdoublets2,
                               features = orig.commands$ScaleData.RNA$features,
                               model.use = orig.commands$ScaleData.RNA$model.use,
                               do.scale = orig.commands$ScaleData.RNA$do.scale,
                               do.center = orig.commands$ScaleData.RNA$do.center,
                               scale.max = orig.commands$ScaleData.RNA$scale.max,
                               block.size = orig.commands$ScaleData.RNA$block.size,
                               min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
    rm(seu_wdoublets2); 
    gc()
    Sys.sleep(10) # Free up memory

    print("Running PCA...#1")
    seu_wdoublets4 <- RunPCA(seu_wdoublets3,
                            features = orig.commands$ScaleData.RNA$features,
                            npcs = length(PCs),
                            rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
                            weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                            verbose=FALSE)
    # rm(seu_wdoublets3); 
    # gc() # Free up memory
    # Sys.sleep(10)

    pca.coord <- seu_wdoublets4@reductions$pca@cell.embeddings[ , PCs]
    cell.names <- rownames(seu_wdoublets4@meta.data)
    nCells <- length(cell.names)
    # rm(seu_wdoublets4); 
    gc() # Free up memory
    Sys.sleep(10)

  } else {
    require(sctransform)

    print("Running SCTransform...")

    ## turn off warings excessive warning logs
    options(warn=-1)
    seu_wdoublets1 <- SCTransform(seu_wdoublets)
    options(warn=0)
    
    rm(seu_wdoublets); 
    gc()
    Sys.sleep(10)

    print("Running PCA...#2")
    seu_wdoublets2 <- RunPCA(seu_wdoublets1, npcs = length(PCs))
    
    rm(seu_wdoublets1)
    gc()
    Sys.sleep(10)
 
    pca.coord <- seu_wdoublets2@reductions$pca@cell.embeddings[ , PCs]
    cell.names <- rownames(seu_wdoublets2@meta.data)
    nCells <- length(cell.names)
    
    rm(seu_wdoublets2); 
    gc()
    Sys.sleep(10)
  }

  # Compute pANN
  print("Computing pANN...")
  pANN <- matrix(-1L, nrow = n_real.cells, ncol = 1, dimnames = list(real.cells, c("pANN")))
  doublet.neighbors.list <- list()
  k <- round(nCells * pK)

  step <- ifelse(is.infinite(batch.size), n_real.cells, batch.size)
  for (i in seq(1, n_real.cells, by = step)) {
    max.i <- min(i + step - 1, n_real.cells)

    dists <- rdist(pca.coord[i:max.i,], pca.coord)
    for (j in seq(1, nrow(dists))) {
      neighbors <- order(dists[j, ])[2:(k + 1)]
      ann <- neighbors[neighbors > n_real.cells]
      pANN[[i+j-1, 1]] <- length(ann)/k
      if(get.neighbor.doublets) doublet.neighbors.list[[i+j-1]] <- (ann - n_real.cells)
      rm(neighbors,ann) gc()
    }

    rm(dists)
    gc()
    
    if(!is.infinite(batch.size))
      print(paste('Completed computing pANN for', round((max.i/n_real.cells)*100, 2) , '% of cells'))
  }

  print("Classifying doublets..")
  classifications <- rep("Singlet",n_real.cells)
  classifications[order(pANN, decreasing=TRUE)[1:nExp]] <- "Doublet"
  object@meta.data[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(object@meta.data),1]
  object@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications


  if(get.neighbor.doublets) {
    names(doublet.neighbors.list) <- real.cells
    Tool(object) <- list(doublet.parents=data.frame(parent1=real.cells1, parent2=real.cells2),
                         neighbor.doublets=doublet.neighbors.list)
  }
  object <- Seurat::LogSeuratCommand(object)
  return(object)
}

    ## Store important pre-processing information

    # orig.commands <- seu@commands

    # ## Pre-process Seurat object
    # if (sct == FALSE) {
    #   print("Creating Seurat object...")
    #   seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

    #   print("Normalizing Seurat object...")
    #   seu_wdoublets <- NormalizeData(seu_wdoublets,
    #                                  normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
    #                                  scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
    #                                  margin = orig.commands$NormalizeData.RNA@params$margin)

    #   print("Finding variable genes...")
    #   seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
    #                                         selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
    #                                         loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
    #                                         clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
    #                                         mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
    #                                         dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
    #                                         num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
    #                                         binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
    #                                         nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
    #                                         mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
    #                                         dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)

    #   print("Scaling data...")
    #   seu_wdoublets <- ScaleData(seu_wdoublets,
    #                              features = orig.commands$ScaleData.RNA$features,
    #                              model.use = orig.commands$ScaleData.RNA$model.use,
    #                              do.scale = orig.commands$ScaleData.RNA$do.scale,
    #                              do.center = orig.commands$ScaleData.RNA$do.center,
    #                              scale.max = orig.commands$ScaleData.RNA$scale.max,
    #                              block.size = orig.commands$ScaleData.RNA$block.size,
    #                              min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)

    #   print("Running PCA...")
    #   seu_wdoublets <- RunPCA(seu_wdoublets,
    #                           features = orig.commands$ScaleData.RNA$features,
    #                           npcs = length(PCs),
    #                           rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
    #                           weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
    #                           verbose=FALSE)
    #   pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
    #   cell.names <- rownames(seu_wdoublets@meta.data)
    #   nCells <- length(cell.names)
    #   rm(seu_wdoublets); gc() # Free up memory
    # }

    # if (sct == TRUE) {
    #   require(sctransform)
    #   print("Creating Seurat object...")
    #   seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

    #   print("Running SCTransform...")
    #   seu_wdoublets <- SCTransform(seu_wdoublets)

    #   print("Running PCA...")
    #   seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
    #   pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
    #   cell.names <- rownames(seu_wdoublets@meta.data)
    #   nCells <- length(cell.names)
    #   rm(seu_wdoublets); gc()
  

    ## Compute PC distance matrix
    # print("Calculating PC distance matrix...")
    # dist.mat <- fields::rdist(pca.coord)

    # ## Compute pANN
    # print("Computing pANN...")
    # pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    # if(!is.null(idents)){
    #   neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = length(levels(doublet_types1))))
    # }
    # rownames(pANN) <- real.cells
    # colnames(pANN) <- "pANN"
    # k <- round(nCells * pK)
    # for (i in 1:n_real.cells) {
    #   neighbors <- order(dist.mat[, i])
    #   neighbors <- neighbors[2:(k + 1)]
    #   pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
    #   if(!is.null(idents)){
    #     for(ct in unique(idents)){
    #       neighbors_that_are_doublets = neighbors[neighbors>n_real.cells]
    #       if(length(neighbors_that_are_doublets) > 0){
    #         neighbor_types[i,] <-
    #           table( doublet_types1[neighbors_that_are_doublets - n_real.cells] ) +
    #           table( doublet_types2[neighbors_that_are_doublets - n_real.cells] )
    #         neighbor_types[i,] <- neighbor_types[i,] / sum( neighbor_types[i,] )
    #       } else {
    #         neighbor_types[i,] <- NA
    #       }
    #     }
    #   }
    # }
    # print("Classifying doublets..")
    # classifications <- rep("Singlet",n_real.cells)
    # classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"
    # object@meta.data[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(object@meta.data), 1]
    # object@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    # if(!is.null(idents)){
    #   colnames(neighbor_types) = levels(doublet_types1)
    #   for(ct in levels(doublet_types1)){
    #     object@meta.data[, paste("DF.doublet.contributors",pN,pK,nExp,ct,sep="_")] <- neighbor_types[,ct]
    #   }
    # }
    # return(object)
