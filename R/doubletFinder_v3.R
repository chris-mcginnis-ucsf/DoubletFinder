doubletFinder_v3 <- function(object, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE, 
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
  
  ## Make merged real-artifical data
  real.cells <- rownames(object@meta.data)
  data <- object@assays$RNA@counts[, real.cells]
  n_real.cells <- length(real.cells)
  n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
  print(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
  real.cells1 <- sample(real.cells, n_doublets, replace = T)
  real.cells2 <- sample(real.cells, n_doublets, replace = T)
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
    seu_wdoublets <- NormalizeData(seu_wdoublets,
                                   normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                   scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                   margin = orig.commands$NormalizeData.RNA@params$margin)
    
    print("Finding variable genes...")
    seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
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
    
    print("Scaling data...")
    seu_wdoublets <- ScaleData(seu_wdoublets,
                               features = orig.commands$ScaleData.RNA$features,
                               model.use = orig.commands$ScaleData.RNA$model.use,
                               do.scale = orig.commands$ScaleData.RNA$do.scale,
                               do.center = orig.commands$ScaleData.RNA$do.center,
                               scale.max = orig.commands$ScaleData.RNA$scale.max,
                               block.size = orig.commands$ScaleData.RNA$block.size,
                               min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
    
    print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets,
                            features = orig.commands$ScaleData.RNA$features,
                            npcs = length(PCs),
                            rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
                            weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                            verbose=FALSE)
    pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
    cell.names <- rownames(seu_wdoublets@meta.data)
    nCells <- length(cell.names)
    rm(seu_wdoublets); gc() # Free up memory
  } else {
    require(sctransform)
    
    print("Running SCTransform...")
    seu_wdoublets <- SCTransform(seu_wdoublets)
    
    print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
    pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
    cell.names <- rownames(seu_wdoublets@meta.data)
    nCells <- length(cell.names)
    rm(seu_wdoublets); gc()
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
    }
    
    rm(dists)
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
