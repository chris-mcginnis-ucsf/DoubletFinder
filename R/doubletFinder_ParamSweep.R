doubletFinder_ParamSweep <- function(seu) {
  require(Seurat); require(fields)
  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
  pN <- seq(0.05,0.3,by=0.05)
  
  sweep.res.list <- list()
  list.ind <- 0
  PCs <- seu@calc.params$RunTSNE$dims.use
  
  ## Down-sample cells to 10000 (when applicable) for computational effiency
  if (length(seu@cell.names) > 10000) {
    real.cells <- seu@cell.names[sample(1:nrow(seu@meta.data), 10000, replace=FALSE)]
    data <- seu@raw.data[ , real.cells]
    n.real.cells <- ncol(data)
  }
  
  if (length(seu@cell.names) <= 10000){
    real.cells <- seu@cell.names
    data <- seu@raw.data
    n.real.cells <- ncol(data)
  }
  
  ## Iterate through pN, computing pANN vectors at varying pK
  for (n in 1:length(pN)) {
    
    ## Make merged real-artifical data
    print(paste("Creating artificial doublets for pN = ", pN[n]*100,"%",sep=""))
    n_doublets <- round(n.real.cells/(1 - pN[n]) - n.real.cells)
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)

    ## Pre-process Seurat object
    print("Creating Seurat object...")    
    seu_wdoublets <- CreateSeuratObject(raw.data = data_wdoublets)
  
    print("Normalizing Seurat object...")
    seu_wdoublets <- NormalizeData(seu_wdoublets, 
                                   normalization.method = seu@calc.params$NormalizeData$normalization.method, 
                                   scale.factor = seu@calc.params$NormalizeData$scale.factor)
  
    print("Finding variable genes...")
    seu_wdoublets <- FindVariableGenes(seu_wdoublets, 
                                       do.plot = FALSE, 
                                       x.low.cutoff = seu@calc.params$FindVariableGenes$x.low.cutoff, 
                                       x.high.cutoff = seu@calc.params$FindVariableGenes$x.high.cutoff, 
                                       y.high.cutoff = seu@calc.params$FindVariableGenes$y.high.cutoff, 
                                       y.cutoff = seu@calc.params$FindVariableGenes$y.cutoff)
  
    print("Scaling data...")
    seu_wdoublets <- ScaleData(seu_wdoublets, display.progress = TRUE)
  
    print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets, 
                            pc.genes = seu_wdoublets@var.genes, 
                            pcs.print = 0, 
                            pcs.compute = length(PCs))
    
    ## Compute PC distance matrix
    print("Calculating PC distance matrix...")
    nCells <- length(seu_wdoublets@cell.names)
    pca.coord <- seu_wdoublets@dr$pca@cell.embeddings[, PCs]
    rm(seu_wdoublets)
    gc()
    dist.mat <- fields::rdist(pca.coord)[,1:n.real.cells]
    
    ## Pre-order PC distance matrix prior to iterating across pK for pANN computations
    print("Defining neighborhoods...")
    for (i in 1:n.real.cells) { 
      dist.mat[,i] <- order(dist.mat[,i]) 
    }
    
    ## Trim PC distance matrix for faster manipulations
    ind <- round(nCells * max(pK))+5
    dist.mat <- dist.mat[1:ind, ]
  
    ## Compute pANN across pK sweep
    print("Computing pANN across all pK...")
    for (k in 1:length(pK)) {
      print(paste("pK = ", pK[k], "...", sep = ""))
      pk.temp <- round(nCells * pK[k])
      pANN <- as.data.frame(matrix(0L, nrow = n.real.cells, ncol = 1))
      colnames(pANN) <- "pANN"
      rownames(pANN) <- real.cells
      list.ind <- list.ind + 1

      for (i in 1:n.real.cells) {
        neighbors <- dist.mat[2:(pk.temp + 1),i]
        pANN$pANN[i] <- length(which(neighbors > n.real.cells))/pk.temp
      }
    
      sweep.res.list[[list.ind]] <- pANN

    }
  }

  ## Assign names to list of results
  name.vec <- NULL
  for (j in 1:length(pN)) { 
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_" )) 
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
  
}





