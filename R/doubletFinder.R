doubletFinder <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE) {
  require(Seurat); require(fields); require(KernSmooth)

  ## Generate new list of doublet classificatons from existing pANN vector to save time
  if (reuse.pANN != FALSE ) {
    pANN.old <- seu@meta.data[ , reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[,ncol(seu@meta.data)+1] <- classifications
    colnames(seu@meta.data)[ncol(seu@meta.data)] <- paste("DF.classifications",pN,pK,nExp,sep="_")
    return(seu)
  }

  if (reuse.pANN == FALSE) {
    ## Make merged real-artifical data
    real.cells <- colnames(seu@data)
    data <- seu@raw.data[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
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
                                       do.plot = FALSE, x.low.cutoff = seu@calc.params$FindVariableGenes$x.low.cutoff,
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
    cell.names <- seu_wdoublets@cell.names
    nCells <- length(cell.names)

    ## Compute PC distance matrix
    print("Calculating PC distance matrix")
    pca.coord <- seu_wdoublets@dr$pca@cell.embeddings[, PCs]
    rm(seu_wdoublets) ## Free-up memory
    gc() ## Free-up memory
    dist.mat <- fields::rdist(pca.coord)
    colnames(dist.mat) <- cell.names
    rownames(dist.mat) <- cell.names

    ## Compute pANN
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      neighbor.names <- rownames(dist.mat)[neighbors]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
    }

    print("Classifying doublets..")
    classifications <- rep("Singlet",n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[,ncol(seu@meta.data)+1] <- pANN[seu@cell.names, 1]
    colnames(seu@meta.data)[ncol(seu@meta.data)] <- paste("pANN",pN,pK,nExp,sep="_")
    seu@meta.data[,ncol(seu@meta.data)+1] <- classifications
    colnames(seu@meta.data)[ncol(seu@meta.data)] <- paste("DF.classifications",pN,pK,nExp,sep="_")
    return(seu)
  }
}
