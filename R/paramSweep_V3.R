paramSweep_v3 <- function(seu) {
  #require(Seurat); require(fields)
  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
  pN <- seq(0.05,0.3,by=0.05)

  sweep.res.list <- list()
  list.ind <- 0
  PCs <- seu@commands$RunTSNE$dims
  orig.commands <- seu@commands

  ## Down-sample cells to 10000 (when applicable) for computational effiency
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 10000, replace=FALSE)]
    data <- seu@assays$RNA@counts[ , real.cells]
    n.real.cells <- ncol(data)
  }

  if (nrow(seu@meta.data) <= 10000){
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts
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
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

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
                            npcs = orig.commands$RunPCA.RNA$npcs,
                            rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
                            weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                            verbose=FALSE)

    ## Compute PC distance matrix
    print("Calculating PC distance matrix...")
    nCells <- nrow(seu_wdoublets@meta.data)
    pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
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





