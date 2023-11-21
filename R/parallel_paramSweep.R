parallel_paramSweep <- function(n, n.real.cells, real.cells, pK, pN, data, orig.commands, PCs, sct)  {

  sweep.res.list = list()
  list.ind = 0

  ## Make merged real-artifical data
  print(paste("Creating artificial doublets for pN = ", pN[n]*100,"%",sep=""))
  n_doublets <- round(n.real.cells/(1 - pN[n]) - n.real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[, real.cells1] + data[, real.cells2])/2
  colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
  data_wdoublets <- cbind(data, doublets)

  ## Pre-process Seurat object
  if (sct == FALSE) {
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
                            npcs = length(PCs),
                            rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
                            weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                            verbose=FALSE)
  }

  if (sct == TRUE) {
    require(sctransform)
    print("Creating Seurat object...")
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

    print("Running SCTransform...")
    seu_wdoublets <- SCTransform(seu_wdoublets)

    print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
  }

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

  return(sweep.res.list)
}
