#' parallel_paramSweep
#'
#' Internal parallelization function for paramSweep.
#'
#'
#' @param n pN iteration counter.
#' @param n.real.cells Number of real cells. Set automatically during
#' paramSweep_v3.
#' @param real.cells Vector of real cell IDs. Set automatically during
#' paramSweep_v3.
#' @param pK The PC neighborhood size used to compute pANN, expressed as a
#' proportion of the merged real-artificial data. No default is set, as pK
#' should be adjusted for each scRNA-seq dataset. Optimal pK values can be
#' determined using mean-variance-normalized bimodality coefficient.
#' @param pN The number of generated artificial doublets, expressed as a
#' proportion of the merged real-artificial data. Default is set to 0.25, based
#' on observation that DoubletFinder performance is largely pN-invariant (see
#' McGinnis, Murrow and Gartner 2019, Cell Systems).
#' @param data Count matrix. Set automatically during paramSweep_v3.
#' @param orig.commands Count matrix. Set automatically during paramSweep_v3.
#' @param PCs Number of statistically-sigificant PCs. Set according to
#' paramSweep_v3 arguments.
#' @param sct Logical representing whether Seurat object was pre-processed
#' using 'sctransform'. Set according to paramSweep_v3 arguments (default = F).
#' @return Parallelization function compatible with mclapply.
#' @author Nathan Skeene, June 2019.
#' @examples
#'
#'
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
  if (!sct) {
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
  } else {
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
