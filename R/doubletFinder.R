#' doubletFinder
#'
#' Core doublet prediction function of the DoubletFinder package. Generates
#' artifical doublets from an existing, pre-processed Seurat object. Real and
#' artificial data are then merged and pre-processed using parameters utilized
#' for the existing Seurat object. PC distance matrix is then computed and used
#' the measure the proportion of artificial nearest neighbors (pANN) for every
#' real cell. pANN is then thresholded according to the number of expected
#' doublets to generate final doublet predictions.
#'
#'
#' @param seu A fully-processed Seurat object (i.e., After NormalizeData,
#' FindVariableGenes, ScaleData, and RunPCA have all been performed).
#' @param PCs Number of statistically-significant principal components (e.g.,
#' as estimated from PC elbow plot)
#' @param pN The number of generated artificial doublets, expressed as a
#' proportion of the merged real-artificial data. Default is set to 0.25, based
#' on observation that DoubletFinder performance is largely pN-invariant (see
#' McGinnis, Murrow and Gartner 2019, Cell Systems).
#' @param pK The PC neighborhood size used to compute pANN, expressed as a
#' proportion of the merged real-artificial data. No default is set, as pK
#' should be adjusted for each scRNA-seq dataset. Optimal pK values can be
#' determined using mean-variance-normalized bimodality coefficient.
#' @param nExp The total number of doublet predictions produced. This value can
#' best be estimated from cell loading densities into the 10X/Drop-Seq device,
#' and adjusted according to the estimated proportion of homotypic doublets.
#' @param reuse.pANN Seurat metadata column name for previously-generated pANN
#' results. Argument should be set to NULL (default) for initial DoubletFinder
#' runs. Enables fast adjusting of doublet predictions for different nExp.
#' @param sct Logical representing whether SCTransform was used during original
#' Seurat object pre-processing (default = FALSE).
#' @param annotations annotations.
#' @return Seurat object with updated metadata including pANN and doublet
#' classifications.
#' @author Chris McGinnis
#' @export
#' @importFrom Seurat GetAssayData Cells CreateSeuratObject NormalizeData
#'   FindVariableFeatures ScaleData RunPCA SCTransform
#' @importFrom SeuratObject LayerData
#' @examples
#'
#' data(pbmc_small)
#' seu <- pbmc_small
#' nExp_poi <- round(0.15*nrow(seu@@meta.data))
#' seu <- doubletFinder(seu, PCs = 1:10,
#' pN = 0.25, pK = 0.01,
#' nExp = nExp_poi,
#' reuse.pANN = NULL, sct=FALSE)
#'
#' ## With homotypic adjustment
#' annotations <- seu$RNA_snn_res.1
#' homotypic.prop <- modelHomotypic(annotations)
#' nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#' seu <- doubletFinder(seu, PCs = 1:10, pN = 0.25,
#'                      pK = 0.01, nExp = nExp_poi.adj,
#'                      reuse.pANN = NULL,
#'                      sct=FALSE)
#'
doubletFinder <- function(seu,
                          PCs,
                          pN = 0.25,
                          pK,
                          nExp,
                          reuse.pANN = NULL,
                          sct = FALSE,
                          annotations = NULL) {

  ## Generate new list of doublet classificatons from existing pANN vector to save time
  if (is.null(reuse.pANN)) {
    pANN.old <- seu@meta.data[ , reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    return(seu)
  } else {
    ## Make merged real-artifical data
    real.cells <- rownames(seu@meta.data)
    if (SeuratObject::Version(seu)>= '5.0') {
         counts <- LayerData(seu, assay = "RNA", layer = "counts")
    } else {
         counts <- GetAssayData(object = seu, assay = "RNA", slot = "counts")
    }
    data <- counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    # Keep track of the types of the simulated doublets
    if(!is.null(annotations)){
      stopifnot(typeof(annotations)=="character")
      stopifnot(length(annotations)==length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    ## Store important pre-processing information
    orig.commands <- seu@commands

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
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc() # Free up memory
    } else {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)

      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc()
    }

    ## Compute PC distance matrix
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)

    ## Compute pANN
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    if(!is.null(annotations)){
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      if(!is.null(annotations)){
        for(ct in unique(annotations)){
          neighbors_that_are_doublets = neighbors[neighbors>n_real.cells]
          if(length(neighbors_that_are_doublets) > 0){
            neighbor_types[i,] <-
              table( doublet_types1[neighbors_that_are_doublets - n_real.cells] ) +
              table( doublet_types2[neighbors_that_are_doublets - n_real.cells] )
            neighbor_types[i,] <- neighbor_types[i,] / sum( neighbor_types[i,] )
          } else {
            neighbor_types[i,] <- NA
          }
        }
      }
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet",n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(seu@meta.data), 1]
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    if(!is.null(annotations)){
      colnames(neighbor_types) = levels(doublet_types1)
      for(ct in levels(doublet_types1)){
        seu@meta.data[, paste("DF.doublet.contributors",pN,pK,nExp,ct,sep="_")] <- neighbor_types[,ct]
      }
    }
    return(seu)
  }
}
