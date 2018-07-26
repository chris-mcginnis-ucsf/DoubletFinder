#' Doublet detection in single-cell RNA sequencing data
#'
#' This function generaetes artificial nearest neighbors from existing single-cell RNA
#' sequencing data. First, real and artificial data are merged. Second, dimension reduction
#' is performed on the merged real-artificial dataset using PCA. Third, the proportion of
#' artificial nearest neighbors is defined for each real cell. Finally, real cells are rank-
#' ordered and predicted doublets are defined via thresholding based on the expected number
#' of doublets.
#'
#' @param seu A fully-processed Seurat object (i.e. after normalization, variable gene definition,
#' scaling, PCA, and tSNE).
#' @param expected.doublets The number of doublets expected to be present in the original data.
#' This value can best be estimated from cell loading densities into the 10X/Drop-Seq device.
#' @param porportion.artificial The proportion (from 0-1) of the merged real-artificial dataset
#' that is artificial. In other words, this argument defines the total number of artificial doublets.
#' Default is set to 25%, based on optimization on PBMCs (see McGinnis, Murrow and Gartner 2018, BioRxiv).
#' @param proportion.NN The proportion (from 0-1) of the merged real-artificial dataset used to define
#' each cell's neighborhood in PC space. Default set to 1%, based on optimization on PBMCs (see McGinnis,
#' Murrow and Gartner 2018, BioRxiv).
#' @return An updated Seurat object with metadata for pANN values and doublet predictions.
#' @export
#' @examples
#' seu <- doubletFinder(seu, expected.doublets = 1000, proportion.artificial = 0.25, proportion.NN = 0.01)


doubletFinder <- function(seu, expected.doublets = 0, proportion.artificial = 0.25, proportion.NN = 0.01) {

  if (expected.doublets == 0) {  stop("Need to set number of expected doublets...")  }

  ## Step 1: Generate artificial doublets from Seurat object input
  print("Creating artificial doublets...")
  data <- seu@raw.data[ , seu@cell.names]
  real.cells <- seu@cell.names
  n_real.cells <- length(real.cells)
  n_doublets <- round(n_real.cells/(1-proportion.artificial)-n_real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[ , real.cells1] + data[ , real.cells2])/2
  colnames(doublets) <- paste("X", 1:n_doublets, sep="")
  data_wdoublets <- cbind(data, doublets)

  ## Step 2: Pre-process real-artificial merged data using Seurat
  print("Creating Seurat object...")
  seu_wdoublets <- Seurat::CreateSeuratObject(raw.data = data_wdoublets)
  print("Normalizing Seurat object...")
  seu_wdoublets <- Seurat::NormalizeData(seu_wdoublets,
                                 normalization.method = seu@calc.params$NormalizeData$normalization.method,
                                 scale.factor = seu@calc.params$NormalizeData$scale.factor)
  print("Finding variable genes...")
  seu_wdoublets <- Seurat::FindVariableGenes(seu_wdoublets, do.plot = FALSE,
                                     x.low.cutoff = seu@calc.params$FindVariableGenes$x.low.cutoff,
                                     x.high.cutoff = seu@calc.params$FindVariableGenes$x.high.cutoff,
                                     y.high.cutoff = seu@calc.params$FindVariableGenes$y.high.cutoff,
                                     y.cutoff = seu@calc.params$FindVariableGenes$y.cutoff)
  print("Scaling data...")
  seu_wdoublets <- Seurat::ScaleData(seu_wdoublets, display.progress = TRUE)
  print("Running PCA...")
  PCs <- seu@calc.params$RunTSNE$dims.use
  if (length(PCs) == 0) {  stop("Need to run tSNE on original Seurat object...")  }
  seu_wdoublets <- Seurat::RunPCA(seu_wdoublets, pc.genes = seu_wdoublets@var.genes, pcs.print = 0,pcs.compute = length(PCs))

  ## Step 3: Record cell names and number
  cell.names <- seu_wdoublets@cell.names
  nCells <- length(cell.names)

  ## Step 4: Create distance matrix from PCA results
  print("Calculating PC distance matrix...")

  pca.coord <- seu_wdoublets@dr$pca@cell.embeddings[, PCs]
  rm(seu_wdoublets)  ## Free up memory
  gc()               ## Free up memory
  dist.mat <- as.matrix(dist(pca.coord))
  dist.mat <- dist.mat[,-grep("X", colnames(dist.mat))]

  ## Step 5: Initialize pANN structure
  pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
  rownames(pANN) <- real.cells
  colnames(pANN) <- "pANN"

  ## Step 6: Calculate pANN
  k <- round(nCells*proportion.NN)
  for (i in 1:n_real.cells) {
    neighbors <- order(dist.mat[,i])
    neighbors <- neighbors[2:(k+1)]
    neighbor.names <- rownames(dist.mat)[neighbors]
    pANN[i,1] <- length(grep("X", neighbor.names)) / k
  }

  ## Step 7: Add results to metadata slot of original Seurat object
  seu <- Seurat::AddMetaData(seu, metadata = pANN, col.name = "pANN")
  predictions <- as.data.frame(rep("Singlet", n_real.cells), ncol = 1, stringsAsFactors = FALSE)
  rownames(predictions) <- real.cells
  doublet.predictions <- rownames(seu@meta.data)[order(seu@meta.data$pANN, decreasing = TRUE)]
  doublet.predictions <- doublet.predictions[1:expected.doublets]
  predictions[doublet.predictions, ] <- "Doublet"
  colnames(predictions) <- "pANNPredictions"
  seu <- Seurat::AddMetaData(seu, metadata = predictions, col.name = "pANNPredictions")
  return(seu)
}
