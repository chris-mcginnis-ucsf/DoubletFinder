#' Doublet detection in single-cell RNA sequencing data
#'
#' This function generaetes artificial nearest neighbors from existing single-cell RNA
#' sequencing data. First, real and artificial data are merged. Second, dimension reduction
#' is performed on the merged real-artificial dataset using PCA. Third, the proportion of
#' artificial nearest neighbors is defined for each real cell. Finally, real cells are rank-
#' ordered and predicted doublets are defined via thresholding based on the expected number
#' of doublets. This function has methods for both Seurat and loom objects. To use the method for loom,
#' please install loomR (from GitHub, not on CRAN yet) and the loom branch of Seurat. See
#' https://satijalab.org/seurat/mca_loom.html for more instructions.
#'
#' @param obj A fully-processed Seurat or loom object (i.e. after normalization, variable gene definition,
#' scaling, PCA, and tSNE).
#' @rdname doubletFinder
#' @export doubletFinder
#'
doubletFinder <- function(obj, ...) {
  UseMethod(generic = "doubletFinder", object = obj)
}

#' @param expected.doublets The number of doublets expected to be present in the original data.
#' This value can best be estimated from cell loading densities into the 10X/Drop-Seq device.
#' @param proportion.artificial The proportion (from 0-1) of the merged real-artificial dataset
#' that is artificial. In other words, this argument defines the total number of artificial doublets.
#' Default is set to 25\%, based on optimization on PBMCs (see McGinnis, Murrow and Gartner 2018, BioRxiv).
#' @param proportion.NN The proportion (from 0-1) of the merged real-artificial dataset used to define
#' each cell's neighborhood in PC space. Default set to 1\%, based on optimization on PBMCs (see McGinnis,
#' Murrow and Gartner 2018, BioRxiv).
#' @return An updated Seurat object with metadata for pANN values and doublet predictions.
#' @describeIn doubletFinder Doublet detection in single-cell RNA sequencing data
#' @export doubletFinder.seurat
#' @method doubletFinder seurat
#' @examples
#' seu <- doubletFinder(obj, expected.doublets = 1000, proportion.artificial = 0.25, proportion.NN = 0.01)
#'
doubletFinder.seurat <- function(obj, expected.doublets = 0, proportion.artificial = 0.25, proportion.NN = 0.01) {

  if (expected.doublets == 0) {  stop("Need to set number of expected doublets...")  }

  ## Step 1: Generate artificial doublets from Seurat object input
  print("Creating artificial doublets...")
  real.cells <- colnames(obj@data)
  data <- obj@raw.data[ , real.cells]
  n_real.cells <- length(real.cells)
  n_doublets <- round(n_real.cells/(1-proportion.artificial)-n_real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[ , real.cells1] + data[ , real.cells2])/2
  colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
  data_wdoublets <- cbind(data, doublets)

  ## Step 2: Pre-process real-artificial merged data using Seurat
  print("Creating Seurat object...")
  seu_wdoublets <- Seurat::CreateSeuratObject(raw.data = data_wdoublets)
  print("Normalizing Seurat object...")
  seu_wdoublets <- Seurat::NormalizeData(seu_wdoublets,
                                 normalization.method = obj@calc.params$NormalizeData$normalization.method,
                                 scale.factor = obj@calc.params$NormalizeData$scale.factor)
  print("Finding variable genes...")
  seu_wdoublets <- Seurat::FindVariableGenes(seu_wdoublets, do.plot = FALSE,
                                     x.low.cutoff = obj@calc.params$FindVariableGenes$x.low.cutoff,
                                     x.high.cutoff = obj@calc.params$FindVariableGenes$x.high.cutoff,
                                     y.high.cutoff = obj@calc.params$FindVariableGenes$y.high.cutoff,
                                     y.cutoff = obj@calc.params$FindVariableGenes$y.cutoff)
  print("Scaling data...")
  seu_wdoublets <- Seurat::ScaleData(seu_wdoublets, display.progress = TRUE)
  print("Running PCA...")
  PCs <- obj@calc.params$RunTSNE$dims.use
  if (length(PCs) == 0) {  stop("Need to run tSNE on original Seurat object...")  }
  seu_wdoublets <- Seurat::RunPCA(seu_wdoublets, pc.genes = seu_wdoublets@var.genes, pcs.print = 0, pcs.compute = length(PCs))

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
  obj <- Seurat::AddMetaData(obj, metadata = pANN, col.name = "pANN")
  predictions <- as.data.frame(rep("Singlet", n_real.cells), ncol = 1, stringsAsFactors = FALSE)
  rownames(predictions) <- real.cells
  doublet.predictions <- rownames(obj@meta.data)[order(obj@meta.data$pANN, decreasing = TRUE)]
  doublet.predictions <- doublet.predictions[1:expected.doublets]
  predictions[doublet.predictions, ] <- "Doublet"
  colnames(predictions) <- "pANNPredictions"
  obj <- Seurat::AddMetaData(obj, metadata = predictions, col.name = "pANNPredictions")
  return(obj)
}

#' @importFrom progress progress_bar
NULL
#'
#' @param doublet.block.size Number of doublets to add at a time, so not to contain the entire
#' doublet matrix in memory in case there are many doublets.
#' @param gene.names Dataset name in loom object for gene names; duplicates are not allowed.
#' @param cell.names Dataset name in loom object for cell names; duplicates are not allowed.
#' @param pcs.compute Number of principal components to compute to determine distance between
#' cells.
#' @param overwrite Whether to overwrite existing temporary file for the data with simulated
#' doublets.
#' @param chunk.size Chunk size to iterate over
#' @param chunk.dims HDF5 chunk dimensions
#' @param \dots Extra arguments to RunPCA.loom, such as \code{online.pca}, which determines
#' whether the online PCA method should be used. Online PCA should only be used for very large
#' datasets.
#' @return (loom) Nothing; loom objects are modified in place. This function will add pANN values
#' and doublet predictions to column attributes.
#' @describeIn doubletFinder Doublet detection in single-cell RNA sequencing data
#' @export doubletFinder.loom
#' @method doubletFinder loom
#' @examples doubletFinder(pbmc_loom, expected.doublets = 1000, proportion.artificial = 0.25,
#'               proportion.NN = 0.01, gene.names = "row_attrs/Accession",
#'               cell.names = "col_attrs/CellID",
#'               overwrite = TRUE, pcs.compute = 10, online.pca = FALSE)

doubletFinder.loom <- function(obj, expected.doublets = 0,
                               proportion.artificial = 0.25, proportion.NN = 0.01,
                               doublet.block.size = 50,
                               gene.names = "row_attrs/gene_names",
                               cell.names = "col_attrs/cell_names",
                               overwrite = FALSE, pcs.compute = 55,
                               chunk.size = 1000, chunk.dims = c(256,256), ...) {
  if (expected.doublets == 0) {  stop("Need to set number of expected doublets...")  }
  ## Step 1: Generate artificial doublets from loom object input
  print("Creating artificial doublets...")
  real.cells <- obj[[cell.names]][]
  genes <- obj[[gene.names]][]
  n_real.cells <- length(real.cells)
  n_doublets <- round(n_real.cells/(1 - proportion.artificial) - n_real.cells)
  nCells <- n_real.cells + n_doublets
  cells_both <- c(real.cells, paste0("X", 1:n_doublets))

  # Make temporary loom file for data with simulated doublets
  if (overwrite && file.exists("loom_wdoublets.loom")) {
    invisible(file.remove("loom_wdoublets.loom"))
  }
  invisible(file.copy(obj$filename, "loom_wdoublets.loom"))
  # Generate doublets
  real_ind1 <- sample(1:n_real.cells, n_doublets, replace = TRUE)
  real_ind2 <- sample(1:n_real.cells, n_doublets, replace = TRUE)
  # If not chunked, this will load about half of the data into memory
  # To conserve memory, get doublets by chunk
  cn <- strsplit(cell.names, split = "/")[[1]][2]
  loom_wdoublets <- connect("loom_wdoublets.loom", "r+")
  # Delete the temporary file and clean up on exit
  on.exit({loom_wdoublets$close_all(); invisible(file.remove("loom_wdoublets.loom"))
    gc(verbose = FALSE)
    hdf5r::h5garbage_collect()})
  if (n_doublets <= doublet.block.size) {
    doublet_mat <-
      (obj[["matrix"]][real_ind1,] + obj[["matrix"]][real_ind2,]) / 2
    doublet_attrs <- extend_col_attrs(n_doublets, loom_wdoublets, cell.names = cn)
    loom_wdoublets$add.cells(matrix.data = doublet_mat,
                             attributes.data = doublet_attrs,
                             layers.data = list(norm_data = doublet_mat, scale_data = doublet_mat),
                             do.transpose = FALSE, display.progress = FALSE)
  } else {
    n_chunks <- floor(n_doublets / doublet.block.size) # the last chunk is added separately
    pb <- progress_bar$new(format = "|:bar| :percent :eta")
    for (i in 1:n_chunks) {
      ind1 <- 1 + (i - 1) * doublet.block.size
      ind2 <- i * doublet.block.size
      doublet_mat <-
        (obj[["matrix"]][real_ind1[ind1:ind2],] + obj[["matrix"]][real_ind2[ind1:ind2],]) / 2
      doublet_attrs <- extend_col_attrs(doublet.block.size, loom_wdoublets, cell.names = cn)
      loom_wdoublets$add.cells(matrix.data = doublet_mat,
                               attributes.data = doublet_attrs,
                               layers.data = list(norm_data = doublet_mat, scale_data = doublet_mat),
                               do.transpose = FALSE, display.progress = FALSE)
      pb$update(i / (n_chunks + 1))
    }
    remainder <- n_doublets %% doublet.block.size
    if (remainder > 0) {
      doublet_mat <-
        (obj[["matrix"]][real_ind1[(ind2 + 1):n_doublets],] + obj[["matrix"]][real_ind2[(ind2 + 1):n_doublets],]) / 2
      doublet_attrs <- extend_col_attrs(remainder, loom_wdoublets, cell.names = cn)
      loom_wdoublets$add.cells(matrix.data = doublet_mat,
                               attributes.data = doublet_attrs,
                               layers.data = list(norm_data = doublet_mat, scale_data = doublet_mat),
                               do.transpose = FALSE, display.progress = FALSE)
      pb$update(1)
    }
    pb$terminate()
  }
  loom_wdoublets[[cell.names]][] <- cells_both
  # Free up memory
  rm(doublet_mat)
  gc(verbose = FALSE)
  hdf5r::h5garbage_collect()
  ## Step 2: Pre-process real-artificial merged data using Seurat
  print("Normalizing Seurat object...")
  Seurat::NormalizeData(loom_wdoublets, overwrite = TRUE,
                scale.factor = as.integer(hdf5r::h5attributes(obj[["layers/norm_data"]])$calc_params[1]),
                chunk.size = chunk.size, chunk.dims = chunk.dims)
  print("Finding variable genes...")
  fvg_params <- as.numeric(hdf5r::h5attributes(obj[["row_attrs/var_genes"]])$calc_params[1:4])
  Seurat::FindVariableGenes(loom_wdoublets, overwrite = TRUE,
                            x.low.cutoff = fvg_params[1],
                            x.high.cutooff = fvg_params[2],
                            y.cutoff = fvg_params[3],
                            y.high.cutoff = fvg_params[4],
                            chunk.size = chunk.size)
  print("Scaling data...")
  Seurat::ScaleData(loom_wdoublets, display.progress = TRUE, overwrite = TRUE,
                    chunk.size = chunk.size, chunk.dims = chunk.dims)
  print("Running PCA...")
  # Get the names of the highly variable genes
  var_gene_inds <- loom_wdoublets[["row_attrs/var_genes"]][]
  var_genes <- genes[var_gene_inds]
  Seurat::RunPCA(loom_wdoublets, pcs.compute = pcs.compute, overwrite = TRUE,
                 pc.genes = var_genes, gene.names = gene.names, cell.names = cell.names,
                 do.print = FALSE,
                 chunk.size = chunk.size, ...)

  ## Step 4: Create distance matrix from PCA results
  print("Calculating PC distance matrix...")
  # This will be loaded into memory; since the number of PCs computed is small,
  # this shouldn't be too bad
  pca.coord <- t(loom_wdoublets$col.attrs$pca_cell_embeddings[,])
  dist.mat <- as.matrix(dist(pca.coord))
  # Add column names
  colnames(dist.mat) <- rownames(dist.mat) <- cells_both
  dist.mat <- dist.mat[,-grep("X", colnames(dist.mat))]

  ## Step 5: Initialize pANN structure
  pANN <- data.frame(pANN = rep(0, n_real.cells),
                     pANNPredictions = "Singlet",
                     stringsAsFactors = FALSE)
  #rownames(pANN) <- real.cells

  ## Step 6: Calculate pANN
  print("Calculating pANN")
  k <- round(nCells * proportion.NN)
  for (i in 1:n_real.cells) {
    neighbors <- order(dist.mat[,i])
    neighbors <- neighbors[2:(k + 1)]
    neighbor.names <- rownames(dist.mat)[neighbors]
    pANN[i,1] <- length(grep("X", neighbor.names)) / k
  }

  ## Step 7: Add results to column attributes of the original loom object
  doublet_inds <- order(pANN$pANN, decreasing = TRUE)[1:expected.doublets]
  pANN$pANNPredictions[doublet_inds] <- "Doublet"
  obj$add.attribute(pANN, MARGIN = 2, overwrite = TRUE)

  invisible(x = obj)
}
