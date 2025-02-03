#' paramSweep
#'
#' Performs pN-pK parameter sweeps on a 10,000-cell subset of a pre-processed
#' Seurat object. Will use all cells if Seurat object contains less than 10,000
#' cells. Results are fed into 'summarizeSweep' and 'find.pK' functions during
#' optimal pK parameter selection workflow. Parameters tested: pN = 0.05-0.3,
#' pK = 0.0005-0.3.
#'
#'
#' @param seu A fully-processed Seurat object (i.e., After NormalizeData,
#' FindVariableGenes, ScaleData, RunPCA, and RunTSNE have all been performed).
#' @param PCs Number of statistically-significant principal components (e.g.,
#' as estimated from PC elbow plot)
#' @param sct Logical representing whether SCTransform was used during original
#' Seurat object pre-processing (default = FALSE).
#' @param num.cores Number of cores to use for parallelization, default=1.
#' @return List of pANN vectors for every pN and pK combination. Output also
#' contains pANN information for artificial doublets.
#' @author Chris McGinnis
#' @importFrom parallel makeCluster stopCluster  mclapply
#' @export
#' @examples
#' data(pbmc_small)
#' seu <- pbmc_small
#' sweep.list <- paramSweep(seu, PCs = 1:10, sct=FALSE)
#' sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
#' bcmvn <- find.pK(sweep.stats)
#'
paramSweep <- function(seu, PCs=1:10, sct = FALSE, num.cores=1) {
  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
  pN <- seq(0.05,0.3,by=0.05)

  ## Remove pK values with too few cells
  min.cells <- round(nrow(seu@meta.data)/(1-0.05) - nrow(seu@meta.data))
  pK.test <- round(pK*min.cells)
  pK <- pK[which(pK.test >= 1)]

  ## Extract pre-processing parameters from original data analysis workflow
  orig.commands <- seu@commands

 if (SeuratObject::Version(seu)>= '5.0') {
      counts <- LayerData(seu, assay = "RNA", layer = "counts")
  } else {
     counts <- GetAssayData(object = seu, assay = "RNA", slot = "counts")
  }

  ## Down-sample cells to 10000 (when applicable) for computational effiency
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 10000, replace=FALSE)]
	# counts <- seu@assays$RNA@layers$counts
	# rownames(counts) <- rownames(seu)
	# colnames(counts) <- colnames(seu)
    data <- counts[ , real.cells]
    n.real.cells <- ncol(data)
  } else {
    real.cells <- rownames(seu@meta.data)
    data <- counts
	# rownames(data) <- rownames(seu)
	# colnames(data) <- real.cells
    n.real.cells <- ncol(data)
  }

  ## Iterate through pN, computing pANN vectors at varying pK
  #no_cores <- detectCores()-1
  if(num.cores>1){
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)),
                        FUN = parallel_paramSweep,
                        n.real.cells,
                        real.cells,
                        pK,
                        pN,
                        data,
                        orig.commands,
                        PCs,
                        sct,mc.cores=num.cores)
    stopCluster(cl)
  }else{
    output2 <- lapply(as.list(1:length(pN)),
                      FUN = parallel_paramSweep,
                      n.real.cells,
                      real.cells,
                      pK,
                      pN,
                      data,
                      orig.commands,
                      PCs,
                      sct)
  }

  ## Write parallelized output into list
  sweep.res.list <- list()
  list.ind <- 0
  for(i in 1:length(output2)){
    for(j in 1:length(output2[[i]])){
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
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
