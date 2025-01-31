#' modelHomotypic
#'
#' Leverages user-provided cell annotations to model the proportion of
#' homotypic doublets. Building on the assumption that literature-supported
#' annotations reflect real transcriptional divergence, homotypic doublet
#' proportions are modeled as the sum of squared annotation frequencies.
#'
#'
#' @param annotations An nCell-length character vector of annotations.
#' @return Numeric proportion of homotypic doublets.
#' @author Chris McGinnis
#' @export
#' @examples
#'
#' ## Initial run, nExp set to Poisson loading estimate (e.g., 913 total
#' doublet predictions)
#' nExp_poi <- round(0.15*length(seu@cell.names))
#' seu <- doubletFinder(seu,
#' pN = 0.25,
#' pK = 0.01,
#' nExp = nExp_poi,
#' reuse.pANN = FALSE)
#'
#' ## With homotypic adjustment
#' homotypic.prop <- modelHomotypic(annotations)
#' nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#' seu <- doubletFinder(seu, pN = 0.25,
#'  pK = 0.01,
#'  nExp = nExp_poi.adj,
#'  reuse.pANN = "pANN_0.25_0.01_913")
#'
modelHomotypic <- function(annotations) {
  anno.freq <- table(annotations)/length(annotations)
  x <- sum(anno.freq^2)
  return(x)
}
