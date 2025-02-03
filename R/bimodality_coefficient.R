#' bimodality_coefficient
#'
#' Internal function to compute bimodality coefficient during BCmvn computation
#' and pK estimation.
#'
#'
#' @param x Guassian kernel density estimation representing pANN distribution
#' @return Bimodality coefficient value
#' @author Chris McGinnis
#' @references Taken from the 'modes' R package (v0.7)

bimodality_coefficient <- function(x) {
  G <- skewness(x)
  sample.excess.kurtosis <- kurtosis(x)
  K <- sample.excess.kurtosis
  n <- length(x)
  B <- ((G^2)+1)/(K+((3*((n-1)^2))/((n-2)*(n-3))))
  return(B)
}
