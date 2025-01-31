#' kurtosis
#'
#' Internal function to compute bimodality coefficient during BCmvn computation
#' and pK estimation.
#'
#'
#' @param x Guassian kernel density estimation representing pANN distribution
#' @return Kurtosis value
#' @author Chris McGinnis
#' @references Taken from the 'modes' R package (v0.7).
#' @examples
#'
#'
kurtosis <- function (x) {
  n <- length(x)
  K <- (1/n)*sum((x-mean(x))^4)/(((1/n)*sum((x-mean(x))^2))^2)-3
  K <- ((n - 1)*((n+1)*K-3*(n-1))/((n-2)*(n-3)))+3
  return(K)
}
