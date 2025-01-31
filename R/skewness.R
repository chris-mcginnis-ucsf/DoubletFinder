#' skewness
#'
#' Internal function to compute skewness during BCmvn computation and pK
#' estimation.
#'
#'
#' @param x Guassian kernel density estimation representing pANN distribution
#' @return Skewness value
#' @author Chris McGinnis
#' @references Taken from the 'modes' R package (v0.7).
#' @examples
#'
#' ## Internal to bimodality_coefficient
#' G <- skewness(x)
#'
skewness <- function(x) {
  n <- length(x)
  S <- (1/n)*sum((x-mean(x))^3)/(((1/n)*sum((x-mean(x))^2))^1.5)
  S <- S*(sqrt(n*(n-1)))/(n-2)
  return(S)
}
