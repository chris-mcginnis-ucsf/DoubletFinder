bimodality_coefficient <- function(x) {
  G <- skewness(x)
  sample.excess.kurtosis <- kurtosis(x)
  K <- sample.excess.kurtosis
  n <- length(x)
  B <- ((G^2)+1)/(K+((3*((n-1)^2))/((n-2)*(n-3))))
  return(B)
}
