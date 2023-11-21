skewness <- function(x) {
  n <- length(x)
  S <- (1/n)*sum((x-mean(x))^3)/(((1/n)*sum((x-mean(x))^2))^1.5)
  S <- S*(sqrt(n*(n-1)))/(n-2)
  return(S)
}
