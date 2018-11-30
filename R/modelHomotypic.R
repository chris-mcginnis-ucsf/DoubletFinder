modelHomotypic <- function(annotations) {
  anno.freq <- table(annotations)/length(annotations)
  x <- sum(anno.freq^s)
  return(x)
}