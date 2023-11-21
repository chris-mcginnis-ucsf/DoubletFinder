find.pK <- function(sweep.stats) {

  ## Implementation for data without ground-truth doublet classifications
  '%ni%' <- Negate('%in%')
  if ("AUC" %ni% colnames(sweep.stats) == TRUE) {
    ## Initialize data structure for results storage
    bc.mvn <- as.data.frame(matrix(0L, nrow=length(unique(sweep.stats$pK)), ncol=5))
    colnames(bc.mvn) <- c("ParamID","pK","MeanBC","VarBC","BCmetric")
    bc.mvn$pK <- unique(sweep.stats$pK)
    bc.mvn$ParamID <- 1:nrow(bc.mvn)

    ## Compute bimodality coefficient mean, variance, and BCmvn across pN-pK sweep results
    x <- 0
    for (i in unique(bc.mvn$pK)) {
      x <- x + 1
      ind <- which(sweep.stats$pK == i)
      bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"])/(sd(sweep.stats[ind, "BCreal"])^2)
    }

    ## Plot for visual validation of BCmvn distribution
    par(mar=rep(1,4))
    x <- plot(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, pch=16, col="#41b6c4", cex=0.75)
    x <- lines(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, col="#41b6c4")
    print(x)

    return(bc.mvn)

  }

  ## Implementation for data with ground-truth doublet classifications (e.g., MULTI-seq, CellHashing, Demuxlet, etc.)
  if ("AUC" %in% colnames(sweep.stats) == TRUE) {
    ## Initialize data structure for results storage
    bc.mvn <- as.data.frame(matrix(0L, nrow=length(unique(sweep.stats$pK)), ncol=6))
    colnames(bc.mvn) <- c("ParamID","pK","MeanAUC","MeanBC","VarBC","BCmetric")
    bc.mvn$pK <- unique(sweep.stats$pK)
    bc.mvn$ParamID <- 1:nrow(bc.mvn)

    ## Compute bimodality coefficient mean, variance, and BCmvn across pN-pK sweep results
    x <- 0
    for (i in unique(bc.mvn$pK)) {
      x <- x + 1
      ind <- which(sweep.stats$pK == i)
      bc.mvn$MeanAUC[x] <- mean(sweep.stats[ind, "AUC"])
      bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"])/(sd(sweep.stats[ind, "BCreal"])^2)
    }

    ## Plot for visual validation of BCmvn distribution
    par(mar=rep(1,4))
    x <- plot(x=bc.mvn$ParamID, y=bc.mvn$MeanAUC, pch=18, col="black", cex=0.75,xlab=NA, ylab = NA)
    x <- lines(x=bc.mvn$ParamID, y=bc.mvn$MeanAUC, col="black", lty=2)
    par(new=TRUE)
    x <- plot(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, pch=16, col="#41b6c4", cex=0.75)
    axis(side=4)
    x <- lines(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, col="#41b6c4")
    print(x)

    return(bc.mvn)

  }
}
