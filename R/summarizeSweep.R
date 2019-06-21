summarizeSweep <- function(sweep.list, GT = FALSE, GT.calls = NULL) {
  require(KernSmooth); require(ROCR); require(modes)
  ## Set pN-pK param sweep ranges
  name.vec <- names(sweep.list)
  name.vec <- unlist(strsplit(name.vec, split="pN_"))
  name.vec <- name.vec[seq(2, length(name.vec), by=2)]
  name.vec <- unlist(strsplit(name.vec, split="_pK_"))
  pN <- as.numeric(unique(name.vec[seq(1, length(name.vec), by=2)]))
  pK <- as.numeric(unique(name.vec[seq(2, length(name.vec), by=2)]))

  ## Initialize data structure w/ or w/o AUC column, depending on whether ground-truth doublet classifications are available
  if (GT == TRUE) {
    sweep.stats <- as.data.frame(matrix(0L, nrow=length(sweep.list), ncol=4))
    colnames(sweep.stats) <- c("pN","pK","AUC","BCreal")
    sweep.stats$pN <- factor(rep(pN, each=length(pK), levels = pN))
    sweep.stats$pK <- factor(rep(pK, length(pN),levels = pK))
  }

  if (GT == FALSE) {
    sweep.stats <- as.data.frame(matrix(0L, nrow=length(sweep.list), ncol=3))
    colnames(sweep.stats) <- c("pN","pK","BCreal")
    sweep.stats$pN <- factor(rep(pN, each=length(pK), levels = pN))
    sweep.stats$pK <- factor(rep(pK, length(pN),levels = pK))
  }

  ## Perform pN-pK parameter sweep summary
  for (i in 1:length(sweep.list)) {
    res.temp <- sweep.list[[i]]

    ## Use gaussian kernel density estimation of pANN vector to compute bimodality coefficient
    gkde <- approxfun(bkde(res.temp$pANN, kernel="normal"))
    x <- seq(from=min(res.temp$pANN), to=max(res.temp$pANN), length.out=nrow(res.temp))
    sweep.stats$BCreal[i] <- bimodality_coefficient(gkde(x))

    if (GT == FALSE) { next }

    ## If ground-truth doublet classifications are available, perform ROC analysis on logistic
    ## regression model trained using pANN vector
    meta <- as.data.frame(matrix(0L, nrow=nrow(res.temp), ncol=2))
    meta[,1] <- GT.calls
    meta[,2] <- res.temp$pANN
    train.ind <- sample(1:nrow(meta), round(nrow(meta)/2), replace=FALSE)
    test.ind <- (1:nrow(meta))[-train.ind]
    colnames(meta) <- c("SinDub","pANN")
    meta$SinDub <- factor(meta$SinDub, levels = c("Doublet","Singlet"))
    model.lm <- glm(SinDub ~ pANN, family="binomial"(link='logit'), data=meta, subset=train.ind)
    prob <- predict(model.lm, newdata=meta[test.ind, ], type="response")
    ROCpred <- ROCR::prediction(predictions=prob, labels=meta$SinDub[test.ind])
    perf.auc <- ROCR::performance(ROCpred, measure="auc")
    sweep.stats$AUC[i] <- perf.auc@y.values[[1]]
  }

  return(sweep.stats)

}
