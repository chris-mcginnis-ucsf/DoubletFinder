\name{find.pK}
\alias{find.pK}
\title{find.pK}
\description{ Computes and visualizes the mean-variance normalized bimodality coefficient (BCmvn) score for each pK value tested during doubletFinder_ParamSweep. Optimal pK for any scRNA-seq data can be manually discerned as maxima in BCmvn distributions. If ground-truth doublet classifications are available, BCmvn is plotted along with mean ROC AUC for each pK.
}
\usage{ bcmvn <- find.pK(sweep.stats) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sweep.stats}{ pN-pK bimodality coefficient dataframe as produced by summarizeSweep.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Dataframe of mean BC, BC variance, and BCmvn scores for each pK value. Includes mean AUC for each pK value if ground-truth doublet classifications are utilized during summarizeSweep.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
}
\examples{
sweep.list <- doubletFinder_ParamSweep(seu)
sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
}

