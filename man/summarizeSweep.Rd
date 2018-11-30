\name{summarizeSweep}
\alias{summarizeSweep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{summarizeSweep}
\description{ Summarizes results from doubletFinder_ParamSweep, computing the bimodality coefficient across pN and pK parameter space. If ground-truth doublet classifications are available, then ROC analysis is performed, enabling optimal DoubletFinder parameter selection.
}
\usage{
sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
sweep.stats <- summarizeSweep(sweep.list, GT = TRUE, GT.calls = classifications)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sweep.list}{ List of pANN vectors across pN-pK space, as produced by doubletFinder_ParamSweep.
  }
  \item{GT}{ Logical set to TRUE when ground-truth doublet classifications are available for ROC analysis. Default set to FALSE.
  }
  \item{GT.calls}{ An nCell-length character vector of ground-truth doublet classifications (e.g., "Singlet" or "Doublet") used to gauge performance of logistic regression models trained using pANN vectors during ROC analysis.
  }

}
\details{
}
\value{
Dataframe with bimodality coefficient values at each pN-pK parameter set. If GT = TRUE, dataframe also includes AUC for each pN-pK parameter set computed during ROC analysis.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
sweep.list <- doubletFinder_ParamSweep(seu)
sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
}

