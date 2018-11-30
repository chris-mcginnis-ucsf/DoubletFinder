\name{doubletFinder_ParamSweep}
\alias{doubletFinder_ParamSweep}
\title{doubletFinder_ParamSweep
}
\description{
Performs pN-pK parameter sweeps on a 10,000-cell subset of a pre-processed Seurat object. Will use all cells if Seurat object contains less than 10,000 cells. Results are fed into 'summarizeSweep' and 'find.pK' functions during optimal pK parameter selection workflow. Parameters tested: pN = 0.05-0.3, pK = 0.0005-0.3.
}
\usage{
sweep.list <- doubletFinder_ParamSweep(seu)
}
\arguments{
  \item{seu}{ A fully-processed Seurat object (i.e., After NormalizeData, FindVariableGenes, ScaleData, RunPCA, and RunTSNE have all been performed).
  }
}
\details{
}
\value{
List of pANN vectors for every pN and pK combination. Output also contains pANN information for artificial doublets.
}
\references{
}
\author{
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
}
\examples{
sweep.list <- doubletFinder_ParamSweep(seu)
sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
}
