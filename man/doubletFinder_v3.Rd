\name{doubletFinder_v3}
\alias{doubletFinder_V3}
\title{doubletFinder_V3}
\description{
Core doublet prediction function of the DoubletFinder package. Generates artifical doublets from an existing, pre-processed Seurat object. Real and artificial data are then merged and pre-processed using parameters utilized for the existing Seurat object. PC distance matrix is then computed and used the measure the proportion of artificial nearest neighbors (pANN) for every real cell. pANN is then thresholded according to the number of expected doublets to generate final doublet predictions.
}
\usage{
doubletFinder_V3(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE)
}
\arguments{
  \item{seu}{ A fully-processed Seurat object (i.e., After NormalizeData, FindVariableGenes, ScaleData, and RunPCA have all been performed).
  }
  \item{PCs}{ Number of statistically-significant principal components (e.g., as estimated from PC elbow plot)
  }
  \item{pN}{ The number of generated artificial doublets, expressed as a proportion of the merged real-artificial data. Default is set to 0.25, based on observation that DoubletFinder performance is largely pN-invariant (see McGinnis, Murrow and Gartner 2019, Cell Systems).
  }
  \item{pK}{ The PC neighborhood size used to compute pANN, expressed as a proportion of the merged real-artificial data. No default is set, as pK should be adjusted for each scRNA-seq dataset. Optimal pK values can be determined using mean-variance-normalized bimodality coefficient.
  }
  \item{nExp}{ The total number of doublet predictions produced. This value can best be estimated from cell loading densities into the 10X/Drop-Seq device, and adjusted according to the estimated proportion of homotypic doublets.
  }
  \item{reuse.pANN}{ Seurat metadata column name for previously-generated pANN results. Argument should be set to FALSE (default) for initial DoubletFinder runs. Enables fast adjusting of doublet predictions for different nExp.
  }
  \item{sct}{ Logical representing whether SCTransform was used during original Seurat object pre-processing (default = FALSE).
  }
}
\details{
}
\value{ Seurat object with updated metadata including pANN and doublet classifications.
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
## Initial run, nExp set to 0.15 Poisson loading estimate (e.g., 1000 total doublet predictions)
nExp_poi <- round(0.15*nrow(seu@meta.data))
seu <- doubletFinder_v3(seu, PCs = 1:10, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct=FALSE)

## With homotypic adjustment
homotypic.prop <- modelHomotypic(annotations)
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seu <- doubletFinder_v3(seu, PCs = 1:10, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_1000", sct=FALSE)
}
