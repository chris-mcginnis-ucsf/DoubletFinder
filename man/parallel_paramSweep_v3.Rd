\name{parallel_paramSweep_v3}
\alias{parallel_paramSweep_v3}
\title{parallel_paramSweep_v3}
\description{
Internal parallelization function for paramSweep_v3.
}
\usage{
**NOT RUN**
parallel_paramSweep_v3(n, n.real.cells, real.cells, pK, pN, data, orig.commands, PCs, sct)
}
\arguments{
  \item{n}{ pN iteration counter.
  }
  \item{n.real.cells}{ Number of real cells. Set automatically during paramSweep_v3.
  }
  \item{real.cells}{ Vector of real cell IDs. Set automatically during paramSweep_v3.
  }
  \item{pK}{ The PC neighborhood size used to compute pANN, expressed as a proportion of the merged real-artificial data. No default is set, as pK should be adjusted for each scRNA-seq dataset. Optimal pK values can be determined using mean-variance-normalized bimodality coefficient.
  }
  \item{pN}{ The number of generated artificial doublets, expressed as a proportion of the merged real-artificial data. Default is set to 0.25, based on observation that DoubletFinder performance is largely pN-invariant (see McGinnis, Murrow and Gartner 2019, Cell Systems).
  }
  \item{data}{ Count matrix. Set automatically during paramSweep_v3.
  }
  \item{orig.commands}{ Count matrix. Set automatically during paramSweep_v3.
  }
  \item{PCs}{ Number of statistically-sigificant PCs. Set according to paramSweep_v3 arguments.
  }
  \item{sct}{ Logical representing whether Seurat object was pre-processed using 'sctransform'. Set according to paramSweep_v3 arguments (default = F).
  }
}
\details{
}
\value{ Parallelization function compatible with mclapply.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{ Implemented by Nathan Skeene, June 2019. }
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
