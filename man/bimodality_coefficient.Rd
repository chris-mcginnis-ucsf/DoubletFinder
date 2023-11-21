\name{bimodality_coefficient}
\alias{bimodality_coefficient}
\title{bimodality_coefficient}
\description{
Internal function to compute bimodality coefficient during BCmvn computation and pK estimation.
}
\usage{
bimodality_coefficient(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Guassian kernel density estimation representing pANN distribution }
}
\details{
}
\value{ Bimodality coefficient value
}
\references{
Taken from the 'modes' R package (v0.7)
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
## Internal to summarizeSweep
gkde <- approxfun(bkde(res.temp$pANN, kernel="normal"))
x <- seq(from=min(res.temp$pANN), to=max(res.temp$pANN), length.out=nrow(res.temp))
sweep.stats$BCreal[i] <- bimodality_coefficient(gkde(x))
}
