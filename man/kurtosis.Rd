\name{kurtosis}
\alias{kurtosis}
\title{kurtosis}
\description{
Internal function to compute bimodality coefficient during BCmvn computation and pK estimation.
}
\usage{
kurtosis(x)
}
\arguments{
  \item{x}{ Guassian kernel density estimation representing pANN distribution }
}
\details{
}
\value{ Kurtosis value
}
\references{
Taken from the 'modes' R package (v0.7).
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
## Internal to bimodality_coefficient
sample.excess.kurtosis <- kurtosis(x)
}
