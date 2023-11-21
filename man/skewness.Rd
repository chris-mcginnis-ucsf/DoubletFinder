\name{skewness}
\alias{skewness}
\title{skewness}
\description{
Internal function to compute skewness during BCmvn computation and pK estimation.
}
\usage{
skewness(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Guassian kernel density estimation representing pANN distribution }
}
\details{
}
\value{ Skewness value
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
G <- skewness(x)
}
