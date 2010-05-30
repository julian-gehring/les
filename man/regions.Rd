\name{regions}

\alias{regions}

\title{
  regions
}

\description{
  Estimates regions of continously elevated \eqn{\Lambda_i}{Lambda}.
}
      
\usage{
  regions(object, limit, minLength = 10, maxGap = 100, verbose = FALSE)
}

\arguments{
  \item{object}{Object of class 'Les' as returned by 'cutoff'}
  \item{limit}{Numeric specifing the threshold \eqn{\Theta}{Theta} for
    \eqn{\Lambda}{Lambda} with regions are defined as
    \eqn{\Lambda_i\geq\Theta}{Lambda>=Theta}. If missing
    \eqn{\Theta}{Theta} will be derived from the estimated number of
    regulated probes as computed by 'cutoff'.}
  \item{minLength}{Integer specifing the minimum number of probes in a
    region to be reported as such.}
  \item{maxGap}{Integer specifing maximum gap in base pairs between two
    neighboring probes in a region. If this value is exceeded the region
    will be split into smaller ones such that maxGap is not exceed by
    any probe pair.}
  \item{verbose}{Logical indicating whether a summary of the estimated
    regions should be printed on screen (default: FALSE).}
}

\details{
  Finds regions above threshold
}

\value{
  Object of class 'Les' with additionally filled slots:
  regions, limit
}

\author{
  Julian Gehring, Clemens Kreutz

  Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}

\seealso{
  Class:
  \code{\linkS4class{Les}}
  
  Functions:
  \code{\link{estimate}}
  \code{\link{cutoff}}
  \code{\link{ci}}
}

\examples{
  \dontrun{
    res <- regions(res, verbose=TRUE)
  }
}

\keyword{htest}
