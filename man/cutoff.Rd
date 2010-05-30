\name{cutoff}

\alias{cutoff}

\title{cutoff}

\description{
  Estimates cutoff for regions
}

\usage{
  cutoff(object, grenander = FALSE, verbose = FALSE)
}

\arguments{
  \item{object}{Object of class 'Les' as returned by 'estimate'.}
  \item{grenander}{Logical indicating whether Grenander correction
    should be used (default: FALSE). Note that this is independent from
    the choice in 'estimate'.}
  \item{verbose}{Logical indicating whether a summary of the estimated
    number of regulated probes should be printed on screen (default: FALSE).}
}

\details{
  This method estimates the number of probes with a significant
  effect R. The estimation is based on the p-value distribution.

  A reasonable estimate for the cutoff value \eqn{\Theta}{Theta} can be
  chosen such that \eqn{\mid\Lambda_i\geq\Theta\mid=R}{|Lambda >= Theta|
  = R}.
}

\value{
  Object of class 'Les' with additionally filled slots:
  nSigProbes, theta	
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
  \code{\link{regions}}
  \code{\link{ci}}
}

\examples{
  \dontrun{
    res <- cutoff(res, verbose=TRUE)
  }
}

\keyword{htest}
