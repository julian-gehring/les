\name{weighting}

\alias{weighting}
\alias{triangWeight}
\alias{rectangWeight}
\alias{gaussWeight}

\title{
  Weighting functions
}

\description{
  Set of functions to compute spatial weights between probes.
}
      
\usage{
  triangWeight(distance, win)
  rectangWeight(distance, win)
  gaussWeight(distance, win) 
}

\arguments{
  \item{distance}{Numeric vector specifing distance of probes from central probe.}
  \item{win}{Integer specifing maximum size of window.}
}

\details{
  The functions 'triangWeight', 'rectangWeight' and 'gaussWeight'
  provide a triangular and rectangular weighting window
  respectively. The weighting function can be specified by the
  'weightingFunction' argument in the 'estimate' method.

  This way it is also possible to use custom weighting functions. In
  general they have to be called the same way as the functions mentioned
  before and have to return a vector of weights of the same length as
  the argument 'distance'. For more details on how to use own weighting
  functions please refer to the vignette of this package.
  
  Please note that the returned weights do not have to be normalized to
  sum up to 1 since this is done at the computation of the weighted CDF.
}
      
\value{
  A numeric vector with weights for each probe in the window.
}

\author{
  Julian Gehring, Clemens Kreutz

  Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}

\seealso{
  Class:
  \code{\linkS4class{Les}}
  
  Methods and functions:
  \code{\link{create}}
  \code{\link{estimate}}
  \code{\link{threshold}}
  \code{\link{regions}}
  \code{\link{ci}}
  \code{\link{exportLambda}}
  \code{\link{weighting}}
}

\examples{
\dontrun{
weight <- triangWeight(dist, win)
weight <- rectangWeight(dist, win)
weight <- gaussWeight(dist, win)

res <- estimate(res, win=20, weighting=rectangWeight)
}
}

\keyword{htest}

