\name{weighting}

\alias{weighting}
\alias{triangWeight}
\alias{rectangWeight}
\alias{gaussWeight}
\alias{epWeight}

\title{Weighting functions}

\description{
  Set of functions to compute spatial weights between probes.
}
      
\usage{
triangWeight(distance, win)
rectangWeight(distance, win)
gaussWeight(distance, win)
epWeight(distance, win)
}

\arguments{
  \item{distance}{Numeric vector specifying the distance of probes from
    the central probe. Negative values refer to probes upstream, positive
    values to probes downstream.}
  \item{win}{Integer specifying maximum size of window.}
}

\details{
  The functions 'triangWeight', 'rectangWeight', 'epWeight' and
  'gaussWeight' provide a triangular, rectangular, Epanechnikov and
  Gaussian weighting window, respectively. The weighting function can be 
  specified by the 'weightingFunction' argument in the 'estimate'
  method.

  This way it is also possible to use custom weighting functions. In
  general they have to be called the same way as the functions mentioned
  before and have to return a vector of weights of the same length as
  the argument 'distance'. For more details on how to use own weighting
  functions please refer to the vignette of this package.
  
  Please note that the returned weights do not have to be normalized
  since this is done at the computation of the weighted cumulative
  density.
}
      
\value{
  A numeric vector with weights for each probe in the window.
}

\author{
  Julian Gehring

  Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}

\seealso{
  Package:
  \code{\link[les]{les-package}}

  Class:
  \code{\linkS4class{Les}}
  
  Methods and functions:
  \code{\link[les]{Les}}
  \code{\link[les]{estimate}}
  \code{\link[les]{threshold}}
  \code{\link[les]{regions}}
  \code{\link[les]{ci}}
  \code{\link[les]{chi2}}
  \code{\link[les]{export}}
  \code{\link[les]{plot}}
}

\examples{
distance <- seq(-50, 50)
win <- 50

weight <- triangWeight(distance, win)
plot(distance, weight, type="l", main="triangWeight")

weight <- rectangWeight(distance, win)
plot(distance, weight, type="l", main="rectangWeight")

weight <- gaussWeight(distance, win)
plot(distance, weight, type="l", main="gaussWeight")

weight <- epWeight(distance, win)
plot(distance, weight, type="l", main="epWeight")

## simple example for a custom weighting function
ownWeighting <- function(distance, win)  {
    weight <- as.integer(abs(distance) < win)
    return(weight)
}
}

\keyword{htest}
\keyword{smooth}

