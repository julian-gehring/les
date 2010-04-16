\name{weighting}
\alias{weighting}
\alias{triangWeight}
\alias{rectangWeight}

\title{Weighting functions
}
\description{
	Weighting functions
}
\usage{
	triangWeight(distance, win)
	rectangWeight(distance, win)
}

\arguments{
  \item{distance}{
	Numeric vector specifing distance of probes from center probe.}
  \item{win}{Integer specifing maximum size of window.}
}
\details{
	Weighting functions
}
\value{
	Numeric vector with weights for each probe in window.
}

\author{
Julian Gehring, Clemens Kreutz

Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}
\note{
	...
}

\seealso{
	Class:
	\code{\linkS4class{Les}}
	
	Functions:
	\code{\link{estimate}}
	\code{\link{cutoff}}
	\code{\link{regions}}
	\code{\link{ci}}
}

\examples{
\dontrun{
	weight <- triangWeight(dist, win)
	weight <- rectangWeight(dist, win)
}
}

\keyword{...}

