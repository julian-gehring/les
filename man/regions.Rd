\name{regions}
\alias{regions}
\title{
regions
}
\description{
	Finds regions above threshold
}
\usage{
	regions(res, ...)
}

\arguments{
  \item{res}{Object of class 'Les'
}
  \item{\dots}{
	...
}
}
\details{
	Finds regions above threshold
}
\value{
	Object of class 'Les'
}

\author{
Julian Gehring, Clemens Kreutz

Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}
\note{
	... notes
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

\keyword{...}
