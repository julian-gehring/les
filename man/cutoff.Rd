\name{cutoff}
\alias{cutoff}
\title{cutoff
}
\description{
	Estimates cutoff for regions
}
\usage{
cutoff(res, ...)
}

\arguments{
  \item{res}{Object of class 'les'
}
  \item{\dots}{...
}
}
\details{
	Estimates cutoff for regions
}
\value{
	Object of class 'Les'
}

\author{
Julian Gehring, Clemens Kreutz

Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}
\note{
	... notes ...
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

\keyword{...}
