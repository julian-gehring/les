\name{ci}
\alias{ci}
\title{ci}
\description{
	Compute confidence intervals with bootstrapping
}
\usage{
	ci(res, ...)
}

\arguments{
  \item{res}{Object of class 'Les'
}
  \item{\dots}{...
}
}
\details{
	Compute confidence intervals with bootstrapping
}
\value{
	Object of class 'Les'
}

\author{
Julian Gehring, Clemens Kreutz

Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}
\note{
	... some notes ...
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
	res <- ci(res, subset=1:100, ci=0.95)
}
}

\keyword{...}
