\name{estimate}
\alias{estimate}
\title{
estimate
}
\description{
	Estimates LESS.
}
\usage{
	estimate(pos, pval, win, weighting = triangWeight, grenander = FALSE, nCores=FALSE)
}

\arguments{
  \item{pos}{Numeric vector containing probe positions.}
  \item{pval}{Numeric vector containing probe p-values. Same length as probe.}
  \item{win}{Integer specifing window size.}
  \item{weighting}{Function specifing weighting (default: triangWeight).}
  \item{grenander}{Logical indicating whether Greander estimate should be used (default: FALSE).}
  \item{nCores}{Logical or numeric indicating the number of cores to use. Requires the 'multicore' package. The value is passed to 'mclapply' as 'n.cores'. If FALSE the package will not be used. If integer it indicates the number of cores to use.}
}
\details{
	Estimates LESS.
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
	\code{\link{regions}}
	\code{\link{ci}}
}

\examples{
\dontrun{
	res <- estimate(pos, pval, win)
}
}

\keyword{...}

