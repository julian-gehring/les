\name{les-internal}
\title{Internal les functions}
\alias{les-internal}
\alias{bootSingle}
\alias{calcSingle}
\alias{fitGSRI}
\alias{fitGSRIboot}
\alias{getCDF}
\alias{mcsapply}

\description{Internal les functions}

\usage{
bootSingle(ind0, pos, pval, win, weighting, nBoot, conf, grenander)
calcSingle(ind0, pos, pval, win, weighting, grenander)
fitGSRI(pval, weight, pos0, nValidProbes, grenander)
fitGSRIboot(pval, index, cweight, nValidProbes, grenander)
getCDF(pval, weight, grenander)
mcsapply(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, mc.cores)
}

\details{
	Internal functions of the 'les' package.
	These are not to be called by the user.
}

\references{
	Modified code from the following packages is used:
	...
}
\author{
	Julian Gehring, Clemens Kreutz
}

\seealso{
	\code{\link{les}}
}
\keyword{internal}

