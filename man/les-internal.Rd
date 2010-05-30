\name{les-internal}

\alias{les-internal}
\alias{calcSingle}
\alias{fitGSRI}
\alias{wcdf}
\alias{mcsapply}
\alias{seFast}
\alias{xvalWeight}

\title{
  Internal les functions
}

\description{Internal les functions}

\usage{
  calcSingle(ind0, pos, pval, win, weighting, grenander, se, nBoot, conf)
  fitGSRI(pval, index, cweight, nValidProbes, grenander, se)
  wcdf(pval, weight, grenander)
  mcsapply(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, mc.cores)
  seFast(x, y, b)
  xvalWeight(distance, win)
  log2ind(log)
  ind2log(ind, n)
  reg2log(reg, x)
  reg2ind(reg, x)
  log2reg(log)
  ind2reg(ind, n)
  scaling(x)
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

