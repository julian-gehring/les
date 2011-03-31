\name{les-internal}

\alias{les-internal}
\alias{calcSingle}
\alias{fitGsri}
\alias{wcdf2}
\alias{mcsapply}
\alias{seFast}
\alias{log2ind}
\alias{ind2log}
\alias{slopeWeight}
\alias{qrSlope}
\alias{itLinReg}
\alias{cdfDuplicates}
\alias{seSlopeWeight}
\alias{setState}
\alias{checkState}
\alias{modGrenander}

\title{
  Internal functions of the 'les' package.
}

\description{Internal les functions}

\usage{
calcSingle(ind0, pos, pval, win, weighting, grenander, se, nBoot, conf,
minProbes, custom)
fitGsri(pval, index, cweight, nValidProbes, grenander, se, custom)
itLinReg(x, y, cweight, nValidProbes, se, custom, noBoot)
gsri(pval, grenander=FALSE, se=TRUE, custom=FALSE)
wcdf2(pval, weight)
mcsapply(X, FUN, ..., mc.cores=NULL, simplify=TRUE)
seFast(x, y, b)
log2ind(log)
ind2log(ind, n)
slopeWeight(x, y, c)
diagSquare(w, n)
qrSlope(x, y, w)
cdfDuplicates(pvalSort, weightSort)
seSlopeWeight(fit, n=length(fit$weights))
setState(state, flag)
checkState(state, flag, text=as.character(flag))
modGrenander(x, y, eta0, unique, check)
getArgs(name, first=NULL, last=NULL, ...)
}

\details{
  Internal functions of the 'les' package.
  These should not be called by the user.
}

\author{
  Julian Gehring, Clemens Kreutz
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

\keyword{internal}

