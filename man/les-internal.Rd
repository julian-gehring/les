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
\alias{cdfCorrect}
\alias{cdfDuplicates}
\alias{grenanderPass}

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
wcdf2(pval, weight, grenander=FALSE)
mcsapply(X, FUN, ..., mc.cores=NULL)
seFast(x, y, b)
log2ind(log)
ind2log(ind, n)
slopeWeight(x, y, c)
diagSquare(w, n)
qrSlope(x, y, w)
cdfCorrect(x, y, q0)
cdfDuplicates(pvalSort, weightSort)
grenanderPass(x, y, un)
}

\details{
  Internal functions of the 'les' package.
  These are not to be called by the user.
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
  \code{\link[les]{create}}
  \code{\link[les]{estimate}}
  \code{\link[les]{threshold}}
  \code{\link[les]{regions}}
  \code{\link[les]{ci}}
  \code{\link[les]{export}}
  \code{\link[les]{plot}}
}

\keyword{internal}

