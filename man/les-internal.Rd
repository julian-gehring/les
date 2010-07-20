\name{les-internal}

\alias{les-internal}
\alias{calcSingle}
\alias{fitGSRI}
\alias{wcdf}
\alias{mcsapply}
\alias{seFast}

\title{
  Internal functions of the 'les' package.
}

\description{Internal les functions}

\usage{
calcSingle(ind0, pos, pval, win, weighting, grenander, se, nBoot, conf, minProbes)
fitGSRI(pval, index, cweight, nValidProbes, grenander, se)
gsri(pval, grenander, se)
wcdf2(pval, weight, grenander=FALSE)
mcsapply(X, FUN, ..., mc.cores=NULL)
seFast(x, y, b)
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

