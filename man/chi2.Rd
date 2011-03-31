\name{chi2,Les-method}

\alias{chi2}
\alias{chi2-methods}
\alias{chi2,Les-method}

\title{chi2}

\description{
  The 'chi2' method can be used to optimize the window size for defined
  regions of interest. It computes \eqn{\chi^2(w)}{chi^2(w)} for window
  sizes \eqn{w}{w} based on the estimates of \eqn{\Lambda}{Lambda} and
  the false-discovery rate.
}

\usage{
chi2(object, winSize, regions, offset, fdr = "lfdr", method,
scaling = les:::scaleNorm, nCores = NULL,  verbose = FALSE, ...)

\S4method{chi2}{Les}(object, winSize, regions, offset, fdr = "lfdr",
method, scaling = les:::scaleNorm, nCores = NULL,  verbose = FALSE, ...)
}

\arguments{
  \item{object}{Object of class 'Les' as returned by 'estimate' or a
    subsequent step.}
  \item{winSize}{Integer vector specifying the window sizes \eqn{w}{w}
    for which \eqn{\chi^2(w)}{chi^2(w)} should be computed. For each value of
    'winSize' and each region a computation as for 'estimate' is
    run. For details please see the 'win' argument in 'estimate'.}
  \item{regions}{Data frame containing the regions of interest. It has to
    contain the columns 'start', 'end' and 'chr' specifying the start
    and end position for each region, as well the chromosome if more than
    one chromosome is present. The structure is related to the 'regions'
    output of the 'regions' method. If missing the data frame from the
    'regions' slot will be used if available. For details please see the
    'regions' method.}
  \item{offset}{Integer or vector of integers specifying the offset for
    the regions given by the 'regions' input argument. If missing the
    regions will be taken as specified. If present start and end of each
    regions will be taken as 'start - offset' and 'end + offset'.}
  \item{fdr}{Character string specifying the fdr method to use for
    \eqn{\chi^2}{chi^2} computation (default: 'lfdr'). Possible values
    are 'lfdr' for local fdr and 'qval' for q-values. For details see
    the 'Details' section below and the 'fdrtool' package.}
  \item{method}{Character string specifying the method used for linear
    regression. It is equivalent to the 'method' argument in the 'estimate'
    method. If missing the value set in the 'estimate' method will be
    used.}
  \item{scaling}{Function specifying the scaling of Lambda and fdr
    (default: les:::scaleNorm). By default both will be scaled to the
    range [0,1].}
  \item{nCores}{Integer indicating the number of cores to use for
    computation. This feature requires the 'multicore' package which
    is only available for certain platforms. The package is used only if
    'library(multicore)' has been called manually by the user before and if
    'nCores' is an integer unequal NULL specifying the number of cores
    to use. The value is passed directly to 'mclapply' as argument
    'n.cores'. For details and benefits please see the 'Details'
    section.}
  \item{verbose}{Logical indicating whether the progress of the
    computation should be printed on screen (default: FALSE).}
  \item{...}{Further arguments passed to subsequent functions.}
}

\details{
  The 'chi2' method can be used to optimize the window size for defined
  regions of interest. It computes \eqn{\chi^2(w)}{chi^2(w)} for each
  window size \eqn{w}{w} based on the estimates of false-discovery rate (fdr) and
  \eqn{\Lambda(w)}{Lambda(w)} with a Leave-One-Out Cross Validation
  (LOOCV). The shape of the \eqn{\chi^2(w)}{chi^2(w)} landscape can
  constrain suitable values for \eqn{w}{w}.
}
      
\value{
  Object of class 'Les' with additionally filled slots:
  winSize, chi2
}

\author{
  Julian Gehring

  Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}

\seealso{
  Package:
  \code{\link[les]{les-package}}
  \code{\link[fdrtool]{fdrtool}}

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
  \code{\link[les]{weighting}}
}

\examples{
\dontrun{
data(statTile)

x <- Les(pos, pval, chr)
x <- estimate(x, 50)
x <- regions(x, limit=0.5)

regions <- x["regions"]
winsize <- seq(10, 100, by=10)
x <- chi2(x, winsize, regions, offset=250, verbose=TRUE)
chi2 <- x["chi2"]
}
}

\keyword{htest}
\keyword{methods}
