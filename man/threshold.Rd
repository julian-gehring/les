\name{threshold,Les-method}

\alias{threshold}
\alias{threshold-methods}
\alias{threshold,Les}
\alias{threshold,Les-method}

\title{threshold}

\description{
  The 'threshold' method estimates a suitable threshold
  \eqn{\Theta}{Theta} from the data. Thresholding provides the ability
  the define distinct Loci of Enhanced Significance (LES) with the
  'regions' method in a later step.
}

\usage{
threshold(object, grenander = FALSE, verbose = FALSE, ...)

\S4method{threshold}{Les}(object, grenander = FALSE, verbose = FALSE,
...)
}

\arguments{
  \item{object}{Object of class 'Les' as returned by 'estimate' or
    'Les'.}
  \item{grenander}{Logical indicating whether the Grenander estimator
    for the cumulative density should be used (default: FALSE). For
    details see below and at the 'GSRI' package.}
  \item{verbose}{Logical indicating whether a summary of the estimated
    number of regulated probes should be printed on screen (default:
    FALSE).}
  \item{...}{Further arguments passed to subsequent functions.}
}

\details{
  This method estimates the number of probes with a significant
  effect \eqn{R}{R}. The estimation is based on the p-value
  distribution. The analysis is based on the 'Gene Set Regulation Index'
  by Bartholome et al., 2009.

  Estimation of the threshold is independent of the computation
  performed by the 'estimate' method.

  A reasonable estimate for the cutoff value \eqn{\Theta}{Theta} can be
  chosen such that \eqn{\mid\Lambda_i\geq\Theta\mid=R}{|Lambda >= Theta|
    = R}.

  The Grenander estimator for the cumulative density results in a
  conservative estimate for the number of significant probes with
  decreased variance.

  A reasonable subsequent step is to call 'regions' to find distinct
  Loci of Enhanced Significance.
}

\value{
  Object of class 'Les' with additionally filled slots:
  \item{nSigProbes}{Estimated number of significant probes.}
  \item{theta}{Estimated threshold \eqn{\Theta}{Theta}.}
}

\author{
  Julian Gehring

  Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}

\references{
  Kilian Bartholome, Clemens Kreutz, and Jens Timmer:
  Estimation of gene induction enables a relevance-based ranking of gene sets,
  Journal of Computational Biology: A Journal of Computational Molecular Cell 
  Biology 16, no. 7 (July 2009): 959-967.
  \url{http://www.liebertonline.com/doi/abs/10.1089/cmb.2008.0226}
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
  \code{\link[les]{export}}
  \code{\link[les]{plot}}
}

\examples{
data(statTile)

x <- Les(pos, pval, chr)
x <- estimate(x, 50)
x <- threshold(x, verbose=TRUE)
}

\keyword{htest}
\keyword{methods}
