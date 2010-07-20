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
    'create'.}
  \item{grenander}{Logical indicating whether the Grenander correction
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

  The Grenander correction for the cumulative density introduces a bias
  but decreases the variability of the estimate. This may be reasonable
  in some cases.

  A reasonable subsequent step is to call 'regions' to find distinct
  regions Loci of Enhanced Significance (LES).
}

\value{
  Object of class 'Les' with additionally filled slots:
  \item{nSigProbes}{Estimated number of significant probes.}
  \item{theta}{Estimated threshold \eqn{\Theta}{Theta}.}
}

\author{
  Julian Gehring, Clemens Kreutz

  Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
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

\examples{
\dontrun{
data(statTile)

x <- create(pos, pval, chr)
x <- estimate(x, win=50)
x <- threshold(x, verbose=TRUE)
}
}

\keyword{htest}
\keyword{methods}