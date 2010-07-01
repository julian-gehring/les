\name{ci}

\alias{ci}

\title{ci}

\description{
  Compute confidence intervals for \eqn{\Lambda}{Lambda} with bootstrapping
}

\usage{
  ci(object, subset, nBoot = 100, conf = 0.95, nCores = NULL)
}

\arguments{
  \item{object}{Object of class 'Les' as returned by 'estimate'.}
  \item{subset}{Vector of logical specifing the probes for which the CIs
    should be computed. If missing CIs will be computed for all probes.}
  \item{nBoot}{Integer specifing the number of bootstrap samples
    (default: 100). For details see 'boot' from the 'boot' package.}
  \item{conf}{Numeric specifing the confidence level (default:
    0.95). For details see 'boot' from the 'boot' package.}
  \item{nCores}{Logical or numeric indicating the number of cores to use
    for computation. This feature requires the 'multicore' package which
    is only available for certain platforms. The value is passed
    directly to 'mclapply' as argument 'n.cores'. If FALSE the package
    will not be used. If an integer it indicates the number of cores to
    use. For details and benefits please see below.}
}

\details{
  Computes confidence intervals by bootstrapping probes in each
  window. Since based on percentiles the resulting CIs are
  unsymmetrical.

  All arguments for computation are taken from 'object' and thereby kept
  the same as for the results from 'estimation'.

  Since bootstrapping is computational demanding and CIs are often only
  wanted for certain regions of interest it may be useful to restrict
  computation with the 'subset' argument. 
}

\value{
  Object of class 'Les' with additionally filled slots:
  ci, subset, nBoot, conf
}

\author{
  Julian Gehring, Clemens Kreutz

  Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
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
    res <- estimate(res, 20)
    res <- ci(res, subset=1:100, conf=0.95)
  }
}

\keyword{htest}
