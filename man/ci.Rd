\name{ci}

\alias{ci}

\title{ci}

\description{
  Computes confidence intervals for \eqn{\Lambda}{Lambda} with bootstrapping.
}

\usage{
ci(object, subset, nBoot = 100, conf = 0.95, nCores = NULL)
}

\arguments{
  \item{object}{Object of class 'Les' as returned by 'estimate' or
    subsequent methods.}
  \item{subset}{Vector of logical specifying the probes for which the CIs
    should be computed. If missing CIs will be computed for all probes.}
  \item{nBoot}{Integer specifying the number of bootstrap samples
    (default: 100). For details see 'boot' from the 'boot' package.}
  \item{conf}{Numeric specifying the confidence level (default:
    0.95). For details see 'boot' from the 'boot' package.}
  \item{nCores}{Integer indicating the number of cores to use for
    computation. This feature requires the 'multicore' package which
    is only available for certain platforms. The package is used only if
    'library(multicore)' has been called manually by the user before and if
    'nCores' is an integer unequal NULL specifying the number of cores
    to use. The value is passed directly to 'mclapply' as argument
    'n.cores'. For details and benefits please see below.}
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

  The 'multicore' package can be used to spread the computation over several
  cores in a simple way. This can be useful on multi-core machines for
  large datasets and intensive CI computations. The 'multicore' package
  is not available on all platforms. To use multicore processing
  'library(multicore)' has to be called beforehand and a number of cores
  to use has to be specified in 'nCores'. For details see the
  documentation of the 'multicore' package.
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
  
  Methods and functions:
  \code{\link{create}}
  \code{\link{estimate}}
  \code{\link{threshold}}
  \code{\link{regions}}
  \code{\link{ci}}
  \code{\link{exportLambda}}
}

\examples{
\dontrun{
data(resTile)

x <- create(pos, pval, chr)
x <- estimate(x, win=25)

x <- ci(x, subset=1:100, conf=0.95)

## multicore computation
## only available of certain platforms
library(multicore)
x <- ci(x, subset=1:100, nCores=2)
}
}

\keyword{htest}
