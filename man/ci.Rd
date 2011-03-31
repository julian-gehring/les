\name{ci,Les-method}

\alias{ci}
\alias{ci-methods}
\alias{ci,Les}
\alias{ci,Les-method}

\title{ci}

\description{
  Computes confidence intervals (CI) for \eqn{\Lambda}{Lambda} using
  bootstrapping.
}

\usage{
ci(object, subset, nBoot = 100, conf = 0.95, nCores = NULL, ...)

\S4method{ci}{Les}(object, subset, nBoot = 100, conf = 0.95, nCores =
NULL, ...)
}

\arguments{
  \item{object}{Object of class 'Les' as returned by 'estimate' or
    'regions'.}
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
    'n.cores'. For details and benefits please see the 'Details'
    section.}
  \item{...}{Further arguments passed to subsequent functions.}
}

\details{
  The 'ci' method computes confidence intervals (CI) by bootstrapping
  probes in each window. Since based on percentiles the resulting CIs
  are asymmetrical.

  All arguments for computation are taken from 'object' and thereby kept
  the same as for the results from 'estimation'.

  Since bootstrapping is computational demanding and CIs are often only
  wanted for certain regions of interest it may be useful to restrict
  computation with the 'subset' argument.

  The 'multicore' package can be used to spread the computation over several
  cores in a simple way. This can be useful on multi-core machines for
  large datasets or computation of confidence intervals for many
  probes. The 'multicore' package is not available on all platforms. To
  use multicore processing 'library(multicore)' has to be called
  beforehand and a number of cores to use has to be specified in
  'nCores'. For details see the documentation of the 'multicore' package.
}

\value{
  Object of class 'Les' with additionally filled slots:
  ci, subset, nBoot, conf
}

\author{
  Julian Gehring

  Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}

\seealso{
  Package:
  \code{\link[les]{les-package}}
  \code{\link[boot]{boot}}

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

\examples{
data(statTile)

x <- Les(pos, pval, chr)
x <- estimate(x, win=50, grenander=FALSE)

x <- ci(x, subset=150:200, conf=0.90, nBoot=50)

plot(x, chr=1, error="ci")

\dontrun{
## multicore computation
## only available on certain platforms
library(multicore)
x <- ci(x, subset=150:200, conf=0.90, nBoot=50, nCores=2)
}
}

\keyword{htest}
\keyword{methods}
