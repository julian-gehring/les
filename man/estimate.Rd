\name{estimate}

\alias{estimate}

\title{
  estimate
}

\description{
  Estimates the Loci of Enhanced Signifiance Index \eqn{\Lambda}{Lambda} for a
  given window.
}

\usage{
  estimate(object, win, weighting = triangWeight, grenander = FALSE, se
  = TRUE, nCores=NULL)
}

\arguments{
  \item{object}{Object of class 'Les' containing experimental data, as
    returned by 'create'.}
  \item{win}{Integer specifing window size.}
  \item{weighting}{Function specifing the shape of the weighting
    window. If not specified the supplied 'triangWeight' function with a
    triangular window will be used. For details on other window
    functions and how to specify own function please see below.}
  \item{grenander}{Logical indicating whether Grenander correction
    should be used (default: FALSE).}
  \item{se}{Logical indicating whether the standard error from the final
    linear model should be computed (default: TRUE). The standard error
    displays the goodness of fit for every probe, but is not needed for
    further computation. If memory or computation time are a critical
    factor computation of the SE can be omitted.}
  \item{nCores}{Logical or numeric indicating the number of cores to use
    for computation. This feature requires the 'multicore' package which
    is only available for certain platforms. The value is passed
    directly to 'mclapply' as argument 'n.cores'. If FALSE the package
    will not be used. If an integer it indicates the number of cores to
    use. For details and benefits please see below.}
}

\details{
  This function estimates \eqn{\Lambda_i}{Lambda(i)} for all probes \eqn{i}{i}.
}
      
\value{
  Object of class 'Les' with additionally filled slots:
  lambda, win, weighting, grenander, nProbes, se (if computed)
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
  \code{\link{cutoff}}
  \code{\link{regions}}
  \code{\link{ci}}
}

\examples{
  \dontrun{
    pos <- 1:100
    pval <- runif(100)
    chr <- rep(1:2, each=50)

    res <- create(pos, pval, chr)

    res <- estimate(res, 20)
  }
}

\keyword{htest}

