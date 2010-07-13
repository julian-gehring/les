\name{estimate}

\alias{estimate}

\title{estimate}

\description{
  The 'estimate' method estimates the index for Loci of Enhanced
  Significance \eqn{\Lambda}{Lambda} with a given weighting window.
}

\usage{
estimate(object, win, weighting = triangWeight, grenander = FALSE, se =
FALSE, nCores=NULL)

}

\arguments{
  \item{object}{Object of class 'Les' containing experimental data, as
    returned by 'create'.}
  \item{win}{Integer specifying window size for the weighting
    function. This value is directely passed to the function specified
    by 'weighting'. For details see the description for the used window
    function 'weighting'.}
  \item{weighting}{Function specifying the shape of the weighting
    window. If not specified the supplied 'triangWeight' function with a
    triangular window will be used. For details on other window
    functions and how to specify own function please see below.}
  \item{grenander}{Logical indicating whether the Grenander correction
    for the cumulative density should be used (default: FALSE). For
    details see below and at the 'GSRI' package.}
  \item{se}{Logical indicating whether the standard error from the final
    linear model should be computed and stored (default: FALSE). The
    standard error displays the goodness of fit for every probe, but is
    not needed for further computation. If computation time is a
    critical factor computation of the SE can be omitted to save some
    time.}
  \item{nCores}{Integer indicating the number of cores to use for
    computation. This feature requires the 'multicore' package which
    is only available for certain platforms. The package is used only if
    'library(multicore)' has been called manually by the user before and if
    'nCores' is an integer unequal NULL specifying the number of cores
    to use. The value is passed directly to 'mclapply' as argument
    'n.cores'. For details and benefits please see below.}
}

\details{
  This function estimates \eqn{\Lambda_i}{Lambda(i)} for all probes
  \eqn{i}{i}. This is normally the first step in the analysis after
  storing the experimental data with 'create'.

  The 'win' argument influences the amount of neighboring probes taken
  into account by the weighting function. The value is passed to the
  function specified in 'weighting'. Larger values result in a smoother
  features. Details on a reasonable choise for this value can be found
  in the references.

  With the 'weighting' argument the applied weighting function can be
  specified from a predifined set or a custom function can be used. In
  the 'les' package the three functions 'triangWeight', 'rectangWeight'
  and 'gaussWeight' are already supplied offering a triangular,
  rectangular or gaussian weighting function respectively. The
  'triangWeight' is in general a good choise. For details on the
  functions itself and how to use custom ones please see the
  documentation of the single functions or the vignette.

  The Grenander correction for the cumulative densityintroduces a bias
  but decreases the variablity of the estimate. This may be reasonable
  in some cases.

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

plot(x)
}
}

\keyword{htest}

