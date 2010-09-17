\name{estimate,Les-method}

\alias{estimate}
\alias{estimate-methods}
\alias{estimate,Les-method}

\title{estimate}

\description{
  The 'estimate' method computes the fraction \eqn{\Lambda}{Lambda} of
  probes with significant effect in the local surrounding of the genome.
}

\usage{
estimate(object, win, weighting = triangWeight, grenander = TRUE,
se = FALSE, minProbes = 3, method = "la", nCores = NULL, verbose =
FALSE, ...)

\S4method{estimate}{Les}(object, win, weighting = triangWeight,
grenander = TRUE, se = FALSE, minProbes = 3, method = "la", nCores =
NULL, verbose = FALSE, ...)
}

\arguments{
  \item{object}{Object of class 'Les' containing experimental data, as
    returned by 'Les'.}
  \item{win}{Integer specifying window size for the weighting
    function. This value is directly passed to the function specified
    by 'weighting'. For details see the description for the used window
    function 'weighting'.}
  \item{weighting}{Function specifying the shape of the weighting
    window. If not specified the supplied 'triangWeight' function with a
    triangular window will be used. For details on other window
    functions and how to specify own functions please see the 'Details'
    section.}
  \item{grenander}{Logical specifying if the Grenander correction for
    the cumulative density should be used (default: TRUE). For details
    see the 'Details' section.}
  \item{se}{Logical indicating whether the standard error (SE) from the
    final linear model should be computed and stored (default: FALSE). The
    standard error displays the goodness of fit for every probe, but is
    not needed for further computation. If computation time is a
    critical factor computation of the SE can be omitted to save some
    time.}
  \item{minProbes}{Integer specifying the minimal number of unqiue
    p-values that must be present for each fit (default: 3). For very
    small number of p-values the cumulative density is not well defined
    and therefore estimation has a high uncertainty. If the number of
    unique p-values is smaller than 'minProbes' no estimation is
    performed for this probe and \eqn{\Lambda=}{Lambda=}NA is returned.}
  \item{method}{Character string specifying the method used for linear
    regression (default: 'la'). Possible options are 'la' for a method
    based on linear algebra or 'qr' for a method based on qr
    decomposition. 'la' will be faster for few probes, 'qr' for many
    probes in a window. The best choice varies between data sets,
    parameters and machines. However this option only influences
    computation time but not the results.}
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
  This function estimates \eqn{\Lambda_i}{Lambda(i)} for all probes
  \eqn{i}{i}. This is normally the first step in the analysis after
  storing the experimental data with 'Les'.

  The 'win' argument influences the number of neighboring probes taken
  into account by the weighting function. The value is passed to the
  function specified in 'weighting'. Larger values result in a smoother
  features. Details on a reasonable choice for this value can be found
  in the references.

  With the 'weighting' argument the applied weighting function can be
  specified from a predefined set or a custom function can be used. In
  the 'les' package the four functions 'triangWeight', 'rectangWeight',
  'epWeight' and 'gaussWeight' are already supplied and offer a triangular,
  rectangular, Epanechnikov and Gaussian weighting function
  respectively. For details on the functions itself and how to use
  custom ones please see the documentation of the single functions or
  the vignette of this package. 

  The Grenander correction for the cumulative density includes the
  general knowledge about the concave shape of the cumulative
  density. This reduces the variance of the estimates and leads to a
  conservative estimation. Please note that using this feature may
  significantly increase computation time.

  The 'multicore' package can be used to spread the computation over several
  cores in a simple way. This can be useful on multi-core machines for
  large datasets. The 'multicore' package is not available on all
  platforms. To use multicore processing 'library(multicore)' has to be
  called beforehand and a number of cores to use has to be specified in
  'nCores'. For details see the documentation of the 'multicore package.

  Please note that calling 'estimate' with an object returned by the
  methods 'ci' and 'regions' will overwrite data stored by these two
  methods. This ensures that no inconsistent data is stored.
}
      
\value{
  Object of class 'Les' with additionally filled slots:
  lambda, win, weighting, grenander, nProbes, se (se only if computed)
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
data(statTile)

x <- Les(pos, pval, chr)
x <- estimate(x, 50)

x
plot(x, chr=1)
}

\keyword{htest}
\keyword{methods}
