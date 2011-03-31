\name{plot,Les-method}

\alias{plot}
\alias{plot-methods}
\alias{plot,Les}
\alias{plot,Les-method}

\title{plot method for class Les}

\description{
  The 'plot' method plots the estimates of the 'les' package along the
  genome. This includeds \eqn{\Lambda}{Lambda} with confidence
  intervals and estimated regions.
}

\usage{
plot(x, y, ...)

\S4method{plot}{Les}(x, y, chr, error="none", region=FALSE, limit=TRUE,
rug=FALSE, xlim, ylim=c(0, 1),  ...)
}

\arguments{
  \item{x}{Object of class 'Les', as returned by 'estimate', 'threshold' or 'ci'.}
  \item{y}{Annotation object, currently not used.}
  \item{chr}{Character or numeric specifying which chromosome to
    plot. Must have a match in 'chr' passed to 'Les'. A value is
    required if the probes are located on more than one chromosome.}
  \item{error}{Character string specifying if error estimates for
    \eqn{\Lambda}{Lambda} should be plotted (default: "none"). Valid
    values are "none" for no error estimates and "ci" for confidence
    intervals computed by the 'ci' method.}
  \item{region}{Logical indicating whether the estimated regions should
    be included in the plot. The 'regions' method must have been called
    beforehand.}
  \item{limit}{Logical specifying whether the estimated threshold
    \eqn{\Theta}{Theta} should be indicated on the y-axis if estimated
    (default: TRUE).}
  \item{rug}{Logical whether the positions of the probes should be
    indicated along the x-axis (default: FALSE). For details see 'rug'.}
  \item{xlim}{Numeric vector with two elements specifying the range on the
    x-axis.}
  \item{ylim}{Numeric vector with two elements specifying the range on the
    y-axis.}
  \item{...}{Optional arguments used in order to customize the plot. See
    the \sQuote{details} section.}
}

\details{
  This method provides high-level plotting for the 'Les' class.

  The \code{plot} method uses a special system in order to customize the
  graphical elements of the figure. It allows to refer to the different
  components with the name of the additional input argument; its value
  is a list containing named graphical parameters for the underlying
  plot function. The following list describes the possible names and
  their contribution.

  \describe{
    \item{\code{plotArgs}}{Arguments for the axis and the labeling, passed to the
      \code{plot} function.}
    \item{\code{borderArgs}}{Arguments for the border lines at
      \eqn{\Lambda}{Lambda} equal to 0 and 1, passed to the
      \code{abline} function.}
    \item{\code{errorArgs}}{Arguments for the confidence interval of
      \eqn{\Lambda}{Lambda}, passed to the \code{plotCI} function of the
      \pkg{gplots} package.}
    \item{\code{probeArgs}}{Arguments for the representation of probes,
      passed to the \code{points} function.}
    \item{\code{limitArgs}}{Arguments for the horizontal line
      representing the threshold \eqn{\Theta}{Theta}, passed to
      the \code{abline} function.}
    \item{\code{sigArgs}}{Arguments for the representation of
      significant probes with \eqn{\Lambda >= \Theta}{Lambda equal or
	greater Theta}, passed to the \code{points} function.}
    \item{\code{rugArgs}}{Arguments for the representation of the probe
      coverage along the genome, passed to the \code{rug} function.}
    \item{\code{regionArgs}}{Arguments  for the representation of
      estimated LES, passed to the \code{rect} function. If 'col' is
      specified as a function it determines the color of each region
      depending on its input \eqn{\Lambda}{Lambda} (default: gray()).
      If 'col' is a vector its elements are used to color the regions
      with recycling.}
  }
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
}

\examples{
data(statTile)

x <- Les(pos, pval, chr)
x <- estimate(x, 50)
x <- threshold(x)
x <- regions(x)

plot(x, chr=1, region=TRUE)
}

\keyword{utilities}
\keyword{methods}
\keyword{hplot}
