\name{plot,Les-method}

\alias{plot}
\alias{plot-methods}
\alias{plot,Les}
\alias{plot,Les-method}

\title{plot method for class Les}

\description{
  The 'plot' method plots \eqn{\Lambda}{Lambda} along the chromosome.
}

\usage{

plot(x, y, ...)

\S4method{plot}{Les}(x, y, chr, region=FALSE, xlim, ylim=c(0, 1),
error="none", probePch=20, probeCol="black", lty=1, sigPch=20,
sigCol="red3", errorCol="azure4", regionCol=gray, rug=FALSE, rugSide=1,
limit=TRUE, main=NULL, ...)
}

\arguments{
  \item{x}{Object of class 'Les', as returned by 'estimate', 'threshold' or 'ci'.}
  \item{y}{Annotation object}
  \item{chr}{Character or numeric specifying which chromosome to
    plot. Must have a match in 'chr' passed to 'create'. A value is
    required if the probes are located on more than one chromosome.}
  \item{region}{Logical indicating whether the estimated regions should
    be included in the plot. The 'regions' method must have been called
    beforehand.}
  \item{xlim}{Vector with two elements specifying the range on the
    x-axis.}
  \item{ylim}{Vector with two elements specifying the range on the
    y-axis.}
  \item{error}{Character string specifying if error estimates for
    \eqn{\Lambda}{Lambda} should be plotted (default: "none"). Valid
    values are "none" for no error estimates and "ci" for confidence
    intervals computed by the 'ci' method.}
  \item{probePch}{Symbol probes should be plotted with (default:
    20). For details see 'par'.}
  \item{probeCol}{Color probes should be plotted with (default:
    "black"). For details see 'par'.}
  \item{lty}{Line type for lines connecting the
    \eqn{\Lambda_i}{Lambda_i} of probes (default: 1). For details see
    'par'}.
  \item{sigPch}{Symbol probes above the threshold should be plotted
    with (default: 20). For details see 'par'.}
  \item{sigCol}{Color probes above the threshold should be plotted
    with (default: "red3"). For details see 'par'.}
  \item{errorCol}{Color error estimates should be plotted with (default:
    "azure4").}
  \item{regionCol}{Function or vector of color used plotting of regions
    (default: gray). By default regions will be colored in grayscale
    with darker colors corresponding to more significant regions. If a
    function it has to return the color given 1-ri (see the 'regions'
    method for details) of each region as an input argument. If a vector
    of colors regions will be filled according to these colors. The vector
    will be recycled if needed. Note that the ordering if by default
    according to 'ri' from the 'regions' slot in the object, not the
    location of the region.}
  \item{rug}{Logical whether the positions of the probes should be
    indicated along the x-axis (default: FALSE). For details see 'rug'.}
  \item{rugSide}{Integer specifying the side at which probe positions
    should be plotted (default: 1). For details see 'rug'}.
  \item{limit}{Logical specifying whether the estimated threshold
    \eqn{\Theta}{Theta} should be indicated on the y-axis (default:
    TRUE) if estimated.}
  \item{main}{Title of the plot (default: no title). For details see 'main'}.
  \item{...}{Further arguments not used so far.}
}

\details{
  This method provides high-level plotting for the 'Les' class.}
      
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
x <- estimate(x, win=25)
x <- threshold(x)
x <- regions(x)

plot(x, chr=1, region=TRUE)
}
}

\keyword{utilities}
\keyword{methods}
\keyword{hplot}