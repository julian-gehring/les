\name{exportLambda}

\alias{exportLambda}

\title{
  exportLambda
}

\description{
  Exports the estimates \eqn{\Lambda_i}{Lambda} of a region of interest
  in a '.wig' file. The exported file can be read by others software or
  loaded into various genome browsers.
}

\usage{
exportLambda(object, file, chr, range, description = "Lambda", precision
= 4)
}

\arguments{
  \item{object}{Object of class 'Les' containing experimental data, as
    returned by 'estimate'.}
  \item{file}{Character string specifing the file path and name for
    export.}
  \item{chr}{Character string specifing results from which chromosome
    should be exported. 'chr' must have a match in the 'chr' argument
    specified in 'create'.}
  \item{range}{Numeric vector specifing the probe range in position
    which should be exported. If missing all probes of the chromsome
    will be exported.}
  \item{description}{Character string with description for the exported
    track (default: 'Lambda'). This will be used as description by
    several programs and genome browsers.}
  \item{precision}{Integer specifing the number of digits the result
    should be rounded to for export (default: 4).}
}

\details{
  This function is useful to export the estimated Lambda to external
  programs and genome browsers.

  The '.wig' format allows the export of large amount of data and is
  compatible most genome browsers. For details about the file format see
  \url{http://genome.ucsc.edu/goldenPath/help/wiggle.html}.
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
x <- threshold(x)
x <- regions(x)

exportLambda(x, "chr1", "out.wig")
}
}

\keyword{htest}

