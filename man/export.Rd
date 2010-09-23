\name{export,Les-method}

\alias{export}
\alias{export-methods}
\alias{export,Les}
\alias{export,Les-method}

\title{export}

\description{
  Exports the results into files for interaction with other
  software. Estimated regions can be exported in 'bed' and 'gff' format,
  \eqn{\Lambda}{Lambda} in 'wig' format.
}

\usage{
export(object, file, format="bed", chr, range, description = "Lambda",
strand=".", group="les", precision=4, ...)

\S4method{export}{Les}(object, file, format="bed", chr, range,
description = "Lambda", strand=".", group="les", precision=4, ...)
}

\arguments{
  \item{object}{Object of class 'Les' containing experimental data, as
    returned by 'estimate' or 'regions'.}
  \item{file}{Character string specifying the file path and name for
    export.}
  \item{format}{Character string with the export format (default:
    'bed'). Possible values are 'bed' and 'gff' for export of the
    estimated regions and 'wig' for export of \eqn{\Lambda}{Lambda} for
    the probes. If set to 'bed' or 'gff' the method 'regions' has to be run
    beforehand. Partial matches are also possible. For details on the
    formats please see the 'Details' section.}
  \item{chr}{Character string specifying the chromosome from which
    results should be exported. This value must be set if exporting to
    'wig' format, for other formats is optional. 'chr' must have exactly
    one match in the 'chr' argument specified in 'Les'.}
  \item{range}{Numeric vector specifying the range of probe positions
    which should be exported. If missing all probes of the chromosome
    will be exported. This value has only an effect if 'format' is set
    to 'wig'.}
  \item{description}{Character string with description for the exported
    track (default: 'Lambda'). This will be used as description by
    several programs and genome browsers.}
  \item{strand}{Character string with strand specification for 'gff'
    export (default: '.'). Possible values are '+', '-' or '.'.}
  \item{group}{Character string with group specifications of the
    resulting regions in 'gff' format (default: 'les').}
  \item{precision}{Integer specifying the number of digits
    \eqn{\Lambda}{Lambda} should be rounded to for export to 'wig'
    format (default: 4).}
  \item{...}{Further arguments passed to subsequent functions.}
}

\details{
  This function is useful to export the estimated Lambda to external
  programs and genome browsers.

  The 'bed', 'gff' and 'wig' format provide widely used standard formats
  and are compatible with most genome browsers and related software. For
  details about the file formats see
  \url{http://genome.ucsc.edu/FAQ/FAQformat.html}.
  
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
\dontrun{
data(statTile)

x <- Les(pos, pval, chr)
x <- estimate(x, 50)
x <- threshold(x)
x <- regions(x)

export(x, file="test.bed")
export(x, file="test.gff", format="gff")
export(x, file="test.wig", format="wig", chr=1)
}
}

\keyword{utilities}
\keyword{methods}
\keyword{IO}
