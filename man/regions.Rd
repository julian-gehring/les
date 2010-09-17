\name{regions,Les-method}

\alias{regions}
\alias{regions-methods}
\alias{regions,Les}
\alias{regions,Les-method}

\title{regions}

\description{
  Estimates distinct regions of continuously elevated
  \eqn{\Lambda_i}{Lambda} above a threshold \eqn{\Theta}{Theta}.
}
      
\usage{
regions(object, limit=NULL, minLength=2, maxGap=Inf, verbose = FALSE,
...)

\S4method{regions}{Les}(object, limit=NULL, minLength=2, maxGap=Inf,
verbose = FALSE, ...)
}

\arguments{
  \item{object}{Object of class 'Les' as returned by 'threshold'.}
  \item{limit}{Numeric specifying the threshold \eqn{\Theta}{Theta} for
    \eqn{\Lambda}{Lambda} with regions are defined as
    \eqn{\Lambda_i\geq\Theta}{Lambda>=Theta} (default: NULL). If 'NULL'
    \eqn{\Theta}{Theta} will be derived from the estimated number of
    regulated probes as computed by 'threshold'.}
  \item{minLength}{Integer specifying the minimum number of probes in a
    region (default: 2).}
  \item{maxGap}{Integer specifying maximum gap in base pairs between two
    neighboring probes in a region (default: Inf). If this value is
    exceeded the region will be split into smaller ones such that one
    region does not have neighboring probes with larger distance than
    'maxGap'. If not specified arbitrary large gap sizes are allowed.}
  \item{verbose}{Logical indicating whether a summary of the estimated
    regions should be printed on screen (default: FALSE).}
  \item{...}{Further arguments passed to subsequent functions.}
}

\details{
  This method finds distinct regions in \eqn{\Lambda}{Lambda} by
  thresholding with \eqn{\Theta}{Theta}. The regions have to meet the
  following criteria:

  (1) For all probes in the region
  \eqn{\Lambda_i\geq\Theta}{Lambda>=Theta} has to hold.
  (2) Each region has to contain at least as many probes as specified in
  'minLength'.
  (3) The gap between to neighboring probes has to be smaller or equal
  to 'maxGap'.

  Along with the boundaries of the regions the number of regulated
  probes within each region is estimated. This is used to sort the
  regions in order to get a list of top regions. The resulting data
  frame containing the regions can be accessed with the '[' method.
 }

\value{
  Object of class 'Les' with additionally filled slots:
  regions, limit

  'regions' is a data frame with variables:
  \item{chr}{Chromosome the regions is located on.}
  \item{start}{Position of the beginning of the region.}
  \item{end}{Position of the end of the region.}
  \item{size}{Extend of the region in base pairs (measured in the same
    units as input 'pos', normally base pairs).}
  \item{nProbes}{Number of probes in the region.}
  \item{ri}{Regulation Index of the region indicating the fraction of
    regulated probes in the region.}
  \item{se}{Standard error of the estimation of 'ri'.}
  \item{rs}{Regulation score defined as 'ri'/'se'. The data frame is
    sorted according to this variable.}
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
x <- estimate(x, win=50)
x <- threshold(x, grenander=TRUE, verbose=TRUE)
x <- regions(x, verbose=TRUE)

print(x["regions"])
}

\keyword{htest}
\keyword{methods}
