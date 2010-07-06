\name{create}

\alias{create}

\title{create}

\description{
  Creates an object of class LES and stores experimental tiling
  microarray data. This is the initial step for all further analysis.
}

\usage{
create(pos, pval, chr)
}

\arguments{
  \item{pos}{Integer vector containing the probe positions of one strand.}
  \item{pval}{Numeric vector containing the p-values corresponding to
    'pos'. It must be of the same length as 'pos'.}
  \item{chr}{Vector specifing the chromosome each probe in located
    on. If missing all probes are located on one chromosome (default: 'chr0'). If
    specified it must be of the same length as 'pos'.}
}

\details{
  This function gathers all data neccessary for subsequent computation,
  checks for valid inputs and stores it in an object of class 'Les'.

  The data is checked for the following criteria:
  'pos' and 'chr' must not contain any NAs. 'pval' may contain NAs but
  such probes (including corresponding 'pos' and 'chr') are discarded
  for subsequent computation since they contain no useable
  information. Please note that in such a case fewer probes are stored
  in the resulting object then were passed to 'create'. In case of
  duplicate probe positions on one chromosome a warning is printed. This
  normally indicates that probes from both strands are present in the
  data. For this analysis it is recommended to analyse only data from
  one strand at a time.
  
  After storing the experimental data with 'create' in an object the
  variables containing the original data can be deleted from the
  workspace. All further steps of 'les' will get their data from this
  obejct. This can be useful in cases when memory useage is a critical
  factor.
}
      
\value{
  Object of class 'Les' with filled slots:
  pos, pval, chr
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
data(statTile)
x <- create(pos, pval, chr)
}
}

\keyword{utilities}

