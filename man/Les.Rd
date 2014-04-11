\name{Les,ANY,numeric-method}

\alias{Les}
\alias{Les-methods}
\alias{Les,ANY,numeric-method}

\title{Les}

\description{
  Constructs an object of class 'Les' and stores experimental tiling
  microarray data. This is the initial step for all further analysis
  with the 'les' package.
}

\usage{
Les(pos, pval, chr)
}

\arguments{
  \item{pos}{Integer vector containing the probe positions.}
  \item{pval}{Numeric vector containing the p-values corresponding to
    'pos'. It must be of the same length as 'pos'.}
  \item{chr}{Vector specifying the chromosome each probe in located
    on. If missing all probes are located on one chromosome (default: 'chr0'). If
    specified it must be of the same length as 'pos'. Internally it will
    be stored as a factor.}
}

\details{
  This method gathers all data necessary for subsequent analysis,
  checks for valid inputs and stores it in an object of class 'Les'.

  The data is checked for the following criteria:
  'pos' and 'chr' must not contain any NAs. 'pval' may contain NAs but
  such probes (including corresponding 'pos' and 'chr') are discarded
  for subsequent computation since they contain no usable
  information. Please note that in such a case fewer probes are stored
  in the resulting object then were passed to 'Les'. In case of
  duplicate probe positions on one chromosome a warning is printed. This
  normally indicates that probes from both strands are present in the
  data.
  
  After storing the experimental data with 'Les' in an object the
  variables containing the original data can be deleted from the
  workspace. All further steps of 'les' will get their data from this
  object. This can be useful in cases when memory usage is a critical
  factor.
}
      
\value{
  Object of class 'Les' with filled slots:
  pos, pval, chr
}

\author{
  Julian Gehring
  
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
data(spikeInStat)

x <- Les(pos, pval)
x
}

\keyword{utilities}
\keyword{methods}

