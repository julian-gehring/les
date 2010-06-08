\name{create}

\alias{create}

\title{
  create
}

\description{
  Creates an object of class LES and stores experimental data in it.
}

\usage{
  create(pos, pval, chr)
}

\arguments{
  \item{pos}{Integer vector containing the probe positions.}
  \item{pval}{Numeric vector containing the p-values corresponding to
    'pos'. Must be of the same length as 'pos'.}
  \item{chr}{Vector specifing the chromosome each probe in located
    on. If missing all probes are located on one chromosome ('chr0'). If
    specified must be of the same length as 'pos'.}
}

\details{
  This function gathers all data neccessary for subsequent computation,
  checks for validity and stores it in an object of class 'Les'.

  The data is checked for certain criteria:
  'pos' and 'chr' must not contain any NAs. 'pval' may contain NAs but
  such probes (including corresponding 'pos' and 'chr') are discarded
  and therefore not taken into account for subsequent computation since
  they contain no useable information. Please note that in such a case
  fewer probes are stored in the resulting object then were given into
  'create'.
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
  \code{\link{estimate}}
  \code{\link{cutoff}}
  \code{\link{regions}}
  \code{\link{ci}}
}

\examples{
## Generate a small data set
set.seed(0)
pos <- 1:1000
pval <- runif(1000)
pval[201:250] <- rexp(50, 20)
chr <- rep(1:2, each=500)

res <- create(pos, pval, chr)
}

\keyword{utilities}

