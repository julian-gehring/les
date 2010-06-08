\name{optimalKernel}

\alias{optimalKernel}

\title{
  optimal
}

\description{
  ... to be added ... 
}

\usage{
  optimalKernel(object, winSize, regions, verbose=FALSE)
}

\arguments{
  \item{object}{Object of class 'Les' containing experimental data, as
    returned by 'create'. ...}
  \item{winSize}{Integer specifing window size. ...}
  \item{regions}{Function specifing the shape of the weighting
    window. If not specified the supplied 'triangWeight' function with a
    triangular window will be used. For details on other window
    functions and how to specify own function please see below. ...}
  \item{verbose}{Logical indicating whether Grenander correction
    should be used (default: FALSE). ...}
}

\details{
  This function estimates \eqn{\Lambda_i}{Lambda(i)} for all probes \eqn{i}{i}.
}
      
\value{
  Object of class 'Les' with additionally filled slots:
  chi2
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
  \code{\link{cutoff}}
  \code{\link{regions}}
  \code{\link{ci}}
}

\examples{
  \dontrun{
    pos <- 1:100
    pval <- runif(100)
    chr <- rep(1:2, each=50)

    res <- create(pos, pval, chr)

    res <- estimate(res, 20)
  }
}

\keyword{htest}

