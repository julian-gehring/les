\name{les-package}
\alias{les-package}
\alias{les}
\docType{package}
\title{Finds loci of enhanced significance for differential expression in tiling array data}
\description{
Finds loci of enhanced significance in tiling array data using the LES score estimating the amount of regulation.
}
\details{
\tabular{ll}{
Package: \tab les\cr
Type: \tab Package\cr
Version: \tab 0.5.0\cr
Date: \tab 2010-04-15\cr
License: \tab What license is it under?\cr
License: \tab GPL-2\cr
Depends: \tab methods\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Julian Gehring, Clemens Kreutz

Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}
\references{
...
}

\keyword{package}

\seealso{
	Class:
	\code{\linkS4class{Les}}
	
	Functions:
	\code{\link{estimate}}
	\code{\link{cutoff}}
	\code{\link{regions}}
	\code{\link{ci}}
}

\examples{
pos <- 1:300
pval <- runif(300)
pval[101:150] <- runif(50, 0, 0.05)
pval[201:250] <- runif(50, 0, 0.05)
win <- 50

anno <- list(pos=matrix(c(101,150, 201,250), nrow=2),
             strand=c(1, 1),
             name=c("gene1", "RNA5"),
             col=c("lightblue", "lightgreen")
             )
}
