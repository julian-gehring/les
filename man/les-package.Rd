\name{les-package}

\alias{les-package}
\alias{les}

\docType{package}

\title{Loci of Enhanced Significance for Tiling Microarray Data}

\description{
  The 'les' package estimates the Loci of Enhanced Significance (LES) in tiling microarray data. These are regions of regulation such as found in differential transcription, chip-CHIP or methylation experiments. The package provides a universal framework suitable of identifying regulatory effects in various data sets.
}

\details{
  \tabular{ll}{
    Package: \tab les\cr
    Type: \tab Package\cr
    Version: \tab 0.9.0\cr
    Date: \tab 2010-07-08\cr
    License: \tab GPL-3\cr
    Depends: \tab methods\cr
  }

  The 'les' package provides a universal framework to detect regions of
  regulation in tiling microarray experiments.

  It is universal in the sense that the user is free to choose any
  statistical test to estimate the effect of differential regulation for
  each probe. Given the p-values for each probe and the corresponding
  positions of the probes 'les' uses a sliding window approach to
  estimate the fraction of regulated probes in the local surrounding of
  each probe. The approach is related to computing a locally weighted
  false discovery rate and yields a interpretable statistical feature
  \eqn{\Lambda}{Lambda}.

  Resulting regions can be scored according to their overall
  effect. Methods for high-level plotting and export of the results to
  other software and genome browsers are provided.
}

\author{
  Julian Gehring, Clemens Kreutz

  Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}

\references{
  ... to be added ...
}

\keyword{package}

\seealso{
  Class:
  \code{\linkS4class{Les}}
  
  Methods and functions:
  \code{\link{create}}
  \code{\link{estimate}}
  \code{\link{threshold}}
  \code{\link{regions}}
  \code{\link{ci}}
  \code{\link{export}}
}

\examples{
\dontrun{
data(statTile)

x <- create(pos, pval, chr)
x <- estimate(x, win=25)
x <- threshold(x)
x <- regions(x)
x <- ci(x)

plot(x)
export(x)
}
}
