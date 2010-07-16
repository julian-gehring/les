\name{les-package}

\alias{les-package}
\alias{les}

\docType{package}

\title{Loci of Enhanced Significance for Tiling Microarray Data}

\description{
  The 'les' package estimates Loci of Enhanced Significance (LES) in
  tiling microarray data. These are regions of regulation such as found
  in differential transcription, ChiP-chip or methylation
  experiments. The package provides a universal framework suitable of
  identifying regulatory effects in various data sets.
}

\details{
  The 'les' package provides a universal framework for detecting regions
  of regulation in tiling microarray experiments.

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

  The 'les' package is published under the GPL-3 license.
}

\author{
  Julian Gehring, Clemens Kreutz

  Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}

\references{
  in preparation
}

\seealso{
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
x <- estimate(x, 25)
x <- threshold(x)
x <- regions(x)
x <- ci(x)

plot(x)
export(x)
}
}

\keyword{package}
\keyword{htest}
