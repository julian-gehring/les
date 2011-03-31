\name{les-package}

\alias{les-package}
\alias{les}

\docType{package}

\title{Identifying Differential Effects in Tiling Microarray Data}

\description{
  The 'les' package estimates Loci of Enhanced Significance (LES) in
  tiling microarray data. These are regions of regulation such as found
  in differential transcription, CHiP-chip, or DNA modification
  analysis. The package provides a universal framework suitable for
  identifying differential effects in tiling microarray data sets, and
  is independent of the underlying statistics at the level of single
  probes.
}

\details{
  The 'les' package provides a universal framework for detecting
  differential effects in tiling microarray experiments.

  It is universal in the sense that one is free to choose any
  statistical test to estimate the effect of differential effect for
  each probe on the tiling microarray. Provided with the p-values for each probe
  and the corresponding positions of the probes, 'les' uses a sliding
  window approach to estimate the fraction of regulated probes in the
  local surrounding of each probe. The approach is related to computing
  a spatially resolved and weighted false discovery rate, and yields a
  interpretable statistical feature \eqn{\Lambda}{Lambda}.

  Resulting regions can be scored according to their overall
  effect. Methods for high-level plotting and export of the results to
  other software and genome browsers are provided.

  The 'les' package is published under the GPL-3 license.
}

\author{
  Julian Gehring, Clemens Kreutz, Jens Timmer

  Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}

\references{
  in preparation

  This package is based on:
  
  Kilian Bartholome, Clemens Kreutz, and Jens Timmer:
  Estimation of gene induction enables a relevance-based ranking of gene sets,
  Journal of Computational Biology: A Journal of Computational Molecular Cell 
  Biology 16, no. 7 (July 2009): 959-967.
  \url{http://www.liebertonline.com/doi/abs/10.1089/cmb.2008.0226}
}

\seealso{
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
x <- estimate(x, 50)
x <- threshold(x)
x <- regions(x)
x <- ci(x, subset=150:200, conf=0.90, nBoot=50)

## plot data of chromosome 1
plot(x, chr=1, region=TRUE)
plot(x, chr=1, region=TRUE, error="ci")

\dontrun{
## export data of chromosome 1
export(x, file="les_out.bed", chr=1)
export(x, file="les_out.wig", format="wig", chr=1)
}
}

\keyword{package}
\keyword{htest}
