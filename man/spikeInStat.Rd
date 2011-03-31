\name{spikeInStat}

\docType{data}

\alias{pos}
\alias{pval}

\title{Spike-in ChIP-chip data set}

\description{
  This data set is part of a quality control study for tiling
  microarrays (Johnson et al., 2008), in which spike-ins were used to
  assess the influence of microarray platforms, preparation procedures,
  and analysis algorithms on the accuracy and reproducibility of
  ChIP-chip experiments. Here, the expression intensities of one region
  from the 'undiluted' data set investigated with Affymetrix arrays is
  selected, consisting of 452 probes and two conditions with three
  replicates each. The data has been normalized using quantile
  normalization and probe positions remapped to a common reference.
}

\usage{data(spikeInStat)}

\format{
  \describe{
    \item{pos}{Vector with probe positions.}
    \item{pval}{Vector with p-values assessing the differential effect
      between the conditions with a modified t-test. For details, please
      see the vignette of this package.}
  }
}

\source{GEO accession IDs: GSM248996, ..., GSM249001}

\references{
  Johnson, D. S. et al. (2008).
  Systematic evaluation of variability in ChIP-chip experiments using
  predefined DNA targets.
  Genome Research, 18(3):393-403.
}

\keyword{datasets}

