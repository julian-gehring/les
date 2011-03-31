\name{spikeInData}

\docType{data}

\alias{exprs}
\alias{reference}

\title{Spike-in ChIP-chip data set}

\description{
  This data set is part of a quality control study for tiling
  microarrays (Johnson et al., 2008), in which spike-ins were used to
  assess microarray platforms, preparation procedures, and analysis
  algorithms on the accuracy and reproducibility of ChIP-chip
  experiments. Here, the expression intensities of one region from the
  'undiluted' data set investigated with Affymetrix arrays is selected,
  consisting of 452 probes and two conditions with three replicates
  each. The data has been normalized using quantile normalization and
  probe positions remapped to a common reference.
}

\usage{data(spikeInData)}

\format{
  \describe{
    \item{exprs}{Matrix with expression intensities, with rows
      representing probes and columns arrays. The names of the rows and
      columns contain the probe position and the treatment,
      respectively.}
    \item{pval}{Vector with p-values assessing the differential effect
      between the conditions with a modified t-test. For details, please
      see the vignette of this package.}
    \item{reference}{Data frame containing the location of the
      spike-in.}
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

