\name{Les-class}

\Rdversion{1.1}
\docType{class}

\alias{Les-class}
\alias{[<-,Les-method}
\alias{[,Les-method}
\alias{show,Les-method}
\alias{summary,Les-method}

\title{Class 'Les'}

\description{
  The class 'Les' is used by the package 'les'.
}

\section{Objects from the Class}{
  Objects of class 'Les' are created by calling the function 'create'.
}
       
\section{Slots}{
  \describe{
    \item{\code{pos}:}{Integer vector with probe positions. Provided by
      the user through the 'create' function.}
    \item{\code{chr}:}{Factor with chromosomes for each probe. Provided
      by the user through the 'create' function.}
    \item{\code{lambda}:}{Numeric vector with Lambda for each probe}
    \item{\code{se}:}{Numeric vector with standard error from linear
      fitting for each probe.}
    \item{\code{nProbes}:}{Integer with number of probes used in fit for
      each probe}
    \item{\code{pval}:}{Numeric vector with p-values for each
      probe. Provided by the user through the 'create' function.}
    \item{\code{win}:}{Integer with the window size. Set by the user in
      the 'estimate' function.}
    \item{\code{weighting}:}{Function with the weighting function used
      for probe weighting according to position. Set by the user in the
      'estimate' function.}
    \item{\code{grenander}:}{Logical indicating usage of the Grenander
      correction. Set by the user in the 'estimate' function.}
    \item{\code{ci}:}{Data frame with the confidence intervals for each
      probe specified in 'subset'.}
    \item{\code{nBoot}:}{Integer with the number of bootstraps to be
      drawn. Set by the user in the 'ci' function.}
    \item{\code{conf}:}{Numeric with the confidence level. Set by the
      user in the 'ci' function.}
    \item{\code{subset}:}{Integer vector with indices of bootstrap
      subset. Set by the user in the 'ci' function.}
    \item{\code{theta}:}{Numeric with the threshold value Theta.}
    \item{\code{nSigProbes}:}{Integer with the number of estimated
      significant probes.xs}
    \item{\code{regions}:}{Data frame with estimated Loci of Enhanced
      Significance.}
    \item{\code{limit}:}{Numeric with the threshold value for estimation of
      'regions'. Set by the user in the 'regions' function.}
    \item{\code{nChr}:}{Integer with the number of chromosomes.}
    \item{\code{maxGap}:}{Numeric specifying the largest gap allowed in
      one region. Set by the user in the 'regions' function.}
    \item{\code{minLength}:}{Integer specifying the minimal number of
      probes in one region. Set by the user in the 'regions'
      function.}
    \item{\code{minProbes}:}{Integer specifying the minimal number of
      unique p-values allowed for each fit. Set by the user in the
      'estimate' function.}
    \item{\code{method}:}{Character specifying the method used for
      linear regression.}
  }
}

\author{
  Julian Gehring, Clemens Kreutz

  Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}

\seealso{
  Package:
  \code{\link[les]{les-package}}

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
showClass("Les")
}
}

\keyword{classes}
\keyword{methods}
