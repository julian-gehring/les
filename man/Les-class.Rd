\name{Les-class}
\Rdversion{1.1}
\docType{class}
\alias{Les-class}
\alias{[<-,Les-method}
\alias{[,Les-method}
\alias{ci,Les-method}
\alias{cutoff,Les-method}
\alias{plot,Les-method}
\alias{regions,Les-method}
\alias{show,Les-method}
\alias{summary,Les-method}

\title{Class "Les" ~~~ }
\description{Class 'Les' used by the 'les' package}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Les", ...)}.
	 ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{pos}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{q}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{se}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nUsedProbes}:}{Object of class \code{"integer"} ~~ }
    \item{\code{pval}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{win}:}{Object of class \code{"integer"} ~~ }
    \item{\code{weighting}:}{Object of class \code{"function"} ~~ }
    \item{\code{grenander}:}{Object of class \code{"logical"} ~~ }
    \item{\code{ci}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{nBoot}:}{Object of class \code{"integer"} ~~ }
    \item{\code{conf}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{subset}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{cutoff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nSigProbes}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{regions}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{limit}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{[<-}{\code{signature(x = "Les")}: ... }
    \item{[}{\code{signature(x = "Les")}: ... }
    \item{ci}{\code{signature(res = "Les")}: ... }
    \item{cutoff}{\code{signature(res = "Les")}: ... }
    \item{plot}{\code{signature(x = "Les")}: ... }
    \item{regions}{\code{signature(res = "Les")}: ... }
    \item{show}{\code{signature(object = "Les")}: ... }
    \item{summary}{\code{signature(object = "Les")}: ... }
	 }
}
\references{ ... }
\author{
Julian Gehring, Clemens Kreutz

Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}
\note{ ~~further notes~~ }

\seealso{
	objects to See Also as \code{\link{estimate}}
}
\examples{
showClass("Les")
}
\keyword{classes}
\keyword{methods}
