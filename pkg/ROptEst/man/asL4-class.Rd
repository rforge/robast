\name{asL4-class}
\docType{class}
\alias{asL4-class}

\title{Asymptotic mean power 4 error}
\description{Class of asymptotic mean power 4 error.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("asL4", ...)}.
  More frequently they are created via the generating function 
  \code{asL4}.   
}
\section{Slots}{
  \describe{
    \item{\code{type}}{Object of class \code{"character"}: 
      \dQuote{asymptotic mean square error}. }
    \item{\code{biastype}}{Object of class \code{"BiasType"}: 
      symmetric, one-sided or asymmetric }
    \item{\code{normtype}}{Object of class \code{"NormType"}: 
      norm in which a multivariate parameter is considered}
  }
}
\section{Extends}{
Class \code{"asGRisk"}, directly.\cr
Class \code{"asRiskwithBias"}, by class \code{"asGRisk"}.\cr
Class \code{"asRisk"}, by class \code{"asRiskwithBias"}.\cr
Class \code{"RiskType"}, by class \code{"asGRisk"}.
}
\section{Methods}{
No methods defined with class "asL4" in the signature.
}
\references{
  Ruckdeschel, P. and Rieder, H. (2004) Optimal Influence Curves for
  General Loss Functions. Statistics & Decisions \emph{22}, 201-223.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{asGRisk-class}}, \code{\link{asMSE}}, \code{\link{asMSE-class}}, \code{\link{asL1-class}}, \code{\link{asL4}}}
\examples{
new("asMSE")
}
\concept{asymptotic risk}
\concept{risk}
\keyword{classes}
