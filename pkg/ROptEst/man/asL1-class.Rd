\name{asL1-class}
\docType{class}
\alias{asL1-class}

\title{Asymptotic mean absolute error}
\description{Class of asymptotic mean absolute error.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("asL1", ...)}.
  More frequently they are created via the generating function 
  \code{asL1}.   
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
No methods defined with class "asL1" in the signature.
}
\references{
  Ruckdeschel, P. and Rieder, H. (2004) Optimal Influence Curves for
  General Loss Functions. Statistics & Decisions \emph{22}, 201-223.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{asGRisk-class}}, \code{\link{asMSE}}, \code{\link{asMSE-class}}, \code{\link{asL4-class}}, \code{\link{asL1}}}
\examples{
new("asMSE")
}
\concept{mean absolute error risk}
\concept{risk}
\keyword{classes}
