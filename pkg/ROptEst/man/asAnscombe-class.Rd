\name{asAnscombe-class}
\docType{class}
\alias{asAnscombe-class}
\alias{eff}
\alias{eff,asAnscombe-method}
\alias{show,asAnscombe-method}

\title{Asymptotic Anscombe risk}
\description{Class of asymptotic Anscombe risk which is
  the ARE (asymptotic relative efficiency) in the ideal model
  obtained by an optimal bias robust IC . }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("asAnscombe", ...)}.
  More frequently they are created via the generating function 
  \code{asAnscombe}. 
}
\section{Slots}{
  \describe{
    \item{\code{type}}{Object of class \code{"character"}: 
      \dQuote{optimal bias robust IC (OBRI) for given ARE (asymptotic relative efficiency)}. }
    \item{\code{eff}}{Object of class \code{"numeric"}: 
      given ARE (asymptotic relative efficiency) to be attained in the ideal model. }
    \item{\code{biastype}}{Object of class \code{"BiasType"}: 
      symmetric, one-sided or asymmetric }
  }
}
\section{Extends}{
Class \code{"asRiskwithBias"}, directly.\cr
Class \code{"asRisk"}, by class \code{"asRiskwithBias"}.
Class \code{"RiskType"}, by class \code{"asRisk"}.
}
\section{Methods}{
  \describe{
    \item{eff}{\code{signature(object = "asAnscombe")}: 
      accessor function for slot \code{eff}. }
    \item{show}{\code{signature(object = "asAnscombe")}}
  }
}
\references{Hampel et al. (1986) \emph{Robust Statistics}. 
  The Approach Based on Influence Functions. New York: Wiley.
  
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@fraunhofer.itwm.de}}
%\note{}
\seealso{\code{\link{asRisk-class}}, \code{\link{asAnscombe}}}
\examples{
new("asAnscombe")
}
\concept{Hampel risk}
\concept{risk}
\keyword{classes}
