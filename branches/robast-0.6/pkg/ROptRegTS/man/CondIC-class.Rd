\name{CondIC-class}
\docType{class}
\alias{CondIC-class}
\alias{CallL2Fam<-,CondIC-method}
\alias{checkIC,CondIC,missing-method}
\alias{checkIC,CondIC,L2RegTypeFamily-method}
\alias{show,CondIC-method}

\title{Conditionally centered partial influence curve}
\description{Class of conditionally centered partial influence curves.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("CondIC", ...)}.
  More frequently they are created via the generating function 
  \code{CondIC}. 
}
\section{Slots}{
  \describe{
    \item{\code{CallL2Fam}:}{Object of class \code{"call"}:
      creates an object of class \code{"L2RegTypeFamily"}. }
    \item{\code{name}:}{Object of class \code{"character"}: name }
    \item{\code{Curve}:}{Object of class \code{"EuclRandVariable"}: curve.}
    \item{\code{Risks}:}{Object of class \code{"list"}: 
      list of risks; cf. \code{RiskType-class}. }
    \item{\code{Infos}:}{Object of class \code{"matrix"}
      with two columns named \code{method} and \code{message}:
      additional informations. }
  }
}
\section{Extends}{
Class \code{"IC"}, directly.\cr
Class \code{"InfluenceCurve"}, by class \code{"IC"}.
}
\section{Methods}{
  \describe{
    \item{CallL2Fam<-}{\code{signature(object = "IC")}: 
      replacement function for slot \code{CallL2Fam}. }

    \item{checkIC}{\code{signature(IC = "CondIC", L2Fam = "missing")}: 
      check conditional centering and Fisher consistency of \code{IC} 
      assuming the L2-differentiable regression-type family 
      which can be created via the slot \code{CallL2Fam} of \code{IC}. }

    \item{checkIC}{\code{signature(IC = "CondIC", L2Fam = "L2RegTypeFamily")}: 
      check conditional centering and Fisher consistency of \code{IC} 
      assuming the L2-differentiable regression-type family \code{L2Fam}. }

    \item{show}{\code{signature(object = "CondIC")}}
  }
}
\references{
  Hampel et al. (1986) \emph{Robust Statistics}. 
  The Approach Based on Influence Functions. New York: Wiley.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{InfluenceCurve-class}, \code{IC-class}}
\examples{
new("CondIC")
}
\concept{partial influence curve}
\concept{influence curve}
\keyword{classes}
