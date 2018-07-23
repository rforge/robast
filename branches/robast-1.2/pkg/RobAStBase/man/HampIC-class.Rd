\name{HampIC-class}
\docType{class}
\alias{HampIC-class}
\alias{lowerCase}
\alias{lowerCase,HampIC-method}
\alias{neighborRadius}
\alias{neighborRadius,ANY-method}
\alias{neighborRadius,HampIC-method}
\alias{neighborRadius<-}
\alias{neighborRadius<-,HampIC-method}
\alias{stand}
\alias{stand,HampIC-method}
\alias{weight,HampIC-method}
\alias{biastype,HampIC-method}
\alias{normtype,HampIC-method}

\title{Influence curve of Hampel type}
\description{Class of (partial) influence curves of Hampel (= total variation or contamination) type;
used as common mother class for classes \code{ContIC} and \code{TotalVarIC}.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("HampIC", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{CallL2Fam}}{ object of class \code{"call"}:
      creates an object of the underlying L2-differentiable
      parametric family. }

    \item{\code{name}}{ object of class \code{"character"} }

    \item{\code{Curve}}{ object of class \code{"EuclRandVarList"}}

    \item{\code{modifyIC}}{ Object of class \code{"OptionalFunction"}:
      function of two arguments, which are an L2 parametric family 
      and an optional influence curve. Returns an object of 
      class \code{"IC"}. This slot is mainly used for internal
      computations! }

    \item{\code{Risks}}{ object of class \code{"list"}:
      list of risks; cf. \code{\link[distrMod]{RiskType-class}}. }

    \item{\code{Infos}}{ object of class \code{"matrix"}
      with two columns named \code{method} and \code{message}:
      additional informations. }

    \item{\code{stand}}{ object of class \code{"matrix"}:
      standardizing matrix. }

    \item{\code{weight}}{ object of class \code{"RobWeight"}:
      weight function }

    \item{\code{biastype}}{ object of class \code{"BiasType"}:
      bias type (symmetric/onsided/asymmetric) }

    \item{\code{normtype}}{ object of class \code{"NormType"}:
      norm type (Euclidean, information/self-standardized)}

    \item{\code{lowerCase}}{ object of class \code{"OptionalNumeric"}:
      optional constant for lower case solution. }

    \item{\code{neighborRadius}}{ object of class \code{"numeric"}:
      radius of the corresponding (unconditional) contamination
      neighborhood. }
  }
}
\section{Extends}{
Class \code{"IC"}, directly.\cr
Class \code{"InfluenceCurve"}, by class \code{"IC"}.
}
\section{Methods}{
  \describe{

    \item{stand}{\code{signature(object = "HampIC")}:
      accessor function for slot \code{stand}. }

    \item{weight}{\code{signature(object = "HampIC")}:
      accessor function for slot \code{weight}. }

    \item{biastype}{\code{signature(object = "HampIC")}:
      accessor function for slot \code{biastype}. }
    \item{normtype}{\code{signature(object = "HampIC")}:
      accessor function for slot \code{normtype}. }
    \item{lowerCase}{\code{signature(object = "HampIC")}:
      accessor function for slot \code{lowerCase}. }

    \item{neighborRadius}{\code{signature(object = "HampIC")}:
      accessor function for slot \code{neighborRadius}. }

    \item{neighborRadius<-}{\code{signature(object = "HampIC")}:
      replacement function for slot \code{neighborRadius}. }

    \item{neighborRadius}{\code{signature(object = "ANY")}:
      returns \code{NULL}. }
  }
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Hampributions to the Asymptotic Theory of Robustness}.
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{IC-class}}}
\examples{
IC1 <- new("HampIC")
plot(IC1)
}
\concept{influence curve}
\keyword{classes}
