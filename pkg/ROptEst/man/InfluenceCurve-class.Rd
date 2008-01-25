\name{InfluenceCurve-class}
\docType{class}
\alias{InfluenceCurve-class}
\alias{addInfo<-}
\alias{addInfo<-,InfluenceCurve-method}
\alias{addRisk<-}
\alias{addRisk<-,InfluenceCurve-method}
\alias{Curve}
\alias{Curve,InfluenceCurve-method}
\alias{Domain,InfluenceCurve-method}
\alias{Infos}
\alias{Infos,InfluenceCurve-method}
\alias{Infos<-}
\alias{Infos<-,InfluenceCurve-method}
\alias{Map,InfluenceCurve-method}
\alias{name,InfluenceCurve-method}
\alias{name<-,InfluenceCurve-method}
\alias{Range,InfluenceCurve-method}
\alias{Risks}
\alias{Risks,InfluenceCurve-method}
\alias{Risks<-}
\alias{Risks<-,InfluenceCurve-method}
\alias{show,InfluenceCurve-method}

\title{Influence curve}
\description{Class of influence curves (functions).}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("InfluenceCurve", ...)}.
  More frequently they are created via the generating function 
  \code{InfluenceCurve}. 
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{ object of class \code{"character"} }
    \item{\code{Curve}:}{ object of class \code{"EuclRandVarList"} }
    \item{\code{Risks}:}{ object of class \code{"list"}: 
      list of risks; cf. \code{\link{RiskType-class}}. }
    \item{\code{Infos}:}{ object of class \code{"matrix"}
      with two columns named \code{method} and \code{message}:
      additional informations. }
  }
}
\section{Methods}{
  \describe{
    \item{name}{\code{signature(object = "InfluenceCurve")}: 
      accessor function for slot \code{name}. }

    \item{name<-}{\code{signature(object = "InfluenceCurve")}: 
      replacement function for slot \code{name}. }

    \item{Curve}{\code{signature(object = "InfluenceCurve")}: 
      accessor function for slot \code{Curve}. }

    \item{Map}{\code{signature(object = "InfluenceCurve")}: 
      accessor function for slot \code{Map} of slot \code{Curve}. }

    \item{Domain}{\code{signature(object = "InfluenceCurve")}:
      accessor function for slot \code{Domain} of slot \code{Curve}. }

    \item{Range}{\code{signature(object = "InfluenceCurve")}:
      accessor function for slot \code{Range} of slot \code{Curve}. }

    \item{Infos}{\code{signature(object = "InfluenceCurve")}: 
      accessor function for slot \code{Infos}. }

    \item{Infos<-}{\code{signature(object = "InfluenceCurve")}: 
      replacement function for slot \code{Infos}. }

    \item{addInfo<-}{\code{signature(object = "InfluenceCurve")}: 
      function to add an information to slot \code{Infos}. }

    \item{Risks}{\code{signature(object = "InfluenceCurve")}: 
      accessor function for slot \code{Risks}. }

    \item{Risks<-}{\code{signature(object = "InfluenceCurve")}: 
      replacement function for slot \code{Risks}. }

    \item{addRisk<-}{\code{signature(object = "InfluenceCurve")}: 
      function to add a risk to slot \code{Risks}. }

    \item{show}{\code{signature(object = "InfluenceCurve")}}
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
\seealso{\code{\link{InfluenceCurve}}, \code{\link{RiskType-class}}}
\examples{
new("InfluenceCurve")
}
\concept{influence curve}
\keyword{classes}
