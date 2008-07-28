\name{ALEstimate-class}
\docType{class}
\alias{ALEstimate-class}
\alias{pIC}
\alias{pIC,ALEstimate-method}

\title{ALEstimate-class.}
\description{Class of asymptotically linear estimates.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ALEstimate", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Object of class \code{"character"}:
      name of the estimator. }
    \item{\code{estimate}:}{Object of class \code{"ANY"}:
      estimate. }
    \item{\code{pIC}:}{Optional object of class \code{InfluenceCurve}:
      influence curve. }
    \item{\code{Infos}:}{ object of class \code{"matrix"}
      with two columns named \code{method} and \code{message}:
      additional informations. }
  }
}
\section{Extends}{
Class \code{"Estimate"}, directly.
}
\section{Methods}{
  \describe{
    \item{pIC}{\code{signature(object = "ALEstimate")}: 
      accessor function for slot \code{pIC}. }
  }
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[distrMod]{Estimate-class}}}
%\examples{}
\concept{estimate}
\keyword{classes}
