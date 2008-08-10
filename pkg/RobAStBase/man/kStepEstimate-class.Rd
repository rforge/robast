\name{kStepEstimate-class}
\docType{class}
\alias{kStepEstimate-class}
\alias{steps}
\alias{steps,kStepEstimate-method}
\alias{show,kStepEstimate-method}

\title{kStepEstimate-class.}
\description{Class of asymptotically linear estimates.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("kStepEstimate", ...)}.
  More frequently they are created via the generating function 
  \code{kStepEstimator}.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Object of class \code{"character"}:
      name of the estimator. }
    \item{\code{estimate}:}{Object of class \code{"ANY"}:
      estimate. }
    \item{\code{samplesize}:}{Object of class \code{"numeric"}:
      sample size. }
    \item{\code{asvar}:}{Optional object of class \code{"matrix"}:
      asymptotic variance. }
    \item{\code{asbias}:}{Optional object of class \code{"numeric"}:
      asymptotic bias. }
    \item{\code{pIC}:}{Optional object of class \code{InfluenceCurve}:
      influence curve. }
    \item{\code{nuis.idx}:}{ object of class \code{"OptionalNumeric"}: 
     indices of \code{estimate} belonging to the nuisance part. }
    \item{\code{steps}:}{Object of class \code{"integer"}: number
      of steps. }
    \item{\code{Infos}:}{ object of class \code{"matrix"}
      with two columns named \code{method} and \code{message}:
      additional informations. }
  }
}
\section{Extends}{
Class \code{"ALEstimate"}, directly.\cr
Class \code{"Estimate"}, by class \code{"ALEstimate"}
}
\section{Methods}{
  \describe{
    \item{steps}{\code{signature(object = "kStepEstimate")}: 
      accessor function for slot \code{steps}. }

    \item{show}{\code{signature(object = "kStepEstimate")} }
  }
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{ALEstimate-class}}}
%\examples{}
\concept{estimate}
\keyword{classes}
