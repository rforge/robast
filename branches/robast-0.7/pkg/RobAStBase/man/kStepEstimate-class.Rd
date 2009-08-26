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
    \item{\code{name}}{Object of class \code{"character"}:
      name of the estimator. }
    \item{\code{estimate}}{Object of class \code{"ANY"}:
      estimate. }
    \item{\code{estimate.call}}{Object of class \code{"call"}:
      call by which estimate was produced.}
    \item{\code{samplesize}}{ object of class \code{"numeric"} ---
      the samplesize (only complete cases are counted)
      at which the estimate was evaluated. }
    \item{\code{completecases}:}{ object of class \code{"logical"} ---
         complete cases at which the estimate was evaluated. }
    \item{\code{asvar}}{ object of class \code{"OptionalNumericOrMatrix"}
      which may contain the asymptotic (co)variance of the estimator. }
    \item{\code{asbias}}{Optional object of class \code{"numeric"}:
      asymptotic bias. }
    \item{\code{pIC}}{Optional object of class \code{InfluenceCurve}:
      influence curve. }
    \item{\code{nuis.idx}}{ object of class \code{"OptionalNumeric"}: 
     indices of \code{estimate} belonging to the nuisance part. }
    \item{\code{fixed}}{ object of class \code{"OptionalNumeric"}:
     the fixed and known part of the parameter.}
    \item{\code{steps}}{Object of class \code{"integer"}: number
      of steps. }
    \item{\code{Infos}}{ object of class \code{"matrix"}
      with two columns named \code{method} and \code{message}:
      additional informations. }
    \item{\code{trafo}}{ object of class \code{"list"}:
     a list with components \code{fct} and \code{mat} (see below). }
    \item{\code{untransformed.estimate}:}{Object of class \code{"ANY"}:
      untransformed estimate.}
    \item{\code{untransformed.asvar}:}{ object of class \code{"OptionalNumericOrMatrix"}
      which may contain the asymptotic (co)variance of the untransformed
      estimator. }
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

    \item{show}{\code{signature(object = "kStepEstimate")}: a show method; }
  }
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{ALEstimate-class}}}
%\examples{}
\concept{estimate}
\keyword{classes}
