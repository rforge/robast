\name{LDEstimate-class}
\docType{class}
\alias{LDEstimate-class}
\alias{dispersion}
\alias{dispersion,LDEstimate-method}
\alias{location,LDEstimate-method}
\alias{show,LDEstimate-method}

\title{LDEstimate-class.}
\description{Class of Location Dispersion estimates.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("LDEstimate", ...)}.
  More frequently they are created via the generating function
  \code{LDEstimator}.
}
\section{Slots}{
  \describe{
    \item{\code{name}}{Object of class \code{"character"}:
      name of the estimator. }
    \item{\code{estimate}}{Object of class \code{"ANY"}:
      estimate.}
    \item{\code{estimate.call}}{Object of class \code{"call"}:
      call by which estimate was produced.}
    \item{\code{dispersion}}{Object of class \code{"numeric"}:
      the value of the fitted dispersion.}
    \item{\code{location}}{Object of class \code{"numeric"}:
      the value of the fitted location.}
    \item{\code{Infos}}{ object of class \code{"matrix"}
      with two columns named \code{method} and \code{message}:
      additional informations. }
    \item{\code{asvar}}{ object of class \code{"OptionalMatrix"}
      which may contain the asymptotic (co)variance of the estimator. }
    \item{\code{samplesize}}{ object of class \code{"numeric"} ---
      the samplesize at which the estimate was evaluated. }
    \item{\code{nuis.idx}}{ object of class \code{"OptionalNumeric"}: 
     indices of \code{estimate} belonging to the nuisance part}
    \item{\code{fixed}}{ object of class \code{"OptionalNumeric"}:
     the fixed and known part of the parameter. }
    \item{\code{trafo}}{ object of class \code{"list"}:
     a list with components \code{fct} and \code{mat} (see below). }
    \item{\code{untransformed.estimate}}{Object of class \code{"ANY"}:
      untransformed estimate.}
    \item{\code{untransformed.asvar}}{ object of class \code{"OptionalNumericOrMatrix"}
      which may contain the asymptotic (co)variance of the untransformed 
      estimator. }
    \item{\code{completecases}}{ object of class \code{"logical"} ---
         complete cases at which the estimate was evaluated. }
  }
}
\section{Extends}{
Class \code{"Estimate"}, directly.
}
\section{Methods}{
  \describe{
    \item{dispersion}{\code{signature(object = "LDEstimate")}:
      accessor function for slot \code{dispersion}. }

    \item{location}{\code{signature(object = "LDEstimate")}:
      accessor function for slot \code{location}. }

    \item{show}{\code{signature(object = "LDEstimate")}}

  }
}

%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{Estimate-class}}, \code{\link{LDEstimator}},
         \code{\link{MCEstimator}}}
\examples{
## (empirical) Data
x <- rgamma(50, scale = 0.5, shape = 3)

## parametric family of probability measures
G <- GammaFamily(scale = 1, shape = 2)

(S <- medQn(x, G))
dispersion(S)
location(S)
}
\concept{estimate}
\keyword{classes}
