\name{RegTypeFamily-class}
\docType{class}
\alias{RegTypeFamily-class}
\alias{ErrorDistr}
\alias{ErrorDistr,RegTypeFamily-method}
\alias{ErrorSymm}
\alias{ErrorSymm,RegTypeFamily-method}
\alias{RegDistr}
\alias{RegDistr,RegTypeFamily-method}
\alias{Regressor}
\alias{Regressor,RegTypeFamily-method}
\alias{RegSymm}
\alias{RegSymm,RegTypeFamily-method}
\alias{show,RegTypeFamily-method}

\title{Parametric regression-type family}
\description{Class for parametric regression-type families.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("RegTypeFamily", ...)}.
  More frequently they are created via the generating function 
  \code{RegTypeFamily}. 
}
\section{Slots}{
  \describe{
    \item{\code{ErrorDistr}:}{ object of class \code{"Distribution"}: 
      error distribution. }
    \item{\code{ErrorSymm}:}{ object of class \code{"DistributionSymmetry"}: 
      symmetry of the error distribution.}
    \item{\code{RegDistr}:}{ object of class \code{"Distribution"}: 
      regressor distribution. }
    \item{\code{RegSymm}:}{ object of class \code{"DistributionSymmetry"}: 
      symmetry of the regressor distribution.}
    \item{\code{Regressor}:}{ object of class \code{"EuclRandVariable"}:
      regressor. }
    \item{\code{param}:}{ object of class \code{"ParamFamParameter"}:
      parameter of the family. }
    \item{\code{props}:}{ object of class \code{"character"}:
      properties of the family. }
    \item{\code{name}:}{ object of class \code{"character"}:
      name of the family. }
    \item{\code{distribution}:}{ object of class \code{"CondDistribution"}:
      distribution given the regressor. }
  }
}
\section{Extends}{
Class \code{"ParamFamily"}, directly.\cr
Class \code{"ProbFamily"}, by class \code{"ParamFamily"}.
}
\section{Methods}{
  \describe{
    \item{ErrorDistr}{\code{signature(object = "RegTypeFamily")}: 
      accessor function for slot \code{ErrorDistr}. }

    \item{ErrorSymm}{\code{signature(object = "RegTypeFamily")}: 
      accessor function for slot \code{ErrorSymm}. }

    \item{RegDistr}{\code{signature(object = "RegTypeFamily")}: 
      accessor function for slot \code{RegDistr}. }

    \item{Regressor}{\code{signature(object = "RegTypeFamily")}: 
      accessor function for slot \code{Regressor}. }

    \item{RegSymm}{\code{signature(object = "RegTypeFamily")}: 
      accessor function for slot \code{RegSymm}. }

    \item{show}{\code{signature(object = "RegTypeFamily")}}
  }
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{ParamFamily-class}}
\examples{
new("RegTypeFamily")
}
\keyword{classes}
\keyword{models}
