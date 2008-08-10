\name{MEstimate-class}
\docType{class}
\alias{MEstimate-class}
\alias{Mroot}
\alias{Mroot,MEstimate-method}
\alias{show,MEstimate-method}

\title{MEstimate-class.}
\description{Class of asymptotically linear estimates.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("MEstimate", ...)}.
  More frequently they are created via the generating function 
  \code{locMEstimator}.
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
    \item{\code{Mroot}:}{Object of class \code{"numeric"}: value of
      the M equation at the estimate. }
    \item{\code{Infos}:}{ object of class \code{"matrix"}
      with two columns named \code{method} and \code{message}:
      additional informations. }
  }
}
\section{Extends}{
Class \code{"ALEstimate"}, directly.\cr
Class \code{"Estimate"}, by class \code{"ALEstimate"}.
}
\section{Methods}{
  \describe{
    \item{Mroot}{\code{signature(object = "MEstimate")}: 
      accessor function for slot \code{Mroot}. }

    \item{show}{\code{signature(object = "MEstimate")} }
  }
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{ALEstimate-class}}}
\examples{
## prototype
new("MEstimate")
}
\concept{estimate}
\keyword{classes}
