\name{CondNeighborhood-class}
\docType{class}
\alias{CondNeighborhood-class}
\alias{radiusCurve}
\alias{radiusCurve,CondNeighborhood-method}
\alias{show,CondNeighborhood-method}

\title{Conditional neighborhood}
\description{Class of conditonal (error-free-variables) neighborhoods.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      type of the neighborhood. }
    \item{\code{radius}:}{Object of class \code{"numeric"}:
      neighborhood radius. }
    \item{\code{radiusCurve}:}{Object of class \code{"function"}:
      radius curve.}
  }
}
\section{Extends}{
Class \code{"Neighborhood"}, directly.
}
\section{Methods}{
  \describe{
    \item{radiusCurve}{\code{signature(object = "CondNeighborhood")}:
      accessor function for slot \code{radiusCurve}. }
    \item{show}{\code{signature(object = "CondNeighborhood")}}
  }
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{Neighborhood-class}}
%\examples{}
\concept{conditional neighborhood}
\concept{neighborhood}
\keyword{classes}
\keyword{models}
\keyword{robust}
