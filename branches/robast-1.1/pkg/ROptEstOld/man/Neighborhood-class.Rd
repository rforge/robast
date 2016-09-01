\name{Neighborhood-class}
\docType{class}
\alias{Neighborhood-class}
\alias{radius}
\alias{radius,Neighborhood-method}
\alias{show,Neighborhood-method}
\alias{type,Neighborhood-method}

\title{Neighborhood}
\description{Class of neighborhoods of families of probability measures.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}: 
      type of the neighborhood. }
    \item{\code{radius}:}{Object of class \code{"numeric"}:
      neighborhood radius. }
  }
}
\section{Methods}{
  \describe{
    \item{type}{\code{signature(object = "Neighborhood")}: 
      accessor function for slot \code{type}. }
    \item{radius}{\code{signature(object = "Neighborhood")}: 
      accessor function for slot \code{radius}. }
    \item{show}{\code{signature(object = "Neighborhood")}}
  }
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{ProbFamily-class}}}
%\examples{}
\concept{neighborhood}
\keyword{classes}
\keyword{models}
