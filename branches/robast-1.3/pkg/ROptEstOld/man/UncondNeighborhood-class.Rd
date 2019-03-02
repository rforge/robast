\name{UncondNeighborhood-class}
\docType{class}

\alias{UncondNeighborhood-class}
\title{Unconditional neighborhood}
\description{Class of unconditonal (errors-in-variables) neighborhoods.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      type of the neighborhood. }
    \item{\code{radius}:}{Object of class \code{"numeric"}:
      neighborhood radius. }
  }
}
\section{Extends}{
Class \code{"Neighborhood"}, directly.
}
%\section{Methods}{}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{Neighborhood-class}}}
%\examples{}
\concept{unconditional neighborhood}
\concept{neighborhood}
\keyword{classes}
\keyword{models}
