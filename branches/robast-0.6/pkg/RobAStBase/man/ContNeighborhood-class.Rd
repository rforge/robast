\name{ContNeighborhood-class}
\docType{class}
\alias{ContNeighborhood-class}

\title{Contamination Neighborhood}
\description{Class of (unconditional) contamination neighborhoods.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ContNeighborhood", ...)}.
  More frequently they are created via the generating function 
  \code{ContNeighborhood}. 
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      \dQuote{(uncond.) convex contamination neighborhood}. }

    \item{\code{radius}:}{Object of class \code{"numeric"}: 
      neighborhood radius. }
  }
}
\section{Extends}{
Class \code{"UncondNeighborhood"}, directly.\cr
Class \code{"Neighborhood"}, by class \code{"UncondNeighborhood"}.
}
\section{Methods}{
No methods defined with class "ContNeighborhood" in the signature.
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{ContNeighborhood}}, \code{\link{UncondNeighborhood-class}}}
\examples{
new("ContNeighborhood")
}
\concept{contamination neighborhood}
\concept{neighborhood}
\keyword{classes}
\keyword{models}
