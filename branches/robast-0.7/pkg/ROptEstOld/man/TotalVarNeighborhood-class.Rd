\name{TotalVarNeighborhood-class}
\docType{class}
\alias{TotalVarNeighborhood-class}

\title{Total variation neighborhood}
\description{Class of (unconditional) total variation neighborhoods.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("TotalVarNeighborhood", ...)}.
  More frequently they are created via the generating function 
  \code{TotalVarNeighborhood}. 
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      \dQuote{(uncond.) total variation neighborhood}. }
    \item{\code{radius}:}{Object of class \code{"numeric"}:
      neighborhood radius. }
  }
}
\section{Extends}{
Class \code{"UncondNeighborhood"}, directly.\cr
Class \code{"Neighborhood"}, by class \code{"UncondNeighborhood"}.
}
\section{Methods}{
No methods defined with class "TotalVarNeighborhood" in the signature.
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{TotalVarNeighborhood}}, \code{\link{UncondNeighborhood-class}}}
\examples{
new("TotalVarNeighborhood")
}
\concept{total variation neighborhood}
\concept{neighborhood}
\keyword{classes}
\keyword{models}
