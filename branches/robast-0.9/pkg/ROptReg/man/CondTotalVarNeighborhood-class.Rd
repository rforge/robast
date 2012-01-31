\name{CondTotalVarNeighborhood-class}
\docType{class}

\alias{CondTotalVarNeighborhood-class}
\title{Conditional total variation neighborhood}
\description{Class of conditional (error-free-variables) 
  total variation neighborhoods.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("CondTotalVarNeighborhood", ...)}.
  More frequently they are created via the generating function 
  \code{CondTotalVarNeighborhood}. 
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}: 
      \dQuote{conditional total variation neighborhood}. }
    \item{\code{radius}:}{Object of class \code{"numeric"}:
      neighborhood radius. }
    \item{\code{radiusCurve}:}{Object of class \code{"function"}:
      radius curve }
  }
}
\section{Extends}{
Class \code{"CondNeighborhood"}, directly.\cr
Class \code{"Neighborhood"}, by class \code{"CondNeighborhood"}.
}
%\section{Methods}{
%No methods defined with class "CondTotalVarNeighborhood" in the signature.
%}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{CondTotalVarNeighborhood}}, \code{\link{CondNeighborhood-class}}}
\examples{
new("CondTotalVarNeighborhood")
}
\concept{conditional total variation neighborhood}
\concept{conditional neighborhood}
\concept{neighborhood}
\keyword{classes}
\keyword{models}
\keyword{robust}
