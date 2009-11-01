\name{CondContNeighborhood-class}
\docType{class}
\alias{CondContNeighborhood-class}

\title{Conditional contamination neighborhood}
\description{Class of conditional (error-free-variables) 
  convex contamination neighborhoods.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("CondContNeighborhood", ...)}.
  More frequently they are created via the generating function 
  \code{CondContNeighborhood}. 
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}: 
      \dQuote{conditional convex contamination neighborhood}. }
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
%No methods defined with class "CondContNeighborhood" in the signature.
%}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{CondContNeighborhood}}, \code{\link{CondNeighborhood-class}}}
\examples{
new("CondContNeighborhood")
}
\concept{conditional contamination neighborhood}
\concept{conditional neighborhood}
\concept{neighborhood}
\keyword{classes}
\keyword{models}
\keyword{robust}
