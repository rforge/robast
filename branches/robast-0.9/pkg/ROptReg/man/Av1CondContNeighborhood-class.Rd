\name{Av1CondContNeighborhood-class}
\docType{class}
\alias{Av1CondContNeighborhood-class}

\title{Average conditional contamination neighborhood}
\description{Class of average conditional contamination 
  neighborhoods (exponent \code{== 1}); i.e. only radius 
  curves \eqn{\varepsilon}{varepsilon} with 
  \eqn{\|\varepsilon\|_1\le 1}{||varepsilon||_1 <= 1}.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("Av1CondContNeighborhood", ...)}.
  More frequently they are created via the generating function 
  \code{Av1CondContNeighborhood}. 
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      \dQuote{average conditional convex contamination neighborhood}. }
    \item{\code{radius}:}{Object of class \code{"numeric"}:
      neighborhood radius. }
    \item{\code{radiusCurve}:}{Object of class \code{"function"}:
      radius curve with L1 norm \code{<= 1}.}
    \item{\code{exponent}:}{ equal to \code{1}. }
  }
}
\section{Extends}{
Class \code{"Av1CondNeighborhood"}, directly.\cr
Class \code{"AvCondNeighborhood"}, by class \code{"Av1CondNeighborhood"}.\cr
Class \code{"CondNeighborhood"}, by class \code{"Av1CondNeighborhood"}.\cr
Class \code{"Neighborhood"}, by class \code{"Av1CondNeighborhood"}.
}
%\section{Methods}{
%No methods defined with class "Av1CondContNeighborhood" in the signature.
%}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{Av1CondNeighborhood-class}}}
\examples{
new("Av1CondContNeighborhood")
}
\concept{average conditional neighborhood}
\concept{conditional neighborhood}
\concept{neighborhood}
\keyword{classes}
\keyword{models}
\keyword{robust}
