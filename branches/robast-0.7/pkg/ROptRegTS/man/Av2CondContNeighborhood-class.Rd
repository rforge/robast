\name{Av2CondContNeighborhood-class}
\docType{class}
\alias{Av2CondContNeighborhood-class}

\title{Average square conditional contamination neighborhood}
\description{Class of average square conditional contamination 
  neighborhoods (exponent \code{== 2}); i.e. only radius 
  curves \eqn{\varepsilon}{varepsilon} with 
  \eqn{\|\varepsilon\|_2\le 1}{||varepsilon||_2 <= 1}.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("Av2CondContNeighborhood", ...)}.
  More frequently they are created via the generating function 
  \code{Av2CondContNeighborhood}. 
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      \dQuote{average square conditional convex contamination neighborhood}. }
    \item{\code{radius}:}{Object of class \code{"numeric"}:
      neighborhood radius. }
    \item{\code{radiusCurve}:}{Object of class \code{"function"}:
      radius curve with L2 norm \code{<= 1}.}
    \item{\code{exponent}:}{ equal to \code{2}. }
  }
}
\section{Extends}{
Class \code{"Av2CondNeighborhood"}, directly.\cr
Class \code{"AvCondNeighborhood"}, by class \code{"Av2CondNeighborhood"}.\cr
Class \code{"CondNeighborhood"}, by class \code{"Av2CondNeighborhood"}.\cr
Class \code{"Neighborhood"}, by class \code{"Av2CondNeighborhood"}.
}
%\section{Methods}{
%No methods defined with class "Av2CondContNeighborhood" in the signature.
%}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{Av2CondNeighborhood-class}}}
\examples{
new("Av2CondContNeighborhood")
}
\concept{average conditional neighborhood}
\concept{conditional neighborhood}
\concept{neighborhood}
\keyword{classes}
\keyword{models}
\keyword{robust}
