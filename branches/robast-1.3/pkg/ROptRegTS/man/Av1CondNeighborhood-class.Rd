\name{Av1CondNeighborhood-class}
\docType{class}
\alias{Av1CondNeighborhood-class}

\title{Average conditional neighborhood}
\description{Class of average conditional neighborhoods 
  (exponent \code{== 1}); i.e. only radius curves 
  \eqn{\varepsilon}{varepsilon} with 
  \eqn{\|\varepsilon\|_1\le 1}{||varepsilon||_1 <= 1}.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      type of the neighborhood. }
    \item{\code{radius}:}{Object of class \code{"numeric"}:
      neighborhood radius. }
    \item{\code{radiusCurve}:}{Object of class \code{"function"}:
      radius curve with L1 norm \code{<= 1}.}
    \item{\code{exponent}:}{ equal to \code{1}. }
  }
}
\section{Extends}{
Class \code{"AvCondNeighborhood"}, directly.\cr
Class \code{"CondNeighborhood"}, by class \code{"AvCondNeighborhood"}.\cr
Class \code{"Neighborhood"}, by class \code{"AvCondNeighborhood"}.
}
%\section{Methods}{
%No methods defined with class "Av1CondNeighborhood" in the signature.
%}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{AvCondNeighborhood-class}}}
%\examples{}
\concept{average conditional neighborhood}
\concept{conditional neighborhood}
\concept{neighborhood}
\keyword{classes}
\keyword{models}
\keyword{robust}
