\name{AvCondNeighborhood-class}
\docType{class}
\alias{AvCondNeighborhood-class}
\alias{show,AvCondNeighborhood-method}

\title{Average conditional neighborhood}
\description{Class of average conditional neighborhoods; i.e. 
  only radius curves \eqn{\varepsilon}{varepsilon} with 
  \eqn{\|\varepsilon\|_\alpha\le 1}{||varepsilon||_alpha <= 1}
  for given \code{exponent} \eqn{\alpha}{alpha}.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      type of the neighborhood. }
    \item{\code{radius}:}{Object of class \code{"numeric"}:
      neighborhood radius. }
    \item{\code{radiusCurve}:}{Object of class \code{"function"}:
      radius curve.}
    \item{\code{exponent}:}{Object of class \code{"numeric"}:
      positive integer or \code{Inf}. }
  }
}
\section{Extends}{
Class \code{"CondNeighborhood"}, directly.\cr
Class \code{"Neighborhood"}, by class \code{"CondNeighborhood"}.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "AvCondNeighborhood")}}
  }
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{CondNeighborhood-class}}}
%\examples{}
\concept{average conditional neighborhood}
\concept{conditional neighborhood}
\concept{neighborhood}
\keyword{classes}
\keyword{models}
\keyword{robust}
