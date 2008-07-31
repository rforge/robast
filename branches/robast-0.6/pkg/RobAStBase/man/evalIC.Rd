\name{evalIC}
\alias{evalIC}

\title{Generic function for evaluating ICs}
\description{
  Generic function for evaluating ICs.
}
\usage{
evalIC(IC, x)
}
\arguments{
  \item{IC}{ object of class \code{"IC"} }
  \item{x}{ numeric vector or matrix }
}
\details{ The list of random variables contained in the slot \code{Curve} 
  is evaluated at \code{x}. }
\value{ In case \code{x} is numeric a vector and in case \code{x} 
  is matrix a matrix is returned. } 
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{IC-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
