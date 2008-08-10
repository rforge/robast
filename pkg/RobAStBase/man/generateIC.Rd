\name{generateIC}
\alias{generateIC}

\title{Generic function for the generation of influence curves}
\description{
  This function is rarely called directly. It is used 
  by other functions to create objects of class \code{"IC"}.
}
\usage{
generateIC(neighbor, L2Fam, ...)
}
\arguments{
  \item{neighbor}{ Object of class \code{"Neighborhood"}. }
  \item{L2Fam}{ L2-differentiable family of probability measures. }
  \item{\dots}{ additional parameters }
}
%\details{}
\value{Object of class \code{"IC"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{IC-class}}, \code{\link{ContIC-class}}, \code{\link{TotalVarIC-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
