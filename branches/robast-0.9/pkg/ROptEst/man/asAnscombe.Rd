\name{asAnscombe}
\alias{asAnscombe}

\title{Generating function for asAnscombe-class}
\description{
  Generates an object of class \code{"asAnscombe"}.
}
\usage{asAnscombe(eff = .95, biastype = symmetricBias(), normtype = NormType())}
\arguments{
  \item{eff}{ value in (0,1]: ARE in the ideal model }
  \item{biastype}{ a bias type of class \code{BiasType}}
  \item{normtype}{ a norm type of class \code{NormType}}
}
%\details{}
\value{Object of class \code{asAnscombe}}
\references{
  Hampel et al. (1986) \emph{Robust Statistics}. 
  The Approach Based on Influence Functions. New York: Wiley.
  
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@fraunhofer.itwm.de}}
%\note{}
\seealso{\code{\link{asAnscombe-class}}}
\examples{
asAnscombe()

## The function is currently defined as
function(eff = .95, biastype = symmetricBias(), normtype = NormType()){ 
    new("asAnscombe", eff = eff, biastype = biastype, normtype = normtype) }
}
\concept{Hampel risk}
\concept{risk}
\keyword{robust}
