\name{asHampel}
\alias{asHampel}

\title{Generating function for asHampel-class}
\description{
  Generates an object of class \code{"asHampel"}.
}
\usage{asHampel(bound = Inf)}
\arguments{
  \item{bound}{ positive real: bias bound }
}
%\details{}
\value{Object of class \code{asHampel}}
\references{
  Hampel et al. (1986) \emph{Robust Statistics}. 
  The Approach Based on Influence Functions. New York: Wiley.
  
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{asHampel-class}}}
\examples{
asHampel()

## The function is currently defined as
function(bound = Inf){ new("asHampel", bound = bound) }
}
\concept{Hampel risk}
\concept{risk}
\keyword{robust}
