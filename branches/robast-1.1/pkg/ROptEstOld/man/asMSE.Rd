\name{asMSE}
\alias{asMSE}

\title{Generating function for asMSE-class}
\description{
  Generates an object of class \code{"asMSE"}.
}
\usage{asMSE()}
%\details{}
\value{Object of class \code{"asMSE"}}
\references{ 
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{asMSE-class}}}
\examples{
asMSE()

## The function is currently defined as
function(){ new("asMSE") }
}
\concept{asymptotic mean square error}
\concept{risk}
\keyword{robust}
