\name{asBias}
\alias{asBias}

\title{Generating function for asBias-class}
\description{
  Generates an object of class \code{"asBias"}.
}
\usage{asBias()}
%\details{}
\value{Object of class \code{"asBias"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{asBias-class}}}
\examples{
asBias()

## The function is currently defined as
function(){ new("asBias") }
}
\concept{asymptotic bias}
\concept{risk}
\keyword{robust}
