\name{asUnOvShoot}
\alias{asUnOvShoot}

\title{Generating function for asUnOvShoot-class}
\description{
  Generates an object of class \code{"asUnOvShoot"}.
}
\usage{asUnOvShoot(width = 1.960)}
\arguments{
  \item{width}{ positive real: half the width of given confidence interval. }
}
%\details{}
\value{Object of class \code{"asUnOvShoot"}}
\references{
  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{asUnOvShoot-class}}}
\examples{
asUnOvShoot()

## The function is currently defined as
function(width = 1.960){ new("asUnOvShoot", width = width) }
}
\concept{asymptotic risk}
\concept{risk}
\keyword{robust}
