\name{asL4}
\alias{asL4}

\title{Generating function for asL4-class}
\description{
  Generates an object of class \code{"asL4"}.
}
\usage{asL4(biastype = symmetricBias(), normtype = NormType())}
\arguments{
  \item{biastype}{a bias type of class \code{BiasType}}
  \item{normtype}{ a norm type of class \code{NormType}}
}

%\details{}
\value{Object of class \code{"asL4"}}
\references{ 
  Ruckdeschel, P. and Rieder, H. (2004) Optimal Influence Curves for
  General Loss Functions. Statistics & Decisions \emph{22}, 201-223.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{asL4-class}}, \code{\link{asMSE}}, \code{\link{asL1}}}
\examples{
asL4()

## The function is currently defined as
function(biastype = symmetricBias(), normtype = NormType()){ 
         new("asL4", biastype = biastype, normtype = normtype) }
}
\concept{asymptotic mean square error}
\concept{risk}
\keyword{robust}
