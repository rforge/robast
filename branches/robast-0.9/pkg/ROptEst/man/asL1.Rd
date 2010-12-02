\name{asL1}
\alias{asL1}

\title{Generating function for asMSE-class}
\description{
  Generates an object of class \code{"asMSE"}.
}
\usage{asL1(biastype = symmetricBias(), normtype = NormType())}
\arguments{
  \item{biastype}{a bias type of class \code{BiasType}}
  \item{normtype}{ a norm type of class \code{NormType}}
}

%\details{}
\value{Object of class \code{"asMSE"}}
\references{ 
  Ruckdeschel, P. and Rieder, H. (2004) Optimal Influence Curves for
  General Loss Functions. Statistics & Decisions \emph{22}, 201-223.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{asL1-class}}, \code{\link{asMSE}}, \code{\link{asL4}}}
\examples{
asL1()

## The function is currently defined as
function(biastype = symmetricBias(), normtype = NormType()){ 
         new("asL1", biastype = biastype, normtype = normtype) }
}
\concept{asymptotic mean square error}
\concept{risk}
\keyword{robust}
