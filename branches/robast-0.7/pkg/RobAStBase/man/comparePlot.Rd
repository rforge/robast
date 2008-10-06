\name{comparePlot-methods}
\docType{methods}
\alias{comparePlot}
\alias{comparePlot-methods}
\alias{comparePlot,IC,IC-method}

\title{Compare - Plots}
\description{
  Plots 2-4 influence curves to the same model.
}
\usage{
comparePlot(obj1, obj2, ... )
\S4method{comparePlot}{IC,IC}(obj1, obj2, obj3 = NULL, obj4 = NULL, ... )
}
\arguments{
  \item{obj1}{ object of class \code{"InfluenceCurve"} }
  \item{obj2}{ object of class \code{"InfluenceCurve"} to be compared with \code{obj1}}
  \item{obj3}{ optional: object of class \code{"InfluenceCurve"} to be compared with \code{obj1}}
  \item{obj4}{ optional: object of class \code{"InfluenceCurve"} to be compared with \code{obj1}}
  \item{\dots}{further arguments to be passed to \code{plot}}
}
%\value{}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link[distrMod]{L2ParamFamily-class}}, \code{\link{IC-class}}, \code{\link[base]{plot}}}
\examples{
if(require(ROptEst)){

N0 <- NormLocationScaleFamily(mean=0, sd=1) 
N0.Rob1 <- InfRobModel(center = N0, neighbor = ContNeighborhood(radius = 0.5))

IC1 <- optIC(model = N0, risk = asCov())
IC2 <- optIC(model = N0.Rob1, risk = asMSE())

comparePlot(IC1,IC2)
}
}
\keyword{robust}
