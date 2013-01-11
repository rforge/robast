\name{comparePlot-methods}
\docType{methods}
\alias{comparePlot}
\alias{comparePlot-methods}
\alias{comparePlot,IC,IC-method}

\title{Compare - Plots}
\description{
  Plots 2-4 influence curves to the same model.
}
\details{
S4-Method \code{comparePlot} for signature \code{IC,IC} has been enhanced compared to
its original definition in \pkg{RobAStBase} so that if
argument \code{MBRB} is \code{NA}, it is filled automatically by a call
to \code{optIC} which computes the MBR-IC on the fly.}
}
\examples{
if(require(ROptEst)){

N0 <- NormLocationScaleFamily(mean=0, sd=1)
N0.Rob1 <- InfRobModel(center = N0,
           neighbor = ContNeighborhood(radius = 0.5))

IC1 <- optIC(model = N0, risk = asCov())
IC2 <- optIC(model = N0.Rob1, risk = asMSE())

comparePlot(IC1,IC2, withMBR=TRUE, MBRB=FALSE)
}
}
\keyword{robust}
