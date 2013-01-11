\name{plot-methods}
\docType{methods}
\alias{plot}
\alias{plot-methods}
\alias{plot,IC,missing-method}
\title{ Methods for Function plot in Package `ROptEst' }
\description{plot-methods}
\details{
S4-Method \code{plot} for for signature \code{IC,missing} has been enhanced
compared to its original definition in \pkg{RobAStBase} so that if
argument \code{MBRB} is \code{NA}, it is filled automatically by a call
to \code{optIC} which computes the MBR-IC on the fly.}
}
\examples{
N <- NormLocationScaleFamily(mean=0, sd=1)
IC <- optIC(model = N, risk = asCov())
plot(IC2, main = TRUE, panel.first= grid(),
     col = "blue", cex.main = 2, cex.inner = 0.6,
     withMBR=TRUE, MBRB=FALSE)
}
\keyword{methods}
\keyword{distribution}

