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
to \code{optIC} which computes the MBR-IC on the fly. To this end, there
is an additional argument \code{n.MBR} defaulting to 10000
to determine the number of evaluation points.
points.
}
\examples{
N <- NormLocationScaleFamily(mean=0, sd=1)
IC <- optIC(model = N, risk = asCov())
## Don't run to reduce check time on CRAN
\donttest{
plot(IC, main = TRUE, panel.first= grid(),
     col = "blue", cex.main = 2, cex.inner = 0.6,
     withMBR=TRUE)
}
}
\keyword{methods}
\keyword{distribution}

