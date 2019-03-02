\name{asvarQBCC}
\alias{asvarQBCC}

\title{ Function to compute asymptotic variance of QuantileBCC estimator}
\description{
  Function \code{asvarQBCC} computes the asymptotic (co)variance of
  a QuantileBCC estimator at a Weibull model.
}
\usage{
asvarQBCC( model, p1 = 1/3, p2 = 2/3)
}
\arguments{
  \item{model}{an object of class \code{"ScaleShapeUnion"}. }
  \item{p1,p2}{levels of the quantiles; maximal breakdown point is achieved
               for \eqn{p1=p2-p1=1-p2=1/3} which is the default.}
}
\details{
All terms are analytic.
}
\value{
  A 2x2 matrix; the covariance. }


%\references{  }
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{QuantileBCCEstimator}} }
\examples{
GP <- WeibullFamily(scale=1,shape=0.7)
asvarQBCC(GP)
asvarQBCC(GP, p1=1/4, p2= 5/8)
}
\keyword{asymptotic variance}
