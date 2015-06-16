\name{asvarPickands}
\alias{asvarPickands}

\title{ Function to compute asymptotic variance of Pickands estimator}
\description{
  Function \code{asvarPickands} computes the asymptotic (co)variance of
  a Pickands estimator at a GPD or GEVD model.
}
\usage{
asvarPickands( model, alpha=2)
}
\arguments{
  \item{model}{an object of class \code{"ScaleShapeUnion"}. }
  \item{alpha}{ numeric > 1; determines the variant of the Pickands-Estimator
   based on matching the empirical \eqn{a_1=1-1/\alpha}{a1=1-1/alpha} and
   \eqn{a_1=1-1/\alpha^2}{a1=1-1/alpha^2} quantiles against the
   population counter parts. The ``classical'' Pickands Estimator is
   obtained for \code{alpha=2} (GPD) resp. for \code{alpha=1/log(2)} (GEVD). }
}
\details{
All terms are analytic.
}
\value{
  A 2x2 matrix; the covariance. }

\references{
Ruckdeschel, P. and Horbenko, N. (2011): Optimally-Robust Estimators in Generalized
Pareto Models. ArXiv 1005.1476. To appear at \emph{Statistics}.
DOI: 10.1080/02331888.2011.628022. \cr

}

%\references{  }
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{PickandsEstimator}} }
\examples{
GP <- GParetoFamily(scale=1,shape=0.7)
asvarPickands(GP)
asvarPickands(GP,alpha=2.3)
}
\keyword{asymptotic variance}
