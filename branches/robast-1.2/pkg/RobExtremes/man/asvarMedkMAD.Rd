\name{asvarMedkMAD}
\alias{asvarMedkMAD}

\title{ Function to compute asymptotic variance of MedkMAD estimator}
\description{
  Function \code{asvarMedkMAD} computes the asymptotic (co)variance of
  a MedkMAD estimator at a Scale-Shape model.
}
\usage{
asvarMedkMAD( model, k=1)
}
\arguments{
  \item{model}{an object of class \code{"ScaleShapeUnion"}. }
  \item{k}{numeric (>0); additional parameter for \code{\link{kMAD}}. }
}
\details{
For the Generalized Pareto Family all terms are analytic; in case
of the general scale-shape model, numerical integration is used.
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
\seealso{\code{\link{LDEstimator}} }
\examples{
GP <- GParetoFamily(scale=1,shape=0.7)
asvarMedkMAD(GP,k=1)

## for didactical purposes turn GP into a non-GPD
setClass("noGP",contains="L2ScaleShapeUnion")
GP2 <- GP
class(GP2) <- "noGP"
asvarMedkMAD(GP2,k=1) ### uses numerical integration
}
\keyword{asymptotic variance}
