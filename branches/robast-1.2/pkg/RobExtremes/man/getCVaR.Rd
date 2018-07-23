\name{getCVaR}
\alias{getVaR}
\alias{getCVaR}
\alias{getEL}
\alias{print.riskMeasure}

\title{Risk Measures for Scale-Shape Families}
\description{
  Functions to compute Value-at-Risk (VaR), Conditional Value-at-Risk (CVaR)
  and Expected Loss (EL) at data from scale-shape families.
}
\usage{
getVaR(data, model, level, rob=TRUE)
getCVaR(data, model, level, rob=TRUE)
getEL(data, model, N0, rob=TRUE)
\method{print}{riskMeasure}(x, level=NULL, ...)
}
\arguments{
  \item{data}{data at which to compute the risk measure. }
  \item{model}{an object of class \code{"L2ScaleShapeFamily"}.
               The parametric family at which to evaluate the risk measure. }
  \item{level}{real: probability needed for VaR and CVaR. }
  \item{N0}{real: expected frequency for expected loss. }
  \item{rob}{logical; if \code{TRUE} (default) the RMXE-parametric estimator is
             used; otherwise the MLE. }
  \item{x}{an object of (S3-)class \code{"riskmeasure"}. }
  \item{...}{further arguments for \code{print}. }
 }
\value{The risk measures \code{getVaR}, \code{getCVaR}, \code{getEL} return
an (S3) object of class \code{"riskMeasure"}, i.e., a numeric vector
of length 2 with components \code{"Risk"} and \code{"varofRisk"}
containing the respective risk measure
and a corresponding (asymptotic) standard error for the risk
measure. To the return class \code{"riskMeasure"},
there is a particular \code{print}-method; if the corresponding argument
\code{level} is \code{NULL} (default) the corresponding standard error
is printed together with the risk measure; otherwise a corresponding
CLT-based confidence interval for the risk meausre is produced.}
\references{
P. Ruckdeschel, N. Horbenko (2013): Optimally-Robust Estimators in Generalized
Pareto Models. Statistics 47(4), 762--791.
N. Horbenko, P. Ruckdeschel, T. Bae (2011): Robust Estimation of Operational Risk.
Journal of Operational Risk 6(2), 3--30.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
\seealso{\code{\link{GParetoFamily}}, \code{\link{GEVFamily}}, \code{\link{WeibullFamily}}, \code{\link{GammaFamily}}}
\examples{
  \donttest{ # to reduce checking time
  set.seed(123)
  GPD <- GParetoFamily(loc=20480, scale=7e4, shape=0.3)
  data <- r(GPD)(500)
  getCVaR(data,GPD,0.99)
  getVaR(data,GPD,0.99)
  getEL(data,GPD,5)
  getVaR(data,GPD,0.99, rob=FALSE)
  getEL(data,GPD,5, rob=FALSE)
  getCVaR(data,GPD,0.99, rob=FALSE)
  }
}
\concept{estimator}
\keyword{risk measure}
