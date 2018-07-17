\name{showdown}
\Rdversion{1.1}
\alias{showdown}
\title{Estimator Showdown by Monte-Carlo Study.}
\description{
  The function \code{showdown} can be used to perform Monte-Carlo studies 
  comparing a competitor with rmx estimators in case of normal location and scale. 
  In addition, maximum likelihood (ML) estimators (mean and sd) and median and 
  MAD are computed. The comparison is based on the empirical MSE.
}
\usage{
showdown(n, M, eps, contD, seed = 123, estfun, estMean, estSd,
         eps.lower = 0, eps.upper = 0.05, steps = 3L, fsCor = TRUE, 
         plot1 = FALSE, plot2 = FALSE, plot3 = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{integer; sample size, should be at least 3.}
  \item{M}{integer; Monte-Carlo replications.}
  \item{eps}{amount of contamination in [0, 0.5].}
  \item{contD}{object of class \code{"UnivariateDistribution"}; contaminating distribution.}
  \item{seed}{random seed.}
  \item{estfun}{function to compute location and scale estimator; see details below.}
  \item{estMean}{function to compute location estimator; see details below.}
  \item{estSd}{function to compute scale estimator; see details below.}
  \item{eps.lower}{used by rmx estimator.}
  \item{eps.upper}{used by rmx estimator.}
  \item{steps}{integer; steps used for estimator construction.}
  \item{fsCor}{logical; use finite-sample correction.}
  \item{plot1}{logical; plot cdf of ideal and real distribution.}
  \item{plot2}{logical; plot 20 (or M if M < 20) randomly selected samples.}
  \item{plot3}{logical; generate boxplots of the results.}
}
\details{
Normal location and scale with mean = 0 and sd = 1 is used as ideal model (without
restriction due to equivariance).

Since there is no estimator which yields reliable results if 50 percent or more of the
observations are contaminated, we use a modification where we re-simulate all samples
including at least 50 percent contaminated data.

If \code{estfun} is specified it has to compute and return a location and scale estimate
(vector of length 2). One can also specify the location and scale estimator separately 
by using \code{estMean} and \code{estSd} where \code{estMean} computes and returns
the location estimate and \code{estSd} the scale estimate.

We use funtion \code{\link{rowRoblox}} for the computation of the rmx estimator.
}
\value{Data.frame including empirical MSE (standardized by sample size n) and
relMSE with respect to the rmx estimator.
}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2008) The Costs of not Knowing
  the Radius. Statistical Methods and Applications \emph{17}(1) 13-40.
  Extended version: \url{http://r-kurs.de/RRlong.pdf}
  
  M. Kohl, P. Ruckdeschel, and H. Rieder (2010). Infinitesimally Robust Estimation 
  in General Smoothly Parametrized Models. \emph{Statistical Methods and Application}, 
  \bold{19}(3):333-354. 
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{rowRoblox}}}
\examples{
library(MASS)
## compare with Huber's Proposal 2
showdown(n = 20, M = 100, eps = 0.02, contD = Norm(mean = 3, sd = 3), 
         estfun = function(x){ unlist(hubers(x)) },
         plot1 = TRUE, plot2 = TRUE, plot3 = TRUE)

## compare with Huber M estimator with MAD scale
showdown(n = 20, M = 100, eps = 0.02, contD = Norm(mean = 3, sd = 3), 
         estfun = function(x){ unlist(huber(x)) },
         plot1 = TRUE, plot2 = TRUE, plot3 = TRUE)
}
\concept{Monte-Carlo study}
\keyword{robust}
