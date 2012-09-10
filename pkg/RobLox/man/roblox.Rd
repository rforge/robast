\name{roblox}
\alias{roblox}
\title{Optimally robust estimator for location and/or scale}
\description{
  The function \code{roblox} computes the optimally robust estimator
  and corresponding IC for normal location und/or scale and 
  (convex) contamination neighborhoods. The definition of 
  these estimators can be found in Rieder (1994) or Kohl (2005),
  respectively.
}
\usage{
roblox(x, mean, sd, eps, eps.lower, eps.upper, initial.est, k = 1L, 
       fsCor = TRUE, returnIC = FALSE, mad0 = 1e-4, na.rm = TRUE)
}
\arguments{
  \item{x}{ vector \code{x} of data values, may also be a matrix or data.frame
        with one row, respectively one column/(numeric) variable. }
  \item{mean}{ specified mean. }
  \item{sd}{ specified standard deviation which has to be positive. }
  \item{eps}{ positive real (0 < \code{eps} <= 0.5): amount of gross errors. 
        See details below. }
  \item{eps.lower}{ positive real (0 <= \code{eps.lower} <= \code{eps.upper}): 
        lower bound for the amount of gross errors. See details below. }
  \item{eps.upper}{ positive real (\code{eps.lower} <= \code{eps.upper} <= 0.5): 
        upper bound for the amount of gross errors. See details below. }
  \item{initial.est}{ initial estimate for \code{mean} and/or \code{sd}. If missing 
        median and/or MAD are used. }
  \item{k}{ positive integer. k-step is used to compute the optimally robust estimator. }
  \item{fsCor}{ logical: perform finite-sample correction. See function \code{\link{finiteSampleCorrection}}. }
  \item{returnIC}{ logical: should IC be returned. See details below. }
  \item{mad0}{ scale estimate used if computed MAD is equal to zero}
  \item{na.rm}{logical: if  \code{TRUE}, the estimator is evaluated at \code{complete.cases(x)}.}
}
\details{
  Computes the optimally robust estimator for location with scale specified,
  scale with location specified, or both if neither is specified. The computation
  uses a k-step construction with an appropriate initial estimate for location
  or scale or location and scale, respectively. Valid candidates are e.g. 
  median and/or MAD (default) as well as Kolmogorov(-Smirnov) or von Mises minimum 
  distance estimators; cf. Rieder (1994) and Kohl (2005).

  If the amount of gross errors (contamination) is known, it can be 
  specified by \code{eps}. The radius of the corresponding infinitesimal 
  contamination neighborhood is obtained by multiplying \code{eps} 
  by the square root of the sample size. 

  If the amount of gross errors (contamination) is unknown, try to find a 
  rough estimate for the amount of gross errors, such that it lies 
  between \code{eps.lower} and \code{eps.upper}.

  In case \code{eps.lower} is specified and \code{eps.upper} is missing, 
  \code{eps.upper} is set to 0.5. In case \code{eps.upper} is specified and
  \code{eps.lower} is missing, \code{eps.lower} is set to 0.

  If neither \code{eps} nor \code{eps.lower} and/or \code{eps.upper} is 
  specified, \code{eps.lower} and \code{eps.upper} are set to 0 and 0.5, 
  respectively.

  If \code{eps} is missing, the radius-minimax estimator in sense of 
  Rieder et al. (2008), respectively Section 2.2 of Kohl (2005) is returned.
  
  In case of location, respectively scale one additionally has to specify
  \code{sd}, respectively \code{mean} where \code{sd} and \code{mean} have
  to be a single number.

  For sample size <= 2, median and/or MAD are used for estimation.
  
  If \code{eps = 0}, mean and/or sd are computed. In this situation it's better
  to use function \code{\link[distrMod]{MLEstimator}}.
}
\value{Object of class \code{"kStepEstimate"}. }
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2008) The Costs of not Knowing
  the Radius. Statistical Methods and Applications \emph{17}(1) 13-40.
  Extended version: \url{http://www.stamats.de/RRlong.pdf}
  
  M. Kohl, P. Ruckdeschel, and H. Rieder (2010). Infinitesimally Robust Estimation 
  in General Smoothly Parametrized Models. \emph{Statistical Methods and Application}, 
  \bold{19}(3):333-354. 
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[RobAStBase]{ContIC-class}}, \code{\link{rlOptIC}}, 
         \code{\link{rsOptIC}}, \code{\link{rlsOptIC.AL}},
         \code{\link[RobAStBase]{kStepEstimate-class}},
         \code{\link[ROptEst]{roptest}} }
\examples{
ind <- rbinom(100, size=1, prob=0.05) 
x <- rnorm(100, mean=ind*3, sd=(1-ind) + ind*9)

## amount of gross errors known
res1 <- roblox(x, eps = 0.05, returnIC = TRUE)
estimate(res1)
confint(res1)
confint(res1, method = symmetricBias())
pIC(res1)
checkIC(pIC(res1))
Risks(pIC(res1))
Infos(pIC(res1))
plot(pIC(res1))
infoPlot(pIC(res1))

## amount of gross errors unknown
res2 <- roblox(x, eps.lower = 0.01, eps.upper = 0.1, returnIC = TRUE)
estimate(res2)
confint(res2)
confint(res2, method = symmetricBias())
pIC(res2)
checkIC(pIC(res2))
Risks(pIC(res2))
Infos(pIC(res2))
plot(pIC(res2))
infoPlot(pIC(res2))

## estimator comparison
# classical optimal (non-robust)
c(mean(x), sd(x))

# most robust
c(median(x), mad(x))

# optimally robust (amount of gross errors known)
estimate(res1)

# optimally robust (amount of gross errors unknown)
estimate(res2)

# Kolmogorov(-Smirnov) minimum distance estimator (robust)
(ks.est <- MDEstimator(x, ParamFamily = NormLocationScaleFamily()))

# optimally robust (amount of gross errors known)
roblox(x, eps = 0.05, initial.est = estimate(ks.est))

# Cramer von Mises minimum distance estimator (robust)
(CvM.est <- MDEstimator(x, ParamFamily = NormLocationScaleFamily(), distance = CvMDist))

# optimally robust (amount of gross errors known)
roblox(x, eps = 0.05, initial.est = estimate(CvM.est))
}
\concept{normal location}
\concept{normal scale}
\concept{normal location and scale}
\concept{influence curve}
\keyword{robust}
