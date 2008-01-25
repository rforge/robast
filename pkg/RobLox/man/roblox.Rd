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
roblox(x, mean, sd, eps, eps.lower, eps.upper, initial.est = "ksMD", 
       tol = 1e-06, A.loc.start = 1, a.sc.start = 0, A.sc.start = 0.5, 
       bUp = 1000, itmax = 100, returnIC = FALSE)
}
\arguments{
  \item{x}{ vector \code{x} of data values }
  \item{mean}{ specified mean.}
  \item{sd}{ specified standard deviation.}
  \item{eps}{ positive real (0 < \code{eps} <= 0.5): amount of gross errors. 
        See details below. }
  \item{eps.lower}{ positive real (0 <= \code{eps.lower} <= \code{eps.upper}): 
        lower bound for the amount of gross errors. See details below }
  \item{eps.upper}{ positive real (\code{eps.lower} <= \code{eps.upper} <= 0.5): 
        upper bound for the amount of gross errors. See details below }
  \item{initial.est}{ indicates which initial estimate should be used. Must be 
        one of \code{"ksMD"} or \code{"med"}. In case \code{"ksMD"} the 
        Kolmogorov(-Smirnov) minimum distance estimator and in case \code{"med"}
        median and/or mad are computed, respectively. }
  \item{tol}{ the desired accuracy (convergence tolerance). }
  \item{A.loc.start}{ positive real: starting value for 
    the standardizing constant of the location part. }
  \item{a.sc.start}{ real: starting value for centering
    constant of the scale part. }
  \item{A.sc.start}{ positive real: starting value for 
    the standardizing constant of the scale part. }
  \item{bUp}{ positive real: the upper end point of the 
    interval to be searched for the clipping bound b. }
  \item{itmax}{ the maximum number of iterations. }
  \item{returnIC}{ logical: should IC be returned. See details below. }
}
\details{
  Computes the optimally robust estimator for location with scale specified,
  scale with location specified, or both if neither is specified.
  
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
  Rieder et al. (2001), respectively Section 2.2 of Kohl (2005) is returned.
}
\value{
  list of location and scale estimates
  \item{mean }{Description of 'comp1'}
  \item{sd }{Description of 'comp2'}
  if 'returnIC' is 'TRUE' the list also contains
  \item{optIC}{ object of class \code{"ContIC"}; optimally robust IC }
}
\references{ 
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2001) The Costs of not Knowing
  the Radius. Submitted. Appeared as discussion paper Nr. 81. 
  SFB 373 (Quantification and Simulation of Economic Processes),
  Humboldt University, Berlin; also available under
  \url{www.uni-bayreuth.de/departments/math/org/mathe7/RIEDER/pubs/RR.pdf}

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[ROptEst]{ContIC-class}}, \code{\link{rlOptIC}}, 
         \code{\link{rsOptIC}}, \code{\link{rlsOptIC.AL}}}
\examples{
ind <- rbinom(100, size=1, prob=0.05) 
x <- rnorm(100, mean=ind*3, sd=(1-ind) + ind*9)

## amount of gross errors known
res1 <- roblox(x, eps = 0.05, returnIC = TRUE)
res1$mean
res1$sd
res1$optIC
checkIC(res1$optIC)
Risks(res1$optIC)
Infos(res1$optIC)
plot(res1$optIC)
infoPlot(res1$optIC)

## amount of gross errors unknown
res2 <- roblox(x, eps.lower = 0.01, eps.upper = 0.1, returnIC = TRUE)
res2$mean
res2$sd
res2$optIC
checkIC(res2$optIC)
Risks(res2$optIC)
Infos(res2$optIC)
plot(res2$optIC)
infoPlot(res2$optIC)

## estimator comparison
# classical optimal (non-robust)
c(mean(x), sd(x))
# most robust
c(median(x), mad(x))
# Kolmogorov(-Smirnov) minimum distance estimator (robust)
ksEstimator(x, distribution = Norm())
# optimally robust (amount of gross errors known)
c(res1$mean, res1$sd)
# optimally robust (amount of gross errors unknown)
c(res2$mean, res2$sd)

}
\concept{normal location}
\concept{normal scale}
\concept{normal location and scale}
\concept{influence curve}
\keyword{robust}
