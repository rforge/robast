\name{roptest}
\alias{roptest}
\title{ Optimally robust estimation }
\description{
  Function to compute optimally robust estimates for L2-differentiable
  parametric families via k-step construction.
}
\usage{
roptest(x, L2Fam, eps, eps.lower, eps.upper, initial.est, 
        neighbor = ContNeighborhood(), risk = asMSE(), steps = 1, 
        distance = CvMDist, interval, par, ...)
}
\arguments{
  \item{x}{ sample }
  \item{L2Fam}{ object of class \code{"L2ParamFamily"} }
  \item{eps}{ positive real (0 < \code{eps} <= 0.5): amount of gross errors. 
          See details below. }
  \item{eps.lower}{ positive real (0 <= \code{eps.lower} <= \code{eps.upper}): 
        lower bound for the amount of gross errors. See details below. }
  \item{eps.upper}{ positive real (\code{eps.lower} <= \code{eps.upper} <= 0.5): 
        upper bound for the amount of gross errors. See details below. }
  \item{initial.est}{ initial estimate for unknown parameter. If missing 
        minimum distance estimator is computed. }
  \item{neighbor}{ object of class \code{"UncondNeighborhood"} }
  \item{risk}{ object of class \code{"RiskType"} }
  \item{steps}{ positive integer: number of steps used for k-steps construction }
  \item{distance}{ distance function }
  \item{interval}{ parameter interval for univariate parameters }
  \item{par}{ initial parameter value for multivariate parameters.
        If missing, the parameters of \code{L2Fam} are used. }
  \item{\dots}{ further arguments }
}
\details{
  Computes the optimally robust estimator for a given L2 differentiable
  parametric family. The computation uses a k-step construction with an 
  appropriate initial estimate. Valid candidates are e.g. Kolmogorov(-Smirnov) 
  or von Mises minimum distance estimators (default); cf. Rieder (1994) and Kohl (2005).

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
  Rieder et al. (2001, 2008), respectively Section 2.2 of Kohl (2005) is returned.
}
\value{Object of class \code{"kStepEstimate"}. }
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2008) The Costs of not Knowing
  the Radius. Statistical Methods and Applications \emph{17}(1) 13-40.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2001) The Costs of not Knowing
  the Radius. Submitted. Appeared as discussion paper Nr. 81. 
  SFB 373 (Quantification and Simulation of Economic Processes),
  Humboldt University, Berlin; also available under
  \url{www.uni-bayreuth.de/departments/math/org/mathe7/RIEDER/pubs/RR.pdf}
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{ \code{\link[RobLox]{roblox}}, 
          \code{\link[distrMod]{L2ParamFamily-class}}
          \code{\link[RobAStBase]{UncondNeighborhood-class}},
          \code{\link[distrMod]{RiskType-class}} }
\examples{
#############################
## 1. Binomial data
#############################
## generate a sample of contaminated data
ind <- rbinom(100, size=1, prob=0.05) 
x <- rbinom(100, size=25, prob=(1-ind)*0.25 + ind*0.9)

## ML-estimate
MLEstimator(x, BinomFamily(size = 25), interval = c(0, 1))

## compute optimally robust estimator (known contamination)
roptest(x, BinomFamily(size = 25), eps = 0.05, interval = c(0, 1), steps = 3)

## compute optimally robust estimator (unknown contamination)
roptest(x, BinomFamily(size = 25), eps.lower = 0, eps.upper = 0.1, interval = c(0, 1), steps = 3)


#############################
## 2. Poisson data
#############################
## Example: Rutherford-Geiger (1910); cf. Feller~(1968), Section VI.7 (a)
x <- c(rep(0, 57), rep(1, 203), rep(2, 383), rep(3, 525), rep(4, 532), 
       rep(5, 408), rep(6, 273), rep(7, 139), rep(8, 45), rep(9, 27), 
       rep(10, 10), rep(11, 4), rep(12, 0), rep(13, 1), rep(14, 1))

## ML-estimate
MLEstimator(x, PoisFamily(), interval = c(0, 10))

## compute optimally robust estimator (unknown contamination)
roptest(x, PoisFamily(), eps.upper = 0.05, interval = c(0, 10), steps = 3)

#############################
## 3. Normal (Gaussian) location and scale
#############################
## Generate a contaminated sample
ind <- rbinom(100, size=1, prob=0.05) 
x <- rnorm(100, mean=0, sd=(1-ind) + ind*9)

## ML-estimate
MLEstimator(x, NormLocationScaleFamily())

## compute optimally robust estimator (known contamination)
## takes some time
roptest(x, NormLocationScaleFamily(), eps = 0.05, steps = 3)

## compute optimally robust estimator (unknown contamination)
## takes some time
roptest(x, NormLocationScaleFamily(), eps.upper = 0.1, steps = 3)
}
\concept{k-step construction}
\concept{optimally robust estimation}
\keyword{robust}
