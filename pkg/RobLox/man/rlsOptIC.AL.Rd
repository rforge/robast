\name{rlsOptIC.AL}
\alias{rlsOptIC.AL}
\title{Computation of the optimally robust IC for AL estimators}
\description{
  The function \code{rlsOptIC.AL} computes the optimally robust IC for 
  AL estimators in case of normal location with unknown scale and 
  (convex) contamination neighborhoods. The definition of 
  these estimators can be found in Section 8.2 of Kohl (2005).
}
\usage{
rlsOptIC.AL(r, mean = 0, sd = 1, A.loc.start = 1, a.sc.start = 0, 
            A.sc.start = 0.5, bUp = 1000, delta = 1e-6, itmax = 100, 
            check = FALSE, computeIC = TRUE)
}
\arguments{
  \item{r}{ non-negative real: neighborhood radius. }
  \item{mean}{ specified mean.}
  \item{sd}{ specified standard deviation.}
  \item{A.loc.start}{ positive real: starting value for 
    the standardizing constant of the location part. }
  \item{a.sc.start}{ real: starting value for centering
    constant of the scale part. }
  \item{A.sc.start}{ positive real: starting value for 
    the standardizing constant of the scale part. }
  \item{bUp}{ positive real: the upper end point of the 
    interval to be searched for the clipping bound b. }
  \item{delta}{ the desired accuracy (convergence tolerance). }
  \item{itmax}{ the maximum number of iterations. }
  \item{check}{ logical: should constraints be checked. }
  \item{computeIC}{ logical: should IC be computed. See details below. }
}
\details{The Lagrange multipliers contained in the expression
  of the optimally robust IC can be accessed via the
  accessor functions \code{cent}, \code{clip} and \code{stand}.
  If 'computeIC' is 'FALSE' only the Lagrange multipliers 'A', 'a', 
  and 'b' contained in the optimally robust IC are computed.
}
\value{
  If 'computeIC' is 'TRUE' an object of class \code{"ContIC"} is returned, 
  otherwise a list of Lagrane multipliers
  \item{A}{ standardizing matrix }
  \item{a}{ centering vector }
  \item{b}{ optimal clipping bound }
}
\references{ 
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[ROptEst]{ContIC-class}}, \code{\link{roblox}}}
\examples{
IC1 <- rlsOptIC.AL(r = 0.1, check = TRUE)
distrExOptions("ErelativeTolerance" = 1e-12)
checkIC(IC1)
distrExOptions("ErelativeTolerance" = .Machine$double.eps^0.25) # default
Risks(IC1)
cent(IC1)
clip(IC1)
stand(IC1)
plot(IC1)
infoPlot(IC1)

## one-step estimation
## see also: ?roblox
## 1. data: random sample
ind <- rbinom(100, size=1, prob=0.05) 
x <- rnorm(100, mean=0, sd=(1-ind) + ind*9)
mean(x)
sd(x)
median(x)
mad(x)

## 2. Kolmogorov(-Smirnov) minimum distance estimator
## -> we use it as initial estimate for one-step construction
(est0 <- ksEstimator(x=x, Norm()))

## 3. one-step estimation: radius known
IC1 <- rlsOptIC.AL(r = 0.5, mean = est0$mean, sd = est0$sd)
(est1 <- oneStepEstimator(x, IC1, est0))

## 4. one-step estimation: radius unknown
## take least favorable radius r = 0.579
## cf. Table 8.1 in Kohl(2005)
IC2 <- rlsOptIC.AL(r = 0.579, mean = est0$mean, sd = est0$sd)
(est2 <- oneStepEstimator(x, IC2, est0))
}
\concept{normal location and scale}
\concept{influence curve}
\keyword{robust}
