\name{rgsOptIC.ALc}
\alias{rgsOptIC.ALc}

\title{Computation of the optimally robust IC for AL estimators}
\description{
  The function \code{rgsOptIC.ALc} computes the optimally robust 
  conditionally centered IC for AL estimators in case of linear 
  regression with unknown scale and average conditional 
  (convex) contamination neighborhoods where the regressor is 
  random; confer Subsubsection 7.2.1.2 of Kohl (2005).
}
\usage{
rgsOptIC.ALc(r, K, theta, scale = 1, A.rg.start, a.sc.start, A.sc.start = 0.5, 
             bUp = 1000, delta = 1e-06, itmax = 50, check = FALSE)
}
\arguments{
  \item{r}{ non-negative real: neighborhood radius. }
  \item{K}{ object of class \code{"DiscreteDistribution"} or
    object of class \code{"DisreteMVDistribution"}. }
  \item{theta}{ specified regression parameter.}
  \item{scale}{ specified error scale.}
  \item{A.rg.start}{ positive definite and symmetric matrix:
    starting value for the standardizing matrix of the
    regression part. }
  \item{a.sc.start}{ real vector: starting values for centering
    function of the scale part. }
  \item{A.sc.start}{ positive real: starting value for 
    the standardizing constant of the scale part. }
  \item{bUp}{ positive real: the upper end point of the
    interval to be searched for b. }
  \item{delta}{ the desired accuracy (convergence tolerance). }
  \item{itmax}{ the maximum number of iterations. }
  \item{check}{ logical. Should constraints be checked. }
}
\details{
  If \code{theta} is missing, it is set to 0.
  If \code{A.rg.start} is missing, the inverse of the 
  second moment matrix of \code{K} is used. In case 
  \code{a.sc.start} is missing, it is set to a null 
  vector with length of the support of \code{K}. }
\value{Object of class \code{"Av1CondContIC"}}
\references{ 
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{Av1CondContIC-class}}
\examples{
K <- DiscreteDistribution(1:5) # = Unif({1,2,3,4,5})
IC1 <- rgsOptIC.ALc(r = 0.1, K = K)
checkIC(IC1)
Risks(IC1)
}
\concept{regression and scale}
\concept{influence curve}
\keyword{robust}
