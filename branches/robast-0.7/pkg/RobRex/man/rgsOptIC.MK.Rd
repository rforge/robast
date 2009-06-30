\name{rgsOptIC.MK}
\alias{rgsOptIC.MK}

\title{Computation of the optimally robust IC for MK estimators}
\description{
  The function \code{rgsOptIC.MK} computes the optimally robust IC
  for MK estimators in case of linear regression with unknown
  scale and (convex) contamination neighborhoods where the
  regressor is random; confer Subsubsection 7.2.2.1 of Kohl (2005).
}
\usage{
rgsOptIC.MK(r, K, ggLo = 0.5, ggUp = 1, a1.start = -0.25, a3.start = 0.25, 
             B.start, bUp = 1000, delta = 1e-06, itmax = 1000, check = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{ non-negative real: neighborhood radius. }
  \item{K}{ object of class \code{"Distribution"}. }
  \item{ggLo}{ positive real: the lower end point of the interval 
    to be searched for \eqn{\gamma}{gamma}. }
  \item{ggUp}{ positive real: the upper end point of the interval 
    to be searched for \eqn{\gamma}{gamma}. }
  \item{a1.start}{ real: starting value for 
    Lagrange multiplier \eqn{\alpha_1}{alpha_1}. }
  \item{a3.start}{ real: starting value for 
    Lagrange multiplier \eqn{\alpha_3}{alpha_3}. }
  \item{B.start}{ symmetric matrix: starting value for 
    Lagrange multiplier B. }
  \item{bUp}{ positive real: the upper end point of the
    interval to be searched for b. }
  \item{delta}{ the desired accuracy (convergence tolerance). }
  \item{itmax}{ the maximum number of iterations. }
  \item{check}{ logical. Should constraints be checked. }
}
%\details{}
\value{Object of class \code{"IC"}}
\references{ 
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{IC-class}}
\examples{
## code takes some time
\dontrun{
K <- DiscreteDistribution(1:5) # = Unif({1,2,3,4,5})
IC1 <- rgsOptIC.MK(r = 0.1, K = K)
checkIC(IC1)
Risks(IC1)
}
}
\concept{regression and scale}
\concept{influence curve}
\keyword{robust}
