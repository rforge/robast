\name{rgsOptIC.Mc}
\alias{rgsOptIC.Mc}

\title{Computation of the optimally robust IC for Mc estimators}
\description{
  The function \code{rgsOptIC.Mc} computes the optimally robust 
  conditionally centered IC for Mc estimators in case of linear 
  regression with unknown scale and average conditional (convex) 
  contamination neighborhoods where the regressor is random; 
  confer Subsubsection 7.2.2.2 of Kohl (2005).
}
\usage{
rgsOptIC.Mc(r, K, ggLo = 0.5, ggUp = 1, a1.x.start, a3.start = 0.25, 
             bUp = 1000, delta = 1e-05, itmax = 1000, check = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{ non-negative real: neighborhood radius. }
  \item{K}{ object of class \code{"DiscreteDistribution"} }
  \item{ggLo}{ positive real: the lower end point of the interval 
    to be searched for \eqn{\gamma}{gamma}. }
  \item{ggUp}{ positive real: the upper end point of the interval 
    to be searched for \eqn{\gamma}{gamma}. }
  \item{a1.x.start}{ real: starting value for 
    the Lagrange multiplier function \eqn{\alpha_1(x)}{alpha_1(x)}. }
  \item{a3.start}{ real: starting value for 
    Lagrange multiplier \eqn{\alpha_3}{alpha_3}. }
  \item{bUp}{ positive real: the upper end point of the
    interval to be searched for b. }
  \item{delta}{ the desired accuracy (convergence tolerance). }
  \item{itmax}{ the maximum number of iterations. }
  \item{check}{ logical. Should constraints be checked. }
}
%\details{}
\value{Object of class \code{"CondIC"}}
\references{ 
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{CondIC-class}}
\examples{
## code takes some time
\dontrun{
K <- DiscreteDistribution(1:5) # = Unif({1,2,3,4,5})
IC1 <- rgsOptIC.Mc(r = 0.1, K = K)
checkIC(IC1)
Risks(IC1)
}
}
\concept{regression and scale}
\concept{influence curve}
\keyword{robust}
