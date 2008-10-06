\name{rgsOptIC.Ms}
\alias{rgsOptIC.Ms}

\title{Computation of the optimally robust IC for Ms estimators}
\description{
  The function \code{rgsOptIC.Ms} computes the optimally robust 
  conditionally centered IC for Ms estimators in case of linear 
  regression with unknown scale and average conditional (convex) 
  contamination neighborhoods where the regressor is random;
  confer Subsection 7.3.2 of Kohl (2005).
}
\usage{
rgsOptIC.Ms(r, K, a1.x.start, a3.start = 0.25, b.sc.start = 1.5, 
             bUp = 1000, ggLo = 0.5, ggUp = 1, delta = 1e-06, 
             itmax = 1000, check = FALSE)
}
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
  \item{b.sc.start}{ positive real: starting value for 
    the clipping bound \eqn{b_{\rm sc}}{b_sc}. }
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
K <- DiscreteDistribution(1:5) # = Unif({1,2,3,4,5})
IC1 <- rgsOptIC.Ms(r = 0.1, K = K)
checkIC(IC1)
Risks(IC1)
}
\concept{regression and scale}
\concept{influence curve}
\keyword{robust}
