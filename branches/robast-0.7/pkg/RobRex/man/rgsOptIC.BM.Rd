\name{rgsOptIC.BM}
\alias{rgsOptIC.BM}

\title{Computation of the optimally robust IC for BM estimators}
\description{
  The function \code{rgsOptIC.BM} computes the optimally robust IC
  for BM estimators in case of linear regression with unknown
  scale and (convex) contamination neighborhoods where the
  regressor is random. These estimators were proposed 
  by Bednarski and Mueller (2001); confer also 
  Subsection 7.3.3 of Kohl (2005).
}
\usage{
rgsOptIC.BM(r, K, b.rg.start = 2.5, b.sc.0.x.start, delta = 1e-06, 
            MAX = 100, itmax = 1000)
}
\arguments{
  \item{r}{ non-negative real: neighborhood radius. }
  \item{K}{ object of class \code{"DiscreteDistribution"} }
  \item{b.rg.start}{ positive real: starting value for \eqn{b_{\rm rg}}{b_rg}. }
  \item{b.sc.0.x.start}{ positive real: starting value for \eqn{b_{{\rm sc},0,x}}{b_sc,0,x}. }
  \item{delta}{ the desired accuracy (convergence tolerance). }
  \item{itmax}{ the maximum number of iterations. }
  \item{MAX}{ if \eqn{b_{\rm loc}}{b_loc} or \eqn{b_{{\rm sc},0}}{b_sc,0} 
    are beyond the admitted values, \code{MAX} is returned. }
}
\details{
  The computation of the optimally robust IC for BM estimators
  is based on \code{optim} where \code{MAX} is used to 
  control the constraints on \eqn{b_{\rm rg}}{b_rg} 
  and \eqn{b_{{\rm sc},0,x}}{b_sc,0,x}.
}
\value{Object of class \code{"CondIC"}}
\references{ 
  Bednarski, T and Mueller, C.H. (2001) Optimal bounded influence
  regression and scale M-estimators in the context of experimental
  design. Statistics, \bold{35}(4): 349--369.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{CondIC-class}}
\examples{
K <- DiscreteDistribution(1:5) # = Unif({1,2,3,4,5})
IC1 <- rgsOptIC.BM(r = 0.1, K = K)
checkIC(IC1)
Risks(IC1)
}
\concept{regression and scale}
\concept{influence curve}
\keyword{robust}
