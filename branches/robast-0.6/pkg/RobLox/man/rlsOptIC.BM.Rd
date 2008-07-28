\name{rlsOptIC.BM}
\alias{rlsOptIC.BM}
\title{Computation of the optimally robust IC for BM estimators}
\description{
  The function \code{rlsOptIC.BM} computes the optimally robust IC for 
  BM estimators in case of normal location with unknown scale and 
  (convex) contamination neighborhoods. These estimators were proposed 
  by Bednarski and Mueller (2001). A definition of these 
  estimators can also be found in Section 8.4 of Kohl (2005).
}
\usage{
rlsOptIC.BM(r, bL.start = 2, bS.start = 1.5, delta = 1e-06, MAX = 100)
}
\arguments{
  \item{r}{ non-negative real: neighborhood radius. }
  \item{bL.start}{ positive real: starting value for \eqn{b_{\rm loc}}{b_loc}. }
  \item{bS.start}{ positive real: starting value for \eqn{b_{{\rm sc},0}}{b_sc,0}. }
  \item{delta}{ the desired accuracy (convergence tolerance). }
  \item{MAX}{ if \eqn{b_{\rm loc}}{b_loc} or \eqn{b_{{\rm sc},0}}{b_sc,0} 
    are beyond the admitted values, \code{MAX} is returned. }
}
\details{
  The computation of the optimally robust IC for BM estimators
  is based on \code{optim} where \code{MAX} is used to 
  control the constraints on \eqn{b_{\rm loc}}{b_loc} 
  and \eqn{b_{{\rm sc},0}}{b_sc,0}. The optimal values of the  
  tuning constants \eqn{b_{\rm loc}}{b_loc}, \eqn{b_{{\rm sc},0}}{b_sc,0}, 
  \eqn{\alpha}{alpha} and \eqn{\gamma}{gamma} can be read off 
  from the slot \code{Infos} of the resulting IC.
}
\value{Object of class \code{"IC"}}
\references{ 
  Bednarski, T and Mueller, C.H. (2001) Optimal bounded influence
  regression and scale M-estimators in the context of experimental
  design. Statistics, \bold{35}(4): 349--369.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[RobAStBase]{IC-class}}}
\examples{
IC1 <- rlsOptIC.BM(r = 0.1)
checkIC(IC1)
Risks(IC1)
Infos(IC1)
plot(IC1)
infoPlot(IC1)
}
\concept{normal location and scale}
\concept{influence curve}
\keyword{robust}
