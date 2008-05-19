\name{rlsOptIC.Hu2}
\alias{rlsOptIC.Hu2}
\title{Computation of the optimally robust IC for Hu2 estimators}
\description{
  The function \code{rlsOptIC.Hu2} computes the optimally robust IC for 
  Hu2 estimators in case of normal location with unknown scale and 
  (convex) contamination neighborhoods. These estimators were
  proposed in Example 6.4.1 of Huber (1981). A definition of these 
  estimators can also be found in Subsection 8.5.1 of Kohl (2005).
}
\usage{
rlsOptIC.Hu2(r, k.start = 1.5, c.start = 1.5, delta = 1e-06, MAX = 100)
}
\arguments{
  \item{r}{ non-negative real: neighborhood radius. }
  \item{k.start}{ positive real: starting value for k. }
  \item{c.start}{ positive real: starting value for c. }
  \item{delta}{ the desired accuracy (convergence tolerance). }
  \item{MAX}{ if k1 or k2 are beyond the admitted values, 
    \code{MAX} is returned. }
}
\details{
  The computation of the optimally robust IC for Hu2 estimators
  is based on \code{optim} where \code{MAX} is used to 
  control the constraints on k and c. The optimal values of 
  the tuning constants k and c can be read off 
  from the slot \code{Infos} of the resulting IC.
}
\value{Object of class \code{"IC"}}
\references{ 
  Huber, P.J. (1981) \emph{Robust Statistics}. New York: Wiley.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[ROptEst]{IC-class}}}
\examples{
IC1 <- rlsOptIC.Hu2(r = 0.1)
checkIC(IC1)
Risks(IC1)
Infos(IC1)
plot(IC1)
infoPlot(IC1)
}
\concept{normal location and scale}
\concept{influence curve}
\keyword{robust}
