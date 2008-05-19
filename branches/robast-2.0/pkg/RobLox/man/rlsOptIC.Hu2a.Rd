\name{rlsOptIC.Hu2a}
\alias{rlsOptIC.Hu2a}
\title{Computation of the optimally robust IC for Hu2a estimators}
\description{
  The function \code{rlsOptIC.Hu2a} computes the optimally robust IC for 
  Hu2a estimators in case of normal location with unknown scale and 
  (convex) contamination neighborhoods. These estimators are a
  simple modification of Huber (1964), Proposal 2 where we, in addition,
  admit a clipping from below. The definition of 
  these estimators can be found in Subsection 8.5.1 of Kohl (2005).
}
\usage{
rlsOptIC.Hu2a(r, k1.start = 0.25, k2.start = 2.5, delta = 1e-06, MAX = 100)
}
\arguments{
  \item{r}{ non-negative real: neighborhood radius. }
  \item{k1.start}{ positive real: starting value for k1. }
  \item{k2.start}{ positive real: starting value for k2. }
  \item{delta}{ the desired accuracy (convergence tolerance). }
  \item{MAX}{ if k1 or k2 are beyond the admitted values, 
    \code{MAX} is returned. }
}
\details{
  The computation of the optimally robust IC for Hu2a estimators
  is based on \code{optim} where \code{MAX} is used to 
  control the constraints on k1 and k2. The optimal values of 
  the tuning constants k1 and k2 can be read off 
  from the slot \code{Infos} of the resulting IC.
}
\value{Object of class \code{"IC"}}
\references{ 
  Huber, P.J. (1964) Robust estimation of a location parameter. 
  Ann. Math. Stat. \bold{35}: 73--101.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[ROptEst]{IC-class}}}
\examples{
IC1 <- rlsOptIC.Hu2a(r = 0.1)
checkIC(IC1)
Risks(IC1)
Infos(IC1)
plot(IC1)
infoPlot(IC1)
}
\concept{normal location and scale}
\concept{influence curve}
\keyword{robust}
