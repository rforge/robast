\name{rlsOptIC.Hu3}
\alias{rlsOptIC.Hu3}
\title{Computation of the optimally robust IC for Hu3 estimators}
\description{
  The function \code{rlsOptIC.Hu3} computes the optimally robust IC for 
  Hu3 estimators in case of normal location with unknown scale and 
  (convex) contamination neighborhoods. The definition of 
  these estimators can be found in Subsection 8.5.1 of Kohl (2005).
}
\usage{
rlsOptIC.Hu3(r, k.start = 1, c1.start = 0.1, c2.start = 0.5, 
             delta = 1e-06, MAX = 100)
}
\arguments{
  \item{r}{ non-negative real: neighborhood radius. }
  \item{k.start}{ positive real: starting value for k. }
  \item{c1.start}{ positive real: starting value for c1. }
  \item{c2.start}{ positive real: starting value for c2. }
  \item{delta}{ the desired accuracy (convergence tolerance). }
  \item{MAX}{ if k or c1 or c2 are beyond the admitted values, 
    \code{MAX} is returned. }
}
\details{
  The computation of the optimally robust IC for Hu2 estimators
  is based on \code{optim} where \code{MAX} is used to
  control the constraints on k, c1 and c2. The optimal values of
  the tuning constants k, c1 and c2 can be read off
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
\seealso{\code{\link[RobAStBase]{IC-class}}}
\examples{
IC1 <- rlsOptIC.Hu3(r = 0.1)
checkIC(IC1)
Risks(IC1)
Infos(IC1)
plot(IC1)
infoPlot(IC1)
}
\concept{normal location and scale}
\concept{influence curve}
\keyword{robust}
