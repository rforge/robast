\name{rlsOptIC.Ha4}
\alias{rlsOptIC.Ha4}
\title{Computation of the optimally robust IC for Ha4 estimators}
\description{
  The function \code{rlsOptIC.Ha4} computes the optimally robust IC for 
  Ha4 estimators in case of normal location with unknown scale and 
  (convex) contamination neighborhoods. The definition of 
  these estimators can be found in Subsection 8.5.2 of Kohl (2005).
}
\usage{
rlsOptIC.Ha4(r, a.start = 0.25, b.start = 2.5, c.start = 5, 
             k.start = 1, delta = 1e-06, MAX = 100)
}
\arguments{
  \item{r}{ non-negative real: neighborhood radius. }
  \item{a.start}{ positive real: starting value for a. }
  \item{b.start}{ positive real: starting value for b. }
  \item{c.start}{ positive real: starting value for c. }
  \item{k.start}{ positive real: starting value for k. }
  \item{delta}{ the desired accuracy (convergence tolerance). }
  \item{MAX}{ if a or b or c or k are beyond the admitted values, 
    \code{MAX} is returned. }
}
\details{
  The computation of the optimally robust IC for Ha4 estimators
  is based on \code{optim} where \code{MAX} is used to 
  control the constraints on a, b, c and k. The optimal values of 
  the tuning constants a, b, c and k can be read off 
  from the slot \code{Infos} of the resulting IC.
}
\value{Object of class \code{"IC"}}
\references{ 
  Marazzi, A. (1993) \emph{Algorithms, routines, and S functions 
  for robust statistics}. Wadsworth and Brooks / Cole.
  
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[RobAStBase]{IC-class}}}
\examples{
IC1 <- rlsOptIC.Ha4(r = 0.1)
checkIC(IC1)
Risks(IC1)
Infos(IC1)
plot(IC1)
infoPlot(IC1)
}
\concept{normal location and scale}
\concept{influence curve}
\keyword{robust}
