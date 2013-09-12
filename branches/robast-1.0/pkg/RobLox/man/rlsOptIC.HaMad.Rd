\name{rlsOptIC.HaMad}
\alias{rlsOptIC.HaMad}
\title{Computation of the optimally robust IC for HuMad estimators}
\description{
  The function \code{rlsOptIC.HuMad} computes the optimally robust IC for 
  HuMad estimators in case of normal location with unknown scale and 
  (convex) contamination neighborhoods. These estimators were 
  considered in Andrews et al. (1972). A definition of these estimators 
  can also be found in Subsection 8.5.2 of Kohl (2005).
}
\usage{
rlsOptIC.HaMad(r, a.start = 0.25, b.start = 2.5, c.start = 5, 
               delta = 1e-06, MAX = 100)
}
\arguments{
  \item{r}{ non-negative real: neighborhood radius. }
  \item{a.start}{ positive real: starting value for a. }
  \item{b.start}{ positive real: starting value for b. }
  \item{c.start}{ positive real: starting value for c. }
  \item{delta}{ the desired accuracy (convergence tolerance). }
  \item{MAX}{ if a or b or c are beyond the admitted values, 
    \code{MAX} is returned. }
}
\details{
  The computation of the optimally robust IC for HaMad estimators
  is based on \code{optim} where \code{MAX} is used to 
  control the constraints on a, b and c. The optimal values of 
  the tuning constants a, b, and c can be read off 
  from the slot \code{Infos} of the resulting IC.
}
\value{Object of class \code{"IC"}}
\references{ 
  Andrews, D.F., Bickel, P.J., Hampel, F.R., Huber, P.J.,
  Rogers, W.H. and Tukey, J.W. (1972) \emph{Robust estimates of location}. 
  Princeton University Press.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[RobAStBase]{IC-class}}}
\examples{
IC1 <- rlsOptIC.HaMad(r = 0.1)
checkIC(IC1)
Risks(IC1)
Infos(IC1)
plot(IC1)
infoPlot(IC1)
}
\concept{normal location and scale}
\concept{influence curve}
\keyword{robust}
