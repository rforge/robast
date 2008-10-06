\name{rlsOptIC.An2}
\alias{rlsOptIC.An2}
\title{Computation of the optimally robust IC for An2 estimators}
\description{
  The function \code{rlsOptIC.An2} computes the optimally robust IC for 
  An2 estimators in case of normal location with unknown scale and 
  (convex) contamination neighborhoods. The definition of 
  these estimators can be found in Subsection 8.5.3 of Kohl (2005).
}
\usage{
rlsOptIC.An2(r, a.start = 1.5, k.start = 1.5, delta = 1e-06, MAX = 100)
}
\arguments{
  \item{r}{ non-negative real: neighborhood radius. }
  \item{a.start}{ positive real: starting value for a. }
  \item{k.start}{ positive real: starting value for k. }
  \item{delta}{ the desired accuracy (convergence tolerance). }
  \item{MAX}{ if a or k are beyond the admitted values, 
    \code{MAX} is returned. }
}
\details{
  The computation of the optimally robust IC for An2 estimators
  is based on \code{optim} where \code{MAX} is used to 
  control the constraints on a and k. The optimal values of the 
  tuning constants a and k can be read off from the slot 
  \code{Infos} of the resulting IC.
}
%\details{}
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
IC1 <- rlsOptIC.An2(r = 0.1)
checkIC(IC1)
Risks(IC1)
Infos(IC1)
plot(IC1)
infoPlot(IC1)
}
\concept{normal location and scale}
\concept{influence curve}
\keyword{robust}
