\name{rlsOptIC.HuMad}
\alias{rlsOptIC.HuMad}
\title{Computation of the optimally robust IC for HuMad estimators}
\description{
  The function \code{rlsOptIC.HuMad} computes the optimally robust IC for 
  HuMad estimators in case of normal location with unknown scale and 
  (convex) contamination neighborhoods. These estimators were 
  proposed by Andrews et al. (1972), p. 12. A definition of these 
  estimators can also be found in Subsection 8.5.1 of Kohl (2005).
}
\usage{
rlsOptIC.HuMad(r, kUp = 2.5, delta = 1e-06)
}
\arguments{
  \item{r}{ non-negative real: neighborhood radius. }
  \item{kUp}{ positive real: the upper end point of the interval 
    to be searched for k. }
  \item{delta}{ the desired accuracy (convergence tolerance). }
}
\details{The optimal value of the tuning constant k can be read off 
  from the slot \code{Infos} of the resulting IC.}
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
IC1 <- rlsOptIC.HuMad(r = 0.1)
checkIC(IC1)
Risks(IC1)
Infos(IC1)
plot(IC1)
infoPlot(IC1)
}
\concept{normal location and scale}
\concept{influence curve}
\keyword{robust}
