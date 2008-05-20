\name{rlsOptIC.Hu1}
\alias{rlsOptIC.Hu1}
\title{Computation of the optimally robust IC for Hu1 estimators}
\description{
  The function \code{rlsOptIC.Hu1} computes the optimally robust IC for 
  Hu1 estimators in case of normal location with unknown scale and 
  (convex) contamination neighborhoods. These estimators were 
  proposed by Huber (1964), Proposal 2. A definition of these 
  estimators can also be found in Subsection 8.5.1 of Kohl (2005).
}
\usage{
rlsOptIC.Hu1(r, kUp = 2.5, delta = 1e-06)
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
  Huber, P.J. (1964) Robust estimation of a location parameter. 
  Ann. Math. Stat. \bold{35}: 73--101.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[ROptEst]{IC-class}}}
\examples{
IC1 <- rlsOptIC.Hu1(r = 0.1)
checkIC(IC1)
Risks(IC1)
Infos(IC1)
plot(IC1)
infoPlot(IC1)
}
\concept{normal location and scale}
\concept{influence curve}
\keyword{robust}
