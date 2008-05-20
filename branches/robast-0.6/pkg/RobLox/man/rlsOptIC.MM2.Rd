\name{rlsOptIC.MM2}
\alias{rlsOptIC.MM2}
\title{Computation of the optimally robust IC for MM2 estimators}
\description{
  The function \code{rlsOptIC.MM2} computes the optimally robust IC for 
  MM2 estimators in case of normal location with unknown scale and 
  (convex) contamination neighborhoods. These estimators are based
  on a proposal of Fraiman et al. (2001), p. 206. A definition of 
  these estimators can also be found in Section 8.6 of Kohl (2005).
}
\usage{
rlsOptIC.MM2(r, c.start = 1.5, d.start = 2, delta = 1e-06, MAX = 100)
}
\arguments{
  \item{r}{ non-negative real: neighborhood radius. }
  \item{c.start}{ positive real: starting value for c. }
  \item{d.start}{ positive real: starting value for d. }
  \item{delta}{ the desired accuracy (convergence tolerance). }
  \item{MAX}{ if a or k are beyond the admitted values, 
    \code{MAX} is returned. }
}
\details{
  The computation of the optimally robust IC for MM2 estimators
  is based on \code{optim} where \code{MAX} is used to 
  control the constraints on c and d. The optimal values of 
  the tuning constants c and d can be read off from the slot 
  \code{Infos} of the resulting IC.
}
%\details{}
\value{Object of class \code{"IC"}}
\references{ 
  Fraiman, R., Yohai, V.J. and Zamar, R.H. (2001) Optimal robust
  M-estimates of location. Ann. Stat. \bold{29}(1): 194--223.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[ROptEst]{IC-class}}}
\examples{
IC1 <- rlsOptIC.MM2(r = 0.1)
checkIC(IC1)
Risks(IC1)
Infos(IC1)
plot(IC1)
infoPlot(IC1)
}
\concept{normal location and scale}
\concept{influence curve}
\keyword{robust}
