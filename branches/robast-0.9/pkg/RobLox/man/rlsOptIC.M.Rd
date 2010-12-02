\name{rlsOptIC.M}
\alias{rlsOptIC.M}
\title{Computation of the optimally robust IC for M estimators}
\description{
  The function \code{rlsOptIC.M} computes the optimally robust IC for 
  M estimators in case of normal location with unknown scale and 
  (convex) contamination neighborhoods. The definition of 
  these estimators can be found in Section 8.3 of Kohl (2005).
}
\usage{
rlsOptIC.M(r, ggLo = 0.5, ggUp = 1.5, a1.start = 0.75, a3.start = 0.25, 
           bUp = 1000, delta = 1e-05, itmax = 100, check = FALSE)
}
\arguments{
  \item{r}{ non-negative real: neighborhood radius. }
  \item{ggLo}{ non-negative real: the lower end point of the interval to be searched 
    for \eqn{\gamma}{gamma}. }
  \item{ggUp}{ positive real: the upper end point of the interval to be searched 
    for \eqn{\gamma}{gamma}. }
  \item{a1.start}{ real: starting value for \eqn{\alpha_1}{alpha_1}. }
  \item{a3.start}{ real: starting value for \eqn{\alpha_3}{alpha_3}. }
  \item{bUp}{ positive real: upper bound used in the 
    computation of the optimal clipping bound b. }
  \item{delta}{ the desired accuracy (convergence tolerance). }
  \item{itmax}{ the maximum number of iterations. }
  \item{check}{ logical. Should constraints be checked. }
}
\details{The optimal values of the tuning constants \eqn{\alpha_1}{alpha_1}, 
  \eqn{\alpha_3}{alpha_3}, b and \eqn{\gamma}{gamma} can be read off 
  from the slot \code{Infos} of the resulting IC.}
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
IC1 <- rlsOptIC.M(r = 0.1, check = TRUE)
distrExOptions("ErelativeTolerance" = 1e-12)
checkIC(IC1, NormLocationScaleFamily())
distrExOptions("ErelativeTolerance" = .Machine$double.eps^0.25)
Risks(IC1)
Infos(IC1)
plot(IC1)
infoPlot(IC1)
}
\concept{normal location and scale}
\concept{influence curve}
\keyword{robust}
