\name{cniperCont}
\alias{cniperCont}
\alias{cniperCont-methods}
\alias{cniperCont,IC,IC,L2ParamFamily,ContNeighborhood,asMSE-method}
\alias{cniperPoint}
\alias{cniperPoint-methods}
\alias{cniperPoint,L2ParamFamily,ContNeighborhood,asMSE-method}
\alias{cniperPointPlot}
\alias{cniperPointPlot-methods}
\alias{cniperPointPlot,L2ParamFamily,ContNeighborhood,asMSE-method}
\title{ Generic Functions for Computation and Plot of Cniper Contamination 
        and Cniper Points. }
\description{
  These generic functions and their methods can be used to determine cniper 
  contamination as well as cniper points. That is, under which (Dirac) 
  contamination is the risk of one procedure larger than the risk of some
  other procedure.
}
\usage{
cniperCont(IC1, IC2, L2Fam, neighbor, risk, ...)
\S4method{cniperCont}{IC,IC,L2ParamFamily,ContNeighborhood,asMSE}(IC1, 
      IC2, L2Fam, neighbor, risk, lower, upper, n = 101)

cniperPoint(L2Fam, neighbor, risk, ...)
\S4method{cniperPoint}{L2ParamFamily,ContNeighborhood,asMSE}(L2Fam, 
      neighbor, risk, lower, upper)

cniperPointPlot(L2Fam, neighbor, risk, ...)
\S4method{cniperPointPlot}{L2ParamFamily,ContNeighborhood,asMSE}(L2Fam, 
      neighbor, risk, lower, upper, n = 101)
}
\arguments{
  \item{IC1}{ object of class \code{IC} }
  \item{IC2}{ object of class \code{IC} }
  \item{L2Fam}{ object of class \code{L2ParamFamily} }
  \item{neighbor}{ object of class \code{Neighborhood} }
  \item{risk}{ object of class \code{RiskType} }
  \item{\dots}{ additional parameters. }
  \item{lower, upper}{ the lower and upper end points of the 
          contamination interval. }
  \item{n}{ number of points between \code{lower} and \code{upper}}
}
\details{
  In case of \code{cniperCont} the difference between the risks of two ICs 
  is plotted.

  The function \code{cniperPoint} can be used to determine cniper
  points. That is, points such that the optimally robust estimator
  has smaller minimax risk than the classical optimal estimator under 
  contamination with Dirac measures at the cniper points. 

  As such points might be difficult to find, we provide the
  function \code{cniperPointPlot} which can be used to obtain a plot
  of the risk difference.

  For more details about cniper contamination and cniper points we refer 
  to Section~3.5 of Kohl et al. (2008) as well as Ruckdeschel (2004) and 
  the Introduction of Kohl (2005).
}
\value{invisible() resp. cniper point is returned.}
\references{
  Kohl, M. and Ruckdeschel, H. and Rieder, H. (2008). Infinitesimally 
  Robust Estimation in General Smoothly Parametrized Models. Unpublished
  Manuscript.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Ruckdeschel, P. (2004). Higher Order Asymptotics for the MSE of M-Estimators
  on Shrinking Neighborhoods. Unpublished Manuscript.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
## cniper contamination
P <- PoisFamily(lambda = 4)
RobP1 <- InfRobModel(center = P, neighbor = ContNeighborhood(radius = 0.1))
IC1 <- optIC(model=RobP1, risk=asMSE())
RobP2 <- InfRobModel(center = P, neighbor = ContNeighborhood(radius = 1))
IC2 <- optIC(model=RobP2, risk=asMSE())
cniperCont(IC1 = IC1, IC2 = IC2, L2Fam = P, 
           neighbor = ContNeighborhood(radius = 0.5), 
           risk = asMSE(),
           lower = 0, upper = 8, n = 101)

## cniper point plot
cniperPointPlot(P, neighbor = ContNeighborhood(radius = 0.5), 
                risk = asMSE(), lower = 0, upper = 10)

## cniper point
cniperPoint(P, neighbor = ContNeighborhood(radius = 0.5), 
            risk = asMSE(), lower = 0, upper = 4)
cniperPoint(P, neighbor = ContNeighborhood(radius = 0.5), 
            risk = asMSE(), lower = 4, upper = 8)
}
\concept{cniper contamination}
\concept{cniper point}
\keyword{robust}
