\name{getRadius}
\alias{getRadius}

\title{Computation of the Optimal Radius for Given Clipping Bound}
\description{
  The usual robust optimality problem for given asGRisk searches the optimal
  clipping height b of a Hampel-type IC to given radius of the neighborhood.
  Instead, again for given asGRisk  and for given Hampel-Type IC with
  given clipping height b we may determine the radius of the neighborhood
  for which it is optimal in the sense of the first sentence.
}
\usage{
getRadius(IC, risk, neighbor, L2Fam)
}
\arguments{
  \item{IC}{ an IC of class \code{"HampIC"}.}
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{L2Fam}{ object of class \code{"L2FamParameter"}. Can be missing;
                in this case it is constructed from slot \code{CallL2Fam} of
                \code{IC}. }
}
\value{The optimal radius is computed.}
\references{
  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Ruckdeschel, P. and Rieder, H. (2004) Optimal Influence Curves for
  General Loss Functions. Statistics & Decisions \emph{22}, 201-223.

  Ruckdeschel, P. (2005) Optimally One-Sided Bounded Influence Curves.
  Mathematical Methods in Statistics \emph{14}(1), 105-131.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}.
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link[RobAStBase]{ContIC-class}}, \code{\link[RobAStBase]{TotalVarIC-class}}}
\examples{
N <- NormLocationFamily(mean=0, sd=1)
nb <- ContNeighborhood(); ri <- asMSE()
radIC <- radiusMinimaxIC(L2Fam=N, neighbor=nb, risk=ri, loRad=0.1, upRad=0.5)
getRadius(radIC, L2Fam=N, neighbor=nb, risk=ri)

## taken from script NormalScaleModel.R in folder scripts
N0 <- NormScaleFamily(mean=0, sd=1)
(N0.IC7 <- radiusMinimaxIC(L2Fam=N0, neighbor=nb, risk=ri, loRad=0, upRad=Inf))
##
getRadius(N0.IC7, risk=asMSE(), neighbor=nb, L2Fam=N0)
getRadius(N0.IC7, risk=asL4(), neighbor=nb, L2Fam=N0)
}

\concept{influence curve}
\keyword{robust}
