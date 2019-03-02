\name{getInfRad}
\alias{getInfRad}
\alias{getInfRad-methods}
\alias{getInfRad,numeric,UnivariateDistribution,asMSE,ContNeighborhood-method}
\alias{getInfRad,numeric,UnivariateDistribution,asMSE,TotalVarNeighborhood-method}
\alias{getInfRad,numeric,UnivariateDistribution,asL1,ContNeighborhood-method}
\alias{getInfRad,numeric,UnivariateDistribution,asL1,TotalVarNeighborhood-method}
\alias{getInfRad,numeric,UnivariateDistribution,asL4,ContNeighborhood-method}
\alias{getInfRad,numeric,UnivariateDistribution,asL4,TotalVarNeighborhood-method}
\alias{getInfRad,numeric,EuclRandVariable,asMSE,UncondNeighborhood-method}
\alias{getInfRad,numeric,UnivariateDistribution,asUnOvShoot,UncondNeighborhood-method}
\alias{getInfRad,numeric,UnivariateDistribution,asSemivar,ContNeighborhood-method}

\title{Generic Function for the Computation of the Optimal Radius for Given Clipping Bound}
\description{
  The usual robust optimality problem for given asGRisk searches the optimal
  clipping height b of a Hampel-type IC to given radius of the neighborhood.
  Instead, again for given asGRisk  and for given Hampel-Type IC with
  given clipping height b we may determine the radius of the neighborhood
  for which it is optimal in the sense of the first sentence. This
  radius is determined by \code{getInfRad}. This function is rarely called
  directly. It is used withing \code{\link{getRadius}}.
}
\usage{
getInfRad(clip, L2deriv, risk, neighbor, ...)

\S4method{getInfRad}{numeric,UnivariateDistribution,asMSE,ContNeighborhood}(
          clip, L2deriv, risk, neighbor, biastype, cent, symm, trafo)

\S4method{getInfRad}{numeric,UnivariateDistribution,asMSE,TotalVarNeighborhood}(
          clip, L2deriv, risk, neighbor, biastype, cent, symm, trafo)

\S4method{getInfRad}{numeric,UnivariateDistribution,asL1,ContNeighborhood}(
          clip, L2deriv, risk, neighbor, biastype, cent, symm, trafo)

\S4method{getInfRad}{numeric,UnivariateDistribution,asL1,TotalVarNeighborhood}(
          clip, L2deriv, risk, neighbor, biastype, cent, symm, trafo)

\S4method{getInfRad}{numeric,UnivariateDistribution,asL4,ContNeighborhood}(
          clip, L2deriv, risk, neighbor, biastype, cent, symm, trafo)

\S4method{getInfRad}{numeric,UnivariateDistribution,asL4,TotalVarNeighborhood}(
          clip, L2deriv, risk, neighbor, biastype, cent, symm, trafo)

\S4method{getInfRad}{numeric,EuclRandVariable,asMSE,UncondNeighborhood}(
          clip, L2deriv, risk, neighbor, biastype, Distr, stand, cent, trafo, ...)

\S4method{getInfRad}{numeric,UnivariateDistribution,asUnOvShoot,UncondNeighborhood}(
          clip, L2deriv, risk, neighbor, biastype, cent, symm, trafo)

\S4method{getInfRad}{numeric,UnivariateDistribution,asSemivar,ContNeighborhood}(
          clip, L2deriv, risk, neighbor, biastype, cent, symm, trafo)
}
\arguments{
  \item{clip}{ positive real: clipping bound }
  \item{L2deriv}{ L2-derivative of some L2-differentiable family 
    of probability measures. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters. }
  \item{biastype}{ object of class \code{"BiasType"} }
  \item{cent}{ optimal centering constant. }
  \item{stand}{ standardizing matrix. }
  \item{Distr}{ object of class \code{"Distribution"}. }
  \item{symm}{ logical: indicating symmetry of \code{L2deriv}. }
  \item{trafo}{ matrix: transformation of the parameter. }
}
%\details{}
\value{The optimal clipping bound is computed.}
\section{Methods}{
\describe{
  \item{clip = "numeric", L2deriv = "UnivariateDistribution", 
        risk = "asMSE", neighbor = "ContNeighborhood"}{ 
    optimal clipping bound for asymtotic mean square error. }


  \item{clip = "numeric", L2deriv = "UnivariateDistribution", 
        risk = "asMSE", neighbor = "TotalVarNeighborhood"}{ 
    optimal clipping bound for asymtotic mean square error. }

  \item{clip = "numeric", L2deriv = "EuclRandVariable", 
        risk = "asMSE", neighbor = "UncondNeighborhood"}{
    optimal clipping bound for asymtotic mean square error. }

  \item{clip = "numeric", L2deriv = "UnivariateDistribution", 
        risk = "asL1", neighbor = "ContNeighborhood"}{ 
    optimal clipping bound for asymtotic mean absolute error. }

  \item{clip = "numeric", L2deriv = "UnivariateDistribution", 
        risk = "asL1", neighbor = "TotalVarNeighborhood"}{ 
    optimal clipping bound for asymtotic mean absolute error. }

  \item{clip = "numeric", L2deriv = "UnivariateDistribution", 
        risk = "asL4", neighbor = "ContNeighborhood"}{ 
    optimal clipping bound for asymtotic mean power 4 error. }

  \item{clip = "numeric", L2deriv = "UnivariateDistribution", 
        risk = "asL4", neighbor = "TotalVarNeighborhood"}{ 
    optimal clipping bound for asymtotic mean power 4 error. }

  \item{clip = "numeric", L2deriv = "UnivariateDistribution", 
        risk = "asUnOvShoot", neighbor = "UncondNeighborhood"}{ 
    optimal clipping bound for asymtotic under-/overshoot risk. }

  \item{clip = "numeric", L2deriv = "UnivariateDistribution", 
        risk = "asSemivar", neighbor = "ContNeighborhood"}{ 
    optimal clipping bound for asymtotic semivariance.}
}}
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
%\examples{}
\concept{influence curve}
\keyword{robust}
