\name{getInfClip}
\alias{getInfClip}
\alias{getInfClip-methods}
\alias{getInfClip,numeric,UnivariateDistribution,asMSE,ContNeighborhood-method}
\alias{getInfClip,numeric,UnivariateDistribution,asMSE,TotalVarNeighborhood-method}
\alias{getInfClip,numeric,EuclRandVariable,asMSE,UncondNeighborhood-method}
\alias{getInfClip,numeric,UnivariateDistribution,asUnOvShoot,UncondNeighborhood-method}
\alias{getInfClip,numeric,UnivariateDistribution,asSemivar,ContNeighborhood-method}

\title{Generic Function for the Computation of the Optimal Clipping Bound}
\description{
  Generic function for the computation of the optimal clipping bound
  in case of infinitesimal robust models. This function is rarely called 
  directly. It is used to compute optimally robust ICs.
}
\usage{
getInfClip(clip, L2deriv, risk, neighbor, ...)

\S4method{getInfClip}{numeric,UnivariateDistribution,asMSE,ContNeighborhood}(clip, L2deriv, 
                                                risk, neighbor, biastype, cent, symm, trafo)

\S4method{getInfClip}{numeric,UnivariateDistribution,asMSE,TotalVarNeighborhood}(clip, L2deriv, 
                                                    risk, neighbor, biastype, cent, symm, trafo)

\S4method{getInfClip}{numeric,EuclRandVariable,asMSE,UncondNeighborhood}(clip, L2deriv, risk,
                                              neighbor, biastype, Distr, stand, cent, trafo)

\S4method{getInfClip}{numeric,UnivariateDistribution,asUnOvShoot,UncondNeighborhood}(clip, L2deriv, 
                                                        risk, neighbor, biastype, cent, symm, trafo)

\S4method{getInfClip}{numeric,UnivariateDistribution,asSemivar,ContNeighborhood}(clip, L2deriv, 
                                                              risk, neighbor, cent, symm, trafo)
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
        risk = "asUnOvShoot", neighbor = "UncondNeighborhood"}{ 
    optimal clipping bound for asymtotic under-/overshoot risk. }

  \item{clip = "numeric", L2deriv = "UnivariateDistribution", 
        risk = "asSemivar", neighbor = "ContNeighborhood"}{ 
    optimal clipping bound for asymtotic semivariance.}
}}
\references{
  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Ruckdeschel, P. (2005) Optimally One-Sided Bounded Influence Curves.
  Mathematical Methods in Statistics \emph{14}(1), 105-131.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},
        Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link[RobAStBase]{ContIC-class}}, \code{\link[RobAStBase]{TotalVarIC-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
