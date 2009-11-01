\name{getInfGamma}
\alias{getInfGamma}
\alias{getInfGamma-methods}
\alias{getInfGamma,UnivariateDistribution,asMSE,ContNeighborhood,BiasType-method}
\alias{getInfGamma,UnivariateDistribution,asGRisk,TotalVarNeighborhood,BiasType-method}
\alias{getInfGamma,RealRandVariable,asMSE,ContNeighborhood,BiasType-method}
\alias{getInfGamma,RealRandVariable,asMSE,TotalVarNeighborhood,BiasType-method}
\alias{getInfGamma,UnivariateDistribution,asUnOvShoot,ContNeighborhood,BiasType-method}
\alias{getInfGamma,UnivariateDistribution,asMSE,ContNeighborhood,onesidedBias-method}
\alias{getInfGamma,UnivariateDistribution,asMSE,ContNeighborhood,asymmetricBias-method}

\title{Generic Function for the Computation of the Optimal Clipping Bound}
\description{
  Generic function for the computation of the optimal clipping bound.
  This function is rarely called directly. It is called by \code{getInfClip} 
  to compute optimally robust ICs.
}
\usage{
getInfGamma(L2deriv, risk, neighbor, biastype, ...)

\S4method{getInfGamma}{UnivariateDistribution,asMSE,ContNeighborhood,BiasType}(L2deriv, 
     risk, neighbor, biastype, cent, clip)

\S4method{getInfGamma}{UnivariateDistribution,asGRisk,TotalVarNeighborhood,BiasType}(L2deriv, 
     risk, neighbor, biastype, cent, clip)

\S4method{getInfGamma}{RealRandVariable,asMSE,ContNeighborhood,BiasType}(L2deriv, 
     risk, neighbor, biastype, Distr, stand, cent, clip, power = 1L)

\S4method{getInfGamma}{RealRandVariable,asMSE,TotalVarNeighborhood,BiasType}(L2deriv,
     risk, neighbor, biastype, Distr, stand, cent, clip, power = 1L)

\S4method{getInfGamma}{UnivariateDistribution,asUnOvShoot,ContNeighborhood,BiasType}(L2deriv,
     risk, neighbor, biastype, cent, clip)

\S4method{getInfGamma}{UnivariateDistribution,asMSE,ContNeighborhood,onesidedBias}(L2deriv, 
     risk, neighbor, biastype, cent, clip)

\S4method{getInfGamma}{UnivariateDistribution,asMSE,ContNeighborhood,asymmetricBias}(L2deriv, 
    risk, neighbor, biastype, cent, clip)
}
\arguments{
  \item{L2deriv}{ L2-derivative of some L2-differentiable family 
    of probability measures. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{biastype}{ object of class \code{"BiasType"} }
  \item{\dots}{ additional parameters }
  \item{cent}{ optimal centering constant. }
  \item{clip}{ optimal clipping bound. }
  \item{stand}{ standardizing matrix. }
  \item{Distr}{ object of class \code{"Distribution"}. }
  \item{power}{ exponent for the integrand; by default \code{1}, but
   may also be \code{2}, for optimization in \code{getLagrangeMultByOptim}. }
}
\details{
  The function is used in case of asymptotic G-risks; confer
  Ruckdeschel and Rieder (2004).
}
%\value{}
\section{Methods}{
\describe{
  \item{L2deriv = "UnivariateDistribution", risk = "asMSE", 
        neighbor = "ContNeighborhood", 
        biastype = "BiasType"}{ used by \code{getInfClip} for symmetric bias. }

  \item{L2deriv = "UnivariateDistribution", risk = "asGRisk", 
        neighbor = "TotalVarNeighborhood", 
        biastype = "BiasType"}{ used by \code{getInfClip} for symmetric bias. }

  \item{L2deriv = "RealRandVariable", risk = "asMSE", 
        neighbor = "ContNeighborhood", 
        biastype = "BiasType"}{ used by \code{getInfClip} for symmetric bias. }

  \item{L2deriv = "RealRandVariable", risk = "asMSE",
        neighbor = "TotalVarNeighborhood",
        biastype = "BiasType"}{ used by \code{getInfClip} for symmetric bias. }

  \item{L2deriv = "UnivariateDistribution", risk = "asUnOvShoot", 
        neighbor = "ContNeighborhood", 
        biastype = "BiasType"}{ used by \code{getInfClip} for symmetric bias. }

  \item{L2deriv = "UnivariateDistribution", risk = "asMSE", 
        neighbor = "ContNeighborhood", 
        biastype = "onesidedBias"}{ used by \code{getInfClip} for onesided bias. }

  \item{L2deriv = "UnivariateDistribution", risk = "asMSE", 
        neighbor = "ContNeighborhood", 
        biastype = "asymmetricBias"}{ used by \code{getInfClip} for asymmetric bias. }
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
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},
        Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link[distrMod]{asGRisk-class}}, \code{\link[distrMod]{asMSE-class}},
        \code{\link[distrMod]{asUnOvShoot-class}}, \code{\link[RobAStBase]{ContIC-class}}, 
        \code{\link[RobAStBase]{TotalVarIC-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
