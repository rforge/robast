\name{getInfGamma}
\alias{getInfGamma}
\alias{getInfGamma-methods}
\alias{getInfGamma,UnivariateDistribution,asMSE,ContNeighborhood-method}
\alias{getInfGamma,UnivariateDistribution,asGRisk,TotalVarNeighborhood-method}
\alias{getInfGamma,RealRandVariable,asMSE,ContNeighborhood-method}
\alias{getInfGamma,UnivariateDistribution,asUnOvShoot,ContNeighborhood-method}

\title{Generic Function for the Computation of the Optimal Clipping Bound}
\description{
  Generic function for the computation of the optimal clipping bound.
  This function is rarely called directly. It is called by \code{getInfClip} 
  to compute optimally robust ICs.
}
\usage{
getInfGamma(L2deriv, risk, neighbor, ...)

\S4method{getInfGamma}{UnivariateDistribution,asMSE,ContNeighborhood}(L2deriv, risk, neighbor, cent, clip)

\S4method{getInfGamma}{UnivariateDistribution,asGRisk,TotalVarNeighborhood}(L2deriv, risk, neighbor, cent, clip)

\S4method{getInfGamma}{RealRandVariable,asMSE,ContNeighborhood}(L2deriv, risk, neighbor, Distr, stand, cent, clip)

\S4method{getInfGamma}{UnivariateDistribution,asUnOvShoot,ContNeighborhood}(L2deriv, risk, neighbor, cent, clip)
}
\arguments{
  \item{L2deriv}{ L2-derivative of some L2-differentiable family 
    of probability measures. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters }
  \item{cent}{ optimal centering constant. }
  \item{clip}{ optimal clipping bound. }
  \item{stand}{ standardizing matrix. }
  \item{Distr}{ object of class \code{"Distribution"}. }
}
\details{
  The function is used in case of asymptotic G-risks; confer
  Ruckdeschel and Rieder (2004).
}
%\value{}
\section{Methods}{
\describe{
  \item{L2deriv = "UnivariateDistribution", risk = "asMSE", 
        neighbor = "ContNeighborhood"}{ used by \code{getInfClip}. }

  \item{L2deriv = "UnivariateDistribution", risk = "asGRisk", 
        neighbor = "TotalVarNeighborhood"}{ used by \code{getInfClip}. }

  \item{L2deriv = "RealRandVariable", risk = "asMSE", 
        neighbor = "ContNeighborhood"}{ used by \code{getInfClip}. }

  \item{L2deriv = "UnivariateDistribution", risk = "asUnOvShoot", 
        neighbor = "ContNeighborhood"}{ used by \code{getInfClip}. }
}}
\references{
  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Ruckdeschel, P. and Rieder, H. (2004) Optimal Influence Curves for
  General Loss Functions. Statistics & Decisions (submitted).

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{asGRisk-class}}, \code{\link{asMSE-class}},
        \code{\link{asUnOvShoot-class}}, \code{\link{ContIC-class}}, 
        \code{\link{TotalVarIC-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
