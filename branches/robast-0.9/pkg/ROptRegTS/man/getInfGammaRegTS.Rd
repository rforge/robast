\name{getInfGammaRegTS}
\alias{getInfGammaRegTS}
\alias{getInfGammaRegTS-methods}
\alias{getInfGammaRegTS,UnivariateDistribution,UnivariateDistribution,asMSE,ContNeighborhood-method}
\alias{getInfGammaRegTS,UnivariateDistribution,UnivariateDistribution,asMSE,Av1CondContNeighborhood-method}
\alias{getInfGammaRegTS,UnivariateDistribution,UnivariateDistribution,asMSE,Av1CondTotalVarNeighborhood-method}
\alias{getInfGammaRegTS,UnivariateDistribution,MultivariateDistribution,asMSE,ContNeighborhood-method}
\alias{getInfGammaRegTS,UnivariateDistribution,MultivariateDistribution,asMSE,Av1CondContNeighborhood-method}
\alias{getInfGammaRegTS,UnivariateDistribution,MultivariateDistribution,asMSE,Av1CondTotalVarNeighborhood-method}
\alias{getInfGammaRegTS,UnivariateDistribution,Distribution,asMSE,Av2CondContNeighborhood-method}
\alias{getInfGammaRegTS,RealRandVariable,Distribution,asMSE,ContNeighborhood-method}
\alias{getInfGammaRegTS,RealRandVariable,Distribution,asMSE,Av1CondContNeighborhood-method}
\alias{getInfGammaRegTS,UnivariateDistribution,UnivariateDistribution,asUnOvShoot,ContNeighborhood-method}
\alias{getInfGammaRegTS,UnivariateDistribution,UnivariateDistribution,asUnOvShoot,TotalVarNeighborhood-method}
\alias{getInfGammaRegTS,UnivariateDistribution,numeric,asUnOvShoot,CondContNeighborhood-method}
\alias{getInfGammaRegTS,UnivariateDistribution,numeric,asUnOvShoot,CondTotalVarNeighborhood-method}

\title{ Generic Function for the Computation of the Optimal Clipping Bound }
\description{
  Generic function for the computation of the optimal clipping bound.
  This function is rarely called directly. It is called by \code{getInfClipRegTS} 
  to compute optimally robust ICs.
}
\usage{
getInfGammaRegTS(ErrorL2deriv, Regressor, risk, neighbor, ...)

\S4method{getInfGammaRegTS}{UnivariateDistribution,UnivariateDistribution,asMSE,ContNeighborhood}(
             ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent, clip)

\S4method{getInfGammaRegTS}{UnivariateDistribution,UnivariateDistribution,asMSE,Av1CondContNeighborhood}(
             ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent, clip)

\S4method{getInfGammaRegTS}{UnivariateDistribution,UnivariateDistribution,asMSE,Av1CondTotalVarNeighborhood}(
             ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent, clip)

\S4method{getInfGammaRegTS}{UnivariateDistribution,MultivariateDistribution,asMSE,ContNeighborhood}(
             ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent, clip)

\S4method{getInfGammaRegTS}{UnivariateDistribution,MultivariateDistribution,asMSE,Av1CondContNeighborhood}(
             ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent, clip)

\S4method{getInfGammaRegTS}{UnivariateDistribution,MultivariateDistribution,asMSE,Av1CondTotalVarNeighborhood}(
             ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent, clip)

\S4method{getInfGammaRegTS}{UnivariateDistribution,Distribution,asMSE,Av2CondContNeighborhood}(
             ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent, clip)

\S4method{getInfGammaRegTS}{RealRandVariable,Distribution,asMSE,ContNeighborhood}(
             ErrorL2deriv, Regressor, risk, neighbor, ErrorDistr, stand, cent,
             clip)

\S4method{getInfGammaRegTS}{RealRandVariable,Distribution,asMSE,Av1CondContNeighborhood}(
             ErrorL2deriv, Regressor, risk, neighbor, ErrorDistr, stand, cent,
             clip)

\S4method{getInfGammaRegTS}{UnivariateDistribution,UnivariateDistribution,asUnOvShoot,ContNeighborhood}(
             ErrorL2deriv, Regressor, risk, neighbor, cent, clip)

\S4method{getInfGammaRegTS}{UnivariateDistribution,UnivariateDistribution,asUnOvShoot,TotalVarNeighborhood}(
             ErrorL2deriv, Regressor, risk, neighbor, cent, clip)

\S4method{getInfGammaRegTS}{UnivariateDistribution,numeric,asUnOvShoot,CondContNeighborhood}(
             ErrorL2deriv, Regressor, risk, neighbor, clip)

\S4method{getInfGammaRegTS}{UnivariateDistribution,numeric,asUnOvShoot,CondTotalVarNeighborhood}(
             ErrorL2deriv, Regressor, risk, neighbor, clip)
}
\arguments{
  \item{ErrorL2deriv}{ L2-derivative of \code{ErrorDistr}. }
  \item{Regressor}{ regressor. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters. }
  \item{clip}{ optimal clipping bound. }
  \item{cent}{ optimal centering constant/function. }
  \item{stand}{ standardizing matrix. }
  \item{z.comp}{ which components of the centering constant/function 
    have to be computed. }
  \item{ErrorDistr}{ error distribution. }
}
\details{
  The function is used in case of asymptotic G-risks; confer
  Ruckdeschel and Rieder (2004).
}
%\value{}
\section{Methods}{
\describe{
  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", risk = "asMSE", neighbor = "ContNeighborhood"}{ 
    used by \code{getInfClipRegTS}. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", risk = "asMSE", neighbor = "Av1CondContNeighborhood"}{ 
    used by \code{getInfClipRegTS}. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", risk = "asMSE", neighbor = "Av1CondTotalVarNeighborhood"}{ 
    used by \code{getInfClipRegTS}. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "MultivariateDistribution", risk = "asMSE", neighbor = "ContNeighborhood"}{ 
    used by \code{getInfClipRegTS}. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "MultivariateDistribution", risk = "asMSE", neighbor = "Av1CondContNeighborhood"}{ 
    used by \code{getInfClipRegTS}. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "MultivariateDistribution", risk = "asMSE", neighbor = "Av1CondTotalVarNeighborhood"}{ 
    used by \code{getInfClipRegTS}. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", risk = "asMSE", neighbor = "Av2CondContNeighborhood"}{ 
    used by \code{getInfClipRegTS}. }

  \item{ErrorL2deriv = "RealRandVariable", Regressor = "Distribution", risk = "asMSE", neighbor = "ContNeighborhood"}{ 
    used by \code{getInfClipRegTS}. }

  \item{ErrorL2deriv = "RealRandVariable", Regressor = "Distribution", risk = "asMSE", neighbor = "Av1CondContNeighborhood"}{ 
    used by \code{getInfClipRegTS}. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", risk = "asUnOvShoot", neighbor = "ContNeighborhood"}{ 
    used by \code{getInfClipRegTS}. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", risk = "asUnOvShoot", neighbor = "TotalVarNeighborhood"}{ 
    used by \code{getInfClipRegTS}. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "numeric", risk = "asUnOvShoot", neighbor = "CondContNeighborhood"}{ 
    used by \code{getInfClipRegTS}. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "numeric", risk = "asUnOvShoot", neighbor = "CondTotalVarNeighborhood"}{ 
    used by \code{getInfClipRegTS}. }
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
\seealso{\code{\link[ROptEstOld]{asMSE-class}}, \code{\link[ROptEstOld]{asUnOvShoot-class}}, 
        \code{\link[ROptEstOld]{ContIC-class}}, \code{\link{Av1CondContIC-class}}, 
        \code{\link{Av2CondContIC-class}}, \code{\link{Av1CondTotalVarIC-class}},
        \code{\link{CondContIC-class}}, \code{\link{CondTotalVarIC-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}

