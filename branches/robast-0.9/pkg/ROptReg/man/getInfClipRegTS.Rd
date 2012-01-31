\name{getInfClipRegTS}
\alias{getInfClipRegTS}
\alias{getInfClipRegTS,numeric,UnivariateDistribution,Distribution,asMSE,Neighborhood-method}
\alias{getInfClipRegTS,numeric,UnivariateDistribution,Distribution,asMSE,Av1CondTotalVarNeighborhood-method}
\alias{getInfClipRegTS,numeric,EuclRandVariable,Distribution,asMSE,Neighborhood-method}
\alias{getInfClipRegTS,numeric,UnivariateDistribution,UnivariateDistribution,asUnOvShoot,UncondNeighborhood-method}
\alias{getInfClipRegTS,numeric,UnivariateDistribution,numeric,asUnOvShoot,CondNeighborhood-method}

\title{Generic Function for the Computation of the Optimal Clipping Bound}
\description{
  Generic function for the computation of the optimal clipping bound/function.
  This function is rarely called directly. It is used to 
  compute optimally robust ICs in case infinitesimal models.
}
\usage{
getInfClipRegTS(clip, ErrorL2deriv, Regressor, risk, neighbor, ...)

\S4method{getInfClipRegTS}{numeric,UnivariateDistribution,Distribution,asMSE,Neighborhood}(clip, 
                ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent)

\S4method{getInfClipRegTS}{numeric,UnivariateDistribution,Distribution,asMSE,Av1CondTotalVarNeighborhood}(clip, 
                ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent)

\S4method{getInfClipRegTS}{numeric,EuclRandVariable,Distribution,asMSE,Neighborhood}(clip, ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorDistr, stand, cent, trafo)

\S4method{getInfClipRegTS}{numeric,UnivariateDistribution,UnivariateDistribution,asUnOvShoot,UncondNeighborhood}(clip, 
                ErrorL2deriv, Regressor, risk, neighbor, z.comp, cent)

\S4method{getInfClipRegTS}{numeric,UnivariateDistribution,numeric,asUnOvShoot,CondNeighborhood}(clip, 
                ErrorL2deriv, Regressor, risk, neighbor)
}
\arguments{
  \item{clip}{ optimal clipping bound. }
  \item{ErrorL2deriv}{ L2-derivative of \code{ErrorDistr}. }
  \item{Regressor}{ regressor. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters. }
  \item{cent}{ optimal centering constant/function. }
  \item{stand}{ standardizing matrix. }
  \item{z.comp}{ which components of the centering constant/function 
    have to be computed. }
  \item{ErrorDistr}{ error distribution. }
  \item{trafo}{ matrix: transformation of the parameter. }  
}
%\details{}
\value{The optimal clipping bound/function is computed.}
\section{Methods}{
\describe{
  \item{clip = "numeric", ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", risk = "asMSE", neighbor = "Neighborhood"}{ 
    optimal clipping bound for asymtotic mean square error. }

  \item{clip = "numeric", ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", risk = "asMSE", neighbor = "Av1CondTotalVarNeighborhood"}{ 
    optimal clipping bound for asymtotic mean square error. }

  \item{clip = "numeric", ErrorL2deriv = "EuclRandVariable", Regressor = "Distribution", risk = "asMSE", neighbor = "Neighborhood"}{ 
    optimal clipping bound for asymtotic mean square error. }

  \item{clip = "numeric", ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", risk = "asUnOvShoot", neighbor = "UncondNeighborhood"}{ 
    optimal clipping bound for asymtotic under-/overshoot risk. }

  \item{clip = "numeric", ErrorL2deriv = "UnivariateDistribution", Regressor = "numeric", risk = "asUnOvShoot", neighbor = "CondNeighborhood"}{ 
    optimal clipping function for asymtotic under-/overshoot risk. }
}}
\references{
  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[ROptEstOld]{ContIC-class}}, \code{\link[ROptEstOld]{TotalVarIC-class}}, 
        \code{\link{Av1CondContIC-class}}, \code{\link{Av2CondContIC-class}}, 
        \code{\link{Av1CondTotalVarIC-class}}, \code{\link{CondContIC-class}},
        \code{\link{CondTotalVarIC-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}

