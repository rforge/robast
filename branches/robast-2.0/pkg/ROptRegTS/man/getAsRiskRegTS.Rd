\name{getAsRiskRegTS}
\alias{getAsRiskRegTS}
\alias{getAsRiskRegTS-methods}
\alias{getAsRiskRegTS,asMSE,UnivariateDistribution,Distribution,Neighborhood-method}
\alias{getAsRiskRegTS,asMSE,UnivariateDistribution,Distribution,Av2CondContNeighborhood-method}
\alias{getAsRiskRegTS,asMSE,EuclRandVariable,Distribution,Neighborhood-method}
\alias{getAsRiskRegTS,asBias,UnivariateDistribution,UnivariateDistribution,ContNeighborhood-method}
\alias{getAsRiskRegTS,asBias,UnivariateDistribution,UnivariateDistribution,Av1CondContNeighborhood-method}
\alias{getAsRiskRegTS,asBias,UnivariateDistribution,UnivariateDistribution,Av1CondTotalVarNeighborhood-method}
\alias{getAsRiskRegTS,asBias,UnivariateDistribution,MultivariateDistribution,ContNeighborhood-method}
\alias{getAsRiskRegTS,asBias,UnivariateDistribution,MultivariateDistribution,Av1CondContNeighborhood-method}
\alias{getAsRiskRegTS,asBias,UnivariateDistribution,MultivariateDistribution,Av1CondTotalVarNeighborhood-method}
\alias{getAsRiskRegTS,asBias,UnivariateDistribution,Distribution,Av2CondContNeighborhood-method}
\alias{getAsRiskRegTS,asBias,RealRandVariable,Distribution,ContNeighborhood-method}
\alias{getAsRiskRegTS,asBias,RealRandVariable,Distribution,Av1CondContNeighborhood-method}
\alias{getAsRiskRegTS,asUnOvShoot,UnivariateDistribution,UnivariateDistribution,UncondNeighborhood-method}
\alias{getAsRiskRegTS,asUnOvShoot,UnivariateDistribution,UnivariateDistribution,CondNeighborhood-method}

\title{Generic Function for Computation of Asymptotic Risks in case of Regression-Type Models}
\description{
  Generic function for the computation of asymptotic risks in case
  of regression-type models. This function is rarely called directly. 
  It is used by other functions.
}
\usage{
getAsRiskRegTS(risk, ErrorL2deriv, Regressor, neighbor, ...)

\S4method{getAsRiskRegTS}{asMSE,UnivariateDistribution,Distribution,Neighborhood}(risk, ErrorL2deriv, 
                Regressor, neighbor, clip, cent, stand, trafo)

\S4method{getAsRiskRegTS}{asMSE,UnivariateDistribution,Distribution,Av2CondContNeighborhood}(risk, ErrorL2deriv, 
                Regressor, neighbor, clip, cent, stand, trafo)

\S4method{getAsRiskRegTS}{asMSE,EuclRandVariable,Distribution,Neighborhood}(risk, ErrorL2deriv, 
                Regressor, neighbor, clip, cent, stand, trafo)

\S4method{getAsRiskRegTS}{asBias,UnivariateDistribution,UnivariateDistribution,ContNeighborhood}(risk, ErrorL2deriv, 
                Regressor, neighbor, ErrorL2derivDistrSymm, trafo, maxiter, tol)

\S4method{getAsRiskRegTS}{asBias,UnivariateDistribution,UnivariateDistribution,Av1CondContNeighborhood}(risk, 
                ErrorL2deriv, Regressor, neighbor, ErrorL2derivDistrSymm, trafo, maxiter, tol)

\S4method{getAsRiskRegTS}{asBias,UnivariateDistribution,UnivariateDistribution,Av1CondTotalVarNeighborhood}(risk, 
                ErrorL2deriv, Regressor, neighbor, ErrorL2derivDistrSymm, trafo, maxiter, tol)

\S4method{getAsRiskRegTS}{asBias,UnivariateDistribution,MultivariateDistribution,ContNeighborhood}(risk, 
                ErrorL2deriv, Regressor, neighbor, ErrorL2derivDistrSymm, trafo, maxiter, tol)

\S4method{getAsRiskRegTS}{asBias,UnivariateDistribution,MultivariateDistribution,Av1CondContNeighborhood}(risk, 
                ErrorL2deriv, Regressor, neighbor, ErrorL2derivDistrSymm, trafo, maxiter, tol)

\S4method{getAsRiskRegTS}{asBias,UnivariateDistribution,MultivariateDistribution,Av1CondTotalVarNeighborhood}(risk, 
                ErrorL2deriv, Regressor, neighbor, ErrorL2derivDistrSymm, trafo, maxiter, tol)

\S4method{getAsRiskRegTS}{asBias,UnivariateDistribution,Distribution,Av2CondContNeighborhood}(risk, 
                ErrorL2deriv, Regressor, neighbor, ErrorL2derivDistrSymm, trafo, maxiter, tol)

\S4method{getAsRiskRegTS}{asBias,RealRandVariable,Distribution,ContNeighborhood}(risk, 
                ErrorL2deriv, Regressor, neighbor, ErrorDistr, trafo, z.start, A.start, maxiter, tol)

\S4method{getAsRiskRegTS}{asBias,RealRandVariable,Distribution,Av1CondContNeighborhood}(risk, 
                ErrorL2deriv, Regressor, neighbor, ErrorDistr, trafo, z.start, A.start, maxiter, tol)

\S4method{getAsRiskRegTS}{asUnOvShoot,UnivariateDistribution,UnivariateDistribution,UncondNeighborhood}(risk, 
                ErrorL2deriv, Regressor, neighbor, clip, cent, stand)

\S4method{getAsRiskRegTS}{asUnOvShoot,UnivariateDistribution,UnivariateDistribution,CondNeighborhood}(risk, 
                ErrorL2deriv, Regressor, neighbor, clip, cent, stand)
}

\arguments{
  \item{risk}{ object of class \code{"asRisk"}. }
  \item{ErrorL2deriv}{ L2-derivative of \code{ErrorDistr}. }
  \item{Regressor}{ regressor. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters. }
  \item{clip}{ optimal clipping bound. }
  \item{cent}{ optimal centering constant/function. }
  \item{stand}{ standardizing matrix. }
  \item{trafo}{ matrix: transformation of the parameter. }
  \item{ErrorDistr}{ error distribution. }
  \item{ErrorL2derivDistrSymm}{ symmetry of \code{ErrorL2derivDistr}. }
  \item{maxiter}{ the maximum number of iterations }
  \item{tol}{ the desired accuracy (convergence tolerance). }
  \item{z.start}{ initial value for the centering constant/function. }
  \item{A.start}{ initial value for the standardizing matrix. }
}
%\details{}
\value{The asymptotic risk is computed.}
\section{Methods}{
\describe{
  \item{risk = "asMSE", ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", neighbor = "Neighborhood"}{
    computes asymptotic mean square error in methods for function \code{getInfRobRegTypeIC}. }

  \item{risk = "asMSE", ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", neighbor = "Av2CondContNeighborhood"}{ 
    computes asymptotic mean square error in methods for function \code{getInfRobRegTypeIC}. }

  \item{risk = "asMSE", ErrorL2deriv = "EuclRandVariable", Regressor = "Distribution", neighbor = "Neighborhood"}{
    computes asymptotic mean square error in methods for function \code{getInfRobRegTypeIC}. }

  \item{risk = "asBias", ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", neighbor = "ContNeighborhood"}{
    computes standardized asymptotic bias in methods for function \code{getInfRobRegTypeIC}. }

  \item{risk = "asBias", ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", neighbor = "Av1CondContNeighborhood"}{ 
    computes standardized asymptotic bias in methods for function \code{getInfRobRegTypeIC}. }

  \item{risk = "asBias", ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", neighbor = "Av1CondTotalVarNeighborhood"}{ 
    computes standardized asymptotic bias in methods for function \code{getInfRobRegTypeIC}. }

  \item{risk = "asBias", ErrorL2deriv = "UnivariateDistribution", Regressor = "MultivariateDistribution", neighbor = "ContNeighborhood"}{
    computes standardized asymptotic bias in methods for function \code{getInfRobRegTypeIC}. }

  \item{risk = "asBias", ErrorL2deriv = "UnivariateDistribution", Regressor = "MultivariateDistribution", neighbor = "Av1CondContNeighborhood"}{ 
    computes standardized asymptotic bias in methods for function \code{getInfRobRegTypeIC}. }

  \item{risk = "asBias", ErrorL2deriv = "UnivariateDistribution", Regressor = "MultivariateDistribution", neighbor = "Av1CondTotalVarNeighborhood"}{ 
    computes standardized asymptotic bias in methods for function \code{getInfRobRegTypeIC}. }

  \item{risk = "asBias", ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", neighbor = "Av2CondContNeighborhood"}{ 
    computes standardized asymptotic bias in methods for function \code{getInfRobRegTypeIC}. }

  \item{risk = "asBias", ErrorL2deriv = "RealRandVariable", Regressor = "Distribution", neighbor = "ContNeighborhood"}{ 
    computes standardized asymptotic bias in methods for function \code{getInfRobRegTypeIC}. }

  \item{risk = "asBias", ErrorL2deriv = "RealRandVariable", Regressor = "Distribution", neighbor = "Av1CondContNeighborhood"}{ 
    computes standardized asymptotic bias in methods for function \code{getInfRobRegTypeIC}. }

  \item{risk = "asUnOvShoot", ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", neighbor = "UncondNeighborhood"}{
    computes asymptotic under-/overshoot risk in methods for function \code{getInfRobRegTypeIC}. }

  \item{risk = "asUnOvShoot", ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", neighbor = "CondNeighborhood"}{
    computes asymptotic under-/overshoot risk in methods for function \code{getInfRobRegTypeIC}. }
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
\seealso{\code{\link[ROptEst]{asRisk-class}}}
%\examples{}
\concept{asymptotic risk}
\concept{risk}
\keyword{}
