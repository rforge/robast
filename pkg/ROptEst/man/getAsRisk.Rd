\name{getAsRisk}
\alias{getAsRisk}
\alias{getAsRisk-methods}
\alias{getAsRisk,asMSE,UnivariateDistribution,Neighborhood,BiasType-method}
\alias{getAsRisk,asMSE,EuclRandVariable,Neighborhood,BiasType-method}
\alias{getAsRisk,asBias,UnivariateDistribution,ContNeighborhood,BiasType-method}
\alias{getAsRisk,asBias,UnivariateDistribution,ContNeighborhood,onesidedBias-method}
\alias{getAsRisk,asBias,UnivariateDistribution,ContNeighborhood,asymmetricBias-method}
\alias{getAsRisk,asBias,UnivariateDistribution,TotalVarNeighborhood,BiasType-method}
\alias{getAsRisk,asBias,RealRandVariable,ContNeighborhood,BiasType-method}
\alias{getAsRisk,asCov,UnivariateDistribution,ContNeighborhood,BiasType-method}
\alias{getAsRisk,asCov,UnivariateDistribution,TotalVarNeighborhood,BiasType-method}
\alias{getAsRisk,asCov,RealRandVariable,ContNeighborhood,BiasType-method}
\alias{getAsRisk,trAsCov,UnivariateDistribution,UncondNeighborhood,BiasType-method}
\alias{getAsRisk,trAsCov,RealRandVariable,ContNeighborhood,BiasType-method}
\alias{getAsRisk,asUnOvShoot,UnivariateDistribution,UncondNeighborhood,BiasType-method}
\alias{getAsRisk,asSemivar,UnivariateDistribution,Neighborhood,onesidedBias-method}

\title{Generic Function for Computation of Asymptotic Risks}
\description{
  Generic function for the computation of asymptotic risks.
  This function is rarely called directly. It is used by 
  other functions.
}
\usage{
getAsRisk(risk, L2deriv, neighbor, biastype, ...)

\S4method{getAsRisk}{asMSE,UnivariateDistribution,Neighborhood,BiasType}(risk, L2deriv, neighbor, biastype, clip, cent, stand, trafo)

\S4method{getAsRisk}{asMSE,EuclRandVariable,Neighborhood,BiasType}(risk, L2deriv, neighbor, biastype, clip, cent, stand, trafo)

\S4method{getAsRisk}{asBias,UnivariateDistribution,ContNeighborhood,BiasType}(risk, L2deriv, neighbor, biastype, trafo)

\S4method{getAsRisk}{asBias,UnivariateDistribution,ContNeighborhood,onesidedBias}(risk, L2deriv, neighbor, biastype, trafo)

\S4method{getAsRisk}{asBias,UnivariateDistribution,ContNeighborhood,asymmetricBias}(risk, L2deriv, neighbor, biastype, trafo)

\S4method{getAsRisk}{asBias,UnivariateDistribution,TotalVarNeighborhood,BiasType}(risk, L2deriv, neighbor, biastype, trafo)

\S4method{getAsRisk}{asBias,RealRandVariable,ContNeighborhood,BiasType}(risk, L2deriv, neighbor, biastype, Distr, L2derivDistrSymm, trafo, 
             z.start, A.start, maxiter, tol)

\S4method{getAsRisk}{asCov,UnivariateDistribution,ContNeighborhood,BiasType}(risk, L2deriv, neighbor, biastype, clip, cent, stand)

\S4method{getAsRisk}{asCov,UnivariateDistribution,TotalVarNeighborhood,BiasType}(risk, L2deriv, neighbor, biastype, clip, cent, stand)

\S4method{getAsRisk}{asCov,RealRandVariable,ContNeighborhood,BiasType}(risk, L2deriv, neighbor, biastype, Distr, clip, cent, stand)

\S4method{getAsRisk}{trAsCov,UnivariateDistribution,UncondNeighborhood,BiasType}(risk, L2deriv, neighbor, biastype, clip, cent, stand)

\S4method{getAsRisk}{trAsCov,RealRandVariable,ContNeighborhood,BiasType}(risk, L2deriv, neighbor, biastype, Distr, clip, cent, stand)

\S4method{getAsRisk}{asUnOvShoot,UnivariateDistribution,UncondNeighborhood,BiasType}(risk, L2deriv, neighbor, biastype, clip, cent, stand, trafo)

\S4method{getAsRisk}{asSemivar,UnivariateDistribution,Neighborhood,onesidedBias}(risk, L2deriv, neighbor, biastype, 
    clip, cent, stand, trafo)
}
\arguments{
  \item{risk}{ object of class \code{"asRisk"}. }
  \item{L2deriv}{ L2-derivative of some L2-differentiable family
    of probability distributions. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{biastype}{ object of class \code{"BiasType"}. }
  \item{\dots}{ additional parameters. }
  \item{clip}{ optimal clipping bound. }
  \item{cent}{ optimal centering constant. }
  \item{stand}{ standardizing matrix. }
  \item{trafo}{ matrix: transformation of the parameter. }
  \item{Distr}{ object of class \code{"Distribution"}. }
  \item{L2derivDistrSymm}{ object of class \code{"DistrSymmList"}. }
  \item{z.start}{ initial value for the centering constant. }
  \item{A.start}{ initial value for the standardizing matrix. }
  \item{maxiter}{ the maximum number of iterations }
  \item{tol}{ the desired accuracy (convergence tolerance).}
}
%\details{}
\value{The asymptotic risk is computed.}
\section{Methods}{
\describe{
  \item{risk = "asMSE", L2deriv = "UnivariateDistribution", neighbor = "Neighborhood", biastype = "BiasType":}{
    computes asymptotic mean square error in methods for
    function \code{getInfRobIC}. }

  \item{risk = "asMSE", L2deriv = "EuclRandVariable", neighbor = "Neighborhood", biastype = "BiasType":}{ 
    computes asymptotic mean square error in methods for
    function \code{getInfRobIC}. }

  \item{risk = "asBias", L2deriv = "UnivariateDistribution", neighbor = "ContNeighborhood", biastype = "BiasType":}{ 
    computes standardized asymptotic bias in methods for function \code{getInfRobIC}. }

  \item{risk = "asBias", L2deriv = "UnivariateDistribution", neighbor = "ContNeighborhood", biastype = "onesidedBias":}{ 
    computes standardized asymptotic bias in methods for function \code{getInfRobIC}. }

  \item{risk = "asBias", L2deriv = "UnivariateDistribution", neighbor = "ContNeighborhood", biastype = "asymmetricBias":}{ 
    computes standardized asymptotic bias in methods for function \code{getInfRobIC}. }

  \item{risk = "asBias", L2deriv = "UnivariateDistribution", neighbor = "TotalVarNeighborhood", biastype = "BiasType":}{ 
    computes standardized asymptotic bias in methods for function \code{getInfRobIC}. }

  \item{risk = "asBias", L2deriv = "RealRandVariable", neighbor = "ContNeighborhood", biastype = "BiasType":}{ 
    computes standardized asymptotic bias in methods for function \code{getInfRobIC}. }

  \item{risk = "asCov", L2deriv = "UnivariateDistribution", neighbor = "ContNeighborhood", biastype = "BiasType":}{ 
    computes asymptotic covariance in methods for function \code{getInfRobIC}. }

  \item{risk = "asCov", L2deriv = "UnivariateDistribution", neighbor = "TotalVarNeighborhood", biastype = "BiasType":}{ 
    computes asymptotic covariance in methods for function \code{getInfRobIC}. }

  \item{risk = "asCov", L2deriv = "RealRandVariable", neighbor = "ContNeighborhood", biastype = "BiasType":}{ 
    computes asymptotic covariance in methods for function \code{getInfRobIC}. }

  \item{risk = "trAsCov", L2deriv = "UnivariateDistribution", neighbor = "UncondNeighborhood", biastype = "BiasType":}{ 
    computes trace of asymptotic covariance in methods 
    for function \code{getInfRobIC}. }

  \item{risk = "trAsCov", L2deriv = "RealRandVariable", neighbor = "ContNeighborhood", biastype = "BiasType":}{ 
    computes trace of asymptotic covariance in methods for 
    function \code{getInfRobIC}. }

  \item{risk = "asUnOvShoot", L2deriv = "UnivariateDistribution", neighbor = "UncondNeighborhood", biastype = "BiasType":}{ 
    computes asymptotic under-/overshoot risk in methods for 
    function \code{getInfRobIC}. }

  \item{risk = "asSemivar", L2deriv = "UnivariateDistribution", neighbor = "Neighborhood", biastype = "onesidedBias":}{
    computes asymptotic semivariance in methods for
    function \code{getInfRobIC}. }
}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Ruckdeschel, P. and Rieder, H. (2004) Optimal Influence Curves for
  General Loss Functions. Statistics & Decisions \emph{22}, 201-223.

  Ruckdeschel, P. (2005) Optimally One-Sided Bounded Influence Curves.
  Mathematical Methods in Statistics \emph{14}(1), 105-131.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[distrMod]{asRisk-class}}}
%\examples{}
\concept{asymptotic risk}
\concept{risk}
\keyword{}
