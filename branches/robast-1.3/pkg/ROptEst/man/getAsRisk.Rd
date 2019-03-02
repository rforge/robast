\name{getAsRisk}
\alias{getAsRisk}
\alias{getAsRisk-methods}
\alias{getAsRisk,asMSE,UnivariateDistribution,Neighborhood,ANY-method}
\alias{getAsRisk,asL1,UnivariateDistribution,Neighborhood,ANY-method}
\alias{getAsRisk,asL4,UnivariateDistribution,Neighborhood,ANY-method}
\alias{getAsRisk,asMSE,EuclRandVariable,Neighborhood,ANY-method}
\alias{getAsRisk,asBias,UnivariateDistribution,ContNeighborhood,ANY-method}
\alias{getAsRisk,asBias,UnivariateDistribution,ContNeighborhood,onesidedBias-method}
\alias{getAsRisk,asBias,UnivariateDistribution,ContNeighborhood,asymmetricBias-method}
\alias{getAsRisk,asBias,UnivariateDistribution,TotalVarNeighborhood,ANY-method}
\alias{getAsRisk,asBias,RealRandVariable,ContNeighborhood,ANY-method}
\alias{getAsRisk,asBias,RealRandVariable,TotalVarNeighborhood,ANY-method}
\alias{getAsRisk,asCov,UnivariateDistribution,ContNeighborhood,ANY-method}
\alias{getAsRisk,asCov,UnivariateDistribution,TotalVarNeighborhood,ANY-method}
\alias{getAsRisk,asCov,RealRandVariable,ContNeighborhood,ANY-method}
\alias{getAsRisk,trAsCov,UnivariateDistribution,UncondNeighborhood,ANY-method}
\alias{getAsRisk,trAsCov,RealRandVariable,ContNeighborhood,ANY-method}
\alias{getAsRisk,asAnscombe,UnivariateDistribution,UncondNeighborhood,ANY-method}
\alias{getAsRisk,asAnscombe,RealRandVariable,ContNeighborhood,ANY-method}
\alias{getAsRisk,asUnOvShoot,UnivariateDistribution,UncondNeighborhood,ANY-method}
\alias{getAsRisk,asSemivar,UnivariateDistribution,Neighborhood,onesidedBias-method}

\title{Generic Function for Computation of Asymptotic Risks}
\description{
  Generic function for the computation of asymptotic risks.
  This function is rarely called directly. It is used by 
  other functions.
}
\usage{
getAsRisk(risk, L2deriv, neighbor, biastype, ...)

\S4method{getAsRisk}{asMSE,UnivariateDistribution,Neighborhood,ANY}(risk,
    L2deriv, neighbor, biastype, normtype = NULL, clip = NULL, cent = NULL,
    stand, trafo, ...)

\S4method{getAsRisk}{asL1,UnivariateDistribution,Neighborhood,ANY}(risk,
    L2deriv, neighbor, biastype, normtype = NULL, clip = NULL, cent = NULL,
    stand, trafo, ...)

\S4method{getAsRisk}{asL4,UnivariateDistribution,Neighborhood,ANY}(risk,
    L2deriv, neighbor, biastype, normtype = NULL, clip = NULL, cent = NULL,
    stand, trafo, ...)

\S4method{getAsRisk}{asMSE,EuclRandVariable,Neighborhood,ANY}(risk,
    L2deriv, neighbor, biastype, normtype = NULL, clip = NULL, cent = NULL,
    stand, trafo, ...)

\S4method{getAsRisk}{asBias,UnivariateDistribution,ContNeighborhood,ANY}(risk,
    L2deriv, neighbor, biastype, normtype = NULL, clip = NULL, cent = NULL,
    stand = NULL, trafo, ...)

\S4method{getAsRisk}{asBias,UnivariateDistribution,ContNeighborhood,onesidedBias}(
    risk, L2deriv, neighbor, biastype, normtype = NULL, clip = NULL, cent = NULL,
    stand = NULL, trafo, ...)

\S4method{getAsRisk}{asBias,UnivariateDistribution,ContNeighborhood,asymmetricBias}(
    risk, L2deriv, neighbor, biastype, normtype = NULL, clip = NULL, cent = NULL,
    stand = NULL, trafo, ...)

\S4method{getAsRisk}{asBias,UnivariateDistribution,TotalVarNeighborhood,ANY}(
    risk, L2deriv, neighbor, biastype, normtype = NULL, clip = NULL, cent = NULL,
    stand = NULL, trafo, ...)

\S4method{getAsRisk}{asBias,RealRandVariable,ContNeighborhood,ANY}(
    risk,L2deriv, neighbor, biastype, normtype = NULL, clip = NULL, cent = NULL,
    stand = NULL, Distr, DistrSymm, L2derivSymm,
    L2derivDistrSymm, Finfo, trafo, z.start, A.start, maxiter, tol,
    warn, verbose = NULL, ...)
\S4method{getAsRisk}{asBias,RealRandVariable,TotalVarNeighborhood,ANY}(
    risk, L2deriv, neighbor, biastype, normtype = NULL, 
    clip = NULL, cent = NULL, stand = NULL, Distr, DistrSymm, L2derivSymm,
    L2derivDistrSymm, Finfo, trafo, z.start, A.start, maxiter, tol,
    warn, verbose = NULL, ...)

\S4method{getAsRisk}{asCov,UnivariateDistribution,ContNeighborhood,ANY}(
    risk, L2deriv, neighbor, biastype, normtype = NULL, clip, cent, stand,
    trafo = NULL, ...)

\S4method{getAsRisk}{asCov,UnivariateDistribution,TotalVarNeighborhood,ANY}(
    risk, L2deriv, neighbor, biastype, normtype = NULL, clip, cent, stand,
    trafo = NULL, ...)

\S4method{getAsRisk}{asCov,RealRandVariable,ContNeighborhood,ANY}(risk,
    L2deriv, neighbor, biastype, normtype = NULL, clip = NULL, cent, stand,
    Distr, trafo = NULL, V.comp =  matrix(TRUE, ncol = nrow(stand),
    nrow = nrow(stand)), w, ...)

\S4method{getAsRisk}{trAsCov,UnivariateDistribution,UncondNeighborhood,ANY}(
    risk, L2deriv, neighbor, biastype, normtype = NULL, clip, cent, stand,
    trafo = NULL, ...)

\S4method{getAsRisk}{trAsCov,RealRandVariable,ContNeighborhood,ANY}(risk,
    L2deriv, neighbor, biastype, normtype, clip, cent, stand, Distr,
    trafo = NULL,  V.comp =  matrix(TRUE, ncol = nrow(stand),
    nrow = nrow(stand)), w, ...)
  
\S4method{getAsRisk}{asAnscombe,UnivariateDistribution,UncondNeighborhood,ANY}(
    risk, L2deriv, neighbor, biastype, normtype = NULL, clip, cent, stand,
    trafo = NULL, FI, ...)

\S4method{getAsRisk}{asAnscombe,RealRandVariable,ContNeighborhood,ANY}(risk,
    L2deriv, neighbor, biastype, normtype, clip, cent, stand, Distr, trafo = NULL, 
    V.comp =  matrix(TRUE, ncol = nrow(stand), nrow = nrow(stand)),
    FI, w, ...)
  

\S4method{getAsRisk}{asUnOvShoot,UnivariateDistribution,UncondNeighborhood,ANY}(
    risk, L2deriv, neighbor, biastype, normtype = NULL, clip, cent, stand,
    trafo, ...)

\S4method{getAsRisk}{asSemivar,UnivariateDistribution,Neighborhood,onesidedBias}(
    risk, L2deriv, neighbor, biastype, normtype = NULL, clip, cent, stand,
    trafo, ...)
}
\arguments{
  \item{risk}{ object of class \code{"asRisk"}. }
  \item{L2deriv}{ L2-derivative of some L2-differentiable family
    of probability distributions. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{biastype}{ object of class \code{"ANY"}. }
  \item{\dots}{ additional parameters; often used to enable flexible calls. }
  \item{clip}{ optimal clipping bound. }
  \item{cent}{ optimal centering constant. }
  \item{stand}{ standardizing matrix. }
  \item{Finfo}{ matrix: the Fisher Information of the parameter. }
  \item{trafo}{ matrix: transformation of the parameter. }
  \item{Distr}{ object of class \code{"Distribution"}. }
  \item{DistrSymm}{ object of class \code{"DistributionSymmetry"}. }
  \item{L2derivSymm}{ object of class \code{"FunSymmList"}. }
  \item{L2derivDistrSymm}{ object of class \code{"DistrSymmList"}. }
  \item{z.start}{ initial value for the centering constant. }
  \item{A.start}{ initial value for the standardizing matrix. }
  \item{maxiter}{ the maximum number of iterations }
  \item{tol}{ the desired accuracy (convergence tolerance).}
  \item{warn}{ logical: print warnings. }
  \item{normtype}{ object of class \code{"NormType"}. }
  \item{V.comp}{ matrix: indication which components of the standardizing
    matrix have to be computed. }
  \item{w}{object of class \code{RobWeight}; current weight}
  \item{FI}{trace of the respective Fisher Information}
  \item{verbose}{logical: if \code{TRUE} some diagnostics are printed out.}
}
\details{ This function is rarely called directly. It is used by 
  other functions/methods. }
\value{The asymptotic risk is computed.}
\section{Methods}{
\describe{
  \item{risk = "asMSE", L2deriv = "UnivariateDistribution",
        neighbor = "Neighborhood", biastype = "ANY":}{
    computes asymptotic mean square error in methods for
    function \code{getInfRobIC}. }

  \item{risk = "asL1", L2deriv = "UnivariateDistribution",
        neighbor = "Neighborhood", biastype = "ANY":}{
    computes asymptotic mean absolute error in methods for
    function \code{getInfRobIC}. }

  \item{risk = "asL4", L2deriv = "UnivariateDistribution",
        neighbor = "Neighborhood", biastype = "ANY":}{
    computes asymptotic mean power 4 error in methods for
    function \code{getInfRobIC}. }

  \item{risk = "asMSE", L2deriv = "EuclRandVariable",
        neighbor = "Neighborhood", biastype = "ANY":}{
    computes asymptotic mean square error in methods for
    function \code{getInfRobIC}. }

  \item{risk = "asBias", L2deriv = "UnivariateDistribution",
        neighbor = "ContNeighborhood", biastype = "ANY":}{
    computes standardized asymptotic bias in methods
    for function \code{getInfRobIC}. }

  \item{risk = "asBias", L2deriv = "UnivariateDistribution",
        neighbor = "ContNeighborhood", biastype = "onesidedBias":}{
    computes standardized asymptotic bias in methods for function
    \code{getInfRobIC}. }

  \item{risk = "asBias", L2deriv = "UnivariateDistribution",
        neighbor = "ContNeighborhood", biastype = "asymmetricBias":}{
    computes standardized asymptotic bias in methods for function
    \code{getInfRobIC}. }

  \item{risk = "asBias", L2deriv = "UnivariateDistribution",
        neighbor = "TotalVarNeighborhood", biastype = "ANY":}{
    computes standardized asymptotic bias in methods for function
    \code{getInfRobIC}. }

  \item{risk = "asBias", L2deriv = "RealRandVariable",
        neighbor = "ContNeighborhood", biastype = "ANY":}{
    computes standardized asymptotic bias in methods for function
    \code{getInfRobIC}. }

  \item{risk = "asCov", L2deriv = "UnivariateDistribution",
        neighbor = "ContNeighborhood", biastype = "ANY":}{
    computes asymptotic covariance in methods for function
    \code{getInfRobIC}. }

  \item{risk = "asCov", L2deriv = "UnivariateDistribution",
        neighbor = "TotalVarNeighborhood", biastype = "ANY":}{
    computes asymptotic covariance in methods for function
    \code{getInfRobIC}. }

  \item{risk = "asCov", L2deriv = "RealRandVariable",
        neighbor = "ContNeighborhood", biastype = "ANY":}{
    computes asymptotic covariance in methods for function
    \code{getInfRobIC}. }

  \item{risk = "trAsCov", L2deriv = "UnivariateDistribution",
        neighbor = "UncondNeighborhood", biastype = "ANY":}{
    computes trace of asymptotic covariance in methods 
    for function \code{getInfRobIC}. }

  \item{risk = "trAsCov", L2deriv = "RealRandVariable",
        neighbor = "ContNeighborhood", biastype = "ANY":}{
    computes trace of asymptotic covariance in methods for 
    function \code{getInfRobIC}. }

  \item{risk = "asAnscombe", L2deriv = "UnivariateDistribution",
        neighbor = "UncondNeighborhood", biastype = "ANY":}{
    computes the ARE in the ideal model in methods 
    for function \code{getInfRobIC}. }

  \item{risk = "asAnscombe", L2deriv = "RealRandVariable",
        neighbor = "ContNeighborhood", biastype = "ANY":}{
    computes the ARE in the ideal model in methods for 
    function \code{getInfRobIC}. }

  \item{risk = "asUnOvShoot", L2deriv = "UnivariateDistribution",
        neighbor = "UncondNeighborhood", biastype = "ANY":}{
    computes asymptotic under-/overshoot risk in methods for 
    function \code{getInfRobIC}. }

  \item{risk = "asSemivar", L2deriv = "UnivariateDistribution",
        neighbor = "Neighborhood", biastype = "onesidedBias":}{
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
\keyword{robust}
