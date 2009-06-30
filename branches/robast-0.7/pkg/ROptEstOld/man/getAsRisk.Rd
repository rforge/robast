\name{getAsRisk}
\alias{getAsRisk}
\alias{getAsRisk-methods}
\alias{getAsRisk,asMSE,UnivariateDistribution,Neighborhood-method}
\alias{getAsRisk,asMSE,EuclRandVariable,Neighborhood-method}
\alias{getAsRisk,asBias,UnivariateDistribution,ContNeighborhood-method}
\alias{getAsRisk,asBias,UnivariateDistribution,TotalVarNeighborhood-method}
\alias{getAsRisk,asBias,RealRandVariable,ContNeighborhood-method}
\alias{getAsRisk,asCov,UnivariateDistribution,ContNeighborhood-method}
\alias{getAsRisk,asCov,UnivariateDistribution,TotalVarNeighborhood-method}
\alias{getAsRisk,asCov,RealRandVariable,ContNeighborhood-method}
\alias{getAsRisk,trAsCov,UnivariateDistribution,UncondNeighborhood-method}
\alias{getAsRisk,trAsCov,RealRandVariable,ContNeighborhood-method}
\alias{getAsRisk,asUnOvShoot,UnivariateDistribution,UncondNeighborhood-method}

\title{Generic Function for Computation of Asymptotic Risks}
\description{
  Generic function for the computation of asymptotic risks.
  This function is rarely called directly. It is used by 
  other functions.
}
\usage{
getAsRisk(risk, L2deriv, neighbor, ...)

\S4method{getAsRisk}{asMSE,UnivariateDistribution,Neighborhood}(risk, L2deriv, neighbor, clip, cent, stand, trafo)

\S4method{getAsRisk}{asMSE,EuclRandVariable,Neighborhood}(risk, L2deriv, neighbor, clip, cent, stand, trafo)

\S4method{getAsRisk}{asBias,UnivariateDistribution,ContNeighborhood}(risk, L2deriv, neighbor, trafo)

\S4method{getAsRisk}{asBias,UnivariateDistribution,TotalVarNeighborhood}(risk, L2deriv, neighbor, trafo)

\S4method{getAsRisk}{asBias,RealRandVariable,ContNeighborhood}(risk, L2deriv, neighbor, Distr, L2derivDistrSymm, trafo, 
             z.start, A.start, maxiter, tol)

\S4method{getAsRisk}{asCov,UnivariateDistribution,ContNeighborhood}(risk, L2deriv, neighbor, clip, cent, stand)

\S4method{getAsRisk}{asCov,UnivariateDistribution,TotalVarNeighborhood}(risk, L2deriv, neighbor, clip, cent, stand)

\S4method{getAsRisk}{asCov,RealRandVariable,ContNeighborhood}(risk, L2deriv, neighbor, Distr, clip, cent, stand)

\S4method{getAsRisk}{trAsCov,UnivariateDistribution,UncondNeighborhood}(risk, L2deriv, neighbor, clip, cent, stand)

\S4method{getAsRisk}{trAsCov,RealRandVariable,ContNeighborhood}(risk, L2deriv, neighbor, Distr, clip, cent, stand)

\S4method{getAsRisk}{asUnOvShoot,UnivariateDistribution,UncondNeighborhood}(risk, L2deriv, neighbor, clip, cent, stand, trafo)
}
\arguments{
  \item{risk}{ object of class \code{"asRisk"}. }
  \item{L2deriv}{ L2-derivative of some L2-differentiable family
    of probability distributions. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
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
  \item{risk = "asMSE", L2deriv = "UnivariateDistribution", neighbor = "Neighborhood":}{
    computes asymptotic mean square error in methods for
    function \code{getInfRobIC}. }

  \item{risk = "asMSE", L2deriv = "EuclRandVariable", neighbor = "Neighborhood":}{ 
    computes asymptotic mean square error in methods for
    function \code{getInfRobIC}. }

  \item{risk = "asBias", L2deriv = "UnivariateDistribution", neighbor = "ContNeighborhood":}{ 
    computes standardized asymptotic bias in methods for function \code{getInfRobIC}. }

  \item{risk = "asBias", L2deriv = "UnivariateDistribution", neighbor = "TotalVarNeighborhood":}{ 
    computes standardized asymptotic bias in methods for function \code{getInfRobIC}. }

  \item{risk = "asBias", L2deriv = "RealRandVariable", neighbor = "ContNeighborhood":}{ 
    computes standardized asymptotic bias in methods for function \code{getInfRobIC}. }

  \item{risk = "asCov", L2deriv = "UnivariateDistribution", neighbor = "ContNeighborhood":}{ 
    computes asymptotic covariance in methods for function \code{getInfRobIC}. }

  \item{risk = "asCov", L2deriv = "UnivariateDistribution", neighbor = "TotalVarNeighborhood":}{ 
    computes asymptotic covariance in methods for function \code{getInfRobIC}. }

  \item{risk = "asCov", L2deriv = "RealRandVariable", neighbor = "ContNeighborhood":}{ 
    computes asymptotic covariance in methods for function \code{getInfRobIC}. }

  \item{risk = "trAsCov", L2deriv = "UnivariateDistribution", neighbor = "UncondNeighborhood":}{ 
    computes trace of asymptotic covariance in methods 
    for function \code{getInfRobIC}. }

  \item{risk = "trAsCov", L2deriv = "RealRandVariable", neighbor = "ContNeighborhood":}{ 
    computes trace of asymptotic covariance in methods for 
    function \code{getInfRobIC}. }

  \item{risk = "asUnOvShoot", L2deriv = "UnivariateDistribution", neighbor = "UncondNeighborhood":}{ 
    computes asymptotic under-/overshoot risk in methods for 
    function \code{getInfRobIC}. }
}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Ruckdeschel, P. and Rieder, H. (2004) Optimal Influence Curves for
  General Loss Functions. Statistics & Decisions (submitted).

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{asRisk-class}}}
%\examples{}
\concept{asymptotic risk}
\concept{risk}
\keyword{robust}
