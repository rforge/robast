\name{getInfRobIC}
\alias{getInfRobIC}
\alias{getInfRobIC-methods}
\alias{getInfRobIC,UnivariateDistribution,asCov,ContNeighborhood-method}
\alias{getInfRobIC,UnivariateDistribution,asCov,TotalVarNeighborhood-method}
\alias{getInfRobIC,RealRandVariable,asCov,ContNeighborhood-method}
\alias{getInfRobIC,UnivariateDistribution,asBias,ContNeighborhood-method}
\alias{getInfRobIC,UnivariateDistribution,asBias,TotalVarNeighborhood-method}
\alias{getInfRobIC,RealRandVariable,asBias,ContNeighborhood-method}
\alias{getInfRobIC,UnivariateDistribution,asHampel,UncondNeighborhood-method}
\alias{getInfRobIC,RealRandVariable,asHampel,ContNeighborhood-method}
\alias{getInfRobIC,UnivariateDistribution,asGRisk,UncondNeighborhood-method}
\alias{getInfRobIC,RealRandVariable,asGRisk,ContNeighborhood-method}
\alias{getInfRobIC,UnivariateDistribution,asUnOvShoot,UncondNeighborhood-method}

\title{ Generic Function for the Computation of Optimally Robust ICs }
\description{
  Generic function for the computation of optimally robust ICs 
  in case of infinitesimal robust models. This function is 
  rarely called directly. 
}
\usage{
getInfRobIC(L2deriv, risk, neighbor, ...)

\S4method{getInfRobIC}{UnivariateDistribution,asCov,ContNeighborhood}(L2deriv, risk, neighbor, Finfo, trafo)

\S4method{getInfRobIC}{UnivariateDistribution,asCov,TotalVarNeighborhood}(L2deriv, risk, neighbor, Finfo, trafo)

\S4method{getInfRobIC}{RealRandVariable,asCov,ContNeighborhood}(L2deriv, risk, neighbor, Distr, Finfo, trafo)

\S4method{getInfRobIC}{UnivariateDistribution,asBias,ContNeighborhood}(L2deriv, risk, neighbor, symm, Finfo, trafo, 
             upper, maxiter, tol, warn)

\S4method{getInfRobIC}{UnivariateDistribution,asBias,TotalVarNeighborhood}(L2deriv, risk, neighbor, symm, Finfo, trafo, 
             upper, maxiter, tol, warn)

\S4method{getInfRobIC}{RealRandVariable,asBias,ContNeighborhood}(L2deriv, risk, neighbor, Distr, DistrSymm, L2derivSymm, 
             L2derivDistrSymm, Finfo, z.start, A.start, trafo, upper, maxiter, tol, warn)

\S4method{getInfRobIC}{UnivariateDistribution,asHampel,UncondNeighborhood}(L2deriv, risk, neighbor, symm, Finfo, trafo, 
             upper, maxiter, tol, warn)

\S4method{getInfRobIC}{RealRandVariable,asHampel,ContNeighborhood}(L2deriv, risk, neighbor, Distr, DistrSymm, L2derivSymm, 
             L2derivDistrSymm, Finfo, trafo, z.start, A.start, upper, maxiter, tol, warn)

\S4method{getInfRobIC}{UnivariateDistribution,asGRisk,UncondNeighborhood}(L2deriv, risk, neighbor, symm, Finfo, trafo, 
             upper, maxiter, tol, warn)

\S4method{getInfRobIC}{RealRandVariable,asGRisk,ContNeighborhood}(L2deriv, risk, neighbor, Distr, DistrSymm, L2derivSymm, 
             L2derivDistrSymm, Finfo, trafo, z.start, A.start, upper, maxiter, tol, warn)

\S4method{getInfRobIC}{UnivariateDistribution,asUnOvShoot,UncondNeighborhood}(L2deriv, risk, neighbor, symm, Finfo, trafo, 
             upper, maxiter, tol, warn)
}
\arguments{
  \item{L2deriv}{ L2-derivative of some L2-differentiable family 
    of probability measures. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters. }
  \item{Distr}{ object of class \code{"Distribution"}. }
  \item{symm}{ logical: indicating symmetry of \code{L2deriv}. }
  \item{DistrSymm}{ object of class \code{"DistributionSymmetry"}. }
  \item{L2derivSymm}{ object of class \code{"FunSymmList"}. }
  \item{L2derivDistrSymm}{ object of class \code{"DistrSymmList"}. }
  \item{Finfo}{ Fisher information matrix. }
  \item{z.start}{ initial value for the centering constant. }
  \item{A.start}{ initial value for the standardizing matrix. }
  \item{trafo}{ matrix: transformation of the parameter. }
  \item{upper}{ upper bound for the optimal clipping bound. }
  \item{maxiter}{ the maximum number of iterations. }
  \item{tol}{ the desired accuracy (convergence tolerance).}
  \item{warn}{ logical: print warnings. }
}
%\details{}
\value{The optimally robust IC is computed.}
\section{Methods}{
\describe{
  \item{L2deriv = "UnivariateDistribution", risk = "asCov", 
        neighbor = "ContNeighborhood"}{ 
    computes the classical optimal influence curve for L2 differentiable 
    parametric families with unknown one-dimensional parameter. }

  \item{L2deriv = "UnivariateDistribution", risk = "asCov", 
        neighbor = "TotalVarNeighborhood"}{
    computes the classical optimal influence curve for L2 differentiable 
    parametric families with unknown one-dimensional parameter. }

  \item{L2deriv = "RealRandVariable", risk = "asCov", 
        neighbor = "ContNeighborhood"}{ 
    computes the classical optimal influence curve for L2 differentiable 
    parametric families with unknown \eqn{k}-dimensional parameter 
    (\eqn{k > 1}) where the underlying distribution is univariate. }

  \item{L2deriv = "UnivariateDistribution", risk = "asBias", 
        neighbor = "ContNeighborhood"}{ 
    computes the bias optimal influence curve for L2 differentiable 
    parametric families with unknown one-dimensional parameter. }

  \item{L2deriv = "UnivariateDistribution", risk = "asBias", 
        neighbor = "TotalVarNeighborhood"}{ 
    computes the bias optimal influence curve for L2 differentiable 
    parametric families with unknown one-dimensional parameter. }

  \item{L2deriv = "RealRandVariable", risk = "asBias", 
        neighbor = "ContNeighborhood"}{ 
    computes the bias optimal influence curve for L2 differentiable 
    parametric families with unknown \eqn{k}-dimensional parameter 
    (\eqn{k > 1}) where the underlying distribution is univariate. }

  \item{L2deriv = "UnivariateDistribution", risk = "asHampel", 
        neighbor = "UncondNeighborhood"}{
    computes the optimally robust influence curve for L2 differentiable 
    parametric families with unknown one-dimensional parameter. }

  \item{L2deriv = "RealRandVariable", risk = "asHampel", 
        neighbor = "ContNeighborhood"}{ 
    computes the optimally robust influence curve for L2 differentiable 
    parametric families with unknown \eqn{k}-dimensional parameter 
    (\eqn{k > 1}) where the underlying distribution is univariate. }

  \item{L2deriv = "UnivariateDistribution", risk = "asGRisk", 
        neighbor = "UncondNeighborhood"}{
    computes the optimally robust influence curve for L2 differentiable 
    parametric families with unknown one-dimensional parameter. }

  \item{L2deriv = "RealRandVariable", risk = "asGRisk", 
        neighbor = "ContNeighborhood"}{ 
    computes the optimally robust influence curve for L2 differentiable 
    parametric families with unknown \eqn{k}-dimensional parameter 
    (\eqn{k > 1}) where the underlying distribution is univariate. }

  \item{L2deriv = "UnivariateDistribution", risk = "asUnOvShoot", 
        neighbor = "UncondNeighborhood"}{
    computes the optimally robust influence curve for one-dimensional
    L2 differentiable parametric families and 
    asymptotic under-/overshoot risk. }
}}
\references{
  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{InfRobModel-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
