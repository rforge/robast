\name{minmaxBias}
\alias{minmaxBias}
\alias{minmaxBias-methods}
\alias{minmaxBias,UnivariateDistribution,ContNeighborhood,BiasType-method}
\alias{minmaxBias,UnivariateDistribution,ContNeighborhood,onesidedBias-method}
\alias{minmaxBias,UnivariateDistribution,ContNeighborhood,asymmetricBias-method}
\alias{minmaxBias,UnivariateDistribution,TotalVarNeighborhood,BiasType-method}
\alias{minmaxBias,RealRandVariable,ContNeighborhood,BiasType-method}
\alias{minmaxBias,RealRandVariable,TotalVarNeighborhood,BiasType-method}

\title{ Generic Function for the Computation of Bias-Optimally Robust ICs }
\description{
  Generic function for the computation of bias-optimally robust ICs 
  in case of infinitesimal robust models. This function is 
  rarely called directly. 
}
\usage{
minmaxBias(L2deriv, neighbor, biastype, ...)

\S4method{minmaxBias}{UnivariateDistribution,ContNeighborhood,BiasType}(L2deriv,
     neighbor, biastype, symm, trafo, maxiter, tol, warn, Finfo, verbose = NULL)

\S4method{minmaxBias}{UnivariateDistribution,ContNeighborhood,asymmetricBias}(
     L2deriv, neighbor, biastype, symm, trafo, maxiter, tol, warn, Finfo, verbose = NULL)

\S4method{minmaxBias}{UnivariateDistribution,ContNeighborhood,onesidedBias}(
     L2deriv, neighbor, biastype, symm, trafo, maxiter, tol, warn, Finfo, verbose = NULL)

\S4method{minmaxBias}{UnivariateDistribution,TotalVarNeighborhood,BiasType}(
     L2deriv, neighbor, biastype, symm, trafo, maxiter, tol, warn, Finfo, verbose = NULL)

\S4method{minmaxBias}{RealRandVariable,ContNeighborhood,BiasType}(L2deriv,
     neighbor, biastype, normtype, Distr, z.start, A.start,  z.comp, A.comp,
     Finfo, trafo, maxiter, tol, verbose = NULL, ...)

\S4method{minmaxBias}{RealRandVariable,TotalVarNeighborhood,BiasType}(L2deriv,
     neighbor, biastype, normtype, Distr, z.start, A.start,  z.comp, A.comp,
     Finfo, trafo, maxiter, tol, verbose = NULL, ...)
}
\arguments{
  \item{L2deriv}{ L2-derivative of some L2-differentiable family 
    of probability measures. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{biastype}{ object of class \code{"BiasType"}. }
  \item{normtype}{ object of class \code{"NormType"}. }
  \item{\dots}{ additional arguments to be passed to \code{E} }
  \item{Distr}{ object of class \code{"Distribution"}. }
  \item{symm}{ logical: indicating symmetry of \code{L2deriv}. }
  \item{z.start}{ initial value for the centering constant. }
  \item{A.start}{ initial value for the standardizing matrix. }
  \item{z.comp}{ \code{logical} indicator which indices need to be computed and which are 0 due to symmetry. }
  \item{A.comp}{ \code{matrix} of \code{logical} indicator which indices need to be computed and which are 0 due to symmetry.}
  \item{trafo}{ matrix: transformation of the parameter. }
  \item{maxiter}{ the maximum number of iterations. }
  \item{tol}{ the desired accuracy (convergence tolerance).}
  \item{warn}{ logical: print warnings. }
  \item{Finfo}{ Fisher information matrix. }
  \item{verbose}{ logical: if \code{TRUE}, some messages are printed }
}
%\details{}
\value{The bias-optimally robust IC is computed.}
\section{Methods}{
\describe{
  \item{L2deriv = "UnivariateDistribution", neighbor = "ContNeighborhood", 
        biastype = "BiasType"}{ 
    computes the bias optimal influence curve for symmetric bias for L2 differentiable 
    parametric families with unknown one-dimensional parameter. }

  \item{L2deriv = "UnivariateDistribution", neighbor = "ContNeighborhood", 
        biastype = "asymmetricBias"}{ 
    computes the bias optimal influence curve for asymmetric bias for L2 differentiable 
    parametric families with unknown one-dimensional parameter. }
    
  \item{L2deriv = "UnivariateDistribution", neighbor = "TotalVarNeighborhood", 
        biastype = "BiasType"}{ 
    computes the bias optimal influence curve for symmetric bias for L2 differentiable 
    parametric families with unknown one-dimensional parameter. }

  \item{L2deriv = "RealRandVariable", neighbor = "ContNeighborhood", 
        biastype = "BiasType"}{ 
    computes the bias optimal influence curve for symmetric bias for L2 differentiable 
    parametric families with unknown \eqn{k}-dimensional parameter 
    (\eqn{k > 1}) where the underlying distribution is univariate. }

  \item{L2deriv = "RealRandVariable", neighbor = "TotalNeighborhood",
        biastype = "BiasType"}{
    computes the bias optimal influence curve for symmetric bias for L2 differentiable
    parametric families in a setting where we are interested in a \eqn{p=1}
    dimensional aspect of an unknown \eqn{k}-dimensional parameter
    (\eqn{k > 1}) where the underlying distribution is univariate. }
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
        Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link[RobAStBase]{InfRobModel-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
