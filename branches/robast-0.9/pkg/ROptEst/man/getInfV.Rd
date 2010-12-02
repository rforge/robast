\name{getInfV}
\alias{getInfV}
\alias{getInfV-methods}
\alias{getInfV,UnivariateDistribution,ContNeighborhood,BiasType-method}
\alias{getInfV,UnivariateDistribution,TotalVarNeighborhood,BiasType-method}
\alias{getInfV,RealRandVariable,ContNeighborhood,BiasType-method}
\alias{getInfV,RealRandVariable,TotalVarNeighborhood,BiasType-method}
\alias{getInfV,UnivariateDistribution,ContNeighborhood,onesidedBias-method}
\alias{getInfV,UnivariateDistribution,ContNeighborhood,asymmetricBias-method}

\title{Generic Function for the Computation of the asymptotic Variance of a Hampel type IC}
\description{
  Generic function for the computation of the optimal clipping bound
  in case of infinitesimal robust models. This function is rarely called 
  directly. It is used to compute optimally robust ICs.
}
\usage{
getInfV(L2deriv, neighbor, biastype, ...)
\S4method{getInfV}{UnivariateDistribution,ContNeighborhood,BiasType}(L2deriv, 
         neighbor, biastype, clip, cent, stand)
\S4method{getInfV}{UnivariateDistribution,TotalVarNeighborhood,BiasType}(L2deriv, 
         neighbor, biastype, clip, cent, stand)
\S4method{getInfV}{RealRandVariable,ContNeighborhood,BiasType}(L2deriv, 
         neighbor, biastype, Distr, V.comp, cent, stand, 
         w)
\S4method{getInfV}{RealRandVariable,TotalVarNeighborhood,BiasType}(L2deriv,
         neighbor, biastype, Distr, V.comp, cent, stand,
         w)
\S4method{getInfV}{UnivariateDistribution,ContNeighborhood,onesidedBias}(L2deriv,
         neighbor, biastype, clip, cent, stand)
\S4method{getInfV}{UnivariateDistribution,ContNeighborhood,asymmetricBias}(L2deriv, 
         neighbor, biastype, clip, cent, stand)
}
\arguments{
  \item{L2deriv}{ L2-derivative of some L2-differentiable family 
    of probability measures. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{biastype}{ object of class \code{"BiasType"} }
  \item{\dots}{ additional parameters. }
  \item{clip}{ positive real: clipping bound }
  \item{cent}{ optimal centering constant. }
  \item{stand}{ standardizing matrix. }
  \item{Distr}{ standardizing matrix. }
  \item{V.comp}{ matrix: indication which components of the standardizing
    matrix have to be computed. }
  \item{w}{object of class \code{RobWeight}; current weight}  
}
%\details{}
\value{The asymptotic variance of an ALE to IC of Hampel type is computed.}
\references{
  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Ruckdeschel, P. (2005) Optimally One-Sided Bounded Influence Curves.
  Mathematical Methods in Statistics \emph{14}(1), 105-131.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link[RobAStBase]{ContIC-class}}, \code{\link[RobAStBase]{TotalVarIC-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
