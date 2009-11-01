\name{getInfCent}
\alias{getInfCent}
\alias{getInfCent-methods}
\alias{getInfCent,UnivariateDistribution,ContNeighborhood,BiasType-method}
\alias{getInfCent,UnivariateDistribution,TotalVarNeighborhood,BiasType-method}
\alias{getInfCent,RealRandVariable,ContNeighborhood,BiasType-method}
\alias{getInfCent,RealRandVariable,TotalVarNeighborhood,BiasType-method}
\alias{getInfCent,UnivariateDistribution,ContNeighborhood,onesidedBias-method}
\alias{getInfCent,UnivariateDistribution,ContNeighborhood,asymmetricBias-method}

\title{Generic Function for the Computation of the Optimal Centering Constant/Lower Clipping Bound}
\description{
  Generic function for the computation of the optimal centering constant
  (contamination neighborhoods) respectively, of the optimal lower clipping
  bound (total variation neighborhood).
  This function is rarely called directly. It is used to 
  compute optimally robust ICs.
}
\usage{
getInfCent(L2deriv, neighbor, biastype, ...)

\S4method{getInfCent}{UnivariateDistribution,ContNeighborhood,BiasType}(L2deriv, 
     neighbor, biastype, clip, cent, tol.z, symm, trafo)

\S4method{getInfCent}{UnivariateDistribution,TotalVarNeighborhood,BiasType}(L2deriv, 
     neighbor, biastype, clip, cent, tol.z, symm, trafo)

\S4method{getInfCent}{RealRandVariable,ContNeighborhood,BiasType}(L2deriv,
     neighbor, biastype, Distr, z.comp, w, tol.z = .Machine$double.eps^.5)

\S4method{getInfCent}{RealRandVariable,TotalVarNeighborhood,BiasType}(L2deriv,
     neighbor, biastype, Distr, z.comp, w, tol.z = .Machine$double.eps^.5)

\S4method{getInfCent}{UnivariateDistribution,ContNeighborhood,onesidedBias}(L2deriv,
     neighbor, biastype, clip, cent, tol.z, symm, trafo)

\S4method{getInfCent}{UnivariateDistribution,ContNeighborhood,asymmetricBias}(L2deriv, 
     neighbor, biastype, clip, cent, tol.z, symm, trafo)
}
\arguments{
  \item{L2deriv}{ L2-derivative of some L2-differentiable family 
    of probability measures. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{biastype}{ object of class \code{"BiasType"} }
  \item{\dots}{ additional parameters. }
  \item{clip}{ optimal clipping bound. }
  \item{cent}{ optimal centering constant. }
  \item{tol.z}{ the desired accuracy (convergence tolerance). }
  \item{symm}{ logical: indicating symmetry of \code{L2deriv}. }
  \item{trafo}{ matrix: transformation of the parameter. }
  \item{Distr}{object of class \code{Distribution}. }
  \item{z.comp}{ logical vector: indication which components of the 
    centering constant have to be computed. }
  \item{w}{object of class \code{RobWeight}; current weight}
}
%\details{}
\value{The optimal centering constant is computed.}
\section{Methods}{
\describe{
  \item{L2deriv = "UnivariateDistribution", neighbor = "ContNeighborhood", biastype = "BiasType"}{ 
    computation of optimal centering constant for symmetric bias. }

  \item{L2deriv = "UnivariateDistribution", neighbor = "TotalVarNeighborhood", biastype = "BiasType"}{ 
    computation of optimal lower clipping bound for symmetric bias. }

  \item{L2deriv = "RealRandVariable", neighbor = "TotalVarNeighborhood", biastype = "BiasType"}{
    computation of optimal centering constant for symmetric bias. }

  \item{L2deriv = "RealRandVariable", neighbor = "ContNeighborhood", biastype = "BiasType"}{
    computation of optimal centering constant for symmetric bias. }

  \item{L2deriv = "UnivariateDistribution", neighbor = "ContNeighborhood", biastype = "onesidedBias"}{
    computation of optimal centering constant for onesided bias. }

  \item{L2deriv = "UnivariateDistribution", neighbor = "ContNeighborhood", biastype = "asymmetricBias"}{ 
    computation of optimal centering constant for asymmetric bias. }
}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Ruckdeschel, P. (2005) Optimally One-Sided Bounded Influence Curves.
  Mathematical Methods in Statistics \emph{14}(1), 105-131.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},
        Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link[RobAStBase]{ContIC-class}}, \code{\link[RobAStBase]{TotalVarIC-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
