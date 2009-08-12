\name{getInfStand}
\alias{getInfStand}
\alias{getInfStand-methods}
\alias{getInfStand,UnivariateDistribution,ContNeighborhood,BiasType-method}
\alias{getInfStand,UnivariateDistribution,TotalVarNeighborhood,BiasType-method}
\alias{getInfStand,RealRandVariable,UncondNeighborhood,BiasType-method}
\alias{getInfStand,UnivariateDistribution,ContNeighborhood,onesidedBias-method}
\alias{getInfStand,UnivariateDistribution,ContNeighborhood,asymmetricBias-method}

\title{ Generic Function for the Computation of the Standardizing Matrix }
\description{
  Generic function for the computation of the standardizing matrix which
  takes care of the Fisher consistency of the corresponding IC. This function 
  is rarely called directly. It is used to compute optimally robust ICs.
}
\usage{
getInfStand(L2deriv, neighbor, biastype, ...)

\S4method{getInfStand}{UnivariateDistribution,ContNeighborhood,BiasType}(L2deriv, 
     neighbor, biastype, clip, cent, trafo)

\S4method{getInfStand}{UnivariateDistribution,TotalVarNeighborhood,BiasType}(L2deriv, 
     neighbor, biastype, clip, cent, trafo)

\S4method{getInfStand}{RealRandVariable,UncondNeighborhood,BiasType}(L2deriv,
     neighbor, biastype, Distr, A.comp, cent, trafo, w)

\S4method{getInfStand}{UnivariateDistribution,ContNeighborhood,BiasType}(L2deriv, 
     neighbor, biastype, clip, cent, trafo)

\S4method{getInfStand}{UnivariateDistribution,ContNeighborhood,BiasType}(L2deriv, 
     neighbor, biastype, clip, cent, trafo)
}
\arguments{
  \item{L2deriv}{ L2-derivative of some L2-differentiable family 
    of probability measures. }
  \item{neighbor}{ object of class \code{"Neighborhood"} }
  \item{biastype}{ object of class \code{"BiasType"} }
  \item{\dots}{ additional parameters }
  \item{clip}{ optimal clipping bound. }
  \item{cent}{ optimal centering constant. }
  \item{Distr}{ object of class \code{"Distribution"}. }
  \item{trafo}{ matrix: transformation of the parameter. }  
  \item{A.comp}{ matrix: indication which components of the standardizing
    matrix have to be computed. }
  \item{w}{object of class \code{RobWeight}; current weight}  
}
%\details{}
\value{The standardizing matrix is computed.}
\section{Methods}{
\describe{
  \item{L2deriv = "UnivariateDistribution", neighbor = "ContNeighborhood", 
        biastype = "BiasType"}{ 
    computes standardizing matrix for symmetric bias. }

  \item{L2deriv = "UnivariateDistribution", neighbor = "TotalVarNeighborhood", 
        biastype = "BiasType"}{ 
    computes standardizing matrix for symmetric bias. }

  \item{L2deriv = "RealRandVariable", neighbor = "UncondNeighborhood",
        biastype = "BiasType"}{ 
    computes standardizing matrix for symmetric bias. }

  \item{L2deriv = "UnivariateDistribution", neighbor = "ContNeighborhood", 
        biastype = "onesidedBias"}{ 
    computes standardizing matrix for onesided bias. }

  \item{L2deriv = "UnivariateDistribution", neighbor = "ContNeighborhood", 
        biastype = "asymmetricBias"}{ 
    computes standardizing matrix for asymmetric bias. }
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
