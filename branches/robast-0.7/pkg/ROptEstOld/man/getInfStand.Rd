\name{getInfStand}
\alias{getInfStand}
\alias{getInfStand-methods}
\alias{getInfStand,UnivariateDistribution,ContNeighborhood-method}
\alias{getInfStand,UnivariateDistribution,TotalVarNeighborhood-method}
\alias{getInfStand,RealRandVariable,ContNeighborhood-method}

\title{ Generic Function for the Computation of the Standardizing Matrix }
\description{
  Generic function for the computation of the standardizing matrix which
  takes care of the Fisher consistency of the corresponding IC. This function 
  is rarely called directly. It is used to compute optimally robust ICs.
}
\usage{
getInfStand(L2deriv, neighbor, ...)

\S4method{getInfStand}{UnivariateDistribution,ContNeighborhood}(L2deriv, neighbor, clip, cent, trafo)

\S4method{getInfStand}{UnivariateDistribution,TotalVarNeighborhood}(L2deriv, neighbor, clip, cent, trafo)

\S4method{getInfStand}{RealRandVariable,ContNeighborhood}(L2deriv, neighbor, Distr, A.comp, stand, clip, cent, trafo)
}
\arguments{
  \item{L2deriv}{ L2-derivative of some L2-differentiable family 
    of probability measures. }
  \item{neighbor}{ object of class \code{"Neighborhood"} }
  \item{\dots}{ additional parameters }
  \item{clip}{ optimal clipping bound. }
  \item{cent}{ optimal centering constant. }
  \item{stand}{ standardizing matrix. }
  \item{Distr}{ object of class \code{"Distribution"}. }
  \item{trafo}{ matrix: transformation of the parameter. }  
  \item{A.comp}{ matrix: indication which components of the standardizing
    matrix have to be computed. }
}
%\details{}
\value{The standardizing matrix is computed.}
\section{Methods}{
\describe{
  \item{L2deriv = "UnivariateDistribution", neighbor = "ContNeighborhood"}{ 
    computes standardizing matrix. }

  \item{L2deriv = "UnivariateDistribution", neighbor = "TotalVarNeighborhood"}{ 
    computes standardizing matrix. }

  \item{L2deriv = "RealRandVariable", neighbor = "ContNeighborhood"}{ 
    computes standardizing matrix. }
}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{ContIC-class}}, \code{\link{TotalVarIC-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
