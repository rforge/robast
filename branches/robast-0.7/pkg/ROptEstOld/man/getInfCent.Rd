\name{getInfCent}
\alias{getInfCent}
\alias{getInfCent-methods}
\alias{getInfCent,UnivariateDistribution,ContNeighborhood-method}
\alias{getInfCent,UnivariateDistribution,TotalVarNeighborhood-method}
\alias{getInfCent,RealRandVariable,ContNeighborhood-method}

\title{Generic Function for the Computation of the Optimal Centering Constant/Lower Clipping Bound}
\description{
  Generic function for the computation of the optimal centering constant
  (contamination neighborhoods) respectively, of the optimal lower clipping
  bound (total variation neighborhood).
  This function is rarely called directly. It is used to 
  compute optimally robust ICs.
}
\usage{
getInfCent(L2deriv, neighbor, ...)

\S4method{getInfCent}{UnivariateDistribution,ContNeighborhood}(L2deriv, neighbor, clip, cent, tol.z, symm, trafo)

\S4method{getInfCent}{UnivariateDistribution,TotalVarNeighborhood}(L2deriv, neighbor, clip, cent, tol.z, symm, trafo)

\S4method{getInfCent}{RealRandVariable,ContNeighborhood}(L2deriv, neighbor, Distr, z.comp, stand, cent, clip)
}
\arguments{
  \item{L2deriv}{ L2-derivative of some L2-differentiable family 
    of probability measures. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters. }
  \item{Distr}{ distribution of L2-differentiable family. }
  \item{clip}{ optimal clipping bound. }
  \item{cent}{ optimal centering constant. }
  \item{stand}{ standardizing matrix. }
  \item{tol.z}{ the desired accuracy (convergence tolerance). }
  \item{symm}{ logical: indicating symmetry of \code{L2deriv}. }
  \item{trafo}{ matrix: transformation of the parameter. }
  \item{z.comp}{ logical vector: indication which components of the 
    centering constant have to be computed. }
}
%\details{}
\value{The optimal centering constant is computed.}
\section{Methods}{
\describe{
  \item{L2deriv = "UnivariateDistribution", neighbor = "ContNeighborhood"}{ 
    computation of optimal centering constant. }

  \item{L2deriv = "UnivariateDistribution", neighbor = "TotalVarNeighborhood"}{ 
    computation of optimal lower clipping bound. }

  \item{L2deriv = "RealRandVariable", neighbor = "ContNeighborhood"}{ 
    computation of optimal centering constant. }
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
