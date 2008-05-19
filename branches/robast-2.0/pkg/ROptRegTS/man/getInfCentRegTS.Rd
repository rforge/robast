\name{getInfCentRegTS}
\alias{getInfCentRegTS}
\alias{getInfCentRegTS-methods}
\alias{getInfCentRegTS,UnivariateDistribution,UnivariateDistribution,ContNeighborhood-method}
\alias{getInfCentRegTS,UnivariateDistribution,UnivariateDistribution,TotalVarNeighborhood-method}
\alias{getInfCentRegTS,UnivariateDistribution,numeric,CondTotalVarNeighborhood-method}
\alias{getInfCentRegTS,UnivariateDistribution,UnivariateDistribution,Av1CondContNeighborhood-method}
\alias{getInfCentRegTS,UnivariateDistribution,UnivariateDistribution,Av1CondTotalVarNeighborhood-method}
\alias{getInfCentRegTS,UnivariateDistribution,MultivariateDistribution,ContNeighborhood-method}
\alias{getInfCentRegTS,UnivariateDistribution,MultivariateDistribution,Av1CondContNeighborhood-method}
\alias{getInfCentRegTS,UnivariateDistribution,MultivariateDistribution,Av1CondTotalVarNeighborhood-method}
\alias{getInfCentRegTS,UnivariateDistribution,Distribution,Av2CondContNeighborhood-method}
\alias{getInfCentRegTS,RealRandVariable,Distribution,ContNeighborhood-method}
\alias{getInfCentRegTS,RealRandVariable,Distribution,Av1CondContNeighborhood-method}

\title{Generic Function for the Computation of the Optimal Centering Constant/Function resp. 
  Lower Clipping Bound/Function}
\description{
  Generic function for the computation of the optimal centering constant/function
  (contamination neighborhoods) respectively, of the optimal lower clipping
  bound/function (total variation neighborhoods).
  This function is rarely called directly. It is used to 
  compute optimally robust ICs.
}
\usage{
getInfCentRegTS(ErrorL2deriv, Regressor, neighbor, ...)

\S4method{getInfCentRegTS}{UnivariateDistribution,UnivariateDistribution,ContNeighborhood}(ErrorL2deriv, 
                Regressor, neighbor, clip, cent, stand, z.comp)

\S4method{getInfCentRegTS}{UnivariateDistribution,UnivariateDistribution,TotalVarNeighborhood}(ErrorL2deriv, 
                Regressor, neighbor, clip, cent, z.comp)

\S4method{getInfCentRegTS}{UnivariateDistribution,numeric,CondTotalVarNeighborhood}(ErrorL2deriv, 
                Regressor, neighbor, clip, cent, z.comp)

\S4method{getInfCentRegTS}{UnivariateDistribution,UnivariateDistribution,Av1CondContNeighborhood}(ErrorL2deriv, 
                Regressor, neighbor, clip, cent, stand, z.comp, x.vec)

\S4method{getInfCentRegTS}{UnivariateDistribution,UnivariateDistribution,Av1CondTotalVarNeighborhood}(ErrorL2deriv, 
                Regressor, neighbor, clip, cent, stand, z.comp, x.vec, tol.z)

\S4method{getInfCentRegTS}{UnivariateDistribution,MultivariateDistribution,ContNeighborhood}(ErrorL2deriv, 
                Regressor, neighbor, clip, cent, stand, z.comp)

\S4method{getInfCentRegTS}{UnivariateDistribution,MultivariateDistribution,Av1CondContNeighborhood}(ErrorL2deriv, 
                Regressor, neighbor, clip, cent, stand, z.comp, x.vec)

\S4method{getInfCentRegTS}{UnivariateDistribution,Distribution,Av2CondContNeighborhood}(ErrorL2deriv, 
                Regressor, neighbor, clip, cent, stand, z.comp, tol.z)

\S4method{getInfCentRegTS}{RealRandVariable,Distribution,ContNeighborhood}(ErrorL2deriv, 
                Regressor, neighbor, ErrorDistr, stand, cent, clip, z.comp)

\S4method{getInfCentRegTS}{RealRandVariable,Distribution,Av1CondContNeighborhood}(ErrorL2deriv, 
                Regressor, neighbor, ErrorDistr, stand, cent, clip, z.comp, x.vec)
}
\arguments{
  \item{ErrorL2deriv}{ L2-derivative of \code{ErrorDistr}. }
  \item{Regressor}{ regressor. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters. }
  \item{clip}{ optimal clipping bound. }
  \item{cent}{ optimal centering constant/function. }
  \item{stand}{ standardizing matrix. }
  \item{z.comp}{ which components of the centering constant/function 
    have to be computed. }
  \item{x.vec}{ (approximated) support of \code{Regressor}. }
  \item{tol.z}{ the desired accuracy (convergence tolerance). }
  \item{ErrorDistr}{ error distribution. }
}
%\details{}
\value{The optimal centering constant/function is computed.}
\section{Methods}{
\describe{
  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", neighbor = "ContNeighborhood"}{ 
    computation of optimal centering constant. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", neighbor = "TotalVarNeighborhood"}{ 
    computation of lower clipping bound. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "numeric", neighbor = "CondTotalVarNeighborhood"}{ 
    computation of lower clipping bound. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", neighbor = "Av1CondContNeighborhood"}{ 
    computation of optimal centering function. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", neighbor = "Av1CondTotalVarNeighborhood"}{ 
    computation of optimal lower clipping function. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "MultivariateDistribution", neighbor = "ContNeighborhood"}{ 
    computation of optimal centering constant. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "MultivariateDistribution", neighbor = "Av1CondContNeighborhood"}{ 
    computation of optimal centering function. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "MultivariateDistribution", neighbor = "Av1CondTotalVarNeighborhood"}{ 
    computation of optimal lower clipping function. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", neighbor = "Av2CondContNeighborhood"}{ 
    computation of optimal centering constant. }

  \item{ErrorL2deriv = "RealRandVariable", Regressor = "Distribution", neighbor = "ContNeighborhood"}{ 
    computation of optimal centering constant. }

  \item{ErrorL2deriv = "RealRandVariable", Regressor = "Distribution", neighbor = "Av1CondContNeighborhood"}{ 
    computation of optimal centering function. }
}}
\references{
  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[ROptEst]{ContIC-class}}, \code{\link{Av1CondContIC-class}}, 
        \code{\link{Av2CondContIC-class}}, \code{\link{Av1CondTotalVarIC-class}}, 
        \code{\link{CondContIC-class}}, \code{\link{CondTotalVarIC-class}}}
%\examples{}
\concept{influence curve}
\keyword{}
