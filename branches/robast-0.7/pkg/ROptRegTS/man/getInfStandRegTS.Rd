\name{getInfStandRegTS}
\alias{getInfStandRegTS}
\alias{getInfStandRegTS-methods}
\alias{getInfStandRegTS,UnivariateDistribution,UnivariateDistribution,ContNeighborhood-method}
\alias{getInfStandRegTS,UnivariateDistribution,UnivariateDistribution,TotalVarNeighborhood-method}
\alias{getInfStandRegTS,UnivariateDistribution,UnivariateDistribution,CondTotalVarNeighborhood-method}
\alias{getInfStandRegTS,UnivariateDistribution,UnivariateDistribution,Av1CondContNeighborhood-method}
\alias{getInfStandRegTS,UnivariateDistribution,UnivariateDistribution,Av1CondTotalVarNeighborhood-method}
\alias{getInfStandRegTS,UnivariateDistribution,MultivariateDistribution,ContNeighborhood-method}
\alias{getInfStandRegTS,UnivariateDistribution,MultivariateDistribution,Av1CondContNeighborhood-method}
\alias{getInfStandRegTS,UnivariateDistribution,MultivariateDistribution,Av1CondTotalVarNeighborhood-method}
\alias{getInfStandRegTS,UnivariateDistribution,Distribution,Av2CondContNeighborhood-method}
\alias{getInfStandRegTS,RealRandVariable,Distribution,ContNeighborhood-method}
\alias{getInfStandRegTS,RealRandVariable,Distribution,Av1CondContNeighborhood-method}

\title{ Generic Function for the Computation of the Standardizing Matrix }
\description{
  Generic function for the computation of the standardizing matrix which
  takes care of the Fisher consistency of the corresponding IC. This function 
  is rarely called directly. It is used to compute optimally robust ICs.
}
\usage{
getInfStandRegTS(ErrorL2deriv, Regressor, neighbor, ...)

\S4method{getInfStandRegTS}{UnivariateDistribution,UnivariateDistribution,ContNeighborhood}(ErrorL2deriv, 
		Regressor, neighbor, z.comp, clip, cent, stand, trafo)

\S4method{getInfStandRegTS}{UnivariateDistribution,UnivariateDistribution,TotalVarNeighborhood}(ErrorL2deriv, 
		Regressor, neighbor, clip, cent)

\S4method{getInfStandRegTS}{UnivariateDistribution,UnivariateDistribution,CondTotalVarNeighborhood}(ErrorL2deriv, 
		Regressor, neighbor, clip, cent)

\S4method{getInfStandRegTS}{UnivariateDistribution,UnivariateDistribution,Av1CondContNeighborhood}(ErrorL2deriv, 
		Regressor, neighbor, z.comp, clip, cent, stand, trafo)

\S4method{getInfStandRegTS}{UnivariateDistribution,UnivariateDistribution,Av1CondTotalVarNeighborhood}(ErrorL2deriv, 
		Regressor, neighbor, z.comp, clip, cent, stand, trafo)

\S4method{getInfStandRegTS}{UnivariateDistribution,MultivariateDistribution,ContNeighborhood}(ErrorL2deriv, 
		Regressor, neighbor, z.comp, clip, cent, stand, trafo)

\S4method{getInfStandRegTS}{UnivariateDistribution,MultivariateDistribution,Av1CondContNeighborhood}(ErrorL2deriv, 
		Regressor, neighbor, z.comp, clip, cent, stand, trafo)

\S4method{getInfStandRegTS}{UnivariateDistribution,MultivariateDistribution,Av1CondTotalVarNeighborhood}(ErrorL2deriv, 
		Regressor, neighbor, z.comp, clip, cent, stand, trafo)

\S4method{getInfStandRegTS}{UnivariateDistribution,Distribution,Av2CondContNeighborhood}(ErrorL2deriv, 
		Regressor, neighbor, z.comp, clip, cent, stand, trafo)

\S4method{getInfStandRegTS}{RealRandVariable,Distribution,ContNeighborhood}(ErrorL2deriv, 
		Regressor, neighbor, ErrorDistr, A.comp, stand, clip, cent, trafo)

\S4method{getInfStandRegTS}{RealRandVariable,Distribution,Av1CondContNeighborhood}(ErrorL2deriv, 
		Regressor, neighbor, ErrorDistr, A.comp, stand, clip, cent, trafo)
}
\arguments{
  \item{ErrorL2deriv}{ L2-derivative of \code{ErrorDistr}. }
  \item{Regressor}{ regressor. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters. }
  \item{ErrorDistr}{ error distribution. }
  \item{clip}{ optimal clipping bound/function. }
  \item{cent}{ optimal centering constant/function. }
  \item{stand}{ standardizing matrix. }
  \item{z.comp}{ which components of the centering constant/function 
    have to be computed. }
  \item{A.comp}{ which components of the standardizing matrix
    have to be computed. }
  \item{trafo}{ matrix: transformation of the parameter. }
}
%\details{}
\value{The standardizing matrix is computed.}
\section{Methods}{
\describe{
  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", neighbor = "ContNeighborhood"}{ 
    computes standardizing matrix. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", neighbor = "TotalVarNeighborhood"}{ 
    computes standardizing constant. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", neighbor = "CondTotalVarNeighborhood"}{ 
    computes standardizing constant. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", neighbor = "Av1CondContNeighborhood"}{ 
    computes standardizing matrix. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", neighbor = "Av1CondTotalVarNeighborhood"}{ 
    computes standardizing matrix. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "MultivariateDistribution", neighbor = "ContNeighborhood"}{ 
    computes standardizing matrix. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "MultivariateDistribution", neighbor = "Av1CondContNeighborhood"}{ 
    computes standardizing matrix. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "MultivariateDistribution", neighbor = "Av1CondTotalVarNeighborhood"}{ 
    computes standardizing matrix. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", neighbor = "Av2CondContNeighborhood"}{ 
    computes standardizing matrix. }

  \item{ErrorL2deriv = "RealRandVariable", Regressor = "Distribution", neighbor = "ContNeighborhood"}{ 
    computes standardizing matrix. }

  \item{ErrorL2deriv = "RealRandVariable", Regressor = "Distribution", neighbor = "Av1CondContNeighborhood"}{ 
    computes standardizing matrix. }
}}
\references{
  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[ROptEstOld]{ContIC-class}}, \code{\link[ROptEstOld]{TotalVarIC-class}}, 
        \code{\link{Av1CondContIC-class}}, \code{\link{Av2CondContIC-class}}, 
        \code{\link{Av1CondTotalVarIC-class}}, \code{\link{CondContIC}},
        \code{\link{CondTotalVarIC}}}
%\examples{}
\concept{influence curve}
\keyword{robust}

