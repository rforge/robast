\name{getModifyIC}
\alias{getModifyIC}
\alias{getModifyIC-methods}
\alias{getModifyIC,L2ParamFamily,Neighborhood,asRisk-method}
\alias{getModifyIC,L2LocationFamily,UncondNeighborhood,asGRisk-method}
\alias{getModifyIC,L2LocationFamily,UncondNeighborhood,fiUnOvShoot-method}
\alias{getModifyIC,L2ScaleFamily,ContNeighborhood,asGRisk-method}
\alias{getModifyIC,L2ScaleFamily,TotalVarNeighborhood,asGRisk-method}
\alias{getModifyIC,L2LocationScaleFamily,ContNeighborhood,asGRisk-method}
\title{Generic Function for the Computation of Functions for Slot modifyIC}

\description{
  This function is used by internal computations and is rarely called directly.
}
\usage{
getModifyIC(L2FamIC, neighbor, risk)
\S4method{getModifyIC}{L2ParamFamily,Neighborhood,asRisk}(L2FamIC, neighbor, risk)
\S4method{getModifyIC}{L2LocationFamily,UncondNeighborhood,asGRisk}(L2FamIC, neighbor, risk)
\S4method{getModifyIC}{L2LocationFamily,UncondNeighborhood,fiUnOvShoot}(L2FamIC, neighbor, risk)
\S4method{getModifyIC}{L2ScaleFamily,ContNeighborhood,asGRisk}(L2FamIC, neighbor, risk)
\S4method{getModifyIC}{L2ScaleFamily,TotalVarNeighborhood,asGRisk}(L2FamIC, neighbor, risk)
\S4method{getModifyIC}{L2LocationScaleFamily,ContNeighborhood,asGRisk}(L2FamIC, neighbor, risk)
}
\arguments{
  \item{L2FamIC}{ object of class \code{L2ParamFamily}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{risk}{ object of class \code{"RiskType"} }
}
\details{ This function is used for internal computations. }
\value{ Function for slot \code{modifyIC} of \code{IC}s }
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{optIC}}, \code{\link[RobAStBase]{IC-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
