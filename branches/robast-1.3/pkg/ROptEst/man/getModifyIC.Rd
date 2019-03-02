\name{getModifyIC}
\alias{getModifyIC}
\alias{scaleUpdateIC}
\alias{getModifyIC-methods}
\alias{scaleUpdateIC-methods}
\alias{scaleUpdateIC,UncondNeighborhood-method}
\alias{scaleUpdateIC,ContNeighborhood-method}
\alias{scaleUpdateIC,TotalVarNeighborhood-method}
\alias{getModifyIC,L2ParamFamily,Neighborhood,asRisk-method}
\alias{getModifyIC,L2LocationFamily,UncondNeighborhood,asGRisk-method}
\alias{getModifyIC,L2LocationFamily,UncondNeighborhood,fiUnOvShoot-method}
\alias{getModifyIC,L2ScaleFamily,UncondNeighborhood,asGRisk-method}
\alias{getModifyIC,L2LocationScaleFamily,UncondNeighborhood,asGRisk-method}
\title{Generic Function for the Computation of Functions for Slot modifyIC}

\description{
  These function is used by internal computations and is rarely called directly.
}
\usage{
getModifyIC(L2FamIC, neighbor, risk,...)
\S4method{getModifyIC}{L2ParamFamily,Neighborhood,asRisk}(L2FamIC,
          neighbor, risk, ...)
\S4method{getModifyIC}{L2LocationFamily,UncondNeighborhood,asGRisk}(L2FamIC,
          neighbor, risk, ...)
\S4method{getModifyIC}{L2LocationFamily,UncondNeighborhood,fiUnOvShoot}(L2FamIC,
          neighbor, risk, ...)
\S4method{getModifyIC}{L2ScaleFamily,UncondNeighborhood,asGRisk}(L2FamIC,
          neighbor, risk, ..., modifyICwarn = NULL)
\S4method{getModifyIC}{L2LocationScaleFamily,UncondNeighborhood,asGRisk}(L2FamIC,
          neighbor, risk, ..., modifyICwarn = NULL)

scaleUpdateIC(neighbor,...)
\S4method{scaleUpdateIC}{UncondNeighborhood}(neighbor, sdneu, sdalt, IC)
\S4method{scaleUpdateIC}{ContNeighborhood}(neighbor, sdneu, sdalt, IC)
\S4method{scaleUpdateIC}{TotalVarNeighborhood}(neighbor, sdneu, sdalt, IC)
}
\arguments{
  \item{L2FamIC}{ object of class \code{L2ParamFamily}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{risk}{ object of class \code{"RiskType"} }
  \item{\dots}{further arguments to be passed over to \code{optIC}.}
  \item{sdneu}{positive numeric of length one; the new scale.}
  \item{sdalt}{positive numeric of length one; the new scale.}
  \item{IC}{a Hampel-IC to be updated.}
  \item{modifyICwarn}{logical: should a (warning) information be added if
     \code{modifyIC} is applied and hence some optimality information could
     no longer be valid? Defaults to \code{NULL} in which case this value
     is taken from \code{RobAStBaseOptions}.}
}
\details{ This function is used for internal computations.
 By setting \code{RobAStBaseOption("all.verbose" = TRUE)} somewhere
 globally, the generated function \code{modifyIC} will generate
 calls to \code{optIC} with argument \code{verbose=TRUE}.
}
\value{\describe{
\item{getmodifyIC}{Function for slot \code{modifyIC} of \code{IC}s}
\item{scaleUpdateIC}{a list to be digested in corresponding methods
of \code{getmodifyIC} by \code{generateIC}}
 }}
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
