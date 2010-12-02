\name{getweight-methods}
\docType{methods}
\alias{getweight-methods}
\alias{getweight}
\alias{getweight,HampelWeight,ContNeighborhood,BiasType-method}
\alias{getweight,HampelWeight,ContNeighborhood,onesidedBias-method}
\alias{getweight,HampelWeight,ContNeighborhood,asymmetricBias-method}
\alias{getweight,BdStWeight,TotalVarNeighborhood,BiasType-method}
\alias{minbiasweight-methods}
\alias{minbiasweight}
\alias{minbiasweight,HampelWeight,ContNeighborhood,BiasType-method}
\alias{minbiasweight,HampelWeight,ContNeighborhood,onesidedBias-method}
\alias{minbiasweight,HampelWeight,ContNeighborhood,asymmetricBias-method}
\alias{minbiasweight,BdStWeight,TotalVarNeighborhood,BiasType-method}

\title{Generating weights}
\description{
  Generates weight functions of Hampel / BdSt type for different bias and norm types.
}
\usage{
getweight(Weight, neighbor, biastype, ...)
minbiasweight(Weight, neighbor, biastype, ...)
\S4method{getweight}{HampelWeight,ContNeighborhood,BiasType}(Weight, neighbor, biastype, normW)
\S4method{minbiasweight}{HampelWeight,ContNeighborhood,BiasType}(Weight, neighbor, biastype, normW)
\S4method{getweight}{HampelWeight,ContNeighborhood,onesidedBias}(Weight, neighbor, biastype, ...)
\S4method{minbiasweight}{HampelWeight,ContNeighborhood,onesidedBias}(Weight, neighbor, biastype, ...)
\S4method{getweight}{HampelWeight,ContNeighborhood,asymmetricBias}(Weight, neighbor, biastype, ...)
\S4method{minbiasweight}{HampelWeight,ContNeighborhood,asymmetricBias}(Weight, neighbor, biastype, ...)
\S4method{getweight}{BdStWeight,TotalVarNeighborhood,BiasType}(Weight, neighbor, biastype, ...)
\S4method{minbiasweight}{BdStWeight,TotalVarNeighborhood,BiasType}(Weight, neighbor, biastype, ...)
}
\arguments{
  \item{Weight}{ Object of class \code{"RobWeight"}. }
  \item{neighbor}{ Object of class \code{"Neighborhood"}. }
  \item{biastype}{ Object of class \code{"BiasType"}. }
  \item{normW}{ Object of class \code{"NormType"} --- only for signature \code{HampelWeight,ContNeighborhood,BiasType}. }
  \item{\dots}{possibly additional (unused) arguments --- like in a call to the less specific methods.} 
}
%\details{}
\value{Object of class \code{"HampelWeight"} resp. \code{"BdStWeight"}}
\references{
  Hampel et al. (1986) \emph{Robust Statistics}. 
  The Approach Based on Influence Functions. New York: Wiley.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\details{These functions generate the weight function in slot \code{weight} in a corresp.
         object of class \code{RobWeight} and descendants.}
\section{Methods}{\describe{
\item{getweight}{\code{signature(Weight = "HampelWeight", neighbor = "ContNeighborhood", 
                       biastype = "BiasType")} with additional argument \code{biastype} of class 
                       \code{"BiasType"}: produces weight slot...}
\item{minbiasweight}{\code{signature(Weight = "HampelWeight", neighbor = "ContNeighborhood", 
                       biastype = "BiasType")} with additional argument \code{biastype} of class 
                       \code{"BiasType"}: produces weight slot...}
\item{getweight}{\code{signature(Weight = "HampelWeight", neighbor = "ContNeighborhood", 
                       biastype = "onesidedBias")}: produces weight slot...}
\item{minbiasweight}{\code{signature(Weight = "HampelWeight", neighbor = "ContNeighborhood", 
                       biastype = "onesidedBias")}: produces weight slot...}
\item{getweight}{\code{signature(Weight = "HampelWeight", neighbor = "ContNeighborhood", 
                       biastype = "asymmetricBias")}: produces weight slot...}
\item{minbiasweight}{\code{signature(Weight = "HampelWeight", neighbor = "ContNeighborhood", 
                       biastype = "asymmetricBias")}: produces weight slot...}
\item{getweight}{\code{signature(Weight = "BdStWeight", neighbor = "TotalVarNeighborhood", 
                       biastype = "BiasType")}: produces weight slot...}
\item{minbiasweight}{\code{signature(Weight = "BdStWeight", neighbor = "TotalVarNeighborhood", 
                       biastype = "BiasType")}: produces weight slot...}
}}

\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{BdStWeight-class}},
         \code{\link{HampelWeight-class}},
         \code{\link{IC-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
