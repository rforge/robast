\name{getBiasIC}
\docType{methods}
\alias{getBiasIC}
\alias{getBiasIC-methods}
\alias{getBiasIC,HampIC,UncondNeighborhood-method}
\alias{getBiasIC,TotalVarIC,UncondNeighborhood-method}

\title{Generic function for the computation of the asymptotic bias for an IC}
\description{
  Generic function for the computation of the asymptotic bias for an IC.
}
\usage{
getBiasIC(IC, neighbor, ...)

\S4method{getBiasIC}{HampIC,UncondNeighborhood}(IC, neighbor, L2Fam, ...)
}
\arguments{
  \item{IC}{ object of class \code{"InfluenceCurve"} }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{L2Fam}{ object of class \code{"L2ParamFamily"}. }
  \item{\dots}{ additional parameters }
}
\details{ This function is rarely called directly. It is used by 
  other functions/methods. }
\value{The bias of the IC is computed.}
\section{Methods}{
\describe{
  \item{IC = "HampIC", neighbor = "UncondNeighborhood"}{
    reads off the as. bias from the risks-slot of the IC. }
  \item{IC = "TotalVarIC", neighbor = "UncondNeighborhood"}{
    reads off the as. bias from the risks-slot of the IC,
    resp. if this is \code{NULL} from the corresponding
    Lagrange Multipliers. }
}}
\references{
  Huber, P.J. (1968) Robust Confidence Limits. Z. Wahrscheinlichkeitstheor.
  Verw. Geb. \bold{10}:269--278.

  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}.
  Bayreuth: Dissertation.

  Ruckdeschel, P. and Kohl, M. (2005) Computation of the Finite Sample Bias
  of M-estimators on Neighborhoods.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
\note{This generic function is still under construction.}
\seealso{\code{\link{getRiskIC-methods}}, \code{\link[RobAStBase]{InfRobModel-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
