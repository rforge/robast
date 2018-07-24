\name{getRiskIC}
\docType{methods}
\alias{getRiskIC}
\alias{getRiskIC-methods}
\alias{getRiskIC,HampIC,asCov,missing,missing-method}
\alias{getRiskIC,HampIC,asCov,missing,L2ParamFamily-method}
\alias{getRiskIC,TotalVarIC,asCov,missing,L2ParamFamily-method}

\title{Generic function for the computation of a risk for an IC}
\description{
  Generic function for the computation of a risk for an IC.
}
\usage{
getRiskIC(IC, risk, neighbor, L2Fam, ...)

\S4method{getRiskIC}{HampIC,asCov,missing,missing}(IC, risk, withCheck= TRUE)

\S4method{getRiskIC}{HampIC,asCov,missing,L2ParamFamily}(IC, risk, L2Fam, withCheck= TRUE)
\S4method{getRiskIC}{TotalVarIC,asCov,missing,L2ParamFamily}(IC, risk, L2Fam, withCheck = TRUE)

}
\arguments{
  \item{IC}{ object of class \code{"InfluenceCurve"} }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}; missing in the methods described here. }
  \item{\dots}{ additional parameters }
  \item{L2Fam}{ object of class \code{"L2ParamFamily"}. }
  \item{withCheck}{logical: should a call to \code{checkIC} be done to
                   check accuracy (defaults to \code{TRUE}; ignored
                   if nothing is computed but simply a slot is read out).}
}
\details{To make sure that the results are valid, it is recommended
  to include an additional check of the IC properties of \code{IC}
  using \code{checkIC}.}
\value{The risk of an IC is computed.}
\section{Methods}{
\describe{
  \item{IC = "HampIC", risk = "asCov", neighbor = "missing", L2Fam = "missing"}{
    asymptotic covariance of \code{IC} read off from corresp. \code{Risks} slot. }

  \item{IC = "HampIC", risk = "asCov", neighbor = "missing", L2Fam = "L2ParamFamily"}{
    asymptotic covariance of \code{IC} under \code{L2Fam} read off from corresp. \code{Risks} slot. }
  \item{IC = "TotalVarIC", risk = "asCov", neighbor = "missing", L2Fam = "L2ParamFamily"}{
   asymptotic covariance of \code{IC} read off from corresp. \code{Risks} slot,
   resp. if this is \code{NULL} calculates it via \code{\link{getInfV}}.}
}}
\references{
  Huber, P.J. (1968) Robust Confidence Limits. Z. Wahrscheinlichkeitstheor.
  Verw. Geb. \bold{10}:269--278.

  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}.
  Bayreuth: Dissertation.

  Ruckdeschel, P. and Kohl, M. (2005) Computation of the Finite Sample Risk
  of M-estimators on Neighborhoods.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
\note{This generic function is still under construction.}
\seealso{\code{\link[ROptEst]{getRiskIC}}, \code{\link[RobAStBase]{InfRobModel-class}}}
\examples{
B <- BinomFamily(size = 25, prob = 0.25)

## classical optimal IC
IC0 <- optIC(model = B, risk = asCov())
getRiskIC(IC0, asCov())
}
\concept{influence curve}
\keyword{robust}
