\name{getRiskIC}
\alias{getRiskIC}
\alias{getRiskIC-methods}
\alias{getRiskIC,IC,asCov,missing,missing-method}
\alias{getRiskIC,IC,asCov,missing,L2ParamFamily-method}
\alias{getRiskIC,IC,trAsCov,missing,missing-method}
\alias{getRiskIC,IC,trAsCov,missing,L2ParamFamily-method}
\alias{getRiskIC,IC,asBias,ContNeighborhood,missing-method}
\alias{getRiskIC,IC,asBias,ContNeighborhood,L2ParamFamily-method}
\alias{getRiskIC,IC,asBias,TotalVarNeighborhood,missing-method}
\alias{getRiskIC,IC,asBias,TotalVarNeighborhood,L2ParamFamily-method}
\alias{getRiskIC,IC,asMSE,UncondNeighborhood,missing-method}
\alias{getRiskIC,IC,asMSE,UncondNeighborhood,L2ParamFamily-method}
\alias{getRiskIC,TotalVarIC,asUnOvShoot,UncondNeighborhood,missing-method}
\alias{getRiskIC,IC,fiUnOvShoot,ContNeighborhood,missing-method}
\alias{getRiskIC,IC,fiUnOvShoot,TotalVarNeighborhood,missing-method}

\title{Generic function for the computation of a risk for an IC}
\description{
  Generic function for the computation of a risk for an IC.
}
\usage{
getRiskIC(IC, risk, neighbor, L2Fam, ...)

\S4method{getRiskIC}{IC,asCov,missing,missing}(IC, risk, tol = .Machine$double.eps^0.25)

\S4method{getRiskIC}{IC,asCov,missing,L2ParamFamily}(IC, risk, L2Fam, tol = .Machine$double.eps^0.25)

\S4method{getRiskIC}{IC,trAsCov,missing,missing}(IC, risk, tol = .Machine$double.eps^0.25)

\S4method{getRiskIC}{IC,trAsCov,missing,L2ParamFamily}(IC, risk, L2Fam, tol = .Machine$double.eps^0.25)

\S4method{getRiskIC}{IC,asBias,ContNeighborhood,missing}(IC, risk, neighbor, tol = .Machine$double.eps^0.25)

\S4method{getRiskIC}{IC,asBias,ContNeighborhood,L2ParamFamily}(IC, risk, neighbor, L2Fam, tol = .Machine$double.eps^0.25)

\S4method{getRiskIC}{IC,asBias,TotalVarNeighborhood,missing}(IC, risk, neighbor, tol = .Machine$double.eps^0.25)

\S4method{getRiskIC}{IC,asBias,TotalVarNeighborhood,L2ParamFamily}(IC, risk, neighbor, L2Fam, tol = .Machine$double.eps^0.25)

\S4method{getRiskIC}{IC,asMSE,UncondNeighborhood,missing}(IC, risk, neighbor, tol = .Machine$double.eps^0.25)

\S4method{getRiskIC}{IC,asMSE,UncondNeighborhood,L2ParamFamily}(IC, risk, neighbor, L2Fam, tol = .Machine$double.eps^0.25)

\S4method{getRiskIC}{TotalVarIC,asUnOvShoot,UncondNeighborhood,missing}(IC, risk, neighbor)

\S4method{getRiskIC}{IC,fiUnOvShoot,ContNeighborhood,missing}(IC, risk, neighbor, sampleSize, Algo = "A", cont = "left")

\S4method{getRiskIC}{IC,fiUnOvShoot,TotalVarNeighborhood,missing}(IC, risk, neighbor, sampleSize, Algo = "A", cont = "left")
}
\arguments{
  \item{IC}{ object of class \code{"InfluenceCurve"} }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{L2Fam}{ object of class \code{"L2ParamFamily"}. }
  \item{\dots}{ additional parameters }
  \item{tol}{ the desired accuracy (convergence tolerance).}  
  \item{sampleSize}{ integer: sample size. }
  \item{Algo}{ "A" or "B". }
  \item{cont}{ "left" or "right". }
}
\details{To make sure that the results are valid, it is recommended
  to include an additional check of the IC properties of \code{IC} 
  using \code{checkIC}.}
\value{The risk of an IC is computed.}
\section{Methods}{
\describe{
  \item{IC = "IC", risk = "asCov", neighbor = "missing", L2Fam = "missing"}{ 
    asymptotic covariance of \code{IC}. }

  \item{IC = "IC", risk = "asCov", neighbor = "missing", L2Fam = "L2ParamFamily"}{ 
    asymptotic covariance of \code{IC} under \code{L2Fam}. }

  \item{IC = "IC", risk = "trAsCov", neighbor = "missing", L2Fam = "missing"}{ 
    asymptotic covariance of \code{IC}. }

  \item{IC = "IC", risk = "trAsCov", neighbor = "missing", L2Fam = "L2ParamFamily"}{ 
    asymptotic covariance of \code{IC} under \code{L2Fam}. }

  \item{IC = "IC", risk = "asBias", neighbor = "ContNeighborhood", L2Fam = "missing"}{ 
    asymptotic bias of \code{IC} under convex contaminations. }

  \item{IC = "IC", risk = "asBias", neighbor = "ContNeighborhood", L2Fam = "L2ParamFamily"}{ 
    asymptotic bias of \code{IC} under convex contaminations and \code{L2Fam}. }

  \item{IC = "IC", risk = "asBias", neighbor = "TotalVarNeighborhood", L2Fam = "missing"}{ 
    asymptotic bias of \code{IC} in case of total variation neighborhoods. }

  \item{IC = "IC", risk = "asBias", neighbor = "TotalVarNeighborhood", L2Fam = "L2ParamFamily"}{ 
    asymptotic bias of \code{IC} under \code{L2Fam} in case of total variation 
    neighborhoods. }

  \item{IC = "IC", risk = "asMSE", neighbor = "UncondNeighborhood", L2Fam = "missing"}{ 
    asymptotic mean square error of \code{IC}. }

  \item{IC = "IC", risk = "asMSE", neighbor = "UncondNeighborhood", L2Fam = "L2ParamFamily"}{ 
    asymptotic mean square error of \code{IC} under \code{L2Fam}. }

  \item{IC = "TotalVarIC", risk = "asUnOvShoot", neighbor = "UncondNeighborhood", L2Fam = "missing"}{ 
    asymptotic under-/overshoot risk of \code{IC}. }

  \item{IC = "IC", risk = "fiUnOvShoot", neighbor = "ContNeighborhood", L2Fam = "missing"}{ 
    finite-sample under-/overshoot risk of \code{IC}. }

  \item{IC = "IC", risk = "fiUnOvShoot", neighbor = "TotalVarNeighborhood", L2Fam = "missing"}{ 
    finite-sample under-/overshoot risk of \code{IC}. }
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
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\note{This generic function is still under construction.}
\seealso{\code{\link{getRiskIC-methods}}, \code{\link{InfRobModel-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
