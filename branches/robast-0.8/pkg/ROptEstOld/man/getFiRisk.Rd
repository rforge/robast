\name{getFiRisk}
\alias{getFiRisk}
\alias{getFiRisk-methods}
\alias{getFiRisk,fiUnOvShoot,Norm,ContNeighborhood-method}
\alias{getFiRisk,fiUnOvShoot,Norm,TotalVarNeighborhood-method}

\title{Generic Function for Computation of Finite-Sample Risks}
\description{
  Generic function for the computation of finite-sample risks.
  This function is rarely called directly. It is used by 
  other functions.
}
\usage{
getFiRisk(risk, Distr, neighbor, ...)

\S4method{getFiRisk}{fiUnOvShoot,Norm,ContNeighborhood}(risk, Distr, neighbor, 
            clip, stand, sampleSize, Algo, cont)

\S4method{getFiRisk}{fiUnOvShoot,Norm,TotalVarNeighborhood}(risk, Distr, neighbor, 
            clip, stand, sampleSize, Algo, cont)
}
\arguments{
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{Distr}{ object of class \code{"Distribution"}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters. }
  \item{clip}{ positive real: clipping bound }
  \item{stand}{ standardizing constant/matrix. }
  \item{sampleSize}{ integer: sample size. }
  \item{Algo}{ "A" or "B". }
  \item{cont}{ "left" or "right". }
}
\details{The computation of the finite-sample under-/overshoot risk
  is based on FFT. For more details we refer to Section 11.3 of Kohl (2005).
}
\value{The finite-sample risk is computed.}
\section{Methods}{
\describe{
  \item{risk = "fiUnOvShoot", Distr = "Norm", neighbor = "ContNeighborhood"}{ 
    computes finite-sample under-/overshoot risk in methods for 
    function \code{getFixRobIC}. }

  \item{risk = "fiUnOvShoot", Distr = "Norm", neighbor = "TotalVarNeighborhood"}{ 
    computes finite-sample under-/overshoot risk in methods for 
    function \code{getFixRobIC}. }
}}
\references{
  Huber, P.J. (1968) Robust Confidence Limits. Z. Wahrscheinlichkeitstheor.
  Verw. Geb. \bold{10}:269--278.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Ruckdeschel, P. and Kohl, M. (2005) Computation of the Finite Sample Risk 
  of M-estimators on Neighborhoods.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{fiRisk-class}}}
%\examples{}
\concept{finite-sample risk}
\concept{risk}
\keyword{robust}
