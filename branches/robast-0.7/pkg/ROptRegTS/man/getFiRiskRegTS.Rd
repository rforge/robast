\name{getFiRiskRegTS}
\alias{getFiRiskRegTS}
\alias{getFiRiskRegTS-methods}
\alias{getFiRiskRegTS,fiUnOvShoot,Norm,UnivariateDistribution,ContNeighborhood-method}
\alias{getFiRiskRegTS,fiUnOvShoot,Norm,UnivariateDistribution,TotalVarNeighborhood-method}
\alias{getFiRiskRegTS,fiUnOvShoot,Norm,UnivariateDistribution,CondContNeighborhood-method}
\alias{getFiRiskRegTS,fiUnOvShoot,Norm,UnivariateDistribution,CondTotalVarNeighborhood-method}

\title{Generic Function for Computation of Finite-Sample Risks}
\description{
  Generic function for the computation of finite-sample risks 
  in regression-type models. This function is rarely called directly. 
  It is used by other functions.
}
\usage{
getFiRiskRegTS(risk, ErrorDistr, Regressor, neighbor, ...)

\S4method{getFiRiskRegTS}{fiUnOvShoot,Norm,UnivariateDistribution,ContNeighborhood}(risk, 
            ErrorDistr, Regressor, neighbor, clip, stand, sampleSize, Algo, cont)

\S4method{getFiRiskRegTS}{fiUnOvShoot,Norm,UnivariateDistribution,TotalVarNeighborhood}(risk, 
            ErrorDistr, Regressor, neighbor, clip, stand, sampleSize, Algo, cont)

\S4method{getFiRiskRegTS}{fiUnOvShoot,Norm,UnivariateDistribution,CondContNeighborhood}(risk, 
            ErrorDistr, Regressor, neighbor, clip, stand, sampleSize, cont)

\S4method{getFiRiskRegTS}{fiUnOvShoot,Norm,UnivariateDistribution,CondTotalVarNeighborhood}(risk, 
            ErrorDistr, Regressor, neighbor, clip, stand, sampleSize, cont)
}
\arguments{
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{ErrorDistr}{ error distribution }
  \item{Regressor}{ regressor }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters. }
  \item{clip}{optimal clipping bound/function.}
  \item{stand}{standardizing matrix.}
  \item{sampleSize}{ integer: sample size. }  
  \item{Algo}{ "A" or "B". }
  \item{cont}{ "left" or "right". }  
}
\details{
  The computation of the finite-sample under-/overshoot risk is
  based on FFT. For more details we refer to Subsections 12.1.3 
  and 12.2.3 of Kohl (2005).
}
\value{The finite-sample risk is computed.}
\section{Methods}{
\describe{
\item{risk = "fiUnOvShoot", ErrorDistr = "Norm", Regressor = "UnivariateDistribution", neighbor = "ContNeighborhood"}{ 
  computes finite-sample under-/overshoot risk in methods for 
  function 'getFixRobRegTypeIC'. }

\item{risk = "fiUnOvShoot", ErrorDistr = "Norm", Regressor = "UnivariateDistribution", neighbor = "TotalVarNeighborhood"}{ 
  computes finite-sample under-/overshoot risk in methods for 
  function 'getFixRobRegTypeIC'. }

\item{risk = "fiUnOvShoot", ErrorDistr = "Norm", Regressor = "UnivariateDistribution", neighbor = "CondContNeighborhood"}{ 
  computes finite-sample under-/overshoot risk in methods for 
  function 'getFixRobRegTypeIC'. }

\item{risk = "fiUnOvShoot", ErrorDistr = "Norm", Regressor = "UnivariateDistribution", neighbor = "CondTotalVarNeighborhood"}{ 
  computes finite-sample under-/overshoot risk in methods for 
  function 'getFixRobRegTypeIC'. }
}}
\references{
  Huber, P.J. (1968) Robust Confidence Limits. Z. Wahrscheinlichkeitstheor.
  Verw. Geb. \bold{10}:269--278.

  Rieder, H. (1989) A finite-sample minimax regression estimator. 
  Statistics \bold{20}(2): 211--221.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{ \code{\link[ROptEst]{fiRisk-class}} }
%\examples{}
\concept{finite-sample risk}
\concept{risk}
\keyword{robust}

