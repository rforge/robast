\name{getFixClipRegTS}
\alias{getFixClipRegTS}
\alias{getFixClipRegTS-methods}
\alias{getFixClipRegTS,numeric,Norm,UnivariateDistribution,fiUnOvShoot,ContNeighborhood-method}
\alias{getFixClipRegTS,numeric,Norm,UnivariateDistribution,fiUnOvShoot,TotalVarNeighborhood-method}
\alias{getFixClipRegTS,numeric,Norm,numeric,fiUnOvShoot,CondContNeighborhood-method}
\alias{getFixClipRegTS,numeric,Norm,numeric,fiUnOvShoot,CondTotalVarNeighborhood-method}

\title{Generic Function for the Computation of the Optimal Clipping Bound}
\description{
  Generic function for the computation of the optimal clipping bound/function.
  This function is rarely called directly. It is used to compute optimally 
  robust ICs in case of fixed robust models.
}
\usage{
getFixClipRegTS(clip, ErrorDistr, Regressor, risk, neighbor, ...)
}
\arguments{
  \item{clip}{ optimal clipping bound. }
  \item{ErrorDistr}{ error distribution. }
  \item{Regressor}{ regressor. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters. }
}
%\details{}
\value{The optimal clipping bound/function is computed.}
\section{Methods}{
\describe{
\item{clip = "numeric", ErrorDistr = "Norm", Regressor = "UnivariateDistribution", risk = "fiUnOvShoot", neighbor = "ContNeighborhood"}{ 
  optimal clipping bound for finite-sample under-/overshoot risk. }

\item{clip = "numeric", ErrorDistr = "Norm", Regressor = "UnivariateDistribution", risk = "fiUnOvShoot", neighbor = "TotalVarNeighborhood"}{ 
  optimal clipping bound for finite-sample under-/overshoot risk. }

\item{clip = "numeric", ErrorDistr = "Norm", Regressor = "numeric", risk = "fiUnOvShoot", neighbor = "CondContNeighborhood"}{ 
  optimal clipping function for finite-sample under-/overshoot risk. }

\item{clip = "numeric", ErrorDistr = "Norm", Regressor = "numeric", risk = "fiUnOvShoot", neighbor = "CondTotalVarNeighborhood"}{ 
  optimal clipping function for finite-sample under-/overshoot risk. }
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
\seealso{\code{\link[ROptEst]{ContIC-class}}, \code{\link[ROptEst]{TotalVarIC-class}}, 
        \code{\link{Av1CondContIC-class}}, \code{\link{Av2CondContIC-class}}, 
        \code{\link{Av1CondTotalVarIC-class}}, \code{\link{CondContIC-class}},
        \code{\link{CondTotalVarIC-class}}}
%\examples{}
\concept{influence curve}
\keyword{}
