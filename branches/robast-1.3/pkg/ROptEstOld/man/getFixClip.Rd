\name{getFixClip}
\alias{getFixClip}
\alias{getFixClip-methods}
\alias{getFixClip,numeric,Norm,fiUnOvShoot,ContNeighborhood-method}
\alias{getFixClip,numeric,Norm,fiUnOvShoot,TotalVarNeighborhood-method}

\title{Generic Function for the Computation of the Optimal Clipping Bound}
\description{
  Generic function for the computation of the optimal clipping bound
  in case of robust models with fixed neighborhoods. This function is 
  rarely called directly. It is used to compute optimally robust ICs.
}
\usage{
getFixClip(clip, Distr, risk, neighbor, ...)

\S4method{getFixClip}{numeric,Norm,fiUnOvShoot,ContNeighborhood}(clip, Distr, risk, neighbor)

\S4method{getFixClip}{numeric,Norm,fiUnOvShoot,TotalVarNeighborhood}(clip, Distr, risk, neighbor)
}
\arguments{
  \item{clip}{ positive real: clipping bound }
  \item{Distr}{ object of class \code{"Distribution"}. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters. }
}
%\details{}
\value{The optimal clipping bound is computed.}
\section{Methods}{
\describe{
  \item{clip = "numeric", Distr = "Norm", risk = "fiUnOvShoot", neighbor = "ContNeighborhood"}{ 
    optimal clipping bound for finite-sample under-/overshoot risk. }

  \item{clip = "numeric", Distr = "Norm", risk = "fiUnOvShoot", neighbor = "TotalVarNeighborhood"}{ 
    optimal clipping bound for finite-sample under-/overshoot risk. }
}}
\references{
  Huber, P.J. (1968) Robust Confidence Limits. Z. Wahrscheinlichkeitstheor.
  Verw. Geb. \bold{10}:269--278.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{ContIC-class}}, \code{\link{TotalVarIC-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
