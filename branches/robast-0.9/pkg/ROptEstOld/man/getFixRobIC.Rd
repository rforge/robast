\name{getFixRobIC}
\alias{getFixRobIC}
\alias{getFixRobIC-methods}
\alias{getFixRobIC,Norm,fiUnOvShoot,UncondNeighborhood-method}

\title{ Generic Function for the Computation of Optimally Robust ICs }
\description{
  Generic function for the computation of optimally robust ICs 
  in case of robust models with fixed neighborhoods. This function is 
  rarely called directly. 
}
\usage{
getFixRobIC(Distr, risk, neighbor, ...)

\S4method{getFixRobIC}{Norm,fiUnOvShoot,UncondNeighborhood}(Distr, risk, neighbor, 
            sampleSize, upper, maxiter, tol, warn, Algo, cont)
}
\arguments{
  \item{Distr}{ object of class \code{"Distribution"}. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters. }
  \item{sampleSize}{ integer: sample size. }  
  \item{upper}{ upper bound for the optimal clipping bound. }
  \item{maxiter}{ the maximum number of iterations. }
  \item{tol}{ the desired accuracy (convergence tolerance).}
  \item{warn}{ logical: print warnings. }
  \item{Algo}{ "A" or "B". }
  \item{cont}{ "left" or "right". }  
}
%\details{}
\value{The optimally robust IC is computed.}
\section{Methods}{
\describe{
  \item{Distr = "Norm", risk = "fiUnOvShoot", neighbor = "UncondNeighborhood"}{ 
    computes the optimally robust influence curve for one-dimensional
    normal location and finite-sample under-/overshoot risk. }
}}
\references{
  Huber, P.J. (1968) Robust Confidence Limits. Z. Wahrscheinlichkeitstheor.
  Verw. Geb. \bold{10}:269--278.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{FixRobModel-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
