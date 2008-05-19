\name{getFixRobRegTypeIC}
\alias{getFixRobRegTypeIC}
\alias{getFixRobRegTypeIC-methods}
\alias{getFixRobRegTypeIC,Norm,UnivariateDistribution,fiUnOvShoot,UncondNeighborhood-method}
\alias{getFixRobRegTypeIC,Norm,UnivariateDistribution,fiUnOvShoot,CondNeighborhood-method}

\title{ Generic Function for the Computation of Optimally Robust Regression-Type ICs }
\description{
  Generic function for the computation of optimally robust regression-type ICs 
  in case of fixed robust models. This function is rarely called directly.
}
\usage{
getFixRobRegTypeIC(ErrorDistr, Regressor, risk, neighbor, ...)

\S4method{getFixRobRegTypeIC}{Norm,UnivariateDistribution,fiUnOvShoot,UncondNeighborhood}(ErrorDistr, 
            Regressor, risk, neighbor, sampleSize, upper, maxiter, tol, warn, Algo, cont)

\S4method{getFixRobRegTypeIC}{Norm,UnivariateDistribution,fiUnOvShoot,CondNeighborhood}(ErrorDistr, 
            Regressor, risk, neighbor, sampleSize, upper, maxiter, tol, warn, Algo, cont)
}
\arguments{
  \item{ErrorDistr}{ error distribution }
  \item{Regressor}{ regressor }
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
\item{ErrorDistr = "Norm", Regressor = "UnivariateDistribution", risk = "fiUnOvShoot", neighbor = "UncondNeighborhood"}{ 
  computes the optimally robust influence curve for one-dimensional normal regression 
  and finite-sample under-/overshoot risk. }

\item{ErrorDistr = "Norm", Regressor = "UnivariateDistribution", risk = "fiUnOvShoot", neighbor = "CondNeighborhood"}{ 
  computes the optimally robust influence curve for one-dimensional normal regression 
  and finite-sample under-/overshoot risk. }
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
\seealso{\code{\link{FixRobRegTypeModel-class}}}
%\examples{}
\concept{influence curve}
\keyword{}
