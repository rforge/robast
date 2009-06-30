\name{optRisk}
\alias{optRisk}
\alias{optRisk-methods}
\alias{optRisk,L2ParamFamily,asCov-method}
\alias{optRisk,InfRobModel,asRisk-method}
\alias{optRisk,FixRobModel,fiUnOvShoot-method}

\title{Generic function for the computation of the minimal risk}
\description{
  Generic function for the computation of the optimal (i.e., minimal) 
  risk for a probability model.
}
\usage{
optRisk(model, risk, ...)

\S4method{optRisk}{InfRobModel,asRisk}(model, risk, z.start = NULL, A.start = NULL, upper = 1e4, 
             maxiter = 50, tol = .Machine$double.eps^0.4, warn = TRUE)

\S4method{optRisk}{FixRobModel,fiUnOvShoot}(model, risk, sampleSize, upper = 1e4, maxiter = 50, 
             tol = .Machine$double.eps^0.4, warn = TRUE, Algo = "A", cont = "left")
}
\arguments{
  \item{model}{ probability model }
  \item{risk}{ object of class \code{RiskType} }
  \item{\dots}{ additional parameters }
  \item{z.start}{ initial value for the centering constant. }
  \item{A.start}{ initial value for the standardizing matrix. }
  \item{upper}{ upper bound for the optimal clipping bound. }
  \item{maxiter}{ the maximum number of iterations }
  \item{tol}{ the desired accuracy (convergence tolerance).}
  \item{warn}{ logical: print warnings. }    
  \item{sampleSize}{ integer: sample size. }
  \item{Algo}{ "A" or "B". }
  \item{cont}{ "left" or "right". }
}
\details{ In case of the finite-sample risk \code{"fiUnOvShoot"} one can choose
  between two algorithms for the computation of this risk where the least favorable
  contamination is assumed to be left or right of some bound. For more details
  we refer to Section 11.3 of Kohl (2005).}
\value{The minimal risk is computed.}
\section{Methods}{
\describe{
  \item{model = "L2ParamFamily", risk = "asCov"}{ 
    asymptotic covariance of L2 differentiable parameteric
    family. }

  \item{model = "InfRobModel", risk = "asRisk"}{ 
    asymptotic risk of a infinitesimal robust model. }

  \item{model = "FixRobModel", risk = "fiUnOvShoot"}{ 
    finite-sample under-/overshoot risk of a robust model
    with fixed neighborhood. }
}}
\references{
  Huber, P.J. (1968) Robust Confidence Limits. Z. Wahrscheinlichkeitstheor.
  Verw. Geb. \bold{10}:269--278.

  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{RiskType-class}}}
\examples{
optRisk(model = NormLocationScaleFamily(), risk = asCov())
}
\concept{risk}
\keyword{robust}
