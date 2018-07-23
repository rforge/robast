\name{optIC-methods}
\docType{methods}
\alias{optIC-methods}
\alias{optIC,L2RegTypeFamily,asCov-method}
\alias{optIC,InfRobRegTypeModel,asRisk-method}
\alias{optIC,InfRobRegTypeModel,asUnOvShoot-method}
\alias{optIC,FixRobRegTypeModel,fiUnOvShoot-method}

\title{Methods for Function optIC in Package `ROptRegTS'}
\description{
  Methods for function \code{optIC} in package \pkg{ROptRegTS}.
}
\usage{
\S4method{optIC}{L2RegTypeFamily,asCov}(model, risk)

\S4method{optIC}{InfRobRegTypeModel,asRisk}(model, risk, z.start = NULL,
                 A.start = NULL, upper = 1e4, maxiter = 50,
                 tol = .Machine$double.eps^0.4, warn = TRUE)

\S4method{optIC}{InfRobRegTypeModel,asUnOvShoot}(model, risk, upper = 1e4,
                 maxiter = 50, tol = .Machine$double.eps^0.4, warn = TRUE)

\S4method{optIC}{FixRobRegTypeModel,fiUnOvShoot}(model, risk, sampleSize,
                 upper = 1e4, maxiter = 50, tol = .Machine$double.eps^0.4,
                 warn = TRUE, Algo = "A", cont = "left")
}
\arguments{
  \item{model}{ probability model. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{z.start}{ initial value for the centering constant. }
  \item{A.start}{ initial value for the standardizing matrix. }
  \item{upper}{ upper bound for the optimal clipping bound. }
  \item{maxiter}{ the maximum number of iterations. }
  \item{tol}{ the desired accuracy (convergence tolerance).}
  \item{warn}{ logical: print warnings. }  
  \item{sampleSize}{ integer: sample size. }
  \item{Algo}{ "A" or "B". }
  \item{cont}{ "left" or "right". }
}
\details{ In case of the finite-sample risk \code{"fiUnOvShoot"} one can choose
  between two algorithms for the computation of this risk where the least favorable
  contamination is assumed to be \dQuote{left} or \dQuote{right} of some boundary 
  curve. For more details we refer to Subsections 12.1.3 and 12.2.3 of Kohl (2005).}
\value{Some optimally robust IC is computed.}
\section{Methods}{
\describe{
  \item{model = "L2RegTypeFamily", risk = "asCov"}{ 
    computes classical optimal influence curve for 
    L2 differentiable regression-type families. }

  \item{model = "InfRobRegTypeModel", risk = "asRisk"}{ 
    computes optimally robust influence curve for 
    robust regression-type models with infinitesimal 
    neighborhoods and various asymptotic risks. }

  \item{model = "InfRobRegTypeModel", risk = "asUnOvShoot"}{ 
    computes optimally robust influence curve for 
    robust regression-type models with infinitesimal 
    neighborhoods and asymptotic under-/overshoot risk. }

  \item{model = "FixRobRegTypeModel", risk = "fiUnOvShoot"}{ 
    computes optimally robust influence curve for 
    robust regression-type models with fixed
    neighborhoods and finite-sample under-/overshoot risk. }
}
}
\references{
  Huber, P.J. (1968) Robust Confidence Limits. Z. Wahrscheinlichkeitstheor.
  Verw. Geb. \bold{10}:269--278.

  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\seealso{\code{\link[ROptEstOld]{optIC}}}
\concept{robust influence curve}
\concept{influence curve}
\keyword{methods}
\keyword{robust}
