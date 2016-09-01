\name{optIC}
\alias{optIC}
\alias{optIC-methods}
\alias{optIC,L2ParamFamily,asCov-method}
\alias{optIC,InfRobModel,asRisk-method}
\alias{optIC,InfRobModel,asUnOvShoot-method}
\alias{optIC,FixRobModel,fiUnOvShoot-method}

\title{Generic function for the computation of optimally robust ICs}
\description{
  Generic function for the computation of optimally robust ICs.
}
\usage{
optIC(model, risk, ...)

\S4method{optIC}{L2ParamFamily,asCov}(model, risk)

\S4method{optIC}{InfRobModel,asRisk}(model, risk, z.start = NULL, A.start = NULL, upper = 1e4, 
             maxiter = 50, tol = .Machine$double.eps^0.4, warn = TRUE)

\S4method{optIC}{InfRobModel,asUnOvShoot}(model, risk, upper = 1e4, maxiter = 50, 
             tol = .Machine$double.eps^0.4, warn = TRUE)

\S4method{optIC}{FixRobModel,fiUnOvShoot}(model, risk, sampleSize, upper = 1e4, maxiter = 50, 
             tol = .Machine$double.eps^0.4, warn = TRUE, Algo = "A", cont = "left")
}
\arguments{
  \item{model}{ probability model. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{\dots}{ additional parameters. }
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
  contamination is assumed to be left or right of some bound. For more details
  we refer to Section 11.3 of Kohl (2005).}
\value{Some optimally robust IC is computed.}
\section{Methods}{
\describe{
  \item{model = "L2ParamFamily", risk = "asCov"}{ computes
    classical optimal influence curve for L2 differentiable 
    parametric families.}

  \item{model = "InfRobModel", risk = "asRisk"}{ 
    computes optimally robust influence curve for 
    robust models with infinitesimal neighborhoods and
    various asymptotic risks. }

  \item{model = "InfRobModel", risk = "asUnOvShoot"}{ 
    computes optimally robust influence curve for 
    robust models with infinitesimal neighborhoods and
    asymptotic under-/overshoot risk. }

  \item{model = "FixRobModel", risk = "fiUnOvShoot"}{ 
    computes optimally robust influence curve for 
    robust models with fixed neighborhoods and
    finite-sample under-/overshoot risk. }
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
\seealso{\code{\link{InfluenceCurve-class}}, \code{\link{RiskType-class}}}
\examples{
B <- BinomFamily(size = 25, prob = 0.25) 

## classical optimal IC
IC0 <- optIC(model = B, risk = asCov())
plot(IC0) # plot IC
checkIC(IC0, B)
}
\concept{robust influence curve}
\concept{influence curve}
\keyword{robust}
