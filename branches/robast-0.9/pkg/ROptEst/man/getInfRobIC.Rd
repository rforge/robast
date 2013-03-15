\name{getInfRobIC}
\alias{getInfRobIC}
\alias{getInfRobIC-methods}
\alias{getInfRobIC,UnivariateDistribution,asCov,ContNeighborhood-method}
\alias{getInfRobIC,UnivariateDistribution,asCov,TotalVarNeighborhood-method}
\alias{getInfRobIC,RealRandVariable,asCov,UncondNeighborhood-method}
\alias{getInfRobIC,UnivariateDistribution,asBias,UncondNeighborhood-method}
\alias{getInfRobIC,RealRandVariable,asBias,UncondNeighborhood-method}
\alias{getInfRobIC,UnivariateDistribution,asHampel,UncondNeighborhood-method}
\alias{getInfRobIC,RealRandVariable,asHampel,UncondNeighborhood-method}
\alias{getInfRobIC,UnivariateDistribution,asAnscombe,UncondNeighborhood-method}
\alias{getInfRobIC,RealRandVariable,asAnscombe,UncondNeighborhood-method}
\alias{getInfRobIC,UnivariateDistribution,asGRisk,UncondNeighborhood-method}
\alias{getInfRobIC,RealRandVariable,asGRisk,UncondNeighborhood-method}
\alias{getInfRobIC,UnivariateDistribution,asUnOvShoot,UncondNeighborhood-method}

\title{ Generic Function for the Computation of Optimally Robust ICs }
\description{
  Generic function for the computation of optimally robust ICs 
  in case of infinitesimal robust models. This function is 
  rarely called directly. 
}
\usage{
getInfRobIC(L2deriv, risk, neighbor, ...)

\S4method{getInfRobIC}{UnivariateDistribution,asCov,ContNeighborhood}(L2deriv,
                       risk, neighbor, Finfo, trafo, verbose = NULL)

\S4method{getInfRobIC}{UnivariateDistribution,asCov,TotalVarNeighborhood}(L2deriv,
                       risk, neighbor, Finfo, trafo, verbose = NULL)

\S4method{getInfRobIC}{RealRandVariable,asCov,UncondNeighborhood}(L2deriv, risk,
                       neighbor, Distr, Finfo, trafo, QuadForm = diag(nrow(trafo)),
                       verbose = NULL)

\S4method{getInfRobIC}{UnivariateDistribution,asBias,UncondNeighborhood}(L2deriv,
                       risk, neighbor, symm, trafo, maxiter, tol, warn, Finfo,
                       verbose = NULL, ...)

\S4method{getInfRobIC}{RealRandVariable,asBias,UncondNeighborhood}(L2deriv, risk,
                       neighbor, Distr, DistrSymm, L2derivSymm,
                       L2derivDistrSymm, z.start, A.start, Finfo, trafo,
                       maxiter, tol, warn, verbose = NULL, ...)

\S4method{getInfRobIC}{UnivariateDistribution,asHampel,UncondNeighborhood}(L2deriv,
                       risk, neighbor, symm, Finfo, trafo, upper = NULL,
                       lower=NULL, maxiter, tol, warn, noLow = FALSE,
                       verbose = NULL, checkBounds = TRUE, ...)

\S4method{getInfRobIC}{RealRandVariable,asHampel,UncondNeighborhood}(L2deriv, risk,
                       neighbor, Distr, DistrSymm, L2derivSymm,
                       L2derivDistrSymm, Finfo, trafo, onesetLM = FALSE,
                       z.start, A.start, upper = NULL, lower=NULL,
                       OptOrIter = "iterate", maxiter, tol, warn,
                       verbose = NULL, checkBounds = TRUE, ...,
                       .withEvalAsVar = TRUE)

\S4method{getInfRobIC}{UnivariateDistribution,asAnscombe,UncondNeighborhood}(
                       L2deriv, risk, neighbor, symm, Finfo, trafo, upper = NULL,
                       lower=NULL, maxiter, tol, warn, noLow = FALSE,
                       verbose = NULL, checkBounds = TRUE, ...)

\S4method{getInfRobIC}{RealRandVariable,asAnscombe,UncondNeighborhood}(L2deriv, 
                       risk, neighbor, Distr, DistrSymm, L2derivSymm,
                       L2derivDistrSymm, Finfo, trafo, onesetLM = FALSE,
                       z.start, A.start, upper = NULL, lower=NULL,
                       OptOrIter = "iterate", maxiter, tol, warn,
                       verbose = NULL, checkBounds = TRUE, ...)

\S4method{getInfRobIC}{UnivariateDistribution,asGRisk,UncondNeighborhood}(L2deriv,
                       risk, neighbor, symm, Finfo, trafo, upper = NULL,
                       lower = NULL, maxiter, tol, warn, noLow = FALSE,
                       verbose = NULL, ...)

\S4method{getInfRobIC}{RealRandVariable,asGRisk,UncondNeighborhood}(L2deriv, risk,
                       neighbor,  Distr, DistrSymm, L2derivSymm,
                       L2derivDistrSymm, Finfo, trafo, onesetLM = FALSE, z.start,
                       A.start, upper = NULL, lower = NULL, OptOrIter = "iterate",
                       maxiter, tol, warn, verbose = NULL, withPICcheck = TRUE,
                       ..., .withEvalAsVar = TRUE)

\S4method{getInfRobIC}{UnivariateDistribution,asUnOvShoot,UncondNeighborhood}(
                       L2deriv, risk, neighbor, symm, Finfo, trafo,
                       upper, lower, maxiter, tol, warn, ...)
}
\arguments{
  \item{L2deriv}{ L2-derivative of some L2-differentiable family 
    of probability measures. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters (mainly for \code{optim}). }
  \item{Distr}{ object of class \code{"Distribution"}. }
  \item{symm}{ logical: indicating symmetry of \code{L2deriv}. }
  \item{DistrSymm}{ object of class \code{"DistributionSymmetry"}. }
  \item{L2derivSymm}{ object of class \code{"FunSymmList"}. }
  \item{L2derivDistrSymm}{ object of class \code{"DistrSymmList"}. }
  \item{Finfo}{ Fisher information matrix. }
  \item{z.start}{ initial value for the centering constant. }
  \item{A.start}{ initial value for the standardizing matrix. }
  \item{trafo}{ matrix: transformation of the parameter. }
  \item{upper}{ upper bound for the optimal clipping bound. }
  \item{lower}{ lower bound for the optimal clipping bound. }
  \item{OptOrIter}{character; which method to be used for determining Lagrange
  multipliers \code{A} and \code{a}: if (partially) matched to \code{"optimize"},
  \code{getLagrangeMultByOptim} is used; otherwise: by default, or if matched to
  \code{"iterate"} or to \code{"doubleiterate"},
  \code{getLagrangeMultByIter} is used. More specifically,
  when using \code{getLagrangeMultByIter}, and if argument \code{risk} is of
  class \code{"asGRisk"}, by default and if matched to \code{"iterate"}
  we use only one (inner) iteration, if matched to \code{"doubleiterate"}
  we use up to \code{Maxiter} (inner) iterations.}
  \item{maxiter}{ the maximum number of iterations. }
  \item{tol}{ the desired accuracy (convergence tolerance).}
  \item{warn}{ logical: print warnings. }
  \item{noLow}{ logical: is lower case to be computed? }
  \item{onesetLM}{ logical: use one set of Lagrange multipliers? }
  \item{QuadForm}{ matrix of (or which may coerced to) class 
                  \code{PosSemDefSymmMatrix} for use of different 
                  (standardizing) norm }
  \item{verbose}{ logical: if \code{TRUE}, some messages are printed }
  \item{checkBounds}{ logical: if \code{TRUE}, minimal and maximal clipping bound are 
      computed to check if a valid bound was specified. }
  \item{withPICcheck}{logical: at the end of the algorithm, shall we check
  how accurately this is a pIC; this will only be done if
  \code{withPICcheck && verbose}.}
  \item{.withEvalAsVar}{logical (of length 1):
     if \code{TRUE}, risks based on covariances are to be
     evaluated (default), otherwise just a call is returned.}

}
%\details{}
\value{The optimally robust IC is computed.}
\section{Methods}{
\describe{
  \item{L2deriv = "UnivariateDistribution", risk = "asCov", 
        neighbor = "ContNeighborhood"}{ 
    computes the classical optimal influence curve for L2 differentiable 
    parametric families with unknown one-dimensional parameter. }

  \item{L2deriv = "UnivariateDistribution", risk = "asCov", 
        neighbor = "TotalVarNeighborhood"}{
    computes the classical optimal influence curve for L2 differentiable 
    parametric families with unknown one-dimensional parameter. }

  \item{L2deriv = "RealRandVariable", risk = "asCov", 
        neighbor = "UncondNeighborhood"}{
    computes the classical optimal influence curve for L2 differentiable 
    parametric families with unknown \eqn{k}-dimensional parameter 
    (\eqn{k > 1}) where the underlying distribution is univariate;
    for total variation neighborhoods only is implemented for the case
    where there is a \eqn{1\times k}{1 x k} transformation \code{trafo}
    matrix. }

  \item{L2deriv = "UnivariateDistribution", risk = "asBias", 
        neighbor = "UncondNeighborhood"}{ 
    computes the bias optimal influence curve for L2 differentiable 
    parametric families with unknown one-dimensional parameter. }

  \item{L2deriv = "RealRandVariable", risk = "asBias", 
        neighbor = "UncondNeighborhood"}{
    computes the bias optimal influence curve for L2 differentiable 
    parametric families with unknown \eqn{k}-dimensional parameter 
    (\eqn{k > 1}) where the underlying distribution is univariate. }

  \item{L2deriv = "UnivariateDistribution", risk = "asHampel", 
        neighbor = "UncondNeighborhood"}{
    computes the optimally robust influence curve for L2 differentiable 
    parametric families with unknown one-dimensional parameter. }

  \item{L2deriv = "RealRandVariable", risk = "asHampel", 
        neighbor = "UncondNeighborhood"}{
    computes the optimally robust influence curve for L2 differentiable 
    parametric families with unknown \eqn{k}-dimensional parameter 
    (\eqn{k > 1}) where the underlying distribution is univariate;
    for total variation neighborhoods only is implemented for the case
    where there is a \eqn{1\times k}{1 x k} transformation \code{trafo}
    matrix. }

  \item{L2deriv = "UnivariateDistribution", risk = "asAnscombe", 
        neighbor = "UncondNeighborhood"}{
    computes the optimally bias-robust influence curve to given ARE in the
    ideal model for L2 differentiable 
    parametric families with unknown one-dimensional parameter. }

  \item{L2deriv = "RealRandVariable", risk = "asAnscombe", 
        neighbor = "UncondNeighborhood"}{
    computes the optimally bias-robust influence curve to given ARE in the
    ideal modelfor L2 differentiable 
    parametric families with unknown \eqn{k}-dimensional parameter 
    (\eqn{k > 1}) where the underlying distribution is univariate;
    for total variation neighborhoods only is implemented for the case
    where there is a \eqn{1\times k}{1 x k} transformation \code{trafo}
    matrix. }

  \item{L2deriv = "UnivariateDistribution", risk = "asGRisk", 
        neighbor = "UncondNeighborhood"}{
    computes the optimally robust influence curve for L2 differentiable 
    parametric families with unknown one-dimensional parameter. }

  \item{L2deriv = "RealRandVariable", risk = "asGRisk", 
        neighbor = "UncondNeighborhood"}{
    computes the optimally robust influence curve for L2 differentiable 
    parametric families with unknown \eqn{k}-dimensional parameter 
    (\eqn{k > 1}) where the underlying distribution is univariate;
    for total variation neighborhoods only is implemented for the case
    where there is a \eqn{1\times k}{1 x k} transformation \code{trafo}
    matrix. }

  \item{L2deriv = "UnivariateDistribution", risk = "asUnOvShoot", 
        neighbor = "UncondNeighborhood"}{
    computes the optimally robust influence curve for one-dimensional
    L2 differentiable parametric families and 
    asymptotic under-/overshoot risk. }
}}
\references{
  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106-115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Ruckdeschel, P. and Rieder, H. (2004) Optimal Influence Curves for
  General Loss Functions. Statistics & Decisions \bold{22}: 201-223.

  Ruckdeschel, P. (2005) Optimally One-Sided Bounded Influence Curves.
  Mathematical Methods in Statistics \emph{14}(1), 105-131.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
        Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link[RobAStBase]{InfRobModel-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
