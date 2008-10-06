\name{getInfRobRegTypeIC}
\alias{getInfRobRegTypeIC}
\alias{getInfRobRegTypeIC-methods}
\alias{getInfRobRegTypeIC,UnivariateDistribution,UnivariateDistribution,asBias,ContNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,UnivariateDistribution,asBias,Av1CondContNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,UnivariateDistribution,asBias,Av1CondTotalVarNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,Distribution,asBias,Av2CondContNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,Distribution,asCov,ContNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,UnivariateDistribution,asCov,TotalVarNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,Distribution,asCov,CondContNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,Distribution,asCov,CondTotalVarNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,Distribution,asCov,Av1CondContNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,Distribution,asCov,Av2CondContNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,Distribution,asCov,Av1CondTotalVarNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,Distribution,asGRisk,ContNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,Distribution,asGRisk,Av1CondContNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,Distribution,asGRisk,Av2CondContNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,Distribution,asGRisk,Av1CondTotalVarNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,MultivariateDistribution,asBias,ContNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,MultivariateDistribution,asBias,Av1CondContNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,MultivariateDistribution,asBias,Av1CondTotalVarNeighborhood-method}
\alias{getInfRobRegTypeIC,RealRandVariable,Distribution,asBias,ContNeighborhood-method}
\alias{getInfRobRegTypeIC,RealRandVariable,Distribution,asBias,Av1CondContNeighborhood-method}
\alias{getInfRobRegTypeIC,RealRandVariable,Distribution,asCov,ContNeighborhood-method}
\alias{getInfRobRegTypeIC,RealRandVariable,Distribution,asCov,Av1CondContNeighborhood-method}
\alias{getInfRobRegTypeIC,RealRandVariable,Distribution,asGRisk,ContNeighborhood-method}
\alias{getInfRobRegTypeIC,RealRandVariable,Distribution,asGRisk,Av1CondContNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,UnivariateDistribution,asUnOvShoot,UncondNeighborhood-method}
\alias{getInfRobRegTypeIC,UnivariateDistribution,UnivariateDistribution,asUnOvShoot,CondNeighborhood-method}

\title{ Generic Function for the Computation of Optimally Robust Regression-Type ICs }
\description{
  Generic function for the computation of optimally robust regression-type ICs 
  in case of infinitesimal robust models. This function is rarely called directly.
}
\usage{
getInfRobRegTypeIC(ErrorL2deriv, Regressor, risk, neighbor, ...)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,UnivariateDistribution,asBias,ContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo, 
                upper, maxiter, tol, warn)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,UnivariateDistribution,asBias,Av1CondContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo, 
                upper, maxiter, tol, warn)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,UnivariateDistribution,asBias,Av1CondTotalVarNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo, 
                upper, maxiter, tol, warn)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,Distribution,asBias,Av2CondContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo, 
                upper, maxiter, tol, warn)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,Distribution,asCov,ContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,UnivariateDistribution,asCov,TotalVarNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,Distribution,asCov,CondContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,Distribution,asCov,CondTotalVarNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,Distribution,asCov,Av1CondContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,Distribution,asCov,Av2CondContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,Distribution,asCov,Av1CondTotalVarNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,Distribution,asGRisk,ContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo, 
                upper, maxiter, tol, warn)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,Distribution,asGRisk,Av1CondContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo, 
                upper, maxiter, tol, warn)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,Distribution,asGRisk,Av2CondContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo, 
                upper, maxiter, tol, warn)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,Distribution,asGRisk,Av1CondTotalVarNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo, 
                upper, maxiter, tol, warn)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,MultivariateDistribution,asBias,ContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo, 
                upper, maxiter, tol, warn)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,MultivariateDistribution,asBias,Av1CondContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo, 
                upper, maxiter, tol, warn)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,MultivariateDistribution,asBias,Av1CondTotalVarNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo, 
                upper, maxiter, tol, warn)

\S4method{getInfRobRegTypeIC}{RealRandVariable,Distribution,asBias,ContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorSymm, RegSymm, ErrorDistr, ErrorL2derivSymm, 
                ErrorL2derivDistrSymm, Finfo, trafo, upper, z.start, A.start, maxiter, tol, warn)

\S4method{getInfRobRegTypeIC}{RealRandVariable,Distribution,asBias,Av1CondContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorSymm, RegSymm, ErrorDistr, ErrorL2derivSymm, 
                ErrorL2derivDistrSymm, Finfo, trafo, upper, z.start, A.start, maxiter, tol, warn)

\S4method{getInfRobRegTypeIC}{RealRandVariable,Distribution,asCov,ContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorDistr, Finfo, trafo)

\S4method{getInfRobRegTypeIC}{RealRandVariable,Distribution,asCov,Av1CondContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorDistr, Finfo, trafo)

\S4method{getInfRobRegTypeIC}{RealRandVariable,Distribution,asGRisk,ContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorSymm, RegSymm, ErrorDistr, ErrorL2derivSymm, 
                ErrorL2derivDistrSymm, Finfo, trafo, upper, z.start, A.start, maxiter, tol, warn)

\S4method{getInfRobRegTypeIC}{RealRandVariable,Distribution,asGRisk,Av1CondContNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorSymm, RegSymm, ErrorDistr, ErrorL2derivSymm, 
                ErrorL2derivDistrSymm, Finfo, trafo, upper, z.start, A.start, maxiter, tol, warn)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,UnivariateDistribution,asUnOvShoot,UncondNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo, 
                upper, maxiter, tol, warn)

\S4method{getInfRobRegTypeIC}{UnivariateDistribution,UnivariateDistribution,asUnOvShoot,CondNeighborhood}(ErrorL2deriv, 
                Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, Finfo, trafo, 
                upper, maxiter, tol, warn)
}
\arguments{
  \item{ErrorL2deriv}{ L2-derivative of \code{ErrorDistr}. }
  \item{Regressor}{ regressor. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters. }
  \item{ErrorSymm}{ symmetry of \code{ErrorDistr}. }
  \item{ErrorL2derivDistrSymm}{ symmetry of \code{ErrorL2derivDistr}. }
  \item{RegSymm}{ symmetry of \code{RegDistr}. }
  \item{ErrorDistr}{ error distribution. }
  \item{ErrorL2derivSymm}{ symmetry of \code{ErrorL2deriv}. }
  \item{Finfo}{ Fisher information matrix. }
  \item{trafo}{ matrix: transformation of the parameter. }
  \item{upper}{ upper bound for the optimal clipping bound. }
  \item{maxiter}{ the maximum number of iterations }
  \item{tol}{ the desired accuracy (convergence tolerance).}
  \item{warn}{ logical: print warnings. }
  \item{z.start}{ initial value for the centering constant/function. }
  \item{A.start}{ initial value for the standardizing matrix. }
}
%\details{}
\value{The optimally robust IC is computed.}
\section{Methods}{
\describe{
  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", risk = "asBias", neighbor = "ContNeighborhood"}{ 
    computes the bias optimal influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", risk = "asBias", neighbor = "Av1CondContNeighborhood"}{ 
    computes the bias optimal influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", risk = "asBias", neighbor = "Av1CondTotalVarNeighborhood"}{ 
    computes the bias optimal influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", risk = "asBias", neighbor = "Av2CondContNeighborhood"}{ 
    computes the bias optimal influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", risk = "asCov", neighbor = "ContNeighborhood"}{ 
    computes the classical optimal influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", risk = "asCov", neighbor = "TotalVarNeighborhood"}{ 
    computes the classical optimal influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", risk = "asCov", neighbor = "CondContNeighborhood"}{ 
    computes the classical optimal influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", risk = "asCov", neighbor = "CondTotalVarNeighborhood"}{ 
    computes the classical optimal influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", risk = "asCov", neighbor = "Av1CondContNeighborhood"}{ 
    computes the classical optimal influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", risk = "asCov", neighbor = "Av2CondContNeighborhood"}{ 
    computes the classical optimal influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", risk = "asCov", neighbor = "Av1CondTotalVarNeighborhood"}{ 
    computes the classical optimal influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", risk = "asGRisk", neighbor = "ContNeighborhood"}{ 
    computes the optimally robust influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", risk = "asGRisk", neighbor = "Av1CondContNeighborhood"}{ 
    computes the optimally robust influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", risk = "asGRisk", neighbor = "Av2CondContNeighborhood"}{ 
    computes the optimally robust influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "Distribution", risk = "asGRisk", neighbor = "Av1CondTotalVarNeighborhood"}{ 
    computes the optimally robust influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "MultivariateDistribution", risk = "asBias", neighbor = "ContNeighborhood"}{ 
    computes the bias optimal influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "MultivariateDistribution", risk = "asBias", neighbor = "Av1CondContNeighborhood"}{ 
    computes the bias optimal influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "MultivariateDistribution", risk = "asBias", neighbor = "Av1CondTotalVarNeighborhood"}{ 
    computes the bias optimal influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "RealRandVariable", Regressor = "Distribution", risk = "asBias", neighbor = "ContNeighborhood"}{ 
    computes the bias optimal influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "RealRandVariable", Regressor = "Distribution", risk = "asCov", neighbor = "ContNeighborhood"}{ 
    computes the classical optimal influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "RealRandVariable", Regressor = "Distribution", risk = "asCov", neighbor = "Av1CondContNeighborhood"}{ 
    computes the classical optimal influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "RealRandVariable", Regressor = "Distribution", risk = "asGRisk", neighbor = "ContNeighborhood"}{ 
    computes the optimally robust influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "RealRandVariable", Regressor = "Distribution", risk = "asGRisk", neighbor = "Av1CondContNeighborhood"}{
    computes the optimally robust influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", risk = "asUnOvShoot", neighbor = "UncondNeighborhood"}{ 
    computes the optimally robust influence curve for L2 differentiable regression-type families. }

  \item{ErrorL2deriv = "UnivariateDistribution", Regressor = "UnivariateDistribution", risk = "asUnOvShoot", neighbor = "CondNeighborhood"}{ 
    computes the optimally robust influence curve for L2 differentiable regression-type families. }
}}
\references{
  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{InfRobRegTypeModel-class}}}
%\examples{}
\concept{influence curve}
\keyword{}
