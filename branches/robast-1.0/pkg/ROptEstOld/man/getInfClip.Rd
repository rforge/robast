\name{getInfClip}
\alias{getInfClip}
\alias{getInfClip-methods}
\alias{getInfClip,numeric,UnivariateDistribution,asMSE,ContNeighborhood-method}
\alias{getInfClip,numeric,UnivariateDistribution,asMSE,TotalVarNeighborhood-method}
\alias{getInfClip,numeric,EuclRandVariable,asMSE,ContNeighborhood-method}
\alias{getInfClip,numeric,UnivariateDistribution,asUnOvShoot,UncondNeighborhood-method}

\title{Generic Function for the Computation of the Optimal Clipping Bound}
\description{
  Generic function for the computation of the optimal clipping bound
  in case of infinitesimal robust models. This function is rarely called 
  directly. It is used to compute optimally robust ICs.
}
\usage{
getInfClip(clip, L2deriv, risk, neighbor, ...)

\S4method{getInfClip}{numeric,UnivariateDistribution,asMSE,ContNeighborhood}(clip, L2deriv, risk, neighbor, cent, symm, trafo)

\S4method{getInfClip}{numeric,UnivariateDistribution,asMSE,TotalVarNeighborhood}(clip, L2deriv, risk, neighbor, cent, symm, trafo)

\S4method{getInfClip}{numeric,EuclRandVariable,asMSE,ContNeighborhood}(clip, L2deriv, risk, neighbor, Distr, stand, cent, trafo)

\S4method{getInfClip}{numeric,UnivariateDistribution,asUnOvShoot,UncondNeighborhood}(clip, L2deriv, risk, neighbor, cent, symm, trafo)
}
\arguments{
  \item{clip}{ positive real: clipping bound }
  \item{L2deriv}{ L2-derivative of some L2-differentiable family 
    of probability measures. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters. }
  \item{cent}{ optimal centering constant. }
  \item{stand}{ standardizing matrix. }
  \item{Distr}{ object of class \code{"Distribution"}. }
  \item{symm}{ logical: indicating symmetry of \code{L2deriv}. }
  \item{trafo}{ matrix: transformation of the parameter. }
}
%\details{}
\value{The optimal clipping bound is computed.}
\section{Methods}{
\describe{
  \item{clip = "numeric", L2deriv = "UnivariateDistribution", 
        risk = "asMSE", neighbor = "ContNeighborhood"}{ 
    optimal clipping bound for asymtotic mean square error. }

  \item{clip = "numeric", L2deriv = "UnivariateDistribution", 
        risk = "asMSE", neighbor = "TotalVarNeighborhood"}{ 
    optimal clipping bound for asymtotic mean square error. }

  \item{clip = "numeric", L2deriv = "EuclRandVariable", 
        risk = "asMSE", neighbor = "ContNeighborhood"}{ 
    optimal clipping bound for asymtotic mean square error. }

  \item{clip = "numeric", L2deriv = "UnivariateDistribution", 
        risk = "asUnOvShoot", neighbor = "UncondNeighborhood"}{ 
    optimal clipping bound for asymtotic under-/overshoot risk. }
}}
\references{
  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{ContIC-class}}, \code{\link{TotalVarIC-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
