\name{locMEstimator}
\alias{locMEstimator}
\alias{locMEstimator-methods}
\alias{locMEstimator,numeric,InfluenceCurve-method}

\title{Generic function for the computation of location M estimators}
\description{
  Generic function for the computation of location M estimators.
}
\usage{
locMEstimator(x, IC, ...)

\S4method{locMEstimator}{numeric,InfluenceCurve}(x, IC, eps = .Machine$double.eps^0.5)
}
\arguments{
  \item{x}{ sample }
  \item{IC}{ object of class \code{"InfluenceCurve"} }
  \item{\dots}{ additional parameters }
  \item{eps}{ the desired accuracy (convergence tolerance). }  
}
%\details{}
\value{Returns a list with component
  \item{loc}{ M estimator of location }
}
\section{Methods}{
\describe{
  \item{x = "numeric", IC = "InfluenceCurve"}{ univariate location. }
}}
\references{
  Huber, P.J. (1964) Robust estimation of a location parameter. 
  Ann. Math. Stat. \bold{35}: 73--101.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{InfluenceCurve-class}}}
%\examples{}
\concept{M estimator}
\concept{estimator}
\keyword{robust}
