\name{oneStepEstimator}
\alias{oneStepEstimator}
\alias{oneStepEstimator-methods}
\alias{oneStepEstimator,numeric,InfluenceCurve,numeric-method}
\alias{oneStepEstimator,numeric,InfluenceCurve,list-method}
\alias{oneStepEstimator,matrix,InfluenceCurve,numeric-method}
\alias{oneStepEstimator,matrix,InfluenceCurve,list-method}

\title{Generic function for the computation of one-step estimators}
\description{
  Generic function for the computation of one-step estimators.
}
\usage{
oneStepEstimator(x, IC, start)
}
\arguments{
  \item{x}{ sample }
  \item{IC}{ object of class \code{"InfluenceCurve"} }
  \item{start}{ initial estimate }
}
\details{
  Given an initial estimation \code{start}, a sample \code{x} 
  and an influence curve \code{IC} the corresponding one-step
  estimator is computed
}
\value{The one-step estimation is computed.}
\section{Methods}{
\describe{
  \item{x = "numeric", IC = "InfluenceCurve", start = "numeric"}{ 
    univariate samples. }
  \item{x = "numeric", IC = "InfluenceCurve", start = "list"}{ 
    univariate samples. }
  \item{x = "matrix", IC = "InfluenceCurve", start = "numeric"}{ 
    multivariate samples. }
  \item{x = "matrix", IC = "InfluenceCurve", start = "list"}{ 
    multivariate samples. }
}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{InfluenceCurve-class}}}
%\examples{}
\concept{one-step estimator}
\concept{estimator}
\keyword{robust}
