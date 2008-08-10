\name{oneStepEstimator}
\alias{oneStepEstimator}
\alias{oneStepEstimator-methods}
\alias{oneStepEstimator,numeric,InfluenceCurve,numeric-method}
\alias{oneStepEstimator,matrix,InfluenceCurve,numeric-method}
\alias{oneStepEstimator,numeric,InfluenceCurve,Estimate-method}
\alias{oneStepEstimator,matrix,InfluenceCurve,Estimate-method}

\title{Generic function for the computation of one-step estimates}
\description{
  Generic function for the computation of one-step estimates.
}
\usage{
oneStepEstimator(x, IC, start, ...)

\S4method{oneStepEstimator}{numeric,InfluenceCurve,numeric}(x, IC, start, 
    useLast = getRobAStBaseOption("kStepUseLast"))
\S4method{oneStepEstimator}{matrix,InfluenceCurve,numeric}(x, IC, start, 
    useLast = getRobAStBaseOption("kStepUseLast"))
\S4method{oneStepEstimator}{numeric,InfluenceCurve,Estimate}(x, IC, start, 
    useLast = getRobAStBaseOption("kStepUseLast"))
\S4method{oneStepEstimator}{matrix,InfluenceCurve,Estimate}(x, IC, start, 
    useLast = getRobAStBaseOption("kStepUseLast"))
}
\arguments{
  \item{x}{ sample }
  \item{IC}{ object of class \code{"InfluenceCurve"} }
  \item{start}{ initial estimate }
  \item{useLast}{ which parameter estimate (initial estimate or
    one-step estimate) shall be used to fill the slots \code{pIC},
    \code{asvar} and \code{asbias} of the return value. }
  \item{...}{ additional arguments }
}
\details{
  Given an initial estimation \code{start}, a sample \code{x} 
  and an influence curve \code{IC} the corresponding one-step
  estimator is computed.

  In case \code{IC} is an object of class \code{"IC"}
  the slots \code{asvar} and \code{asbias} of the return 
  value are filled (based on the initial estimate).

  The default value of argument \code{useLast} is set by the
  global option \code{kStepUseLast} which by default is set to 
  \code{FALSE}. In case of general models \code{useLast} 
  remains unchanged during the computations. However, if 
  slot \code{CallL2Fam} of \code{IC} generates an object of 
  class \code{"L2GroupParamFamily"} the value of \code{useLast} 
  is changed to \code{TRUE}.
  Explicitly setting \code{useLast} to \code{TRUE} should
  be done with care as in this situation the influence curve
  is re-computed using the value of the one-step estimate
  which may take quite a long time depending on the model.

  If \code{useLast} is set to \code{TRUE} and slot \code{modifyIC} 
  of \code{IC} is filled with some function (which can be 
  used to re-compute the IC for a different parameter), the 
  computation of \code{asvar}, \code{asbias} and \code{IC} is 
  based on the one-step estimate.
}
\value{Object of class \code{"kStepEstimate"}}
\section{Methods}{
\describe{
  \item{x = "numeric", IC = "InfluenceCurve", start = "numeric"}{ 
    univariate samples. }
  \item{x = "matrix", IC = "InfluenceCurve", start = "numeric"}{ 
    multivariate samples. }
  \item{x = "matrix", IC = "InfluenceCurve", start = "Estimate"}{ 
    multivariate samples. }
  \item{x = "matrix", IC = "InfluenceCurve", start = "Estimate"}{ 
    multivariate samples. }
}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{InfluenceCurve-class}}, \code{\link{kStepEstimate-class}} }
%\examples{}
\concept{one-step estimator}
\concept{estimator}
\keyword{univar}
\keyword{robust}
