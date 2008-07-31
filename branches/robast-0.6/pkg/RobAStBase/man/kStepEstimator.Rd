\name{kStepEstimator}
\alias{kStepEstimator}
\alias{kStepEstimator-methods}
\alias{kStepEstimator,numeric,IC,numeric-method}
\alias{kStepEstimator,matrix,IC,numeric-method}
\alias{kStepEstimator,numeric,IC,Estimate-method}
\alias{kStepEstimator,matrix,IC,Estimate-method}

\title{Generic function for the computation of k-step estimates}
\description{
  Generic function for the computation of k-step estimates.
}
\usage{
kStepEstimator(x, IC, start, ...)

\S4method{kStepEstimator}{numeric,IC,numeric}(x, IC, start, steps = 1L, 
    useLast = getRobAStBaseOption("kStepUseLast"))
\S4method{kStepEstimator}{matrix,IC,numeric}(x, IC, start, steps = 1L, 
    useLast = getRobAStBaseOption("kStepUseLast"))
\S4method{kStepEstimator}{numeric,IC,Estimate}(x, IC, start, steps = 1L, 
    useLast = getRobAStBaseOption("kStepUseLast"))
\S4method{kStepEstimator}{matrix,IC,Estimate}(x, IC, start, steps = 1L, 
    useLast = getRobAStBaseOption("kStepUseLast"))
}
\arguments{
  \item{x}{ sample }
  \item{IC}{ object of class \code{"IC"} }
  \item{start}{ initial estimate }
  \item{steps}{ integer: number of steps }
  \item{useLast}{ which parameter estimate (initial estimate or
    k-step estimate) shall be used to fill the slots \code{pIC},
    \code{asvar} and \code{asbias} of the return value. }
  \item{...}{ additional parameters }
}
\details{
  Given an initial estimation \code{start}, a sample \code{x} 
  and an influence curve \code{IC} the corresponding k-step
  estimator is computed.

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
  based on the k-step estimate.
}
\value{Object of class \code{"kStepEstimate"}.}
\section{Methods}{
\describe{
  \item{x = "numeric", IC = "IC", start = "numeric"}{ 
    univariate samples. }
  \item{x = "matrix", IC = "IC", start = "numeric"}{ 
    multivariate samples. }
  \item{x = "matrix", IC = "IC", start = "Estimate"}{ 
    multivariate samples. }
  \item{x = "matrix", IC = "IC", start = "Estimate"}{ 
    multivariate samples. }
}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{IC-class}}, \code{\link{kStepEstimate-class}} }
%\examples{}
\concept{k-step estimator}
\concept{estimator}
\keyword{univar}
\keyword{robust}
