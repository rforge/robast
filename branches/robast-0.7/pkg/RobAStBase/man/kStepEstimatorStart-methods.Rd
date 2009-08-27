\name{kStepEstimator.start-methods}
\docType{methods}
\alias{kStepEstimator.start-methods}
\alias{kStepEstimator.start}
\alias{kStepEstimator.start,numeric-method}
\alias{kStepEstimator.start,Estimate-method}
\alias{kStepEstimator.start,function-method}
\title{ Methods for function kStepEstimator.start in Package `RobAStBase' }
\description{kStepEstimator.start-methods; these are called from within
\code{kStepEstimator} to produce a numeric value of for the starting estimator
in the end.}
\usage{
kStepEstimator.start(start, ...)
\S4method{kStepEstimator.start}{numeric}(start, nrvalues, ...)
\S4method{kStepEstimator.start}{Estimate}(start, nrvalues, ...)
\S4method{kStepEstimator.start}{function}(start, x, nrvalues, na.rm, L2Fam, startList)
}
\arguments{
  \item{start}{the start slot of an object of class \code{kStepEstimator}}
  \item{nrvalues}{numeric; dimension \eqn{k} of the original model, i.e.;
  length of the untransformed parameter, or joint length of
  main and nuisance part of the parameter.}
  \item{x}{the data at which the starting estimator is to be evaluated.}
  \item{na.rm}{logical: if  \code{TRUE}, the estimator is evaluated at \code{complete.cases(x)}.}
  \item{startList}{a list of arguments to be given to the call to \code{start}
  if this is a function;}
  \item{L2Fam}{the parametric famliy;}
  \item{\ldots}{further arguments for \code{kStepEstimator.start}.}
}
\value{a numeric vector with the corresponding value of the start estimator
(in \eqn{k} space)}
\section{Methods}{\describe{
\item{kStepEstimator.start}{\code{signature(start = "numeric")}:
returns the unchanged argument \code{start} if it has the correct length;
otherwise throws an error.}
\item{kStepEstimator.start}{\code{signature(start = "Estimate")}:
returns slot  \code{untransformed.estimate} of \code{start} if it is not
\code{NULL}, and else slot \code{estimate} if the latter has dimension
\code{nrvalues}.}
\item{kStepEstimator.start}{\code{signature(start = "function")}:
returns \code{kStepEstimator.start(do.call(start, args=c(list(x,L2Fam),startList)}
where, if \code{na.rm == TRUE},
beforehand \code{x} has been modified to \code{x <- complete.cases(x)}.}
}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

}
\author{
  Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}
  }
\seealso{\code{\link{kStepEstimator}},\code{\link{ALEstimate-class}}}
%\examples{}
\concept{estimate}
\keyword{classes}
