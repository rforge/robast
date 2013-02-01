\name{getStartIC-methods}
\docType{methods}
\alias{getStartIC-methods}
\alias{getStartIC}
\alias{getStartIC,ANY,ANY-method}
\alias{getStartIC,L2ParamFamily,asGRisk-method}
\alias{getStartIC,L2ParamFamily,asBias-method}
\alias{getStartIC,L2ParamFamily,asCov-method}
\alias{getStartIC,L2ParamFamily,trAsCov-method}
\alias{getStartIC,L2ScaleShapeUnion,interpolRisk-method}

\title{Methods for Function getStartIC in Package `ROptEst' }

\description{\code{getStartIC} computes the optimally-robust IC to be used as
argument \code{ICstart} in \code{kStepEstimator}.}

\usage{getStartIC(model, risk, ...)
\S4method{getStartIC}{L2ParamFamily,asRisk}(model, risk, ..., ..debug=FALSE)
}

\arguments{
  \item{model}{normtype of class \code{NormType}}
  \item{risk}{normtype of class \code{NormType}}
  \item{\dots}{ further arguments to be passed to specific methods.}
  \item{..debug}{logical; if \code{TRUE} information for debugging is issued.}
}
\section{Methods}{\describe{
\item{getStartIC}{\code{signature(model = "ANY", risk = "ANY")}:
      issue that this is not yet implemented. }
\item{getStartIC}{\code{signature(model = "L2ParamFamily", risk = "asRisk")}:
      depending on the values of argument \code{eps} (to be passed on through
      the \code{\dots} argument) computes the optimally robust influence
      function on the fly via calls to \code{optIC} or \code{radiusMinimaxIC}. }
\item{getStartIC}{\code{signature(model = "L2ScaleShapeUnion", risk = "interpolRisk")}:
      computes the optimally robust influence function by interpolation
      on a grid (using helper function \code{.getPsi}).}

}}
\value{
An IC of type \code{HampIC}.
}

\details{\code{getStartIC} is used internally in functions \code{robest}
and \code{roptest} to compute the optimally robust influence function
according to the arguments given to them.}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
\seealso{\code{\link{robest}},\code{\link{optIC}}, \code{\link{radiusMinimaxIC}}}
%\examples{}
\concept{asymptotic risk}
\concept{risk}
\keyword{classes}

