\name{getStartIC-methods}
\docType{methods}
\alias{getStartIC-methods}
\alias{getStartIC}
\alias{getStartIC,L2ScaleShapeUnion,interpolRisk-method}
\alias{getStartIC,L2LocScaleShapeUnion,interpolRisk-method}
\alias{getStartIC,ParetoFamily,interpolRisk-method}

\title{Methods for Function getStartIC in Package `RobExtremes' }

\description{\code{getStartIC} computes the optimally-robust IC to be used as
argument \code{ICstart} in \code{kStepEstimator}.}

\usage{getStartIC(model, risk, ...)
\S4method{getStartIC}{L2ScaleShapeUnion,interpolRisk}(model, risk, ...,
   withMakeIC = FALSE, ..debug=FALSE, modifyICwarn = NULL)
\S4method{getStartIC}{L2LocScaleShapeUnion,interpolRisk}(model, risk, ...,
   withMakeIC = FALSE, ..debug=FALSE, modifyICwarn = NULL)
\S4method{getStartIC}{ParetoFamily,interpolRisk}(model, risk, ...,
   withMakeIC = FALSE)
}

\arguments{
  \item{model}{normtype of class \code{NormType}}
  \item{risk}{normtype of class \code{NormType}}
  \item{\dots}{ further arguments to be passed to specific methods.}
  \item{withMakeIC}{logical; if \code{TRUE} the IC is passed through
        \code{makeIC} before return.}
  \item{..debug}{logical; if \code{TRUE} information for debugging is issued.}
  \item{modifyICwarn}{logical: should a (warning) information be added if
     \code{modifyIC} is applied and hence some optimality information could
     no longer be valid? Defaults to \code{NULL} in which case this value
     is taken from \code{RobAStBaseOptions}.}
}
\section{Methods}{\describe{
\item{getStartIC}{\code{signature(model = "L2ScaleShapeUnion", risk = "interpolRisk")}:
      computes the optimally robust influence function by interpolation
      on a grid (using internal helper function \code{.getPsi}).}
\item{getStartIC}{\code{signature(model = "L2LocScaleShapeUnion", risk = "interpolRisk")}:
      computes the optimally robust influence function by interpolation
      on a grid (using internal helper function \code{.getPsi.wL}).}
\item{getStartIC}{\code{signature(model = "ParetoFamily", risk = "interpolRisk")}:
      computes the optimally robust influence function by interpolation
      on a grid (using internal helper function \code{.getPsi.P}).}
}
All of these methods recenter and restandardize the obtained ICs to
warrant centeredness and Fisher consistency.}
\value{
An IC of type \code{HampIC}.
}

\details{\code{getStartIC} is used internally in functions \code{robest}
and \code{roptest} to compute the optimally robust influence function
according to the arguments given to them.}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
\seealso{\code{\link{robest}},\code{\link{optIC}}, \code{\link{radiusMinimaxIC}}}
%\examples{}
\concept{asymptotic risk}
\concept{risk}
\keyword{classes}

