\name{getStartIC-methods}
\docType{methods}
\alias{getStartIC-methods}
\alias{getStartIC}
\alias{getStartIC,ANY,ANY-method}
\alias{getStartIC,L2ParamFamily,asGRisk-method}
\alias{getStartIC,L2ParamFamily,asBias-method}
\alias{getStartIC,L2ParamFamily,asCov-method}
\alias{getStartIC,L2ParamFamily,asAnscombe-method}
\alias{getStartIC,L2ParamFamily,trAsCov-method}
\alias{getStartIC,L2LocationFamily,interpolRisk-method}
\alias{getStartIC,L2ScaleFamily,interpolRisk-method}
\alias{getStartIC,L2LocationScaleFamily,interpolRisk-method}

\title{Methods for Function getStartIC in Package `ROptEst' }

\description{\code{getStartIC} computes the optimally-robust IC to be used as
argument \code{ICstart} in \code{kStepEstimator}.}

\usage{getStartIC(model, risk, ...)
\S4method{getStartIC}{ANY,ANY}(model, risk, ...)
\S4method{getStartIC}{L2ParamFamily,asGRisk}(model, risk, ...,
                      withEvalAsVar = TRUE, withMakeIC = FALSE, ..debug=FALSE,
                      modifyICwarn = NULL, diagnostic = FALSE)
\S4method{getStartIC}{L2ParamFamily,asBias}(model, risk, ..., withMakeIC = FALSE,
        ..debug=FALSE, modifyICwarn = NULL, diagnostic = FALSE)
\S4method{getStartIC}{L2ParamFamily,asCov}(model, risk, ..., withMakeIC = FALSE,
    ..debug=FALSE)
\S4method{getStartIC}{L2ParamFamily,trAsCov}(model, risk, ..., withMakeIC = FALSE,
    ..debug=FALSE)
\S4method{getStartIC}{L2ParamFamily,asAnscombe}(model, risk, ...,
                      withEvalAsVar = TRUE, withMakeIC = FALSE, ..debug=FALSE,
                      modifyICwarn = NULL, diagnostic = FALSE)
\S4method{getStartIC}{L2LocationFamily,interpolRisk}(model, risk, ...)
\S4method{getStartIC}{L2ScaleFamily,interpolRisk}(model, risk, ...)
\S4method{getStartIC}{L2LocationScaleFamily,interpolRisk}(model, risk, ...)
}

\arguments{
  \item{model}{normtype of class \code{NormType}}
  \item{risk}{normtype of class \code{NormType}}
  \item{\dots}{ further arguments to be passed to specific methods.}
  \item{withEvalAsVar}{logical (of length 1):
     if \code{TRUE}, risks based on covariances are to be
     evaluated (default), otherwise just a call is returned.}
  \item{withMakeIC}{logical; if \code{TRUE} the IC is passed through
        \code{makeIC} before return.}
  \item{..debug}{logical; if \code{TRUE} information for debugging is issued.}
  \item{modifyICwarn}{logical: should a (warning) information be added if
     \code{modifyIC} is applied and hence some optimality information could
     no longer be valid? Defaults to \code{NULL} in which case this value
     is taken from \code{RobAStBaseOptions}.}
  \item{diagnostic}{ logical; if \code{TRUE},
    diagnostic information on the performed integrations is gathered and
    shipped out as an attribute \code{diagnostic} of the return value
    of the estimators. }
}
\section{Methods}{\describe{
\item{getStartIC}{\code{signature(model = "ANY", risk = "ANY")}:
      issue that this is not yet implemented. }
\item{getStartIC}{\code{signature(model = "L2ParamFamily", risk = "asGRisk")}:
      depending on the values of argument \code{eps} (to be passed on through
      the \code{\dots} argument) computes the optimally robust influence
      function on the fly via calls to \code{optIC} or \code{radiusMinimaxIC}. }
\item{getStartIC}{\code{signature(model = "L2ParamFamily", risk = "asBias")}:
      computes the most-bias-robust influence function on the fly via
      calls to \code{optIC}. }
\item{getStartIC}{\code{signature(model = "L2ParamFamily", risk = "asCov")}:
      computes the classically optimal influence function on the fly via
      calls to \code{optIC}. }
\item{getStartIC}{\code{signature(model = "L2ParamFamily", risk = "trAsCov")}:
      computes the classically optimal influence function on the fly via
      calls to \code{optIC}. }
}}
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

