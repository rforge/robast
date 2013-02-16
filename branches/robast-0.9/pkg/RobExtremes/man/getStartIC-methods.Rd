\name{getStartIC-methods}
\docType{methods}
\alias{getStartIC-methods}
\alias{getStartIC}
\alias{getStartIC,L2ScaleShapeUnion,interpolRisk-method}

\title{Methods for Function getStartIC in Package `RobExtremes' }

\description{\code{getStartIC} computes the optimally-robust IC to be used as
argument \code{ICstart} in \code{kStepEstimator}.}

\usage{getStartIC(model, risk, ...)
\S4method{getStartIC}{L2ScaleShapeUnion,interpolRisk}(model, risk, ...)
}

\arguments{
  \item{model}{normtype of class \code{NormType}}
  \item{risk}{normtype of class \code{NormType}}
  \item{\dots}{ further arguments to be passed to specific methods.}
}
\section{Methods}{\describe{
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

