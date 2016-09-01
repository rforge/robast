\name{biastype-methods}
\docType{methods}
\alias{biastype,interpolRisk-method}
\alias{biastype}
\title{Methods for Function biastype  in Package `RobAStBase'}
\description{biastype-methods}

\section{Methods}{\describe{

\item{biastype}{\code{signature(object = "interpolrisk")}:
  returns the slot \code{biastype} of an object of class \code{"interpolrisk"}. }
}}
\examples{
myrisk <- MBRRisk(samplesize=100)
biastype(myrisk)
}
\concept{risk}
\keyword{classes}
