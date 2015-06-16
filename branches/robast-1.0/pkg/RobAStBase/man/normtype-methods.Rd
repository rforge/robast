\name{normtype-methods}
\docType{methods}
\alias{normtype,interpolRisk-method}
\alias{normtype}
\title{Methods for Function normtype  in Package `RobAStBase'}
\description{normtype-methods}

\section{Methods}{\describe{

\item{normtype}{\code{signature(object = "interpolrisk")}:
  returns the slot \code{normtype} of an object of class \code{"interpolrisk"}. }
}}
\examples{
myrisk <- MBRRisk(samplesize=100)
normtype(myrisk)
}
\concept{risk}
\keyword{classes}
