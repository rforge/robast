\name{samplesize-methods}
\docType{methods}
\alias{samplesize,interpolRisk-method}
\alias{samplesize}
\alias{samplesize<-,interpolRisk-method}
\alias{samplesize<-}
\title{Methods for Function samplesize  in Package `RobAStBase'}
\description{samplesize-methods}

\section{Methods}{\describe{

\item{samplesize}{\code{signature(object = "interpolrisk")}:
  returns the slot \code{samplesize} of an object of class \code{"interpolrisk"}. }
\item{samplesize<-}{\code{signature(object = "interpolrisk", value = "ANY")}:
  modifies the slot \code{samplesize} of an object of class \code{"interpolrisk"}. }
}}
\examples{
myrisk <- MBRRisk(samplesize=100)
samplesize(myrisk)
samplesize(myrisk) <- 20
}
\concept{risk}
\keyword{classes}
