\name{getRiskFctBV-methods}
\docType{methods}
\alias{getRiskFctBV}
\alias{getRiskFctBV-methods}
\alias{getRiskFctBV,asL1,ANY-method}
\alias{getRiskFctBV,asL4,ANY-method}
\title{Methods for Function getRiskFctBV  in Package `ROptEst'}
\description{getRiskFctBV for a given object of S4 class \code{asGRisk}
       returns a function in bias and variance to compute the asymptotic
       risk.}

\section{Methods}{\describe{
\item{getRiskFctBV}{\code{signature(risk = "asL1", biastype = "ANY")}:
  returns a function with arguments \code{bias} and \code{variance}
  to compute the asymptotic absolute (L1) error for a given ALE at a situation
  where it has bias \code{bias} (including the radius!) and variance \code{variance}. }
\item{getRiskFctBV}{\code{signature(risk = "asL4", biastype = "ANY")}:
  returns a function with arguments \code{bias} and \code{variance}
  to compute the asymptotic L4 error for a given ALE at a situation where
  it has bias \code{bias} (including the radius!) and variance \code{variance}. }
}}
\examples{
myrisk <- asMSE()
getRiskFctBV(myrisk)
}
\concept{risk}
\keyword{classes}


