\name{get.asGRisk.fct-methods}
\docType{methods}
\alias{get.asGRisk.fct-methods}
\alias{get.asGRisk.fct}
\alias{get.asGRisk.fct,asL1-method}
\alias{get.asGRisk.fct,asL4-method}
\alias{get.asGRisk.fct,asMSE-method}
\title{Methods for Function get.asGRisk.fct in Package `ROptEst' }

\description{get.asGRisk.fct-methods to produce a function in r,s,b for computing
 a particular asGRisk}

\usage{get.asGRisk.fct(Risk)
\S4method{get.asGRisk.fct}{asMSE}(Risk)     
\S4method{get.asGRisk.fct}{asL1}(Risk)     
\S4method{get.asGRisk.fct}{asL4}(Risk)     
}

\arguments{
  \item{Risk}{a risk of class \code{"asGRisk"}}
}
\section{Methods}{\describe{
\item{get.asGRisk.fct}{\code{signature(Risk = "asMSE")}: method for asymptotic mean squared error.}
\item{get.asGRisk.fct}{\code{signature(Risk = "asL1")}: method for asymptotic mean absolute error.}
\item{get.asGRisk.fct}{\code{signature(Risk = "asL4")}: method for asymptotic mean power 4 error.}
}}
\value{
\item{get.asGRisk.fct}{a function with arguments \code{r} (radius), \code{s} 
(square root of (trace of) variance), \code{b} bias to compute the respective risk of an IC
with this bias and variance at the respective radius.
}
}

\details{\code{get.asGRisk.fct} is used internally in functions \code{\link{getAsRisk}}
and \code{\link{getReq}}.}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\examples{}
\concept{asymptotic risk}
\concept{risk}
\keyword{classes}

