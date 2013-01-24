\name{internalInterpolate}
\alias{internalInterpolate}
\alias{.getpsi}

\title{Internal / Helper functions for Interpolation of ICs on a grid}

\description{
These functions are internally used helper functions
to produce an object of class \code{HampIC}
from an interpolation grid.}

\usage{
.versionSuff <- function(name)
.getPsi(xi, beta, fct, L2Fam , type)
}
\arguments{
  \item{name}{character; name of the symbol in the \file{sysdata.rda}
     interpolation object without suffix}
  \item{xi}{numeric; shape parameter}
  \item{beta}{numeric; scale parameter}
  \item{fct}{a list of functions}
  \item{L2Fam}{ object of class \code{L2ParamFamily} }
  \item{type}{type of the IC; one of the values  \code{".OMSE"},
        \code{".RMXE"}, \code{".MBRE"}.}
}
\details{
\code{.versionSuff} is a helper function to switch between different return
values of \code{\link{approxfun}} before (a function with body of type
\code{.C("R_approxfun", as.double(x),....}) and after R-2.16
(a function with body of type \code{.approxfun(x, y, v, method, yleft,
yright, f)}); a similar case happens with \code{splinefun};
in our interpolation object, there are values for
both situations with suffices \code{.O} and \code{.N} for before and
after R-2.16, respectively.
\code{.getpsi} reads the respective interpolating function
from an object from \file{sysdata.rda} and generates a respective
\code{HampelIC} object by a call to  \code{generateIC}.
}

\keyword{internal}
\concept{utilities}
\keyword{robust}
