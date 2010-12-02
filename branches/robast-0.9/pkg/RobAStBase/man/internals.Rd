\name{internals_for_RobAStBase}
\alias{internals_for_RobAStBase}
\alias{.eq}
\alias{.getDistr}

\title{Internal / Helper functions of package RobAStBase}

\description{
These functions are used internally by package RobAStBase.}

\usage{
.eq(x,y = 0*x, tol = 1e-7)
.getDistr(L2Fam)
}

\arguments{
  \item{x}{a (numeric) vector}
  \item{y}{a (numeric) vector}
  \item{tol}{numeric --- tolerance}
  \item{L2fam}{object of class \code{L2ParamFamily}}
}

\details{
\code{.eq}checks equality of two vectors up to a given precision;
\code{.getDistr} produces a string with the class of the family and its parameter value;
}


\value{
\item{.eq}{(a vector of) \code{logical}.}
\item{.getDistr}{\code{character}.}
}


\author{
  Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}
  }

\keyword{internal}
\concept{utilities}
