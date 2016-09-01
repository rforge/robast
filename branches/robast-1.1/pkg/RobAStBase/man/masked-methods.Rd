\name{masked-methods}
\docType{methods}
\alias{maskedMethods}
\alias{clip}
\alias{clip-methods}
\alias{start-methods}
\alias{clip,ANY-method}
\alias{start}
\alias{start,ANY-method}
\title{Masked Methods from Packages `stats' and `graphics' in Package `RobAStBase' }
\description{masked methods from packages \pkg{stats} and \pkg{graphics}}
\usage{
clip(x1,...)
\S4method{clip}{ANY}(x1,x2,y1,y2)
start(x,...)
\S4method{start}{ANY}(x,...)
}
\arguments{
  \item{x,...}{see \code{\link[stats]{start}}.}
  \item{x1,x2,y1,y2}{see \code{\link[graphics]{clip}}.}
}
\details{
In order to make accessible the otherwise masked functions
\code{\link[stats]{start}}, \code{\link[graphics]{clip}}, we generate
corresponding S4-methods.
}
\value{
see \code{\link[stats]{start}}, \code{\link[graphics]{clip}}
}
\author{
  Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}
  }

\keyword{methods}

