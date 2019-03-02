\name{internals_for_qqplot}
\alias{internals_for_qqplot}
\alias{.fadeColor}


\title{Internal functions for qqplot of package RobAStBase}

\description{
These functions are used internally by qqplot of package RobAStBase.}

\usage{
.fadeColor(col,x, bg = "white")
}


\arguments{
\item{col}{(a vector of) color(s)}
\item{x}{a (numeric) vector of fading weights}
\item{bg}{a color (of length 1): the background color against which to fade}
}

\details{
\code{.fadecolor} uses function \code{colorRamp} to interpolate between
color \code{col} and \code{bg}, at coordinate given by \code{x}.
}


\value{
\item{.fadecolor}{a vector of (rgb) colors of same length as \code{x} containing
the corresponding colors \code{col} faded towards \code{white} according
to \code{colorRamp}.}
}

\author{
  Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de},
}

\seealso{\code{\link[stats]{ks.test}}, \code{\link[distr]{qqplot}}
,\code{\link[distr:internals-qqplot]{internals_qqplot}}
,\code{\link[distrMod:internals-qqplot]{internals_qqplot}}
, \code{\link[distrMod]{qqplot}}, \code{\link[RobAStBase]{qqplot}}
}
\keyword{internal}
\concept{utilities}
