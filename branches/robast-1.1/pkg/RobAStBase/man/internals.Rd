\name{internals_for_RobAStBase}
\alias{internals_for_RobAStBase}
\alias{.eq}
\alias{.getDistr}
\alias{.msapply}

\title{Internal / Helper functions of package RobAStBase}

\description{
These functions are used internally by package RobAStBase.}

\usage{
.eq(x,y = 0*x, tol = 1e-7)
.getDistr(L2Fam)
.evalListRec(list0)
.msapply(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE)
}
\arguments{
  \item{x}{a (numeric) vector}
  \item{y}{a (numeric) vector}
  \item{tol}{numeric --- tolerance}
  \item{L2fam}{object of class \code{L2ParamFamily}}
  \item{list0}{a list}
  \item{X}{like \code{sapply}: a vector (atomic or list) or an \code{expression} object. Other objects
  (including classed objects) will be coerced by \code{base::as.list}.}
  \item{FUN}{like \code{sapply}: the function to be applied to each element of \code{X}. See also \code{\link[base:lapply]{sapply}}}
  \item{\dots}{like \code{sapply}: optional arguments to \code{FUN}.}
  \item{simplify}{like \code{sapply}: logical or character string;
  should the result be simplified to a vector, matrix or higher dimensional array if possible? See also \code{\link[base:lapply]{sapply}}}
  \item{USE.NAMES}{like \code{sapply}: logical; if \code{TRUE} and if \code{X} is character,
  use \code{X} as names for the result unless it had names already.}
}

\details{
\code{.eq} checks equality of two vectors up to a given precision;
\code{.getDistr} produces a string with the class of the family and its parameter value;
\code{.evalListRec} recursively goes through the entries of a list, evaluating
   each entry.
\code{.msapply} like \code{base::.sapply} but catches \code{NULL}/zero-length arguments \code{X}.
}


\value{
\item{.eq}{(a vector of) \code{logical}.}
\item{.getDistr}{\code{character}.}
\item{.evalListRec}{a list.}
\item{.msapply}{as \code{base::sapply}.}
}


\author{
  Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}
  }

\keyword{internal}
\concept{utilities}
