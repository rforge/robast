\name{internals_for_RobAStBase}
\alias{internals_for_RobAStBase}
\alias{.eq}
\alias{.getDistr}
\alias{.msapply}
\alias{.fixInLiesInSupport}

\title{Internal / Helper functions of package RobAStBase}

\description{
These functions are used internally by package RobAStBase.}

\usage{
.eq(x,y = 0*x, tol = 1e-7)
.getDistr(L2Fam)
.evalListRec(list0)
.msapply(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE)
.fixInLiesInSupport(IC, distr)
.filterEargsWEargList(dots)
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
  \item{IC}{an object of class \code{IC}, i.e., it expects a slot \code{Curve} like an IC.}
  \item{distr}{a distribution}
  \item{dots}{a list, obtained by \code{list(...)}.}
}

\details{
\code{.eq} checks equality of two vectors up to a given precision;
\code{.getDistr} produces a string with the class of the family and its parameter value;
\code{.evalListRec} recursively goes through the entries of a list, evaluating
   each entry.
\code{.msapply} like \code{base::.sapply} but catches \code{NULL}/zero-length arguments \code{X}.
\code{.fixInLiesInSupport} inserts a check into the function(s) in the \code{Map} slot of
the influence curve (IC), whether the arguments at which the IC is to be evaluated lie
in the support of the distribution and accordingly either returns the function value
of the IC, or \code{0}; the check is done via  calling \code{\link[distr]{liesInSupport}}.
\code{.filterEargsWEargList} calls \code{distrEx::.filterEargs} to filter out of \code{dots} 
all relevant arguments for the integrators, \code{integrate}, \code{GLIntegrate},
and \code{distrExIntegrate}; in addition, \code{.filterEargsWEargList} 
checks if an argument "E.argList" is hidden in the \code{dots} argument 
and if so, filters in its entries; in case of collisions with entries filtered
from \code{distrEx::.filterEargs}, it overwrites existing entries. In the 
end it returns a list with the filtered items.
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
