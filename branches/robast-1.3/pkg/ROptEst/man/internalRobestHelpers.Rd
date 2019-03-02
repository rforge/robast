\name{internals_for_ROptEst_robest}
\alias{internals_for_ROptEst_robest}
\alias{.rescalefct}
\alias{.plotRescaledAxis}
\alias{.legendCoord}
\alias{.SelectOrderData}
\alias{.dynScopeEval}

\title{Internal / Helper functions of package ROptEst for function robest}

\description{
These functions are internally used helper functions for \code{\link{robest}},
in package \pkg{ROptEst}.}

\usage{
.dynScopeEval(expr)
.constructArg.list(fun,matchCall, onlyFormal=FALSE, debug =FALSE)
.fix.in.defaults(call.list, fun, withEval=TRUE)
.pretreat(x, na.rm = TRUE)
.check.eps(...)
.isOKsteps(steps)
.isOKfsCor(fsCor)
}
\arguments{
  \item{expr}{an expression. }
  \item{fun}{function, a matched call of which is manipulated. }
  \item{matchCall}{a return value of a call to \code{match.call}. }
  \item{onlyFormal}{logical; shall arguments not explicitely contained in
   the formals of \code{fun} be kept in the matched call? }
  \item{debug}{logical: if switched on, issues information for debugging. }
  \item{call.list}{a list of matched arguments drawn from a call to \code{match.call}
  applied to \code{fun} which is to be supplemented by defaults of
  not-yet-matched formals. }
  \item{withEval}{logical: shall arguments be evaluated? }
  \item{x}{input data \code{x} of \code{robest} or \code{roptest}. }
  \item{na.rm}{logical: if  \code{TRUE}, the estimator is evaluated at
   \code{complete.cases(x)}. }
  \item{\dots}{input from \code{robest} or \code{roptest} from which to conclude
     on radiuses. }
  \item{steps}{number of steps to be used in kStep estimator in \code{robest}. }
  \item{fsCor}{argument \code{fsCor} of \code{robest}. }
}
\details{
\code{.dynScopeEval} marches up the stack of calls to evaluate an expression,
   hence realizes dynamical scoping.

\code{.constructArg.list} takes a function \code{fun} and the return value
of \code{match.call} and, as return value, produces a list of arguments where
the formal arguments of \code{fun} are set to their default values and
with extra item \code{esc}.\cr

If argument \code{onlyFormals} is \code{TRUE} and the formals contain \code{\dots},
the returned list only contains formal arguments of  \code{fun}, filled with
default values from the definition where available, and, in addition,
in element \code{esc} a list with element one of the original matched call
and, as subsequent elements, with the named, evaluated arguments of the
matched call which are no formal arguments of \code{fun}.\cr

If argument \code{onlyFormals} is \code{FALSE} or the formals do not contain
\code{\dots}, the returned list again contains formal arguments of  \code{fun}
filled in with defaults where available, but in addition it contains the arguments
of the matched calls non matched to formals (in particular those passed on through
\code{\dots}). Then element \code{esc} in the returned list with contains
element one of the original matched call coerced to list, i.e., the name of
the called function.

\code{.fix.in.defaults} takes a list of arguments  of \code{fun} taken from a
matched call obtained by \code{match.call} from within a call of \code{fun}
(after coercing to list) and supplements this list by formal arguments of
\code{fun} which are not yet matched but have default arguments (with exactly
these default values). The return value is the prolongated list.

\code{.pretreat},  if \code{is.numeric(x)} is \code{FALSE}, coerces \code{x}
to a numeric matrix (by a call to \code{data.matrix} in case
\code{x} is a data.frame, respectively, by a call to \code{as.matrix} else.
If \code{na.rm} is \code{TRUE}, \code{x} is reduced to \code{na.omit(x)}.
The return value is a list of elements \code{x}, the possibly modified
input \code{x}, and \code{completecases}, the return value of
\code{compeletecases(x)}.

\code{.check.eps} takes its input (possibly empty in part)
and returns a list \code{eps} with elements \code{sqn}, \code{e},
\code{lower}, and \code{upper}. Necessarily the input \code{\dots} must
contain an argument matching to \code{x}, and \code{sqn} is the square root
of either the length of \code{x} (if \code{x} is a vector) or the number of
columns of \code{x} (in case \code{dim(x)==2}). In case \code{\dots} contains
none of the elements \code{eps}, \code{eps.lower}, \code{eps.upper}, elements
\code{lower} and \code{upper} of the return value are set to \code{0} and
\code{0.5}, respectively. Else, if \code{eps} is contained input \code{\dots}
element \code{e} of the return list is set to \code{eps}, and
\code{lower} and \code{upper} are left empty. Otherwise, element \code{e}
of the return list is left empty and  \code{lower} and \code{upper} are filled
with \code{eps.lower} and \code{eps.upper} from input \code{\dots} if available
and else with default values \code{0} and \code{0.5}, respectively.

\code{.isOKsteps} checks whether argument \code{steps} is a valid
argument, i.e., is an integer larger than 0 of length 1 and, accordingly,
returns \code{TRUE} or \code{FALSE}.

\code{.isOKfsCor} checks whether argument \code{fsCor} is a valid
argument, i.e., larger than 0 and of length 1 and, accordingly,
returns \code{TRUE} or \code{FALSE}.
}

\keyword{internal}
\concept{utilities}
\keyword{robust}
