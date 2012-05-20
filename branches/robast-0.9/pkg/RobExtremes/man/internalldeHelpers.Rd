\name{internal_ldehelpers_for_distrMod}
\alias{internal_ldehelpers_for_distrMod}
\alias{.prepend}
\alias{.LDMatch}

\title{Internal helper functions for treating LDEstimators in package distrMod}

\description{
These functions are used internally by function \code{LDEstimator}
in package ``distrMod''.}

\usage{
.prepend(prep, list0, dots = NULL)
.LDMatch(x.0, loc.est.0, disp.est.0, loc.fctal.0, disp.fctal.0,
         ParamFamily.0, loc.est.ctrl.0 = NULL, loc.fctal.ctrl.0=NULL,
         disp.est.ctrl.0 = NULL, disp.fctal.ctrl.0=NULL,
         q.lo.0 =0, q.up.0=Inf, log.q.0 =TRUE, ...)
}


\arguments{
  \item{prep}{a vector; to be prepended as first argument in
              a function; named \code{x} internally. }
  \item{list0}{a list to be appended to \code{prep}; names are taken over. }
  \item{dots}{an optional list (or \code{NULL}) to be appended to \code{prep}
              and \code{list0}; names are taken over. }
  \item{x.0}{a vector (numeric) at which to evaluate the LD-match}
  \item{loc.est.0}{a function expecting \code{x} (a numeric vector)
                   as first argument; location estimator. }
  \item{disp.est.0}{a function expecting \code{x} (a numeric vector)
                    as first argument; dispersion estimator; may only take
                    non-negative values.  }
  \item{loc.fctal.0}{a function expecting a distribution object as first
                     argument; location functional. }
  \item{disp.fctal.0}{a function expecting a distribution object as first
                     argument; dispersion functional; may only take
                    non-negative values. }
  \item{loc.est.ctrl.0}{a list (or \code{NULL}); optional additional arguments
                        for the location estimator. }
  \item{disp.est.ctrl.0}{a list (or \code{NULL}); optional additional arguments
                        for the dispersion estimator.  }
  \item{loc.fctal.ctrl.0}{a list (or \code{NULL}); optional additional arguments
                        for the location functional. }
  \item{disp.fctal.ctrl.0}{a list (or \code{NULL}); optional additional arguments
                        for the dispersion functional.  }
  \item{ParamFamily}{an object of class \code{"ParamFamily"}. The parametric
                     family at which to evaluate the LDEstimator; the respective
                     (main) parameter must contain \code{"scale"}
                     and \code{"shape"}. }
  \item{q.lo.0}{numeric; lower bound for search intervall in shape parameter. }
  \item{q.up.0}{numeric; upper bound for search intervall in shape parameter.  }
  \item{log.q.0}{logical; shall the zero search be done on log-scale? }
  \item{\dots}{further arguments to be passed to location estimator and functional
  and dispersion estimator and functional. }
}

\details{
\code{.prepend} is used to produce the argument list for the calls to
 estimators and functionals. This argument list consists of \code{prep} (first
 argument), named \code{x} internally, the items of \code{list0} and, if
 non-\code{NULL}, the items of  \code{\dots}.
 
\code{.LDMatch} performs the matching of location and dispersion functional
        against empirical values (without any S4-structure).
}

\value{
\item{.prepend}{a named list to be used as arguments in a call. }
\item{.LDMatch}{a named vector with components \code{"shape"} and \code{"scale"}. }
}

\author{
  Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}
}

\seealso{
\code{\link{Estimate-class}}
}
\keyword{internal}
\concept{utilities}
