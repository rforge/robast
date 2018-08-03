\name{LDEstimator}
\alias{LDEstimator}
\alias{medkMAD}
\alias{medkMADhybr}
\alias{medSn}
\alias{medQn}

\title{ Function to compute LD (location-dispersion) estimates }
\description{
  Function \code{LDEstimator} provides a general way to compute
  estimates for a given parametric family of probability measures
  (with a scale and shape parameter) which
  can be obtained by matching location and dispersion functionals
  against empirical counterparts.
}
\usage{
LDEstimator(x, loc.est, disp.est, loc.fctal, disp.fctal, ParamFamily,
            loc.est.ctrl = NULL, loc.fctal.ctrl=NULL,
            disp.est.ctrl = NULL, disp.fctal.ctrl=NULL,
            q.lo =1e-3, q.up=15, log.q =TRUE,
            name, Infos, asvar = NULL, nuis.idx = NULL,
            trafo = NULL, fixed = NULL, asvar.fct  = NULL, na.rm = TRUE,
            ..., .withEvalAsVar = FALSE, vdbg = FALSE)
medkMAD(x, ParamFamily, k=1, q.lo =1e-3, q.up=15, nuis.idx = NULL,
        trafo = NULL, fixed = NULL, asvar.fct = NULL, na.rm = TRUE,
        ..., .withEvalAsVar = FALSE, vdbg = FALSE)
medkMADhybr(x, ParamFamily, k=1, q.lo =1e-3, q.up=15, KK = 20, nuis.idx = NULL,
        trafo = NULL, fixed = NULL,  asvar.fct = NULL, na.rm = TRUE,
        ..., .withEvalAsVar = FALSE)
medSn(x, ParamFamily, q.lo =1e-3, q.up=10, nuis.idx = NULL,
      trafo = NULL, fixed = NULL, asvar.fct = NULL, na.rm = TRUE,
      accuracy = 100, ..., .withEvalAsVar = FALSE)
medQn(x, ParamFamily, q.lo =1e-3, q.up=15, nuis.idx = NULL,
      trafo = NULL, fixed = NULL, asvar.fct = NULL, na.rm = TRUE,
      ..., .withEvalAsVar = FALSE)
}
\arguments{
  \item{x}{ (empirical) data }
  \item{ParamFamily}{an object of class \code{"ParamFamily"}. The parametric
                     family at which to evaluate the LDEstimator; the respective
                     (main) parameter must contain \code{"scale"}
                     and \code{"shape"}. }
  \item{loc.est}{a function expecting \code{x} (a numeric vector)
                   as first argument; location estimator. }
  \item{disp.est}{a function expecting \code{x} (a numeric vector)
                    as first argument; dispersion estimator; may only take
                    non-negative values.  }
  \item{loc.fctal}{a function expecting a distribution object as first
                     argument; location functional. }
  \item{disp.fctal}{a function expecting a distribution object as first
                     argument; dispersion functional; may only take
                    non-negative values. }
  \item{loc.est.ctrl}{a list (or \code{NULL}); optional additional arguments
                        for the location estimator. }
  \item{disp.est.ctrl}{a list (or \code{NULL}); optional additional arguments
                        for the dispersion estimator.  }
  \item{loc.fctal.ctrl}{a list (or \code{NULL}); optional additional arguments
                        for the location functional. }
  \item{disp.fctal.ctrl}{a list (or \code{NULL}); optional additional arguments
                        for the dispersion functional.  }
  \item{k}{numeric; additional parameter for \code{\link{kMAD}}; must be positive
           and of length 1.}
  \item{KK}{numeric; Maximal number of trials with different \code{k} in
   \code{medkMADhybr} . }
  \item{q.lo}{numeric; lower bound for search intervall in shape parameter. }
  \item{q.up}{numeric; upper bound for search intervall in shape parameter.  }
  \item{log.q}{logical; shall the zero search be done on log-scale? }
  \item{name}{ optional name for estimator. }
  \item{Infos}{ character: optional informations about estimator }
  \item{asvar}{ optionally the asymptotic (co)variance of the estimator }
  \item{nuis.idx}{ optionally the indices of the estimate belonging
                  to nuisance parameter}
  \item{fixed}{ optionally (numeric) the fixed part of the parameter}
  \item{trafo}{ an object of class \code{MatrixorFunction} -- a transformation
  for the main parameter}
  \item{asvar.fct}{optionally: a function to determine the corresponding
    asymptotic variance; if given, \code{asvar.fct} takes arguments
    \code{L2Fam}((the parametric model as object of class \code{L2ParamFamily})) 
    and \code{param} (the parameter value as object of class 
    \code{ParamFamParameter}); arguments are called by name; \code{asvar.fct}
     may also process further arguments passed through the \code{\dots} argument}              
  \item{na.rm}{logical: if  \code{TRUE}, the estimator is evaluated at \code{complete.cases(x)}.}
  \item{accuracy}{numeric: argument to be passed on to \code{\link{Sn}}. }
  \item{\dots}{further arguments to be passed to location estimator and functional
  and dispersion estimator and functional. }
  \item{vdbg}{logical; if \code{TRUE}, debugging information is shown.}
  \item{.withEvalAsVar}{logical: shall slot \code{asVar} be evaluated
                   (if \code{asvar.fct} is given) or
                   just the call be returned?}
}
\details{
  The arguments \code{loc.est}, \code{disp.est} (location and dispersion estimators)
  have to be functions with first argument \code{x} (a numeric vector with the
  empirical data) and additional, optional individual arguments to be passed on
  in the respective calls as lists \code{loc.est.ctrl}, \code{disp.est.ctrl},
  and global additional arguments through the \code{\dots} argument.
  Similarly, arguments \code{loc.fctal}, \code{disp.fctal} (location and
  dispersion functionals) have to be functions with first argument an
  object of class \code{UnivariateDistribution}, and additional, optional
  individual arguments to be passed on
  in the respective calls as lists \code{loc.fctal.ctrl}, \code{disp.fctal.ctrl},
  and global additional arguments again through the \code{\dots} argument.
  Uses \code{\link{.LDMatch}} internally.
}
\note{The values for \code{q.lo} and \code{q.up} are a bit delicate and
 have to be found, model by model, by try and error.
 As a rule, \code{medSn} is rather slow, as the evaluation of the \code{Sn}
 functional is quite expensive. So if \code{medSn} is the estimator of choice,
 it pays off, for a given shape-scale family, to evaluate \code{medSn} on a
 grid of shape-values (with scale 1) and then to use an interpolation techniques
 in a particular method to replace the default one for this shape-scale family.
 As an example, we have done so for the GPD family.}
\value{
  An object of S4-class \code{"Estimate"}.
}
\references{
P. Ruckdeschel, N. Horbenko (2011): Yet another breakdown point notion:
EFSBP --illustrated at scale-shape models. ArXiv 1005.1480. To appear at Metrika.
DOI: 10.1007/s00184-011-0366-4.

A. Marazzi and C. Ruffieux (1999): The truncated mean of asymmetric distribution.
Computational Statistics and Data Analysis 32: 79-100 .
}

%\references{  }
\author{Nataliya Horbenko \email{nhorbenko@gmail.com},\cr
        Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{ParamFamily-class}}, \code{\link{ParamFamily}}, 
         \code{\link{Estimate-class}} }
\examples{
## (empirical) Data
set.seed(123)
x <- rgamma(50, scale = 0.5, shape = 3)

## parametric family of probability measures
G <- GammaFamily(scale = 1, shape = 2)

medQn(x = x, ParamFamily = G)
medSn(x = x, ParamFamily = G, q.lo = 0.5, q.up = 4)

\donttest{
## not tested on CRAN because it takes time...
## without speedup for Sn:
LDEstimator(x, loc.est = median, disp.est = Sn, loc.fctal = median,
            disp.fctal = getMethod("Sn","UnivariateDistribution"),
            ParamFamily = G, disp.est.ctrl = list(constant=1))

medkMAD(x = x, ParamFamily = G)
medkMADhybr(x = x, ParamFamily = G)
}
medkMAD(x = x, k=10, ParamFamily = G)

##not at all robust:
LDEstimator(x, loc.est = mean, disp.est = sd,
               loc.fctal = E, disp.fctal = sd,
            ParamFamily = G)
}
\keyword{univar}
