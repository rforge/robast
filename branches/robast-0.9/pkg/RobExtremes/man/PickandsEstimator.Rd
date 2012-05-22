\name{PickandsEstimator}
\alias{PickandsEstimator}
\alias{.PickandsEstimator}

\title{ Function to compute Pickands estimates for the GPD}
\description{
  Function \code{PickandsEstimator} computes Pickands estimator
  (for the GPD) at real data and returns an object of class \code{Estimate}.
}
\usage{
PickandsEstimator(x, alpha=2, ParamFamily=GParetoFamily(),
            name, Infos, asvar = NULL, nuis.idx = NULL,
            trafo = NULL, fixed = NULL, asvar.fct  = NULL, na.rm = TRUE,
            ...)
.PickandsEstimator(x, alpha=2)
}
\arguments{
  \item{x}{ (empirical) data }
  \item{alpha}{ numeric > 1; determines the variant of the Pickands-Estimator
   based on matching the empirical \eqn{a_1=1-1/\alpha}{a1=1-1/alpha} and
   \eqn{a_1=1-1/\alpha^2}{a1=1-1/alpha^2} quantiles against the
   population counter parts. The ``classical'' Pickands Estimator is
   obtained for \code{alpha=2}. }
  \item{ParamFamily}{an object of class \code{"GParetoFamily"}. }
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
  \item{\dots}{not yet used. }
}
\details{
  The actual work is done in \code{.PickandsEstimator}.
  The wrapper \code{PickandsEstimator} pre-treats the data,
  and constructs a respective \code{Estimate} object.
}
\value{
  \item{.PickandsEstimator}{A numeric vector of length \code{2} with components
   named \code{scale} and \code{shape}. }
  \item{PickandsEstimator}{An object of S4-class \code{"Estimate"}. }
}
\references{
P. Ruckdeschel, N. Horbenko (2011): Yet another breakdown point notion:
EFSBP --illustrated at scale-shape models. ArXiv 1005.1480. To appear at Metrika.
DOI: 10.1007/s00184-011-0366-4.

}

%\references{  }
\author{Nataliya Horbenko \email{Nataliya.Horbenko@itwm.fraunhofer.de},\cr
        Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{ParamFamily-class}}, \code{\link{ParamFamily}}, 
         \code{\link{Estimate-class}} }
\examples{
## (empirical) Data
x <- rgpd(50, scale = 0.5, shape = 3)

## parametric family of probability measures
G <- GParetoFamily(scale = 1, shape = 2)

PickandsEstimator(x = x, ParamFamily = G)
}
\keyword{univar}
