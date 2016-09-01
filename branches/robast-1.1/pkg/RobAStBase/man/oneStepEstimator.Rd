\name{oneStepEstimator}
\alias{oneStepEstimator}

\title{Function for the computation of one-step estimates}
\description{
  Function for the computation of one-step estimates.
}
\usage{
oneStepEstimator(x, IC, start = NULL,
      useLast = getRobAStBaseOption("kStepUseLast"),
      withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
      IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
      na.rm = TRUE, startArgList = NULL, ...)
}
\arguments{
  \item{x}{ sample }
  \item{IC}{ object of class \code{"InfluenceCurve"} }
  \item{start}{ initial estimate (for full parameter,i.e. in dimension \eqn{k} respective
   joint length of main and nuisance part of the parameter):
   either a numerical value, or an object of class \code{"Estimate"} or
   a function producing either a numerical value, or an object of class \code{"Estimate"}
   when evaluated at \code{x,...}; if missing or \code{NULL}, we use slot \code{startPar}
   of the L2family \code{L2Fam} from within \code{IC}.}
  \item{useLast}{ which parameter estimate (initial estimate or
    one-step estimate) shall be used to fill the slots \code{pIC},
    \code{asvar} and \code{asbias} of the return value. }
  \item{withUpdateInKer}{if there is a non-trivial trafo in the model with matrix \eqn{D}, shall
  the parameter be updated on \eqn{{\rm ker}(D)}{ker(D)}?}
  \item{IC.UpdateInKer}{if there is a non-trivial trafo in the model with matrix \eqn{D},
  the IC to be used for this; if \code{NULL} the result of \code{getboundedIC(L2Fam,D)} is taken;
  this IC will then be projected onto \eqn{{\rm ker}(D)}{ker(D)}.}
  \item{na.rm}{logical: if  \code{TRUE}, the estimator is evaluated at \code{complete.cases(x)}.}
  \item{startArgList}{a list of arguments to be given to argument \code{start} if the latter
   is a function; this list by default already starts with two unnamed items,
   the sample \code{x}, and the model \code{eval(CallL2Fam(IC))}; in case \code{IC}
   is not of class \code{IC}, the model argument \code{L2Fam} will be set
   to \code{NULL}.}
  \item{...}{ additional arguments }
}
\details{
  Given an initial estimation \code{start}, a sample \code{x} 
  and an influence curve \code{IC} the corresponding one-step
  estimator is computed.

  In case \code{IC} is an object of class \code{"IC"}
  the slots \code{asvar} and \code{asbias} of the return 
  value are filled (based on the initial estimate).

  The default value of argument \code{useLast} is set by the
  global option \code{kStepUseLast} which by default is set to 
  \code{FALSE}. In case of general models \code{useLast} 
  remains unchanged during the computations. However, if 
  slot \code{CallL2Fam} of \code{IC} generates an object of 
  class \code{"L2GroupParamFamily"} the value of \code{useLast} 
  is changed to \code{TRUE}.
  Explicitly setting \code{useLast} to \code{TRUE} should
  be done with care as in this situation the influence curve
  is re-computed using the value of the one-step estimate
  which may take quite a long time depending on the model.

  If \code{useLast} is set to \code{TRUE} and slot \code{modifyIC} 
  of \code{IC} is filled with some function (which can be 
  used to re-compute the IC for a different parameter), the 
  computation of \code{asvar}, \code{asbias} and \code{IC} is 
  based on the one-step estimate.
}
\value{Object of class \code{"kStepEstimate"}}

\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
        Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{InfluenceCurve-class}}, \code{\link{kStepEstimate-class}} }
%\examples{}
\concept{one-step estimator}
\concept{estimator}
\keyword{univar}
\keyword{robust}
