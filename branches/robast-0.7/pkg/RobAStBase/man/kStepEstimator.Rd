\name{kStepEstimator}
\alias{kStepEstimator}

\title{Function for the computation of k-step estimates}
\description{
  Function for the computation of k-step estimates.
}
\usage{
kStepEstimator(x, IC, start = NULL, steps = 1L,
      useLast = getRobAStBaseOption("kStepUseLast"),
      withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
      IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
      withICList = getRobAStBaseOption("withICList"),
      withPICList = getRobAStBaseOption("withPICList"),
      na.rm = TRUE, startArgList = NULL, ...)
}
\arguments{
  \item{x}{ sample }
  \item{IC}{ object of class \code{"IC"} }
  \item{start}{ initial estimate (for full parameter,i.e. in dimension \eqn{k} respective
   joint length of main and nuisance part of the parameter):
   either a numerical value, or an object of class \code{"Estimate"} or
   a function producing either a numerical value, or an object of class \code{"Estimate"}
   when evaluated at \code{x,...}; if missing or \code{NULL}, we use slot \code{startPar}
   of the L2family \code{L2Fam} from within \code{IC}}
  \item{steps}{ integer: number of steps }
  \item{useLast}{ which parameter estimate (initial estimate or
    k-step estimate) shall be used to fill the slots \code{pIC},
    \code{asvar} and \code{asbias} of the return value. }
  \item{withUpdateInKer}{if there is a non-trivial trafo in the model with matrix \eqn{D}, shall
  the parameter be updated on \eqn{{\rm ker}(D)}{ker(D)}?}
  \item{IC.UpdateInKer}{if there is a non-trivial trafo in the model with matrix \eqn{D},
  the IC to be used for this; if \code{NULL} the result of \code{getboundedIC(L2Fam,D)} is taken;
  this IC will then be projected onto \eqn{{\rm ker}(D)}{ker(D)}.}
  \item{na.rm}{logical: if  \code{TRUE}, the estimator is evaluated at \code{complete.cases(x)}.}
  \item{startArgList}{a list of arguments to be given to argument \code{start} if the latter
   is a function; this list by default already starts with two unnamed items,
   the sample \code{x}, and the model \code{eval(CallL2Fam(IC))}.}
  \item{withPICList}{logical: shall slot \code{pICList} of return value
  be filled?}
  \item{withICList}{logical: shall slot \code{ICList} of return value
  be filled?}
  \item{...}{ additional parameters }
}
\details{
  Given an initial estimation \code{start}, a sample \code{x} 
  and an influence curve \code{IC} the corresponding k-step
  estimator is computed.

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
  based on the k-step estimate.
}
\value{Object of class \code{"kStepEstimate"}.}

\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
        Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{IC-class}}, \code{\link{kStepEstimate-class}} }
\examples{
if(require(ROptEst)){
## 1. generate a contaminated sample
ind <- rbinom(100, size=1, prob=0.05)
x <- rnorm(100, mean=0, sd=(1-ind) + ind*9)

## 2. Kolmogorov(-Smirnov) minimum distance estimator
(est0 <- MDEstimator(x=x, NormLocationScaleFamily()))

## 3. k-step estimation: radius known
N1 <- NormLocationScaleFamily(mean=estimate(est0)["mean"], sd=estimate(est0)["sd"])
N1.Rob <- InfRobModel(center = N1, neighbor = ContNeighborhood(radius = 0.5))
IC1 <- optIC(model = N1.Rob, risk = asMSE())
(est1 <- kStepEstimator(x, IC1, est0, steps = 3, withPIC = TRUE))
estimate(est1)
ksteps(est1)
pICList(est1)
start(est1)

## a transformed model
tfct <- function(x){
    nms0 <- c("mean","sd")
    nms  <- "comb"
    fval0 <- x[1]+2*x[2]
    names(fval0) <- nms
    mat0 <- matrix(c(1,2), nrow = 1, dimnames = list(nms,nms0))
    return(list(fval = fval0, mat = mat0))
}

N1.traf <- N1; trafo(N1.traf) <- tfct
N1R.traf <- N1.Rob; trafo(N1R.traf) <- tfct
IC1.traf <- optIC(model = N1R.traf, risk = asMSE())
(est0.traf <- MDEstimator(x, N1.traf))
(est1.traf <- kStepEstimator(x, IC1.traf, est0, steps = 3,
                withIC = TRUE, withPIC = TRUE, withUpdateInKer = FALSE))
(est1a.traf <- kStepEstimator(x, IC1.traf, est0, steps = 3,
                withIC = TRUE, withPIC = TRUE, withUpdateInKer = TRUE))
estimate(est1.traf)
ksteps(est1.traf)
pICList(est1.traf)
startval(est1.traf)

untransformed.estimate(est1.traf)
uksteps(est1.traf)
ICList(est1.traf)
ustartval(est1.traf)

estimate(est1a.traf)
ksteps(est1a.traf)
pICList(est1a.traf)
startval(est1a.traf)

untransformed.estimate(est1a.traf)
uksteps(est1a.traf)
ICList(est1a.traf)
ustartval(est1a.traf)
}
}
\concept{k-step estimator}
\concept{estimator}
\keyword{univar}
\keyword{robust}
