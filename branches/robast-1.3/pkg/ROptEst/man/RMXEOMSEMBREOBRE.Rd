\name{RMXEOMSEMBREOBRE}
\alias{RMXEstimator}
\alias{OMSEstimator}
\alias{OBREstimator}
\alias{MBREstimator}
\title{ Optimally robust estimation: RMXE, OMSE, MBRE, and OBRE }
\description{
  These are wrapper functions to 'roptest' to compute
  optimally robust estimates, more specifically RMXEs, OMSEs, MBREs, and OBREs,
  for L2-differentiable parametric families via k-step construction.
}
\usage{
RMXEstimator(x, L2Fam, fsCor = 1, initial.est, neighbor = ContNeighborhood(),
             steps = 1L, distance = CvMDist, startPar = NULL, verbose = NULL,
             OptOrIter = "iterate", useLast = getRobAStBaseOption("kStepUseLast"),
             withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
             IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
             withICList = getRobAStBaseOption("withICList"),
             withPICList = getRobAStBaseOption("withPICList"), na.rm = TRUE,
             initial.est.ArgList, ..., withLogScale = TRUE, ..withCheck=FALSE,
             withTimings = FALSE, withMDE = NULL, withEvalAsVar = NULL,
             withMakeIC = FALSE, modifyICwarn = NULL, E.argList = NULL,
             diagnostic = FALSE)
OMSEstimator(x, L2Fam, eps=0.5, fsCor = 1, initial.est, neighbor = ContNeighborhood(),
             steps = 1L, distance = CvMDist, startPar = NULL, verbose = NULL,
             OptOrIter = "iterate", useLast = getRobAStBaseOption("kStepUseLast"),
             withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
             IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
             withICList = getRobAStBaseOption("withICList"),
             withPICList = getRobAStBaseOption("withPICList"), na.rm = TRUE,
             initial.est.ArgList, ..., withLogScale = TRUE, ..withCheck=FALSE,
             withTimings = FALSE, withMDE = NULL, withEvalAsVar = NULL,
             withMakeIC = FALSE, modifyICwarn = NULL, E.argList = NULL,
             diagnostic = FALSE)
OBREstimator(x, L2Fam, eff=0.95, fsCor = 1, initial.est, neighbor = ContNeighborhood(),
             steps = 1L, distance = CvMDist, startPar = NULL, verbose = NULL,
             OptOrIter = "iterate", useLast = getRobAStBaseOption("kStepUseLast"),
             withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
             IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
             withICList = getRobAStBaseOption("withICList"),
             withPICList = getRobAStBaseOption("withPICList"), na.rm = TRUE,
             initial.est.ArgList, ..., withLogScale = TRUE, ..withCheck=FALSE,
             withTimings = FALSE, withMDE = NULL, withEvalAsVar = NULL,
             withMakeIC = FALSE, modifyICwarn = NULL, E.argList = NULL,
             diagnostic = FALSE)
MBREstimator(x, L2Fam, fsCor = 1, initial.est, neighbor = ContNeighborhood(),
             steps = 1L, distance = CvMDist, startPar = NULL, verbose = NULL,
             OptOrIter = "iterate", useLast = getRobAStBaseOption("kStepUseLast"),
             withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
             IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
             withICList = getRobAStBaseOption("withICList"),
             withPICList = getRobAStBaseOption("withPICList"), na.rm = TRUE,
             initial.est.ArgList, ..., withLogScale = TRUE, ..withCheck=FALSE,
             withTimings = FALSE, withMDE = NULL, withEvalAsVar = NULL,
             withMakeIC = FALSE, modifyICwarn = NULL, E.argList = NULL,
             diagnostic = FALSE)
}
\arguments{
  \item{x}{ sample }
  \item{L2Fam}{ object of class \code{"L2ParamFamily"} }
  \item{eff}{ positive real (0 <= \code{eff} <= 1): amount of asymptotic
     efficiency loss in the ideal model. See details below. }
  \item{eps}{ positive real (0 < \code{eps} <= 0.5): amount of gross errors.
          See details below. }
  \item{fsCor}{ positive real: factor used to correct the neighborhood radius;
        see details. }
  \item{initial.est}{ initial estimate for unknown parameter. If missing 
        minimum distance estimator is computed. }
  \item{neighbor}{ object of class \code{"UncondNeighborhood"} }
  \item{steps}{ positive integer: number of steps used for k-steps construction }
  \item{distance}{ distance function used in \code{MDEstimator}, which in turn
        is used as (default) starting estimator. }
  \item{startPar}{ initial information used by \code{optimize} resp. \code{optim};
    i.e; if (total) parameter is of length 1, \code{startPar} is 
    a search interval, else it is an initial parameter value; if \code{NULL}
    slot \code{startPar} of \code{ParamFamily} is used to produce it;
    in the multivariate case, \code{startPar} may also be of class \code{Estimate},
    in which case slot \code{untransformed.estimate} is used.}
  \item{verbose}{ logical: if \code{TRUE}, some messages are printed }
  \item{useLast}{ which parameter estimate (initial estimate or
    k-step estimate) shall be used to fill the slots \code{pIC},
    \code{asvar} and \code{asbias} of the return value. }
  \item{OptOrIter}{character; which method to be used for determining Lagrange
  multipliers \code{A} and \code{a}: if (partially) matched to \code{"optimize"},
  \code{getLagrangeMultByOptim} is used; otherwise: by default, or if matched to
  \code{"iterate"} or to \code{"doubleiterate"},
  \code{getLagrangeMultByIter} is used. More specifically,
  when using \code{getLagrangeMultByIter}, and if argument \code{risk} is of
  class \code{"asGRisk"}, by default and if matched to \code{"iterate"}
  we use only one (inner) iteration, if matched to \code{"doubleiterate"}
  we use up to \code{Maxiter} (inner) iterations.}
  \item{withUpdateInKer}{if there is a non-trivial trafo in the model with matrix \eqn{D}, shall
  the parameter be updated on \eqn{{\rm ker}(D)}{ker(D)}?}
  \item{IC.UpdateInKer}{if there is a non-trivial trafo in the model with matrix \eqn{D},
  the IC to be used for this; if \code{NULL} the result of \code{getboundedIC(L2Fam,D)} is taken;
  this IC will then be projected onto \eqn{{\rm ker}(D)}{ker(D)}.}
  \item{withPICList}{logical: shall slot \code{pICList} of return value
  be filled?}
  \item{withICList}{logical: shall slot \code{ICList} of return value
  be filled?}
  \item{na.rm}{logical: if  \code{TRUE}, the estimator is evaluated at \code{complete.cases(x)}.}
  \item{initial.est.ArgList}{a list of arguments to be given to argument \code{start} if the latter
   is a function; this list by default already starts with two unnamed items,
   the sample \code{x}, and the model \code{L2Fam}.}
  \item{\dots}{ further arguments }
  \item{withLogScale}{logical; shall a scale component (if existing and found
   with name \code{scalename}) be computed on log-scale and backtransformed
   afterwards? This avoids crossing 0.}
   \item{..withCheck}{logical: if \code{TRUE}, debugging info is issued.}
  \item{withTimings}{logical: if \code{TRUE}, separate (and aggregate)
   timings for the three steps evaluating the starting value, finding
   the starting influence curve, and evaluating the k-step estimator is
   issued.}
  \item{withMDE}{ logical or \code{NULL}: Shall a minimum distance estimator be used as
       starting estimator---in addition to the function given in slot
       \code{startPar} of the L2 family? If \code{NULL} (default), the content
       of slot \code{.withMDE} in the L2 family is used instead to take
       this decision.}
  \item{withEvalAsVar}{logical or \code{NULL}: if \code{TRUE} (default), tells R
  to evaluate the asymptotic variance or if \code{FALSE} just to produces a call
  to do so. If \code{withEvalAsVar} is \code{NULL} (default),  the content
       of slot \code{.withEvalAsVar} in the L2 family is used instead to take
       this decision.}
  \item{withMakeIC}{logical; if \code{TRUE} the [p]IC is passed through
        \code{makeIC} before return.}
  \item{modifyICwarn}{logical: should a (warning) information be added if
     \code{modifyIC} is applied and hence some optimality information could
     no longer be valid? Defaults to \code{NULL} in which case this value
     is taken from \code{RobAStBaseOptions}.}
  \item{E.argList}{\code{NULL} (default) or a list of arguments to be passed
                   to calls to \code{E} from (a) \code{MDEstimator}
                   (here this additional argument is only used if
                   \code{initial.est} is missing), (b) \code{getStartIC},
                   and (c) \code{kStepEstimator}. Potential clashes with
        arguments of the same name in \code{\dots} are resolved by inserting
        the items of argument list \code{E.argList} as named items, so
        in case of collisions the item of \code{E.argList} overwrites the
        existing one from \code{\dots}.}
  \item{diagnostic}{ logical; if \code{TRUE},
    diagnostic information on the performed integrations is gathered and
    shipped out as an attribute \code{diagnostic} of the return value
    of the estimators. }
}
\details{
  The functions compute optimally robust estimator for a given L2 differentiable
  parametric family; more specifically they are RMXEs, OMSEs, MBREs, and OBREs.
  The computation uses a k-step construction with an
  appropriate initial estimate; cf. also \code{\link[RobAStBase]{kStepEstimator}}. 
  Valid candidates are e.g. Kolmogorov(-Smirnov) or von Mises minimum 
  distance estimators (default); cf. Rieder (1994) and Kohl (2005).

  For OMSE, i.e., the asymptotically linear estimator with minimax mean squared
  error on this neighborhood of given size, the amount of gross errors
  (contamination) is assumed to be known, and is specified by \code{eps}.
  The radius of the corresponding infinitesimal
  contamination neighborhood is obtained by multiplying \code{eps} 
  by the square root of the sample size. 

  If the amount of gross errors (contamination) is unknown, RMXE should be used,
  i.e., the radius-minimax estimator in the sense of
  Rieder et al. (2001, 2008), respectively Section 2.2 of Kohl (2005) is returned.

  The OBRE, i.e., the optimal bias-robust (asymptotically linear) estimator;
  (terminology due to Hampel et al (1985)), expects an efficiency loss (at the
  ideal model) to be specified and then, according to an (asymptotic) Anscombe
  criterion computes the the bias bound achieving this efficiency loss.

  The MBRE, i.e., the most bias-robust (asymptotically linear) estimator;
  (terminology due to Hampel et al (1985)), uses the influence curve with
  minimal possible bias bound, hence minimaxes bias on these neighborhoods
  (in an infinitesimal sense)..

  Finite-sample and higher order results suggest that the asymptotically
  optimal procedure is to liberal. Using \code{fsCor} the radius can be
  modified - as a rule enlarged - to obtain a more conservative estimate.
  In case of normal location and scale there is function 
  \code{\link[RobLox]{finiteSampleCorrection}} which returns a finite-sample 
  corrected (enlarged) radius based on the results of large Monte-Carlo
  studies.

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

  If \code{useLast} is set to \code{TRUE} the computation of \code{asvar}, 
  \code{asbias} and \code{IC} is based on the k-step estimate.

  All these estimators are realized as wrappers to function \code{roptest}.

  Timings for the steps run through in these estimators are available
  in attributes \code{timings}, and for the step of the
  \code{kStepEstimator} in \code{kStepTimings}.

  One may also use the arguments \code{startCtrl}, \code{startICCtrl}, and
  \code{kStepCtrl} of function \code{\link{robest}}. This allows for individual
  settings of \code{E.argList}, \code{withEvalAsVar}, and
  \code{withMakeIC} for the different steps. If any of the three arguments
  \code{startCtrl}, \code{startICCtrl}, and \code{kStepCtrl} is used, the
  respective attributes set in the correspondig argument are used and, if
  colliding with arguments directly passed to the estimator function, the directly
  passed ones are ignored.

  Diagnostics on the involved integrations are available if argument
  \code{diagnostic} is \code{TRUE}. Then there are attributes \code{diagnostic}
  and \code{kStepDiagnostic} attached to the return value, which may be inspected
  and assessed through \code{\link[distrEx:distrExIntegrate]{showDiagnostic}} and
  \code{\link[distrEx:distrExIntegrate]{getDiagnostic}}.

}
\value{Object of class \code{"kStepEstimate"}.  In addition, it has
       an attribute \code{"timings"} where computation time is stored.}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Kohl, M. and Ruckdeschel, P. (2010): R package distrMod: 
  Object-Oriented Implementation of Probability Models. 
  J. Statist. Softw.  \bold{35}(10), 1--27
  
  Kohl, M. and Ruckdeschel, P., and Rieder, H. (2010):
  Infinitesimally Robust Estimation in General Smoothly Parametrized Models.
  \emph{Stat. Methods Appl.}, \bold{19}, 333--354. 

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2008) The Costs of not Knowing
  the Radius. Statistical Methods and Applications \bold{17}(1) 13-40.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2001) The Costs of not Knowing
  the Radius. Appeared as discussion paper Nr. 81. 
  SFB 373 (Quantification and Simulation of Economic Processes),
  Humboldt University, Berlin; also available under
  \url{www.uni-bayreuth.de/departments/math/org/mathe7/RIEDER/pubs/RR.pdf}
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
        Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\examples{
#############################
## 1. Binomial data
#############################
## generate a sample of contaminated data
set.seed(123)
ind <- rbinom(100, size=1, prob=0.05)
x <- rbinom(100, size=25, prob=(1-ind)*0.25 + ind*0.9)

## ML-estimate
MLE.bin <- MLEstimator(x, BinomFamily(size = 25))
## compute optimally robust estimators
OMSE.bin <- OMSEstimator(x, BinomFamily(size = 25), steps = 3)
MBRE.bin <- MBREstimator(x, BinomFamily(size = 25), steps = 3)
estimate(MLE.bin)
estimate(MBRE.bin)
estimate(OMSE.bin)

\donttest{  ## to reduce time load at CRAN tests
RMXE.bin <- RMXEstimator(x, BinomFamily(size = 25), steps = 3)
OBRE.bin <- OBREstimator(x, BinomFamily(size = 25), steps = 3)
estimate(RMXE.bin)
estimate(OBRE.bin)
}
\donttest{  ## to reduce time load at CRAN tests
#############################
## 2. Poisson data
#############################

## Example: Rutherford-Geiger (1910); cf. Feller~(1968), Section VI.7 (a)
x <- c(rep(0, 57), rep(1, 203), rep(2, 383), rep(3, 525), rep(4, 532),
       rep(5, 408), rep(6, 273), rep(7, 139), rep(8, 45), rep(9, 27),
       rep(10, 10), rep(11, 4), rep(12, 0), rep(13, 1), rep(14, 1))

## ML-estimate
MLE.pois <- MLEstimator(x, PoisFamily())
OBRE.pois <- OBREstimator(x, PoisFamily(), steps = 3)
OMSE.pois <- OMSEstimator(x, PoisFamily(), steps = 3)
MBRE.pois <- MBREstimator(x, PoisFamily(), steps = 3)
RMXE.pois <- RMXEstimator(x, PoisFamily(), steps = 3)
estimate(MLE.pois)
estimate(OBRE.pois)
estimate(RMXE.pois)
estimate(MBRE.pois)
estimate(OMSE.pois)
}

\donttest{ ## to reduce time load at CRAN tests
#############################
## 3. Normal (Gaussian) location and scale
#############################
## 24 determinations of copper in wholemeal flour
library(MASS)
data(chem)

MLE.n <- MLEstimator(chem, NormLocationScaleFamily())
MBRE.n <- MBREstimator(chem, NormLocationScaleFamily(), steps = 3)
OMSE.n <- OMSEstimator(chem, NormLocationScaleFamily(), steps = 3)
OBRE.n <- OBREstimator(chem, NormLocationScaleFamily(), steps = 3)
RMXE.n <- RMXEstimator(chem, NormLocationScaleFamily(), steps = 3)

estimate(MLE.n)
estimate(MBRE.n)
estimate(OMSE.n)
estimate(OBRE.n)
estimate(RMXE.n)
}
}

\seealso{ \code{\link{roptest}}, \code{\link{robest}},
          \code{\link[RobLox]{roblox}},
          \code{\link[distrMod]{L2ParamFamily-class}}
          \code{\link[RobAStBase]{UncondNeighborhood-class}},
          \code{\link[distrMod]{RiskType-class}} }
\concept{k-step construction}
\concept{optimally robust estimation}
\keyword{robust}
