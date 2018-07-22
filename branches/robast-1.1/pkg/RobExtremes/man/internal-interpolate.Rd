\name{internal_interpolate_helpers_for_RobExtremes}
\alias{internal_interpolate_helpers_for_RobExtremes}
\alias{.modify.xi.PFam.call}
\alias{.RMXE.xi}
\alias{.OMSE.xi}
\alias{.MBRE.xi}
\alias{.getLMGrid}
\alias{.getPsi}
\alias{.is.na.Psi}
\alias{.svInt}
\alias{.generateInterpGridSn}

\title{Internal helper functions for generating interpolation grids for
       speed up in package RobExtremes}

\description{
These functions are used internally to generate interpolation grids,
for Lagrange multipliers or LDEstimators in package \pkg{RobExtremes},
to be stored in the respective \file{sysdata.rda} file. }

\usage{
.getPsi(param, fct, L2Fam , type, withMakeIC)
.getPsi.wL(param, fct, L2Fam , type, withMakeIC)
.getPsi.P(xi, L2Fam , type, withMakeIC)


.is.na.Psi(param, fct, nam = "shape")

.modify.xi.PFam.call(xi, PFam)

.RMXE.xi(xi, PFam)
.MBRE.xi(xi, PFam)
.OMSE.xi(xi, PFam)

.getLMGrid(xiGrid = getShapeGrid(), PFam = GParetoFamily(scale=1,shape=2),
           optFct = .RMXE.xi, GridFileName="LMGrid.Rdata", withPrint = FALSE,
           len = 13)

.svInt(optF = .RMXE.th, xiGrid = getShapeGrid(700, cutoff.at.0=0.005),
       PFam = GParetoFamily(shape=1,scale=2), radius = 0.5, upper = 1e4,
       lower = 1e-4, OptOrIter = "iterate",  maxiter = 150,
       tol = .Machine$double.eps^0.5, loRad = 0, upRad = Inf, loRad0 = 1e-3,
       loRad.s = 0.2, up.Rad.s = 1, withStartLM = TRUE, len = 13, namFzus = "")

.generateInterpGridSn(xiGrid = getShapeGrid(500, cutoff.at.0=0.005),
                      PFam = GParetoFamily(), withPrint = TRUE)
}

\arguments{
  \item{param}{object of class \code{"ParamFamParameter"}. }
  \item{fct}{list of functions containing the interpolators. }
  \item{L2Fam}{an object of class \code{"L2ParamFamily"}, the parametric family
              at which to evaluate the Lagrange multipliers or LDEstimators;
              in our use case, it is a shape-scale model, hence the respective
              (main) parameter must contain \code{"scale"} and \code{"shape"}. }
  \item{nam}{character; name of the shape parameter. }
  \item{type}{type of the optimality: one of ".OMSE" for maxMSE,
       ".RMXE" for rmx, and  ".MBRE" for MBRE. }
  \item{xi}{numeric of length 1; shape value. }
  \item{PFam}{an object of class \code{"ParamFamily"}, the parametric family
              at which to evaluate the Lagrange multipliers or LDEstimators;
              in our use case, it is a shape-scale model, hence the respective
              (main) parameter must contain \code{"scale"} and \code{"shape"}. }
  \item{xiGrid}{numeric; grid of shape values. }
  \item{optFct, optF}{function with arguments \code{xi} and \code{PFam};
                determines the Lagrange multipliers. }
  \item{GridFileName}{character; if \code{GridFileName!=""}, the pure
            y-grid values are saved under this filename. }
  \item{withPrint}{logical of length 1: shall current shape value be printed out?}
  \item{radius}{ [for OMSE]: positive numeric of length 1: the radius of the
                 neighborhood for which the LM's are to be computed;
                 defaults to 0.5. }
  \item{loRad}{ the lower end point of the interval to be searched
                in the inner optimization (for the least favorable situation
                to the user-guessed radius). }
  \item{upRad}{ the upper end point of the interval to be searched in the
                 inner optimization (for the least favorable situation
                 to the user-guessed radius). }
  \item{loRad.s}{ the lower end point of the interval
                  to be searched in the outer optimization
                  (for the user-guessed radius); if \code{NULL}
                  set to \code{loRad} in the algorithm. }
  \item{upRad.s}{ the upper end point of the interval to be searched in the
                   outer optimization (for the user-guessed radius); if
                   \code{NULL} set to \code{upRad} in the algorithm. }
  \item{upper}{ upper bound for the optimal clipping bound. }
  \item{lower}{ lower bound for the optimal clipping bound. }
  \item{OptOrIter}{character; which method to be used for determining Lagrange
  multipliers \code{A} and \code{a}: if (partially) matched to \code{"optimize"},
  \code{getLagrangeMultByOptim} is used; otherwise: by default, or if matched to
  \code{"iterate"} or to \code{"doubleiterate"},
  \code{getLagrangeMultByIter} is used. More specifically,
  when using \code{getLagrangeMultByIter}, and if argument \code{risk} is of
  class \code{"asGRisk"}, by default and if matched to \code{"iterate"}
  we use only one (inner) iteration, if matched to \code{"doubleiterate"}
  we use up to \code{Maxiter} (inner) iterations. }
  \item{loRad0}{ for numerical reasons: the effective lower bound for the zero search;
   internally set to \code{max(loRad,loRad0)}. }
  \item{withStartLM}{ logical of length 1: shall the LM's of the preceding grid
    value serve as starting value for the next grid value? }
  \item{len}{integer; number of Lagrange multipliers to be calibrated. }
  \item{namFzus}{character; infix for the name of the \file{.csv}-File
                 to which the results are written; used to split the
                 work on xi-grids into chunks.}
  \item{withMakeIC}{logical; if \code{TRUE} the [p]IC is passed through
        \code{makeIC} before return.}
}
\details{
   \code{.getpsi} reads the respective interpolating function
    from an object from \file{sysdata.rda} and generates a respective
   \code{HampelIC} object by a call to  \code{generateIC}.
   \code{.getpsi.wL} does the same thing for the 3-dim model
    \code{GEVFamilyMuUnknown}.
   Last, due to scale equivariance, or the \code{ParetoFamliy}, \code{.getpsi.P}
   reads the LM's for the reference parameter and then generates
   the  respective \code{HampelIC} object by a call to  \code{generateIC}.

   \code{.is.na.Psi} checks whether the shape parameter already lies
   beyond the range for which inter-/extrapolation is admitted
   (and, correspondingly, returns \code{TRUE} if one has to compute the
   IC completely anew.).

  \code{.MBRE.xi} computes the Lagrange multipliers for the MBRE estimator,
  \code{.OMSE.xi} for the OMSE estimator at radius \code{r=0.5},
  and \code{.RMXE.xi} the RMXE estimator.
  
   \code{.svInt} is a short form for \code{ROptEst:::.generateInterpGrid}
   for LM interpolation.

}
\note{These functions are only meant for the developers of package \pkg{RobExtremes}
      (or respective packages).
      They can be used to speed up things by interpolation.
      Our use case is a speed up for further scale-shape families (or enhance
      existing speed-ups) such that the respective grids are stored in
      the \file{sysdata.rda} file of this package and can be used in
      (exported) new particular methods for functional \code{Sn}.
      Special attention has to be paid for R-versions pre and post R-2.16.
      So if interpolation functions are desired for both alternatives, one
      has to run \code{ROptEst:::.recomputeInterpolators} once on each
      version.
      }
\value{
  \item{.getpsi}{an IC. }
  \item{.is.na.Psi}{logical of length 1. }
  \item{.modify.xi.PFam.call}{A call to evaluate the parametric
            family at the new parameter value. }
  \item{.MBRE.xi}{A list with items \code{b} (a number; clipping height),
                  \code{a} (a 2-vector; outer centering),
                  \code{a.w} (a 2-vector; inner centering, in
                   the weight), \code{A} (a 2x2 matrix; outer standardization),
                  \code{A.w} (a 2x2 matrix; inner standardization, in the weight).
                  }
  \item{.OMSE.xi}{as \code{.MBRE.xi}. }
  \item{.RMXE.xi}{as \code{.MBRE.xi}. }
  \item{.getLMGrid}{A list with items \code{grid}, a matrix with the interpolation
                    grid and \code{fct} a function in \code{x} (the shape)
                    and \code{i} deciding on the Lagrange multiplier. }
  \item{.generateInterpGridSn}{  \code{invisible(NULL)}}
  \item{.svInt}{  \code{invisible(NULL)}}
}
\seealso{\code{\link{interpolateSn}}}
\keyword{internal}
\concept{utilities}
