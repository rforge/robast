\name{internal_interpolate_helpers_for_RobExtremes}
\alias{internal_interpolate_helpers_for_RobExtremes}
\alias{.modify.xi.PFam.call}
\alias{.RMXE.xi}
\alias{.OMSE.xi}
\alias{.MBRE.xi}
\alias{.getLMGrid}
\alias{.getPsi}
\alias{.svInt}
\alias{.generateInterpGridSn}

\title{Internal helper functions for generating interpolation grids for
       speed up in package RobExtremes}

\description{
These functions are used internally to generate interpolation grids,
for Lagrange multipliers or LDEstimators in package \pkg{RobExtremes},
to be stored in the respective \file{sysdata.rda} file. }

\usage{
.getPsi(param, fct, L2Fam , type)

.modify.xi.PFam.call(xi, PFam)

.RMXE.xi(xi, PFam)
.MBRE.xi(xi, PFam)
.OMSE.xi(xi, PFam)

.getLMGrid(xiGrid = getShapeGrid(), PFam = GParetoFamily(scale=1,shape=2),
           optFct = .RMXE.xi, GridFileName="LMGrid.Rdata", withPrint = FALSE)

.svInt(optF = .RMXE.xi, xiGrid = getShapeGrid(500, cutoff.at.0=0.005),
       PFam = GParetoFamily(shape=1,scale=2))

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
}
\details{
   \code{.getpsi} reads the respective interpolating function
    from an object from \file{sysdata.rda} and generates a respective
   \code{HampelIC} object by a call to  \code{generateIC}.

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
  \item{.getpsi}{an IC.}
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
\examples{
\dontrun{
### code to produce grid for GPareto:
  .saveInterpGrid(sysRdaFolder = "C:/rtest/RobASt/branches/RobASt-0.9/pkg/RobExtremes/R",
              accuracy = 800)
}
}
\keyword{internal}
\concept{utilities}
