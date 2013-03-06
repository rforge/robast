\name{internal_interpolate_helpers_for_RobExtremes}
\alias{internal_interpolate_helpers_for_RobExtremes}
\alias{.modify.xi.PFam.call}
\alias{.RMXE.xi}
\alias{.OMSE.xi}
\alias{.MBEE.xi}
\alias{.getLMGrid}
\alias{.getpsi}
\alias{.svInt}

\title{Internal helper functions for generating interpolation grids for
       speed up in package RobExtremes}

\description{
These functions are used internally to generate interpolation grids,
for Lagrange multipliers or LDEstimators in package \pkg{RobExtremes},
to be stored in the
respective \file{sysdata.rda} file. }

\usage{
.getPsi(param, fct, L2Fam , type)

.modify.xi.PFam.call(xi, PFam)

.RMXE.xi(xi, PFam)
.MBRE.xi(xi, PFam)
.OMSE.xi(xi, PFam)

.getLMGrid(xiGrid = getShapeGrid(), PFam = GParetoFamily(scale=1,shape=2),
           optFct = .RMXE.xi, GridFileName="LMGrid.Rdata", withSmooth = TRUE,
           withPrint = FALSE, withCall = FALSE)

.saveInterpGrid(xiGrid = getShapeGrid(), PFam = GParetoFamily(),
            sysRdaFolder, sysdataWriteFile = "sysdata.rda", getFun = getSnGrid,
            ..., nameInSysdata = ".SnGrids", GridFileName, withSmooth = TRUE,
            withPrint = TRUE, withCall = FALSE, Y = NULL, elseFun = NULL)

.svInt(optF = .RMXE.xi, xiGrid = getShapeGrid(500, cutoff.at.0=0.005),
                   sysRdafolder, PFam = GParetoFamily(shape=1,scale=2),
                   nam = "GPD")
}

\arguments{
  \item{param}{object of class \code{"ParamFamParameter"}. }
  \item{xi}{numeric of length 1; shape value. }
  \item{PFam}{an object of class \code{"ParamFamily"}, the parametric family
              at which to evaluate the Lagrange multipliers or LDEstimators;
              in our use case, it is a shape-scale model, hence the respective
              (main) parameter must contain \code{"scale"} and \code{"shape"}. }
  \item{withSmooth}{logical of length 1: shall a smoothing spline be used?}
  \item{withPrint}{logical of length 1: shall current shape value be printed out?}
  \item{xiGrid}{numeric; grid of shape values. }
  \item{Y}{in case \code{.MakeGridList}: array or matrix; in case
           \code{.saveInterpGrid} array or \code{NULL}; if non-null,
           contains precomputed y-values, so that call to \code{getFun}
           resp. \code{optFct} can be omitted. }
  \item{optFct, optF}{function with arguments \code{xi} and \code{PFam};
                determines the Lagrange multipliers. }
  \item{GridFileName}{character; if \code{GridFileName!=""}, the pure
            y-grid values are saved under this filename. }
  \item{withCall}{logical of length 1: shall the call be saved, too?}
  \item{xiGrid}{numeric; grid of shape values. }
  \item{sysRdaFolder}{the folder where \pkg{RobExtremes}
       (or the respective package) is being developed; must not be missing. }
  \item{sysdataWriteFile}{filename for the \file{sysdata.rda} on which to write
         the results.}
  \item{getFun}{function with first argument \code{xi}, and second argument
                \code{PFam} and last arguments \code{GridFileName},
                \code{withSmooth}, \code{withPrint}, and \code{withCall};
                produces the y-values for the interpolation grid. }
  \item{\dots}{further arguments to be passed on to \code{getFun}. }
  \item{nameInSysdata, nam}{name under which the list of interpolated grids is stored
                       in file \file{sysdata.rda}. }
  \item{elseFun}{function or \code{NULL}; if \code{Y} is non-null, contains
           function to transform \code{Y} to desired return value. }

  \item{nam}{character; (short) name of the parametric family to be manipulated. }
  \item{grid}{matrix or array; grid to be inserted. }
}
\details{
   \code{.getpsi} reads the respective interpolating function
    from an object from \file{sysdata.rda} and generates a respective
   \code{HampelIC} object by a call to  \code{generateIC}.

  \code{.MBRE.xi} computes the Lagrange multipliers for the MBRE estimator,
  \code{.OMSE.xi} for the OMSE estimator at radius \code{r=0.5},
  and \code{.RMXE.xi} the RMXE estimator.
  
  \code{.saveInterpGrid} is the utility to do the actual computation.
   More specifically, the code first loads the contents of file
   \file{sysdata.rda} into an environment which is particularly created
   for this purpose (if this file exists). It then looks up whether a
   respective entry for family \code{PFam} already exists in
   \file{sysdata.rda}-object \code{nameInSysdata}.
   If this is the case the developer is asked whether he wants to overwrite
   the respective entry, and if so he does so with the results of a respective
   call to \code{getFun}.

   In case there has not been an \file{sysdata.rda}-object \code{nameInSysdata}
   so far, the code creates one and writes the results of a respective
   call to \code{getFun} to it.

   \code{.svInt} is a short form for \code{.saveInterpGrid} for LM interpolation.

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
  \item{.saveInterpGrid}{  \code{invisible(NULL)}}
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
