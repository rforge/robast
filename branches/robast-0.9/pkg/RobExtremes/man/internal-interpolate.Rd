\name{internal_interpolate_helpers_for_RobExtremes}
\alias{internal_interpolate_helpers_for_RobExtremes}
\alias{.modify.xi.PFam.call}
\alias{.RMXE.xi}
\alias{.OMSE.xi}
\alias{.MBEE.xi}
\alias{.getLMGrid}
\alias{.MakeGridList}
\alias{.saveInterpGrid}
\alias{.recomputeInterpolators}
\alias{.svInt}

\title{Internal helper functions for generating interpolation grids for
       speed up in package RobExtremes}

\description{
These functions are used internally to generate interpolation grids,
for Lagrange multipliers or LDEstimators in package \pkg{RobExtremes},
to be stored in the
respective \file{sysdata.rda} file. }

\usage{
.modify.xi.PFam.call(xi, PFam)

.RMXE.xi(xi, PFam)
.MBRE.xi(xi, PFam)
.OMSE.xi(xi, PFam)

.getLMGrid(xiGrid = getShapeGrid(), PFam = GParetoFamily(scale=1,shape=2),
           optFct = .RMXE.xi, GridFileName="LMGrid.Rdata", withSmooth = TRUE,
           withPrint = FALSE, withCall = FALSE)

.MakeGridList(xiGrid, Y, withSmooth = TRUE)

.saveInterpGrid(xiGrid = getShapeGrid(), PFam = GParetoFamily(),
            sysRdaFolder, getFun = getSnGrid, ...,
            nameInSysdata = ".SnGrids", GridFileName, withSmooth = TRUE,
            withPrint = TRUE, withCall = FALSE, Y = NULL, elseFun = NULL)

.svInt(optF = .RMXE.xi, nam = ".RMXE",
       xiGrid = getShapeGrid(500, cutoff.at.0=0.005),
       sysRdafolder, PFam = GParetoFamily(shape=1,scale=2))

.recomputeInterpolators(sysdataFiles, sysRdaFolder = ".", others = FALSE,
                      onlyothers = FALSE, overwrite = TRUE, integrateto = FALSE,
                      onlyCurrent = FALSE, withPrint =TRUE, withSmooth = TRUE,
                      debug = FALSE)
}

\arguments{
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
  \item{getFun}{function with first argument \code{xi}, and second argument
                \code{PFam} and last arguments \code{GridFileName},
                \code{withSmooth}, \code{withPrint}, and \code{withCall};
                produces the y-values for the interpolation grid. }
  \item{\dots}{further arguments to be passed on to \code{getFun}. }
  \item{nameInSysdata, nam}{name under which the list of interpolated grids is stored
                       in file \file{sysdata.rda}. }
  \item{elseFun}{function or \code{NULL}; if \code{Y} is non-null, contains
           function to transform \code{Y} to desired return value. }

  \item{sysdataFiles}{character; filenames of \file{sysdata.rda} files from
                      where to extract the interpolation grids}
  \item{overwrite}{logical; if \code{TRUE} existing interpolation functions
      for the current R version get recomputed.}
  \item{others}{logical; if \code{TRUE} and in the \file{sysdata.rda} files
       to be modified, there are grids not ending to \code{.N} (for R>2.16)
       or \code{.O} (for R<2.16), we also recompute the interpolation
       functions for these grids. }
  \item{onlyothers}{logical; if \code{TRUE}, only the interpolation functions
  for \code{others}-grids (see argument \code{others}) are recomupted. }
 \item{onlyCurrent}{logical; if \code{TRUE} existing interpolation functions
      for the foreign R version are not included in the new \file{sysdata.rda}
      file.}
  \item{integrateto}{logical; if \code{TRUE} and there are grids with
    ending \code{.N}, \code{.O} and ``others'', these are merged. }
  \item{debug}{logical; if \code{TRUE} the \file{sysdata.rda} file is not
     created/overwritten in the end. }
}
\details{
  \code{.MBRE.xi} computes the Lagrange multipliers for the MBRE estimator,
  \code{.OMSE.xi} for the OMSE estimator at radius \code{r=0.5},
  and \code{.RMXE.xi} the RMXE estimator.
  
  \code{.MakeGridList} transforms the return values of the preceding functions
  gathered in matrices in respective grids and can also be used for
  Lagrange multiplier matrices already computed otherwise.
  
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

   \code{.recomputeInterpolators} recomputes the interpolating functions from
     grids in existing \file{sysdata.rda} files -- either to translate them
     to another R version, or to shrink the respective file.
}
\note{These functions are only meant for the developers of package \pkg{RobExtremes}
      (or respective packages).
      They can be used to speed up things by interpolation.
      Our use case is a speed up for further scale-shape families (or enhance
      existing speed-ups) such that the respective grids are stored in
      the \file{sysdata.rda} file of this package and can be used in
      (exported) new particular methods for functional \code{Sn}.
      Special attention has to be paid for R-versions pre and post R-2.16.
      So if interpolation functions are desired for both alternatives,
      
      }
\value{
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
  \item{.MakeGridList}{A list with items \code{grid} and \code{fct} as
                       in the return value of \code{.getLMGrid}. }
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
