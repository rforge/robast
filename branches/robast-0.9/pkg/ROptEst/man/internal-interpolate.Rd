\name{internal_interpolate_helpers}
\alias{internal_interpolate_helpers}
\alias{internalInterpolate}
\alias{.mergeF}
\alias{.copyGrid}
\alias{.renameGridName}
\alias{.MakeGridList}
\alias{.saveInterpGrid}
\alias{.recomputeInterpolators}
\alias{.versionSuff}
\alias{.getLMGrid}
\alias{.RMXE.th}
\alias{.MBRE.th}
\alias{.OMSE.th}

\title{Internal helper functions for generating interpolation grids for
       speed up in package RobExtremes}

\description{
These functions are used internally to generate interpolation grids,
for Lagrange multipliers or LDEstimators in package \pkg{RobExtremes},
to be stored in the
respective \file{sysdata.rda} file. }

\usage{
.versionSuff(name)

.RMXE.th(th, PFam, modifyfct)
.MBRE.th(th, PFam, modifyfct)
.OMSE.th(th, PFam, modifyfct)

.getLMGrid(thGrid, PFam, optFct = .RMXE.th, modifyfct,
           GridFileName="LMGrid.Rdata", withSmooth = TRUE,
           withPrint = FALSE, withCall = FALSE)

.MakeGridList(thGrid, Y, withSmooth = TRUE)

.saveInterpGrid(thGrid, PFam, sysRdaFolder, sysdataWriteFile = "sysdata.rda",
            getFun = .getLMGrid, ..., modifyfct, nameInSysdata, GridFileName,
            withSmooth = TRUE, withPrint = TRUE, withCall = FALSE, Y = NULL,
            elseFun = NULL)

.recomputeInterpolators(sysdataFiles, sysRdaFolder = ".",
                        sysdataWriteFile = "sysdata.rda", excludeGrids = NULL,
                        excludeNams = NULL, others = FALSE,
                      onlyothers = FALSE, translate = TRUE, overwrite = TRUE,
                      integrateto = FALSE, onlyCurrent = FALSE, withPrint =TRUE,
                      withSmooth = TRUE, debug = FALSE)

.renameGridName(gridnam, namOld, namNew, rdafileOld, rdafileNew)
.copyGrid(grid, gridnam, namOld, namNew, rdafileOld, rdafileNew)
.mergeF(file,envir, excludeGrids = NULL, excludeNams = NULL)
}

\arguments{
  \item{name}{Grid name to append a suffix according to the R-version. }
  \item{th}{numeric of length 1; the grid value at which to compute LMs. }
  \item{PFam}{an object of class \code{"ParamFamily"}, the parametric family
              at which to evaluate the Lagrange multipliers or LDEstimators;
              in our use case, it is a shape-scale model, hence the respective
              (main) parameter must contain \code{"scale"} and \code{"shape"}. }
  \item{modifyfct}{function with arguments \code{th} and \code{PFam} to move
       the parametric family to the point of the grid value; returns the
       moved parametric family.}
  \item{withSmooth}{logical of length 1: shall a smoothing spline be used?}
  \item{withPrint}{logical of length 1: shall current grid value be printed out?}
  \item{thGrid}{numeric; grid values. }
  \item{Y}{in case \code{.MakeGridList}: array or matrix; in case
           \code{.saveInterpGrid} array or \code{NULL}; if non-null,
           contains precomputed y-values, so that call to \code{getFun}
           resp. \code{optFct} can be omitted. }
  \item{optFct}{function with arguments \code{theta}, \code{PFam},
                and modifyfct; determines the Lagrange multipliers. }
  \item{GridFileName}{character; if \code{GridFileName!=""}, the pure
            y-grid values are saved under this filename. }
  \item{withCall}{logical of length 1: shall the call be saved, too?}
  \item{thGrid}{numeric; grid values. }
  \item{sysRdaFolder}{the folder where \pkg{RobExtremes}
       (or the respective package) is being developed; must not be missing. }
  \item{getFun}{function with first argument \code{th}, second argument
                \code{PFam} and last arguments \code{GridFileName},
                \code{withSmooth}, \code{withPrint}, and \code{withCall};
                produces the y-values for the interpolation grid. }
  \item{\dots}{further arguments to be passed on to \code{getFun}. }
  \item{nameInSysdata, nam}{name under which the list of interpolated grids is stored
                       in file \file{sysdata.rda}. }
  \item{elseFun}{function or \code{NULL}; if \code{Y} is non-null, contains
           function to transform \code{Y} to desired return value. }

  \item{sysdataFiles}{character; filenames of \file{sysdata.rda} files from
                      where to extract the interpolation grids. }
  \item{sysdataWriteFile}{filename for the \file{sysdata.rda} on which to write
         the results.}
  \item{excludeGrids}{character (or \code{NULL}); grids to be excluded
          from recomputation}
  \item{excludeNams}{character (or \code{NULL}); families to be excluded
          from recomputation}
  \item{overwrite}{logical; if \code{TRUE} foreign grids are translated
  to current R version. }
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
  \item{gridnam}{character; name of the grid to be renamed. }
  \item{namOld}{character; name of the parametric family to be renamed (from). }
  \item{namNew}{character; name of the parametric family to be renamed (to). }
  \item{grid}{matrix or array; grid to be inserted. }
  \item{rdafileOld}{character; filename of the \file{sysdata.rda}-type file
    to be read out.}
  \item{rdafileNew}{character; filename of the \file{sysdata.rda}-type file
    to be written on.}
  \item{file}{character; the name of a file to be read out}
  \item{envir}{an environment}
  \item{name}{character; name of the symbol in the \file{sysdata.rda}
     interpolation object without suffix}
  \item{xi}{numeric; shape parameter}
  \item{beta}{numeric; scale parameter}
  \item{fct}{a list of functions}
  \item{L2Fam}{ object of class \code{L2ParamFamily} }
  \item{type}{type of the IC; one of the values  \code{".OMSE"},
        \code{".RMXE"}, \code{".MBRE"}.}
}
\details{
  \code{.versionSuff}, according to the current R-version, appends a suffix
    ".O" for R<2.16 and ".N" otherwise to argument \code{name}. Needed as
    the return values of \code{splinefun} are incompatible in these two
    situations: i.e., a function with body of type
    \code{.C("R_splinefun", as.double(x),....}) respectively
    a function with body of type \code{.splinefun(....))});
    a similar case happens with \code{approxfun}.

  \code{.MBRE.th} computes the Lagrange multipliers for the MBRE estimator,
  \code{.OMSE.th} for the OMSE estimator at radius \code{r=0.5},
  and \code{.RMXE.th} the RMXE estimator.

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

   \code{.recomputeInterpolators} recomputes the interpolating functions from
     grids in existing \file{sysdata.rda} files -- either to translate them
     to another R version, or to shrink the respective file.

   \code{.renameGridName} is a utility to rename items from a grid. It takes
      grid \code{gridnam} from file \code{rdafileOld} and takes
      the name \code{namOld} of a respective item (i.e., a parametric family),
      renames it to \code{namNew} and writes the result back
      to file \code{rdafileNew}.

   \code{.copyGrid} takes out a respective item \code{namOld} (i.e., a parametric
        family) of grid \code{gridnam} from file \code{rdafileOld}
        copies it to a new grid object onto item \code{namNew}, replaces
        the respective grid-entry by \code{grid},  and saves the result to
        to file \code{rdafileNew}.

   \code{.mergeF} merges the contents of file \code{file} into environment
     \code{envir} in the sense, that if both  \code{file} and \code{envir}
     contain a list object \code{a} also the items of \code{a} are merged,
     where---as for objects themselves--- contents of \code{file} overwrite
     contents of \code{envir}.
}
\note{These functions are only meant for the developers of package
      \pkg{ROptEst} (or respective packages).
      They can be used to speed up things by interpolation.
      Our use case is a speed up for further scale-shape families (or enhance
      existing speed-ups) such that the respective grids are stored in
      a  \file{sysdata.rda} file of an external package \pkg{RobAStRda}
      ---see mail exchange P.Ruckdeschel - U.Ligges on R-devel---
      \url{https://stat.ethz.ch/pipermail/r-devel/2013-February/065794.html}.
      Special attention has to be paid for R-versions pre and post R-2.16
      which is why we use \code{.versionSuff}.
      }
\value{
  \item{.versionSuff}{A character with appended suffix. }
  \item{.MBRE.th}{A list with items \code{b} (a number; clipping height),
                  \code{a} (a 2-vector; outer centering),
                  \code{a.w} (a 2-vector; inner centering, in
                   the weight), \code{A} (a 2x2 matrix; outer standardization),
                  \code{A.w} (a 2x2 matrix; inner standardization, in the weight).
                  }
  \item{.OMSE.th}{as \code{.MBRE.th}. }
  \item{.RMXE.th}{as \code{.MBRE.th}. }
  \item{.getLMGrid}{A list with items \code{grid}, a matrix with the interpolation
                    grid and \code{fct} a function in \code{x} (the shape)
                    and \code{i} deciding on the Lagrange multiplier. }
  \item{.MakeGridList}{A list with items \code{grid} and \code{fct} as
                       in the return value of \code{.getLMGrid}. }
  \item{.saveInterpGrid}{\code{invisible(NULL)}. }
  \item{.recomputeInterpolators}{\code{invisible(NULL)}. }
  \item{.renameGridName}{\code{invisible(NULL)}. }
  \item{.mergeF}{\code{invisible(NULL)}. }
  \item{.copyGrid}{\code{invisible(NULL)}. }
}
\keyword{internal}
\concept{utilities}
