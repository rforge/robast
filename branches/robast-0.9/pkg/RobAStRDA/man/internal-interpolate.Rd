\name{internal_interpolate_helpers}
\alias{internal_interpolate_helpers}
\alias{internalInterpolate}
\alias{.versionSuff}
\alias{.MakeSmoothGridList}
\alias{.readGridFromCSV}
\alias{.generateInterpolators}
\alias{.saveGridToRda}
\alias{.mergeGrid}
\alias{.computeInterpolators}
\alias{.mergeF}
\alias{.copy_smoothGrid}
\alias{.renameGridName}

\title{Internal helper functions for generating interpolation grids for
       speed up in package RobAStRDA}

\description{
These functions are used internally to generate interpolation grids,
for Lagrange multipliers or LDEstimators in package \pkg{RobAStRDA},
to be stored in the respective \file{sysdata.rda} file. }

\usage{
.versionSuff(name)

.MakeSmoothGridList(thGrid, Y, df = NULL, gridRestrForSmooth = NULL)

.readGridFromCSV(fromFileCSV)

.generateInterpolators(Grid, approxOrspline = "spline", extrapol = c(NA,NA))

.saveGridToRda(fromFileCSV, toFileRDA = "sysdata.rda", withMerge = FALSE,
               withPrint = TRUE, withSmooth = TRUE, df = NULL,
               gridRestrForSmooth = NULL)

.mergeGrid(Grid1, Grid2)

.computeInterpolators(sysdataFiles, toFileRDA = "sysdata.rda",
      includeGrids = NULL, includeNams = NULL,
      excludeGrids = NULL, excludeNams = NULL,
      withPrint = TRUE, withSmoothFct = FALSE,
      approxOrspline = "spline", extrapol = c(NA,NA))

.mergeF(file,envir, includeGrids = NULL, includeNams = NULL,
        excludeGrids = NULL, excludeNams = NULL)

.copy_smoothGrid(gridEntry = NULL, rdafileOld, gridnamOld, FamnamOld, rdafileNew,
                 gridnamNew, FamnamNew, withSmooth = FALSE, df = NULL,
                 gridRestrForSmooth = NULL)

.renameGridName(rdafileOld, gridnamOld, FamnamOld, rdafileNew, gridnamNew,
                FamnamNew)

}

\arguments{
  \item{name}{Grid name to append a suffix according to the R-version. }
  \item{thGrid}{numeric; grid values. }
  \item{Y}{in case \code{.MakeGridList}: array or matrix; in case
           \code{.saveInterpGrid} array or \code{NULL}; if non-null,
           contains precomputed y-values, so that call to \code{getFun}
           resp. \code{optFct} can be omitted. }
  \item{df}{argument \code{df} of \code{\link{smooth.spline}}; if \code{NULL}
            (default) it is omitted (and the default of
            \code{\link{smooth.spline}} used); controls the degree to which
            we smooth; can be vectorized; to allow for \code{NULL}-entries
            in some (of the 13) LMs, it can also be a list of length 13,
            some entries being \code{NULL}, some numeric. }
  \item{gridRestrForSmooth}{an expression that can be used as index in
     \code{theta[gridRestrForSmooth]} to restrict the grid-values to
     be smoothed; the excluded grid values are left unchanged. If the argument
     is \code{NULL} no restriction is used. Can be a matrix of same dimension
     as the \code{Y}-grid to allow for column-individual restrictions,
     or a list of same length as number of columns of \code{Y}
     with columnwise restrictions of \code{Y} (and \code{NULL} entries
     are interpreted as no restriction). }
  \item{fromFileCSV}{character; names of the csv files from which the grids are
                     read. }
  \item{Grid, gridEntry}{matrix; grid to be used. }
  \item{approxOrspline}{character; if \code{approxOrspline=="spline"} (default),
    \code{\link{splinefun}} is used for generating the interpolators, otherwise
    we use \code{\link{approxfun}}. }
  \item{extrapol}{numeric of length 2; lower and upper bound, upto which
       extrapolation is done; beyond, the interpolator returns \code{NA};
       if one (or both) entries of \code{extrapol} are \code{NA}, we extrapolate
       beyond limit. }
  \item{toFileRDA}{character; the \file{.rda}-file to which the interpolators
        are saved. }
  \item{withMerge}{logical of length 1: in case a respective grid already
                 exists, shall we merge the new contents into it, or overwrite
                 the object (default)?}
  \item{withPrint}{logical of length 1: shall we issue diagnostic prints?}
  \item{withSmooth}{logical of length 1: shall a smoothing spline be used?}
  \item{Grid1, Grid2}{matrices; grids to be merged. }

  \item{sysdataFiles}{character; filenames of \file{sysdata.rda} files from
                      where to extract the interpolation grids. }
  \item{includeGrids}{character (or \code{NULL}); if non-\code{NULL}, explicit
          restriction to grids used for computation of interpolators. }
  \item{includeNams}{character (or \code{NULL}); if non-\code{NULL}, explicit
          restriction to families used for computation of interpolators. }
  \item{excludeGrids}{character (or \code{NULL}); grids to be excluded
          from computation of interpolators. }
  \item{excludeNams}{character (or \code{NULL}); families to be excluded
          from computation of interpolators. }
  \item{withSmoothFct}{logical of length 1: shall a smoothed grid be used for
                       interpolators? }
  \item{file}{character; the name of a file to be read out. }
  \item{envir}{an environment. }
  \item{gridnamOld}{character; name of the grid to be renamed (from). }
  \item{FamnamOld}{character; name of the parametric family to be renamed (from). }
  \item{rdafileOld}{character; filename of the \file{sysdata.rda}-type file
    to be read out. }
  \item{gridnamNew}{character; name of the grid to be renamed (to). If missing,
          it is set to \code{gridnamOld} internally. }
  \item{FamnamNew}{character; name of the parametric family to be renamed (to).
          If missing, it is set to \code{FamnamOld} internally. }
  \item{rdafileNew}{character; filename of the \file{sysdata.rda}-type file
    to be written on. If missing, it is set to \code{rdafileOld} internally. }
}

\details{
  \code{.versionSuff}, according to the current R-version, appends a suffix
    ".O" for R<2.16 and ".N" otherwise to argument \code{name}. Needed as
    the return values of \code{splinefun} and \code{approxfun} are incompatible
    in these two situations: i.e., a function with body of type
    \code{.C("R_splinefun", as.double(x),....}) respectively
    a function with body of type \code{.splinefun(....))});
    a similar case happens with \code{approxfun}.


  \code{.MakeSmoothGridList} smoothes out a given grid (x is contained in
       \code{thGrid}, y in \code{Y}) by \code{\link{smooth.spline}}; the degree
       of smoothing is controlled by \code{df}.

  \code{.readGridFromCSV} reads in a grid from a csv file together with the
        information given in the corresponding \file{.txt} file. More specifically,
        when generating a \file{.csv} file in e.g.\ \pkg{RobExtremes} by
        \code{.writeGridToCSV}, in addition to writing the interpolation grid
        to the \file{.csv} file, a corresponding \file{.txt} file is created
        containing information on the parametric family and the grid name.
        This information is read in again as well here.

  \code{.generateInterpolators} for a given grid according to argument
    \code{approxOrspline} produces a list of interpolators by calls to either
    \code{splinefun} or \code{approxfun} (i.e., one function per y-column).
    In addition, in case of (x-)values outside the range of the grid, in the
    interpolators, we use linear extrapolation.

  \code{.saveGridToRda} for a vector of files in argument \code{fromFileCSV}
      reads in the grids and writes them to a \file{sysdata.rda} file.
      If \file{sysdata.rda} already exists and corresponding grids already exist
      in \file{sysdata.rda}, argument \code{withMerge} controls whether the
      grids are merged or overwritten. In addition, by argument \code{withSmooth},
      one may generate smoothed out grids as well. So in the end,
      the \file{sysdata.rda} file contains grid objects with names by default
      starting with a dot, say \code{.Sn} or \code{.OMSE}. These are lists
      with parametric families as items (with the names of the parametric
      families as names of the list items). Each family item is again a list
      with (possible) entries \code{grid} (the original grid read out from
      the \file{.csv}-file), \code{gridS} (the smoothed out grid), \code{fct.O},
      the interpolator generated in R<2.16, \code{fct.N}, the interpolator
      generated in R>2.16. \code{.saveGridToRda} only generates the grid entries,
      whereas the interpolators are created in \code{.computeInterpolators}.

   \code{.mergeGrid} merges two grids according to the respective rows of the
     matrices.

   \code{.computeInterpolators} for a vector of input \file{sysdata.rda} files
      contained in argument \code{sysdataFiles} goes through the contained
      grids and produces respective interpolators by calls to
      \code{.generateInterpolators}. If \code{withSmoothFct == TRUE} (and
      a respective \code{gridS} entry exists in the current grid-family
      combination) it uses the \code{gridS} entry, otherwise the \code{grid}
      entry. By means of optional arguments \code{includeGrids},
      \code{includeNams}, \code{excludeGrids}, \code{excludeNams} one may
      explicitely restrict the grid-family combinations for which the
      interpolators are generated.

   \code{.mergeF} merges the contents of file \code{file} into environment
     \code{envir} in the sense, that if both  \code{file} and \code{envir}
     contain a list object \code{a} also the items of \code{a} are merged,
     where---as for objects themselves--- contents of \code{file} overwrite
     contents of \code{envir}. By means of optional arguments \code{includeGrids},
      \code{includeNams}, \code{excludeGrids}, \code{excludeNams} one may
      explicitely restrict the grid-family combinations for which the
      interpolators are generated.

   \code{.copy_smoothGrid} can be used to modify single \code{grid} /
       \code{gridS} entries in an existing grid-family combination in an existing
       \file{sysdata.rda}-file; in particular \code{grid} entries can be replaced
       and \code{gridS} entries can be created (or overwritten). To this end,
       \code{.copy_smoothGrid} takes out a respective item \code{FamnamOld}
       (i.e., a parametric family) of grid \code{gridnamOld} from file
       \code{rdafileOld}. If argument \code{gridEntry} is non-\code{NULL}, it
       replaces the respective \code{grid}-entry by \code{gridEntry}, otherwise
       it uses the one present in combination \code{gridnamOld}-\code{FamnamOld}.
       If argument \code{withSmooth == TRUE}, it also creates a \code{gridS}
       entry, smoothing out the entries of the \code{grid}-entry.
       The results are copies it to a new grid-family combination
       \code{gridnamNew}-\code{FamnamNew} which is saved to file
       \code{rdafileNew}.

   \code{.renameGridName} is a utility to rename grids and items from a grid.
      It takes grid \code{gridnamOld} from file \code{rdafileOld} and takes
      the name \code{FamnamOld} of a respective item (i.e., a parametric family),
      renames the grid-family combination to \code{gridnamNew},
      \code{FamnamNew} and writes the result back to file \code{rdafileNew}.

}
\note{These functions are only meant for the developers of package
      \pkg{ROptEst} (or respective packages).
      They can be used to speed up things by interpolation.
      Our use case is a speed up for further scale-shape families (or enhance
      existing speed-ups) such that the respective grids are stored in
      a  \file{sysdata.rda} file of the present package \pkg{RobAStRDA}
      ---see mail exchange P.Ruckdeschel - U.Ligges on R-devel---
      \url{https://stat.ethz.ch/pipermail/r-devel/2013-February/065794.html}.
      Special attention has to be paid for R-versions pre and post R-2.16
      which is why we use \code{.versionSuff}.
      }
\value{
  \item{.versionSuff}{A character with appended suffix. }
  \item{.MakeSmoothGridList}{the grid, i.e.; a matrix with x and y values .}
  \item{.readGridFromCSV}{ a list with the read-in items, i.e.,
   an item \code{Grid} with the grid, an item \code{namPFam} with the name of
   the parametric family, and \code{namInSysdata}, the name of the read in grid. }
  \item{.generateInterpolators}{ a list with items \code{grid} (the grid,
    i.e.; a matrix with x and y values), and \code{fct} (the list of interpolators,
    i.e.; a list of functions). }
  \item{.saveGridToRda}{\code{invisible(NULL)}. }
  \item{.mergeGrid}{ the merged grid, i.e.; a matrix. }
  \item{.computeInterpolators}{\code{invisible(NULL)}. }
  \item{.mergeF}{\code{invisible(NULL)}. }
  \item{.copy_smoothGrid}{\code{invisible(NULL)}. }
  \item{.renameGridName}{\code{invisible(NULL)}. }
}
\keyword{internal}
\concept{utilities}
