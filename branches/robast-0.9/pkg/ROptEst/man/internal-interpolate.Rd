\name{internal_interpolate_helpers}
\alias{internal_interpolate_helpers}
\alias{internalInterpolate}
\alias{.saveGridToCSV}
\alias{.readGridFromCSV}
\alias{.generateInterpGrid}
\alias{.getLMGrid}
\alias{.RMXE.th}
\alias{.MBRE.th}
\alias{.OMSE.th}

\title{Internal helper functions for generating interpolation grids for
       speed up in package ROptEst}

\description{
These functions are used internally to generate interpolation grids,
for Lagrange multipliers or LDEstimators in package \pkg{RobExtremes},
to be stored in the
respective \file{sysdata.rda} file. }

\usage{

.RMXE.th(th, PFam, modifyfct)
.MBRE.th(th, PFam, modifyfct)
.OMSE.th(th, PFam, modifyfct)

.getLMGrid(thGrid, PFam, optFct = .RMXE.th, modifyfct,
           GridFileName = "LMGrid.Rdata", withPrint = FALSE)


.saveGridToCSV(Grid, toFileCSV, namPFam, nameInSysdata)

.readGridFromCSV <- function(fromFileCSV)

.generateInterpGrid(thGrid, PFam, toFileCSV = "temp.csv",
            getFun = .getLMGrid, ..., modifyfct, nameInSysdata,
            GridFileName, withPrint = TRUE)
}

\arguments{
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
  \item{optFct}{function with arguments \code{theta}, \code{PFam},
                and modifyfct; determines the Lagrange multipliers. }
  \item{GridFileName}{character; if \code{GridFileName!=""}, the pure
            y-grid values are saved under this filename. }
  \item{Grid}{numeric; grid matrix (x- and y-values). }
  \item{toFileCSV}{character; name of the csv file to which the grid is written. }
  \item{namPFam}{character; name of the parametric family for which the grid
                 was generated. }
  \item{nameInSysdata}{character; grid name (e.g., 'OMSE', 'Sn') for which the grid
                 was generated. }
  \item{fromFileCSV}{character; name of the csv file from which the grid is read. }
  \item{getFun}{function with first argument \code{th}, second argument
                \code{PFam} and last arguments \code{GridFileName},
                \code{withPrint}; produces the y-values for the
                interpolation grid. }
  \item{\dots}{further arguments to be passed on to \code{getFun}. }
}
\details{
  \code{.MBRE.th} computes the Lagrange multipliers for the MBRE estimator,
  \code{.OMSE.th} for the OMSE estimator at radius \code{r=0.5},
  and \code{.RMXE.th} the RMXE estimator.

  \code{.getLMGrid} in a large loop computes the Lagrange multipliers for
     optimally robust IFs for each element of a given grid.
     
  \code{.saveGridToCSV} saves a given grid to a csv file, and in addition,
     in a file with same name but with file extension ".txt" writes the
     parametric family and the grid name.

  \code{.readGridFromCSV} reads in a grid from a csv file together with the
    information given in the corresponding ".txt" file.

  \code{.generateInterpGrid} by means of calls to function-argument \code{getFun}
     (e.g. \code{getLMGrid} computes the grid, if desired smoothes it, and
     then saves it to \code{.csv}.
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
  \item{.MBRE.th}{A list with items \code{b} (a number; clipping height),
                  \code{a} (a 2-vector; outer centering),
                  \code{a.w} (a 2-vector; inner centering, in
                   the weight), \code{A} (a 2x2 matrix; outer standardization),
                  \code{A.w} (a 2x2 matrix; inner standardization, in the weight).
                  }
  \item{.OMSE.th}{as \code{.MBRE.th}. }
  \item{.RMXE.th}{as \code{.MBRE.th}. }
  \item{.getLMGrid}{A grid (in form  of a matrix of x and y-values) pasted
                    together by \code{cbind}. }

  \code{.saveGridToCSV}{\code{invisible(NULL)}. }

  \code{.readGridFromCSV}{ a list with the read-in items, i.e.,
   an item \code{Grid} with the grid, an item \code{namPFam} with the name of
   the parametric family, and \code{namInSysdata}, the name of the read in grid. }

  \code{.generateInterpGrid}{\code{invisible(NULL)}. }
}
\keyword{internal}
\concept{utilities}
