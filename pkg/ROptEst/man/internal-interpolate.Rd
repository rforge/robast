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

.RMXE.th(th, PFam, modifyfct, loRad = 0, upRad = Inf, z.start = NULL,
             A.start = NULL, upper = NULL, lower = NULL,
             OptOrIter = "iterate", maxiter = 50,
             tol = .Machine$double.eps^0.4, loRad0 = 1e-3, ...)
.MBRE.th(th, PFam, modifyfct,
             z.start = NULL, A.start = NULL, upper = 1e4,
             lower = 1e-4, OptOrIter = "iterate",
             maxiter = 50, tol = .Machine$double.eps^0.4, ...)
.OMSE.th(th, PFam, modifyfct, radius = 0.5,
             z.start = NULL, A.start = NULL, upper = 1e4,
             lower = 1e-4, OptOrIter = "iterate",
             maxiter = 50, tol = .Machine$double.eps^0.4, ...)

.getLMGrid(thGrid, PFam, optFct = .RMXE.th, modifyfct, radius = 0.5,
           GridFileName = "LMGrid.Rdata", withPrint = FALSE,
           upper = 1e4, lower = 1e-4, OptOrIter = "iterate",
           maxiter = 50, tol = .Machine$double.eps^0.4,
           loRad = 0, upRad = Inf, loRad0 = 1e-3,
           loRad.s = 0.2, up.Rad.s = 1, withStartLM = TRUE)


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
  \item{z.start}{ initial value for the centering constant. }
  \item{A.start}{ initial value for the standardizing matrix. }
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
  \item{maxiter}{ the maximum number of iterations. }
  \item{tol}{ the desired accuracy (convergence tolerance).}
  \item{loRad0}{ for numerical reasons: the effective lower bound for the zero search;
   internally set to \code{max(loRad,loRad0)}. }
  \item{\dots}{ additional parameters. }
  \item{withStartLM}{ logical of length 1: shall the LM's of the preceding grid
    value serve as starting value for the next grid value? }
  \item{withSmooth}{logical of length 1: shall a smoothing spline be used? }
  \item{withPrint}{logical of length 1: shall current grid value be printed out? }
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
  \code{.OMSE.th} for the OMSE estimator at radius \code{radius},
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
