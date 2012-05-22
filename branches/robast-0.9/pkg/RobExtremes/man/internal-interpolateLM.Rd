\name{internal_interpolateLM_helpers_for_RobExtremes}
\alias{internal_interpolateLM_helpers_for_RobExtremes}
\alias{.modify.xi.PFam.call}
\alias{.RMXE.xi}
\alias{.OMSE.xi}
\alias{.MBEE.xi}
\alias{.getLMGrid}
\alias{.MakeGridList}

\title{Internal helper functions for generating interpolation grids for
       Lagrange Multipliers in package RobExtremes}

\description{
These functions are used internally to generate interpolation grids,
for Lagrange multipliers in package \pkg{RobExtremes}, to be stored in the
respective \file{sysdata.rda} file. }

\usage{
.modify.xi.PFam.call(xi, PFam)

.RMXE.xi(xi, PFam)
.MBRE.xi(xi, PFam)
.OMSE.xi(xi, PFam)

.getLMGrid(xiGrid = getShapeGrid(),
                      PFam = GParetoFamily(scale=1, shape=2),
                      optFct = .RMXE.xi,
                      withSmooth = TRUE,
                      withPrint = FALSE)

.MakeGridList(xiGrid, Y, withSmooth = TRUE)
}

\arguments{
  \item{xi}{numeric of length 1; shape value. }
  \item{PFam}{an object of class \code{"ParamFamily"}, the parametric family
              at which to evaluate the Lagrange multipliers; in our use case,
              it is a shape-scale model, hence the respective (main) parameter
              must contain \code{"scale"} and \code{"shape"}. }
  \item{withSmooth}{logical of length 1: shall a smoothing spline be used?}
  \item{withPrint}{logical of length 1: shall current shape value be printed out?}
  \item{xiGrid}{numeric; grid of shape values. }
  \item{Y}{array or matrix; contains precomputed y-values. }
  \item{optFct}{function with arguments \code{xi} and \code{PFam};
                determines the Lagrange multipliers. }
}
\details{
  \code{.MBRE.xi} computes the Lagrange multipliers for the MBRE estimator,
  \code{.OMSE.xi} for the OMSE estimator at radius \code{r=0.5},
  and \code{.RMXE.xi} the RMXE estimator.
  
  \code{.MakeGridList} transforms the return values of the preceding functions
  gathered in matrices in respective grids and can also be used for
  Lagrange multiplier matrices already computed otherwise.
  
  The actual computation can be done with \code{\link{.saveInterpGrid}}.
  
  Instead of \code{.getSnGrid} in our case this will call \code{.getLMGrid}
  to compute the Lagrange multipliers.
}
\note{These functions are only meant for the developers of package \pkg{ROptEst}
      (or respective packages).
      They can be used to speed up things by interpolation.
      Our use case is a speed up for further scale-shape families (or enhance
      existing speed-ups) such that the respective grids are stored in
      the \file{sysdata.rda} file of this package and can be used in
      (exported) new particular methods for functional \code{Sn}.}
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
  \item{.RMXE.xi}{as \code{.MBRE.xi}; in addition \code{r}, the least favorable
                  radius.
                  }
  \item{.getLMGrid}{A list with items \code{grid}, a matrix with the interpolation
                    grid, \code{fct} a function in \code{x} (the shape)
                    and \code{i} deciding on the Lagrange multiplier, and
                    \code{call}, the respective call to \code{.getLMGrid}. }
  \item{.MakeGridList}{A list with items \code{grid} and \code{fct} as
                       in the return value of \code{.getLMGrid}. }
}
\seealso{\code{\link{interpolateSn}}}
\examples{
\dontrun{
### code to produce grid for GPareto:
  .saveInterpGrid(sysRdaFolder = "C:/rtest/distr/branches/distr-2.4/pkg/distrMod/R",
              accuracy = 5000)
}
}
\keyword{internal}
\concept{utilities}
