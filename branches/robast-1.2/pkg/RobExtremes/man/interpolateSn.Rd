\name{interpolateSn}
\alias{interpolateSn}
\alias{getShapeGrid}
\alias{getSnGrid}

\title{ Function to compute LD (location-dispersion) estimates }
\description{
  Function \code{LDEstimator} provides a general way to compute
  estimates for a given parametric family of probability measures
  (with a scale and shape parameter) which
  can be obtained by matching location and dispersion functionals
  against empirical counterparts.
}
\usage{
getShapeGrid(gridsize=1000, centralvalue=0.7,
             withPos=TRUE, cutoff.at.0=1e-4, fac = 2)

getSnGrid(xiGrid = getShapeGrid(), PFam=GParetoFamily(), low=0,
                      upp=1.01, accuracy = 10000, GridFileName="SnGrid.Rdata",
                      withPrint = FALSE)
}
\arguments{
  \item{gridsize}{integer; the size of the grid to be created. }
  \item{centralvalue}{numeric of length 1: the central value of the grid
          (for details see below). }
  \item{withPos}{logical of length 1; are negative values for the shape forbidden? }
  \item{cutoff.at.0}{numeric of length 1: How close may we come to 0? }
  \item{fac}{a scaling factor used for the respective grid values (see below). }
  \item{xiGrid}{numeric; grid of shape values. }
  \item{PFam}{an object of class \code{"ParamFamily"}. The parametric
                     family at which to evaluate the LDEstimator; the respective
                     (main) parameter must contain \code{"scale"}
                     and \code{"shape"}. }
  \item{low}{numeric; argument for \code{\link{Sn}}. }
  \item{upp}{numeric; argument for \code{\link{Sn}}. }
  \item{accuracy}{numeric; argument for \code{\link{Sn}}. }
  \item{GridFileName}{character; if \code{GridFileName!=""}, the pure
  y-grid values are saved under this filename. }
  \item{withPrint}{logical of length 1: shall current shape value be printed out?}
}
\details{
  \code{getShapeGrid} is a helper function to produce an unequally spaced
  grid of shape values xi, with the rationale that we need values close
  to some typical values more often than values at the border. The code
  starts with an equally spaced grid of size \code{gridsize}
  from \code{0.5} to \code{1-0.25/gridsize}. This is reflected at \code{0.5},
  and a grid of respective quantiles of \code{Norm(mean=centralvalue, sd=fac)}
  is produced---with the heuristic rational that most estimators will be
  asymptotically normal around a typical value. If \code{withPos} is TRUE,
  negative values are cut off and replaced by respective higher quantiles of the
  corresponding normal; similarly, values to close to 0 are replaced by values
  between the cutoff value and the next admissible value and again by
  respective higher normal quantiles.
  
  \code{getSnGrid} is a helper function to produce a grid of \code{Sn} values
  for a given grid of shape values and scale equal to 1 in a given
  shape-scale family.
  This result of this function can then be used to speed
  up calls to \code{Sn} (or to \code{medSn}) by providing particular methods
  for \code{Sn}.
  For an example of such a particular method see the body of
  \code{getMethod("Sn", "GPareto")} where object
  \code{sng[["Generalized Pareto Family"]]} is just the result of a call
  \code{getSnGrid(xiGrid = getShapeGrid(), PFam=GParetoFamily())} which
  has been stored in the namespace of package \code{distrMod}.
}
\value{
 \item{getShapeGrid}{a numeric grid of xi-values. }
 \item{getSnGrid}{a grid, i.e.; a matrix with columns
  \code{xi} and \code{Sn}--the respective interpolation grid). }
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
\examples{
## (empirical) Data
getShapeGrid(50)
head(getShapeGrid(withPos=FALSE))

\dontrun{
### code used for the grid stored in the namespace of distrMod:
getSnGrid()
}
}
\keyword{univar}
