\name{internals_for_MBRE_computation}
\alias{internals_for_MBRE_computation}
\alias{.getExtremeCoordIC}

\title{Internal / Helper function of package ROptEst for MBRE calculation}

\description{This function computes the coordinatewise min and max of an IC
numerically.}

\usage{
.getExtremeCoordIC(IC, D, indi, n = 10000)
}
\arguments{
  \item{IC}{object of class \code{IC}}
  \item{D}{a univariate distribution; by means of \code{q(D)},
  we determine the evaluation grid.}
  \item{indi}{integer; the indices of the coordinates at which to compute min
  and max.}
  \item{n}{integer; number of grid points for evaluation.}
}
\value{a matrix with \code{length(indi)} rows and 2 columns \code{min} and \code{max}:
the coordinate-wise min and max of the IC.}
\keyword{internal}
\concept{utilities}
\keyword{robust}
