\name{ismevgpdgevdiag-methods}
\docType{methods}
\alias{gpd.diag}
\alias{gev.diag}
\alias{gpd.prof}
\alias{gev.prof}
\alias{gpd.profxi}
\alias{gev.profxi}
\alias{gpd.diag-methods}
\alias{gev.diag-methods}
\alias{gpd.profxi-methods}
\alias{gev.profxi-methods}
\alias{gpd.diag,GPDEstimate-method}
\alias{gpd.diag,gpd.fit-method}
\alias{gev.diag,GEVEstimate-method}
\alias{gev.diag,gev.fit-method}
\alias{gpd.prof,GPDEstimate-method}
\alias{gpd.prof,gpd.fit-method}
\alias{gev.prof,GEVEstimate-method}
\alias{gev.prof,gev.fit-method}
\alias{gpd.profxi,GPDEstimate-method}
\alias{gpd.profxi,gpd.fit-method}
\alias{gev.profxi,GEVEstimate-method}
\alias{gev.profxi,gev.fit-method}
\title{ Methods for Diagnostic Functions in Package `RobExtremes' }
\description{We provide wrapper to the diagnostic plots
\code{gpd.diag} and \code{gev.diag} of package \pkg{ismev},
as well as to profilers \code{gpd.prof}, \code{gpd.profxi} and \code{gev.prof},
 \code{gev.profxi}. }

\usage{
gpd.diag(z,...)
\S4method{gpd.diag}{gpd.fit}(z)
\S4method{gpd.diag}{GPDEstimate}(z, npy = 365)
gev.diag(z)
\S4method{gev.diag}{gev.fit}(z)
\S4method{gev.diag}{GEVEstimate}(z)
gpd.prof(z,...)
\S4method{gpd.prof}{gpd.fit}(z, m, xlow, xup, npy = 365, conf = 0.95, nint = 100)
\S4method{gpd.prof}{GPDEstimate}(z, m, xlow, xup, npy = 365, conf = 0.95, nint = 100)
gev.prof(z,...)
\S4method{gev.prof}{gev.fit}(z, m, xlow, xup, conf = 0.95, nint = 100)
\S4method{gev.prof}{GEVEstimate}(z, m, xlow, xup, conf = 0.95, nint = 100)
gpd.profxi(z,...)
\S4method{gpd.profxi}{gpd.fit}(z,  xlow, xup, conf = 0.95, nint = 100)
\S4method{gpd.profxi}{GPDEstimate}(z,  xlow, xup, npy = 365, conf = 0.95, nint = 100)
gev.profxi(z,...)
\S4method{gev.profxi}{gev.fit}(z, xlow, xup, conf = 0.95, nint = 100)
\S4method{gev.profxi}{GEVEstimate}(z, xlow, xup, conf = 0.95, nint = 100)
}

\arguments{
  \item{z}{an argument of class \code{gpd.fit}, \code{gev.fit} (recovering
   the original calling convention from package \pkg{ismev} or of class
   \code{GEVFamily} or \code{GParetoFamily}.}

\item{m}{The return level (i.e.\ the profile likelihood is for the
          value that is exceeded with probability \code{1/m}). }

\item{\dots}{further parameters to be passed on the specific methods. }

\item{xlow, xup}{The least and greatest value at which to evaluate the
          profile likelihood. }

\item{npy}{The number of observations per year. }

\item{conf}{The confidence coefficient of the plotted profile confidence
          interval. }

\item{nint}{The number of points at which the profile likelihood is
          evaluated. }

  }
\details{
We provide a coercing of our fits of S4-classes \code{"GPDEstimate"}
and \code{"GEVEstimate"} to the (S3-)classes \code{gpd.fit} and \code{gev.fit}
of package \pkg{ismev} (the latter being cast to an S4 class, internally, in
our package.
}

\value{
For \code{gpd.fit}, \code{gev.fit}
(quoted from package \pkg{ismev}:
     For stationary models four plots are produced; a probability plot,
     a quantile plot, a return level plot and a histogram of data with
     fitted density.

     For non-stationary models two plots are produced; a residual
     probability plot and a residual quantile plot.

For \code{gpd.prof}, \code{gev.prof}
(quoted from package \pkg{ismev}:

     A plot of the profile likelihood is produced, with a horizontal
     line representing a profile confidence interval with confidence
     coefficient \code{conf}.

}
\author{
  Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}
  }

\references{
  ismev: An   Introduction to Statistical Modeling of Extreme Values. R package
  version 1.39. https://CRAN.R-project.org/package=ismev; original S functions
  written by Janet E. Heffernan with R port and R documentation provided by
  Alec G. Stephenson. (2012).

  Coles, S. (2001). \emph{An introduction
  to statistical modeling of extreme values.} London: Springer.
}

\examples{
if(require(ismev)){
  ## from ismev
  data(portpirie)
  data(rain)

  detach(package:ismev)
  ppfit <- ismev::gev.fit(portpirie[,2])
  gev.diag(ppfit)
  ##
  (mlE <- MLEstimator(portpirie[,2], GEVFamilyMuUnknown(withPos=FALSE)))
  gev.diag(mlE)

\donttest{
  ## not tested on CRAN because it takes some time...
  gev.prof(mlE, m = 10, 4.1, 5)
  gev.profxi(mlE, -0.3, 0.3)

  rnfit <- ismev::gpd.fit(rain,10)
  gpd.diag(rnfit)
  ##
  mlE2 <- MLEstimator(rain[rain>10], GParetoFamily(loc=10))
  gpd.diag(mlE2)

  gpd.prof(mlE2, m = 10, 55, 77)
  gpd.profxi(mlE2, -0.02, 0.02)
}
}
}

\keyword{graphics}

