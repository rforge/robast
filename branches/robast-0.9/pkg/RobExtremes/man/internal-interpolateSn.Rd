\name{internal_interpolateSn_helpers_for_distrMod}
\alias{internal_interpolateSn_helpers_for_distrMod}
\alias{.saveInterpGrid}

\title{Internal helper functions for generating interpolation grids for Sn in
       package distrMod}

\description{
This function is used internally to generate interpolation grids,
for example for Sn in package ``distrMod'', to be stored in the
respective \file{sysdata.rda} file. }

\usage{
.saveInterpGrid(xiGrid = getShapeGrid(), PFam = GParetoFamily(),
            sysRdaFolder, getFun = getSnGrid,
            ..., nameInSysdata = ".SnGrids",
            withSmooth = TRUE, withPrint = TRUE, Y = NULL, elseFun = NULL)
}

\arguments{
  \item{xiGrid}{numeric; grid of shape values. }
  \item{PFam}{an object of class \code{"ParamFamily"}. In our use case,
              it is the parametric family at which to evaluate the LDEstimator;
              hence the respective (main) parameter must contain \code{"scale"}
              and \code{"shape"}. }
  \item{sysRdaFolder}{the folder where \pkg{distrMod} (or the respective package)
                      is being developed; must not be missing. }
  \item{getFun}{function with first argument \code{xi}, and second argument
                \code{PFam} and last arguments \code{withSmooth} and
                \code{withPrint}; produces the
                y-value for the interpolation grid; in our case the Sn-value. }
  \item{\dots}{further arguments to be passed on to \code{getFun}. }
  \item{nameInSysdata}{name under which the list of interpolated grids is stored
                       in file \file{sysdata.rda}. }
  \item{withSmooth}{logical of length 1: shall a smoothing spline be used?}
  \item{withPrint}{logical of length 1: shall current shape value be printed out?}
  \item{Y}{array or \code{NULL}; if non-null, contains precomputed y-values, so
           that call to \code{getFun} can be omitted. }
  \item{elseFun}{function or \code{NULL}; if \code{Y} is non-null, contains
           function to transform \code{Y} to desired return value. }
}
\details{
  The code first loads the contents of file \file{sysdata.rda} into an
  environment which is particularly created for this purpose (if this
  file exists). It then looks up whether a respective entry for family
  \code{PFam} already exists in \file{sysdata.rda}-object \code{nameInSysdata}.
  If this is the case the developer is asked whether he wants to overwrite
  the respective entry, and if so he does so with the results of a respective
  call to \code{getFun}.
  
  In case there has not been an \file{sysdata.rda}-object \code{nameInSysdata}
  so far, the code creates one and writes the results of a respective
  call to \code{getFun} to it.
}
\note{This function is only meant for the developers of package \pkg{distrMod}
      (or respective packages).
      It can be used to speed up things by interpolation.
      Our use case is a speed up for further scale-shape families (or enhance
      existing speed-ups) such that the respective grids are stored in
      the \file{sysdata.rda} file of this package and can be used in
      (exported) new particular methods for functional \code{Sn}.}
\value{
  \code{invisible(NULL)}
}
\seealso{\code{\link{interpolateSn}}}
\examples{
\dontrun{
### code to produce grid for GPareto:
RobExtremes:::.saveInterpGrid(sysRdaFolder =
        "C:/rtest/RobASt/branches/robast-0.9/pkg/RobExtremes/R",
              accuracy = 5000)
}
}
\keyword{internal}
\concept{utilities}
