\name{RobRDA-package}
\alias{RobRDA-package}
\alias{RobRDA}
\docType{package}
\title{
RobRDA a sysdata.rda only package
}
\description{
This package only contains sysdata.rda (with corresponding interpolation grids
for speedup); it is currently used in package \pkg{RobExtremes}.
The code to produce its contents is split into two parts: (a) grid construction
and (b) interpolator construction. While the code for (a) can be drawn from
CRAN-packages \pkg{ROptEst} and \pkg{RobExtremes}, more specifically:
see \code{?.RMXE.xi} (\pkg{RobExtremes}) resp. \code{?.RMXE.th} (\pkg{ROptEst}),
as well as the contents of the (system) folder
of package RobExtremes, i.e.,
  \code{dir(file.path(system.file(package="RobExtremes"),"AddMaterial","interpolation"))},
the code for (b) resides in the present package (and does not need to know anything
about the grid construction). As it is not meant for users but rather for
developers, it is not exported to the namespace; still, it is documented,
see \code{?.generateInterpolators}.

The reason to separate the rda file from the actual R packages is to
keep the latter small while we expect the present package to only need
seldom updates.

See also mail exchange P.Ruckdeschel - U.Ligges on R-devel---
\url{https://stat.ethz.ch/pipermail/r-devel/2013-February/065794.html}.
}

\details{
\tabular{ll}{
Package: \tab RobAStRDA \cr
Version: \tab 0.9 \cr
Date: \tab 2013-09-11 \cr
Depends: \tab R (>= 2.14.0), methods \cr
Imports: ROptEst (>= 0.9)
LazyData: \tab yes \cr
LazyLoad: \tab yes \cr
ByteCompile: \tab yes \cr
License: \tab LGPL-3 \cr
URL: \tab http://robast.r-forge.r-project.org/\cr
SVNRevision: \tab 693 \cr
}
}

\keyword{package}

