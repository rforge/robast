\name{RobRDA-package}
\alias{RobRDA-package}
\alias{RobRDA}
\docType{package}
\title{
RobRDA a sysdata.rda only package
}
\description{
This package only contains sysdata.rda (with corresponding interpolation grids
for speedup); it is currently used in package RobExtremes.
The code to produce its contents can be drawn from CRAN-packages
ROptEst and RobExtremes, more specifically: see ?.RMXE.xi (RobExtremes)
resp. ?.RMXE.th (ROptEst), as well as the contents of the (system) folder
of package RobExtremes, i.e.,
  dir(file.path(system.file(package="RobExtremes"),"AddMaterial","interpolation"))

The reason to separate the rda file from the actual R packages is to
keep the latter small while we expect this package to only need seldom updates.
See also mail exchange P.Ruckdeschel - U.Ligges on R-devel---
\url{https://stat.ethz.ch/pipermail/r-devel/2013-February/065794.html}.
}

\details{
\tabular{ll}{
Package: \tab RobAStRDA \cr
Version: \tab 0.9 \cr
Date: \tab 2012-05-17 \cr
Depends: \tab R (>= 2.14.0), methods \cr
Imports: ROptEst (>= 0.9)
LazyData: \tab yes \cr
LazyLoad: \tab yes \cr
ByteCompile: \tab yes \cr
License: \tab LGPL-3 \cr
URL: \tab http://robast.r-forge.r-project.org/\cr
SVNRevision: \tab 450 \cr
}
}

\keyword{package}

