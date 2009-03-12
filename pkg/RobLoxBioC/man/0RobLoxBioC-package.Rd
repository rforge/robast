\name{RobLoxBioC-package}
\alias{RobLoxBioC-package}
\alias{RobLoxBioC}
\docType{package}
\title{
Infinitesimally robust estimators for preprocessing omics data
}
\description{
Functions for the determination of optimally robust influence curves 
and estimators for preprocessing omics data, in particular gene expression data.
}
\details{
\tabular{ll}{
Package: \tab RobLoxBioC\cr
Version: \tab 0.1 \cr
Date: \tab 2009-03-12 \cr
Depends: \tab R (>= 2.8.1), methods, methods, Biobase, affy\cr
LazyLoad: \tab yes\cr
License: \tab LGPL-3\cr
URL: \tab http://robast.r-forge.r-project.org/\cr
}
}
\author{
Matthias Kohl \email{Matthias.Kohl@stamats.de}\cr%

Maintainer: Matthias Kohl  \email{matthias.kohl@stamats.de}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2008) The Costs of not Knowing
  the Radius. Statistical Methods and Applications \emph{17}(1) 13-40.
}
\seealso{
\code{\link[RobLox]{roblox}}, \code{\link[RobLox]{rowRoblox}}
}
\examples{
library(RobLoxBioC)
}
\keyword{package}
