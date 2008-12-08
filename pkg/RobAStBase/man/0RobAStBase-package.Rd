\name{RobAStBase-package}
\alias{RobAStBase-package}
\alias{RobAStBase}
\docType{package}
\title{
Robust Asymptotic Statistics
}
\description{
Base S4-classes and functions for robust asymptotic statistics.
}
\details{
\tabular{ll}{
Package: \tab RobAStBase\cr
Version: \tab 0.1.4\cr
Date: \tab 2008-12-08\cr
Depends: \tab R(>= 2.7.0), methods, distr(>= 2.0), distrEx(>= 2.0),
distrMod(>= 2.0), RandVar(>= 0.6.3)\cr
LazyLoad: \tab yes\cr
License: \tab LGPL-3\cr
URL: \tab http://robast.r-forge.r-project.org/\cr
}
}
\author{
Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de},\cr%
Matthias Kohl \email{Matthias.Kohl@stamats.de}\cr

Maintainer: Matthias Kohl  \email{matthias.kohl@stamats.de}}
\references{
  M. Kohl (2005). Numerical Contributions to the Asymptotic Theory of Robustness.
  Dissertation. University of Bayreuth.
}
\seealso{
\code{\link[distr]{distr-package}}, \code{\link[distrEx]{distrEx-package}},
\code{\link[distrMod]{distrMod-package}}
}
\examples{
library(RobAStBase)

## some L2 differentiable parametric family from package distrMod, e.g.
B <- BinomFamily(size = 25, prob = 0.25) 

## classical optimal IC
IC0 <- optIC(model = B, risk = asCov())
plot(IC0) # plot IC
checkIC(IC0, B)
}
\keyword{package}
