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
Package: \tab RobAStBase \cr
Version: \tab 1.2.0 \cr
Date: \tab 2019-04-05 \cr
Depends: \tab R(>= 3.4), methods, rrcov, distr(>= 2.8.0), distrEx(>= 2.8.0),
    distrMod(>= 2.8.0),RandVar(>= 1.2.0) \cr
Suggests: \tab ROptEst(>= 1.2.0), RUnit(>= 0.4.26) \cr
Imports: \tab startupmsg, graphics, grDevices, stats \cr
ByteCompile: \tab yes \cr
Encoding: \tab  latin1 \cr
License: \tab LGPL-3 \cr
URL: \tab http://robast.r-forge.r-project.org/\cr
VCS/SVNRevision: \tab 1205 \cr
}
}
\author{
Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de},\cr%
Matthias Kohl \email{Matthias.Kohl@stamats.de}\cr
Maintainer: Matthias Kohl  \email{matthias.kohl@stamats.de}}
\references{
  M. Kohl (2005). Numerical Contributions to the Asymptotic Theory of Robustness.
  Dissertation. University of Bayreuth.
}
\seealso{
\code{\link[distr:0distr-package]{distr-package}},
\code{\link[distrEx:0distrEx-package]{distrEx-package}},
\code{\link[distrMod:0distrMod-package]{distrMod-package}}
}
\section{Package versions}{
Note: The first two numbers of package versions do not necessarily reflect
 package-individual development, but rather are chosen for the
 RobAStXXX family as a whole in order to ease updating "depends"
 information.
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
