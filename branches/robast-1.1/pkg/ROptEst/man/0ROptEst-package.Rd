\name{ROptEst-package}
\alias{ROptEst-package}
\alias{ROptEst}
\docType{package}
\title{
Optimally robust estimation
}
\description{
Optimally robust estimation in general smoothly parameterized models 
using S4 classes and methods.
}
\details{
\tabular{ll}{
Package: \tab ROptEst \cr
Version: \tab 1.0 \cr
Date: \tab 2015-05-03 \cr
Depends: \tab R(>= 2.14.0), methods, distr(>= 2.5.2), distrEx(>= 2.5), distrMod(>= 2.5.2), 
RandVar(>= 0.9.2), RobAStBase(>= 1.0) \cr
Suggests: \tab RobLox, MASS \cr
Imports: \tab startupmsg \cr
ByteCompile: \tab yes \cr
License: \tab LGPL-3 \cr
URL: \tab http://robast.r-forge.r-project.org/\cr
SVNRevision: \tab 728 \cr
}
}
\author{
Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de},\cr%
Matthias Kohl \email{Matthias.Kohl@stamats.de}\cr

Maintainer: Matthias Kohl  \email{matthias.kohl@stamats.de}}
\references{
  M. Kohl (2005). Numerical Contributions to the Asymptotic Theory of Robustness.
  Dissertation. University of Bayreuth.

  M. Kohl, P. Ruckdeschel, H. Rieder (2010). Infinitesimally Robust Estimation in 
  General Smoothly Parametrized Models. Statistical Methods and Application 19(3):333-354. 
}
\seealso{
\code{\link[distr:0distr-package]{distr-package}}, 
\code{\link[distrEx:0distrEx-package]{distrEx-package}},
\code{\link[distrMod:0distrMod-package]{distrMod-package}}, 
\code{\link[RandVar:0RandVar-package]{RandVar-package}},
\code{\link[RobAStBase:0RobAStBase-package]{RobAStBase-package}}
}
\section{Package versions}{
Note: The first two numbers of package versions do not necessarily reflect
 package-individual development, but rather are chosen for the
 RobAStXXX family as a whole in order to ease updating "depends"
 information.
}
\examples{
library(ROptEst)

## Example: Rutherford-Geiger (1910); cf. Feller~(1968), Section VI.7 (a)
x <- c(rep(0, 57), rep(1, 203), rep(2, 383), rep(3, 525), rep(4, 532), 
       rep(5, 408), rep(6, 273), rep(7, 139), rep(8, 45), rep(9, 27), 
       rep(10, 10), rep(11, 4), rep(12, 0), rep(13, 1), rep(14, 1))

## ML-estimate from package distrMod
MLest <- MLEstimator(x, PoisFamily())
MLest
## confidence interval based on CLT
confint(MLest)

## compute optimally (w.r.t to MSE) robust estimator (unknown contamination)
robEst <- roptest(x, PoisFamily(), eps.upper = 0.1, steps = 3)
estimate(robEst)
## check influence curve
pIC(robEst)
checkIC(pIC(robEst))
## plot influence curve
plot(pIC(robEst))
## confidence interval based on LAN - neglecting bias
confint(robEst)
## confidence interval based on LAN - including bias
confint(robEst, method = symmetricBias())
}
\keyword{package}
