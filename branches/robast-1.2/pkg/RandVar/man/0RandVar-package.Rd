\name{RandVar-package}
\alias{RandVar-package}
\alias{RandVar}
\docType{package}
\title{
Implementation of Random Variables
}
\description{
Implementation of random variables by means of S4 classes and methods.
}
\details{
\tabular{ll}{
Package: \tab RandVar \cr
Version: \tab 1.2.0 \cr
Date: \tab 2019-04-02 \cr
Depends: \tab R(>= 3.4), methods, distr(>= 2.8.0), distrEx(>= 2.8.0) \cr
Imports: \tab startupmsg \cr
ByteCompile: \tab yes \cr
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
\code{\link[distr:0distr-package]{distr-package}}, \code{\link[distrEx:0distrEx-package]{distrEx-package}}
}
\section{Package versions}{
Note: The first two numbers of package versions do not necessarily reflect
 package-individual development, but rather are chosen for the
 RobAStXXX family as a whole in order to ease updating "depends"
 information.
}

\examples{
library(RandVar)
#vignette("RandVar")
}
\keyword{package}
