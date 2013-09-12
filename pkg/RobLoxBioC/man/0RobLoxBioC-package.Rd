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
Package: \tab RobLoxBioC \cr
Version: \tab 0.9 \cr
Date: \tab 2013-09-12 \cr
Depends: \tab R(>= 2.14.0), methods, Biobase, affy, beadarray, distr, RobLox,
lattice, RColorBrewer \cr
LazyLoad: \tab yes \cr
ByteCompile: \tab yes \cr
License: \tab LGPL-3 \cr
URL: \tab http://robast.r-forge.r-project.org/\cr
SVNRevision: \tab 696 \cr
Encoding: \tab latin1 \cr
}
}
\author{
Matthias Kohl \email{Matthias.Kohl@stamats.de}\cr%

Maintainer: Matthias Kohl  \email{matthias.kohl@stamats.de}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Kohl M. and Deigner H.P. (2010). Preprocessing of gene expression data by optimally 
  robust estimators. \emph{BMC Bioinformatics}, 11:583. 

  M. Kohl, P. Ruckdeschel, and H. Rieder (2010). Infinitesimally Robust Estimation 
  in General Smoothly Parametrized Models. \emph{Statistical Methods and Application}, 
  \bold{19}(3):333-354. 
  
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2008) The Costs of not Knowing
  the Radius. \emph{Statistical Methods and Applications} \bold{17}(1) 13-40.
  Extended version: \url{http://www.stamats.de/RRlong.pdf}
}
\seealso{
\code{\link[RobLox]{roblox}}, \code{\link[RobLox]{rowRoblox}}
}
\section{Package versions}{
Note: The first two numbers of package versions do not necessarily reflect
 package-individual development, but rather are chosen for the
 RobAStXXX family as a whole in order to ease updating "depends"
 information.
}
\examples{
library(RobLoxBioC)
}
\keyword{package}
