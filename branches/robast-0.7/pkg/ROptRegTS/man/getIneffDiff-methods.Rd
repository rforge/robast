\name{getIneffDiff-methods}
\docType{methods}
\alias{getIneffDiff-methods}
\alias{getIneffDiff,numeric,L2RegTypeFamily,Neighborhood,asMSE-method}
\alias{getIneffDiff,numeric,L2RegTypeFamily,Av2CondContNeighborhood,asMSE-method}

\title{ Methods for Function getIneffDiff in Package `ROptRegTS' }
\description{
  Methods for function \code{getIneffDiff} in package \pkg{ROptRegTS}. 
  These methods are rarely called directly. They are used to compute
  the radius minimax IC and the least favorable radius.
}
\section{Methods}{
\describe{
  \item{radius = "numeric", L2Fam = "L2RegTypeFamily", neighbor = "Neighborhood", risk = "asMSE"}{ 
    computes difference of asymptotic MSE--inefficiency for the boundaries of a given radius interval. }

  \item{radius = "numeric", L2Fam = "L2RegTypeFamily", neighbor = "Av2CondContNeighborhood", risk = "asMSE"}{ 
    computes difference of asymptotic MSE--inefficiency for the boundaries of a given radius interval. }
}}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\seealso{\code{\link[ROptEst]{getIneffDiff}}}
\keyword{methods}
\keyword{robust}
