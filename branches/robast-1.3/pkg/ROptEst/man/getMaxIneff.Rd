\name{getMaxIneff}
\alias{getMaxIneff}

\title{getMaxIneff -- computation of the maximal inefficiency of an IC}
\description{
 computes the maximal inefficiency of an IC for the radius range [0,Inf).
}
\usage{getMaxIneff(IC, neighbor, biastype = symmetricBias(), 
                        normtype = NormType(), z.start = NULL, 
                        A.start = NULL, maxiter = 50, 
                        tol = .Machine$double.eps^0.4,
                        warn = TRUE, verbose = NULL, ...)}
\arguments{
  \item{IC}{some IC of class \code{IC}}
  \item{neighbor}{ object of class \code{Neighborhood}; 
                   the neighborhood at which to compute the bias. }
  \item{biastype}{a bias type of class \code{BiasType}}
  \item{normtype}{ a norm type of class \code{NormType}}
  \item{z.start}{ initial value for the centering constant. }
  \item{A.start}{ initial value for the standardizing matrix. }
  \item{maxiter}{ the maximum number of iterations. }
  \item{tol}{ the desired accuracy (convergence tolerance).}
  \item{warn}{ logical: print warnings. }  
  \item{verbose}{ logical: if \code{TRUE}, some messages are printed }
  \item{\dots}{ additional arguments to be passed to \code{E} }

}
%\details{}
\value{The maximal inefficiency, i.e.; a number in [1,Inf).}
\references{
  Hampel et al. (1986) \emph{Robust Statistics}. 
  The Approach Based on Influence Functions. New York: Wiley.
  
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.
  
  Rieder, H., Kohl, M. and Ruckdeschel, P. (2008) The Costs of not Knowing
  the Radius. Statistical Methods and Applications \emph{17}(1) 13-40.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2001) The Costs of not Knowing
  the Radius. Submitted. Appeared as discussion paper Nr. 81. 
  SFB 373 (Quantification and Simulation of Economic Processes),
  Humboldt University, Berlin; also available under
  \url{www.uni-bayreuth.de/departments/math/org/mathe7/RIEDER/pubs/RR.pdf}

}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@fraunhofer.itwm.de}}
%\note{}
\examples{
N0 <- NormLocationFamily(mean=2, sd=3)
## L_2 family + infinitesimal neighborhood
neighbor <- ContNeighborhood(radius = 0.5)
N0.Rob1 <- InfRobModel(center = N0, neighbor = neighbor)
## OBRE solution (ARE 95%)
N0.ICA <- optIC(model = N0.Rob1, risk = asAnscombe(.95))
## OMSE solution radius 0.5
N0.ICM <- optIC(model=N0.Rob1, risk=asMSE())
## RMX solution 
N0.ICR <- radiusMinimaxIC(L2Fam=N0, neighbor=neighbor,risk=asMSE())

getMaxIneff(N0.ICA,neighbor)
getMaxIneff(N0.ICM,neighbor)
getMaxIneff(N0.ICR,neighbor)

## Don't run to reduce check time on CRAN
\donttest{
N0ls <- NormLocationScaleFamily()
ICsc <- makeIC(list(sin,cos),N0ls)
getMaxIneff(ICsc,neighbor)
}
}
\concept{Inefficiency}
\concept{risk}
\keyword{robust}
