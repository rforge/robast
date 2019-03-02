\name{getReq}
\alias{getReq}

\title{getReq -- computation of the radius interval where IC1 is better than IC2.}
\description{
 (tries to) compute a radius interval where IC1 is better than IC2,
 respectively the number of (worst-case) outliers interval where IC1 is
 better than IC2.
}
\usage{getReq(Risk,neighbor,IC1,IC2,n=1,upper=15, radOrOutl=c("radius","Outlier"), ...)}
\arguments{
  \item{Risk}{an object of class \code{"asGRisk"} -- the risk at which IC1 is better than IC2.}
  \item{neighbor}{ object of class \code{"Neighborhood"}; the neighborhood at which to compute the bias. }
  \item{IC1}{some IC of class \code{"IC"}}
  \item{IC2}{some IC of class \code{"IC"}}
  \item{n}{the sample size; by default set to 1; then the radius interval refers to starting radii 
    in the shrinking neighborhood setting of Rieder[94]. Otherwise the radius interval is scaled
    down accordingly.}
  \item{upper}{the upper bound of the radius interval in which to search}  
  \item{radOrOutl}{ a character string specifying whether an interval of
          radii  or a number of outliers is returned; must be one
          of "radius" (default) and "Outlier". }
  \item{\dots}{further arguments to be passed on \code{E()}.}
}
%\details{}
\value{The radius interval (given by its endpoints) where \code{IC1} is better than \code{IC2}
according to the risk. In case \code{IC2} is better than \code{IC1} as to both variance and bias,
the return value is \code{NA}.}
\references{
  Hampel et al. (1986) \emph{Robust Statistics}. 
  The Approach Based on Influence Functions. New York: Wiley.
  
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.
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
## MSE solution
N0.ICM <- optIC(model=N0.Rob1, risk=asMSE())

getReq(asMSE(),neighbor,N0.ICA,N0.ICM,n=1)
getReq(asMSE(),neighbor,N0.ICA,N0.ICM,n=30)

## Don't test to reduce check time on CRAN
\donttest{
## RMX solution
N0.ICR <- radiusMinimaxIC(L2Fam=N0, neighbor=neighbor,risk=asMSE())

getReq(asL1(),neighbor,N0.ICA,N0.ICM,n=30)
getReq(asL4(),neighbor,N0.ICA,N0.ICM,n=30)
getReq(asMSE(),neighbor,N0.ICA,N0.ICR,n=30)
getReq(asL1(),neighbor,N0.ICA,N0.ICR,n=30)
getReq(asL4(),neighbor,N0.ICA,N0.ICR,n=30)
getReq(asMSE(),neighbor,N0.ICM,N0.ICR,n=30)


### when to use MAD and when Qn 
##  for Qn, see C. Croux, P. Rousseeuw (1993). Alternatives to the Median 
##      Absolute Deviation, JASA 88(424):1273-1283
L2M <- NormScaleFamily()
IC.mad <- makeIC(function(x)sign(abs(x)-qnorm(.75)),L2M)
d.qn <- (2^.5*qnorm(5/8))^-1
IC.qn <- makeIC(function(x) d.qn*(1/4 - pnorm(x+1/d.qn) + pnorm(x-1/d.qn)), L2M)
getReq(asMSE(), neighbor, IC.mad, IC.qn)
getReq(asMSE(), neighbor, IC.mad, IC.qn, radOrOutl = "Outlier", n = 30)
# => MAD is better once r > 0.5144 (i.e. for more than 2 outliers for n = 30)
}
}
\concept{Hampel risk}
\concept{risk}
\keyword{robust}
