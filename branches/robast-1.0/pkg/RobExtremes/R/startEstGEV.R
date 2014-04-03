.getXiGrid <- function(){seq(-0.48,5,by=0.5)}


.getBetaXiGEV <- function(x, mu, xiGrid = .getXiGrid(), withPos=TRUE){
  x0 <- x-mu
  s0 <- max(x0)-min(x0)
  crit0 <- Inf
  fu <- function(x,...) .getBetaXiGEV(x,mu,xiGrid = xiGrid,withPos=withPos)
  for(xi in xiGrid){
      funl <- function(sig){
         mygev1 <- GEV(loc=0,scale=sig,shape=xi)
         CvMDist(x0,mygev1)
      }
      intv <- if(xi<0) c(max(1e-5,xi/4*min(x)-mu), max(xi*min(x)-mu,2e-5)) else c(max(1e-5,xi*max(x)-mu),max(s0,2*xi*max(x)-mu))
      sigCvMMD1 <- optimize(funl, interval=intv)$minimum
      mygev <- GEVFamily(loc=0,scale=sigCvMMD1,shape=xi, withPos=withPos,
                         start0Est = fu, ..withWarningGEV=FALSE)
      mde0 <- try(MDEstimator(x0, mygev, distance=CvMDist, startPar=c("scale"=sigCvMMD1,"shape"=xi)),silent=TRUE)
      if(!is(mde0,"try-error")){
          es <- estimate(mde0)
          crit1 <- criterion(mde0)
          if(crit1<crit0&(1+min(es[2]*x0/es[1])>0)){
             mdeb <- mde0
             crit0 <- crit1
          }
      }
  }
  es <- estimate(mdeb)
  names(es) <- c("scale","shape")
  return(es)
}

.getMuBetaXiGEV <- function(x, xiGrid = .getXiGrid(), withPos=TRUE){
  mu <- quantile(x,exp(-1))
  es <- .getBetaXiGEV(x=x, mu=mu, xiGrid=xiGrid, withPos=withPos)
  es0 <- c(mu,es)
  names(es0) <- c("loc","scale","shape")
  return(es0)
}
