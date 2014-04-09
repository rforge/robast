.getXiGrid <- function(){c(0.1,seq(-0.48,5,by=0.5))}


.getBetaXiGEV <- function(x, mu, xiGrid = .getXiGrid(), withPos=TRUE, secLevel = 0.7){

  n <- length(x)
  epsn <- min(floor(secLevel*sqrt(n))+1,n)

  x0 <- x-mu
  s0 <- max(x0)-min(x0)
  crit0 <- Inf
  fu <- function(x,...) .getBetaXiGEV(x,mu,xiGrid = xiGrid,withPos=withPos)
  for(xi in xiGrid){
      funl <- function(sig){
         mygev1 <- GEV(loc=0,scale=sig,shape=xi)
         CvMDist(x0,mygev1)
      }
      intlo <- quantile(-xi*(x-mu),1-epsn/n)
      intv <-  c(max(1e-5,intlo), s0)
      sigCvMMD1 <- optimize(funl, interval=intv)$minimum
      mygev <- GEVFamily(loc=0,scale=sigCvMMD1,shape=xi, withPos=withPos,
                         start0Est = fu, ..withWarningGEV=FALSE)
      mde0 <- try(MDEstimator(x0, mygev, distance=CvMDist, startPar=c("scale"=sigCvMMD1,"shape"=xi)),silent=TRUE)
      es0 <- c(NA,NA)
      if(!is(mde0,"try-error")){
          es <- estimate(mde0)
          crit1 <- criterion(mde0)
          if(quantile(1+es[2]*x0/es[1], epsn/n)>0){
             if(crit1<crit0){
                mdeb <- mde0
                crit0 <- crit1
                es0 <- es
             }
          }else{
             es[1] <- intlo+1e-5
             mygev2 <- GEV(loc=0,scale=es[1],shape=es[2])
             crit1 <- CvMDist(x0,mygev2)
             if(crit1<crit0){
                crit0 <- crit1
                es0 <- es
             }
          }
      }
  }
  names(es) <- c("scale","shape")
  return(es)
}

.getMuBetaXiGEV <- function(x, xiGrid = .getXiGrid(), withPos=TRUE, secLevel = 0.7){
  mu <- quantile(x,exp(-1))
  es <- .getBetaXiGEV(x=x, mu=mu, xiGrid=xiGrid, withPos=withPos, secLevel = secLevel)
  es0 <- c(mu,es)
  names(es0) <- c("loc","scale","shape")
  return(es0)
}
