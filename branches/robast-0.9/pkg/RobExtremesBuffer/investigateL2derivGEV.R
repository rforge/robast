require(RobExtremes)
z <- seq(0.000001,0.999999,by=0.001)
x1 <- seq(0,1,by=0.1)
psf <- function(pf1){
  function(x) (1-pf1(-x)+pf1(x))/2
}
qsf <- function(pf1,qf1){
  u00 <- c(1e-6,1e-5)
  u0 <- sort(c(u00,seq(1e-7,1-1e-7,by=0.0005),1-u00))
  p2 <- psf(pf1)
  li1 <- min(qf1(0),-qf1(1))
  re1 <- max(qf1(1),-qf1(0))
  li <- if(!is.finite(li1)) min(qf1(1e-7),-qf1(1-1e-7)) else li1
  re <- if(!is.finite(re1)) max(qf1(1-1e-7),-qf1(1e-7)) else re1
  fct1 <- function(u){
      fct2 <- function(x) p2(x)-u
      uniroot(fct2,lower=li, upper=re)$root
  }
  q0 <- sapply(u0,fct1)
  q1 <- approxfun(u0,q0)
  q2 <- function(u){
     ni01 <- u<0 | u>1
     i0 <- !ni01& u< 1e-7
     i1 <- !ni01& (1-u)< 1e-7
     ni01.2 <- ni01|i0|i1
     qa <- 0*u
     qa[!ni01.2] <- q1(u[!ni01.2])
     qa[ni01] <- NA
     qa[i0] <- li1
     qa[i1] <- re1
     return(qa)
  }
  return(q2)
}

plotG <- function(G=GEVFamily(shape=0.7,scale=1), pfa,qfa){
  if(missing(pfa)) pfa <- p(G)
  if(missing(qfa)) qfa <- q(G)

  pf0=psf(pfa)
  qf0=qsf(pfa,qfa)
  axi <- function(i,w=q(G)(x1))
             axis(i,at=x1,labels=round(w,2),cex.axis=0.5,las=2)
  tus <- function(){
          axi(1,q(G)(x1))
          axi(2,qf0(x1))
          axi(3,x1)
          box()
        }
  ploti <- function(i){
        plot(z, pf0(L2deriv(G)[[1]]@Map[[i]](q(G)(z))),
              type="l",xlab="x",ylab=if(i==1) "shape" else "scale",axes=F)
        tus()}
  plot2 <- function(){
        plot(z, pf0(L2deriv(G)[[1]]@Map[[1]](q(G)(z)))*pf0(L2deriv(G)[[1]]@Map[[2]](q(G)(z))),
              type="l",xlab="x",ylab="shape x scale",axes=F)
        tus()}
  windows()
  plot(G)
  if(p(G)(0)<=0){
    windows()
    plot(G,log="x")
  }
  windows()
  if(ncol(trafo(G))>1) {par(mfrow=c(3,1))
  ploti(1);
  ploti(2)
  plot2()
  par(mfrow=c(1,1))}else ploti(1)
}
plotG(G=NormLocationScaleFamily())
plotG(GEVFamily(shape=0.5,scale=1))
plotG(GEVFamily(shape=0.7,scale=1))
plotG(GEVFamily(shape=1,scale=1))
plotG(GEVFamily(shape=3,scale=1))
plotG(GEVFamily(shape=5,scale=1))
plotG(GEVFamily(shape=9,scale=1))
