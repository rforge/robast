findLocOptGEV <- function(xi, p1=0.1, p2=0.3){
  G <- GEVFamily(shape=xi,scale=1)
  findf <- function(x, coord){
      L2deriv(G)[[1]]@Map[[coord]](x)
  }
  q1 <- q(G)(p1)
  q2 <- q(G)(p2)
  q0 <- q(G)(0)
  qi <- q(G)(.9999999)
  x11 <- optimize(findf,lower=q0,upper=q1,maximum=FALSE,coord=1)$min
  x12 <- optimize(findf,lower=q0,upper=q1,maximum=TRUE,coord=2)$max
  x1 <- (x11+x12)/2
  x2 <- optimize(findf,lower=q2,upper=qi,maximum=FALSE,coord=2)$min
  P1 <- p(G)(x1)
  P2 <- p(G)(x2)-P1
  P3 <- 1-p(G)(x2)
  return(c(x1,x2,P1,P2,P3))
}

getL2DistrGEV <- function(xi, p1=0.1, p2=0.3, gridN=1000, eps=1e-8, log=TRUE){
  G <- GEVFamily(shape=xi,scale=1)
  findf <- function(x, coord){
      L2deriv(G)[[1]]@Map[[coord]](x)
  }
  q1 <- q(G)(p1)
  q2 <- q(G)(p2)
  q0 <- q(G)(0)
  qi <- q(G)(.99999999)
  x11 <- optimize(findf,lower=q0,upper=q1,maximum=FALSE,coord=1, tol=eps)$min
  x12 <- optimize(findf,lower=q0,upper=q1,maximum=TRUE,coord=2, tol=eps)$max
  x1 <- (x11+x12)/2
  x2 <- optimize(findf,lower=q2,upper=qi,maximum=FALSE,coord=2, tol=eps)$min
  P1 <- p(G)(x1)
  P2 <- p(G)(x2)-P1
  P3 <- 1-p(G)(x2)
  y01 <- findf(q0,1)
  if(is.na(y01)) y01 <- Inf
  y11 <- findf(x1,1)
  if(is.na(y11)) y11 <- -Inf
  y21 <- findf(x2,1)
  y31 <- findf(qi,1)
  if(is.na(y31)) y31 <- Inf
  y02 <- findf(q0,2)
  if(is.na(y02)) y02 <- -Inf
  y12 <- findf(x1,2)
  if(is.na(y12)) y12 <- Inf
  y22 <- findf(x2,2)
  y32 <- findf(qi,2)
  if(is.na(y32)) y32 <- Inf
  m1 <- min(y01,y11,y21,y31)
  m2 <- min(y02,y12,y22,y32)
  M1 <- max(y01,y11,y21,y31)
  M2 <- max(y02,y12,y22,y32)
  u0 <- seq(0,1,length=gridN)
  q0i <- q(G)(u0*P1)
  q1i <- q(G)(P1+u0*P2)
  q2i <- q(G)(P1+P2+u0*P3)
  modifina <- function(x, asc=TRUE){
    l <- length(x)
    if(!asc) x <- rev(x)
    wina <- which(is.na(x))
    wina0 <- wina[wina < l/2]
    wina1 <- wina[wina > l/2]
    x[wina0] <- -Inf
    x[wina1] <- Inf
#    print(x)
    x1 <- sort(unique(x))
#    print(str(x));print(str(x1))
    x <- x1
  return(x)
  }

  f0i1 <- modifina(findf(q0i,1))
  f1i1 <- modifina(findf(q1i,1))
  f2i1 <- modifina(findf(q2i,1))
  f0i2 <- modifina(findf(q0i,2))
  f1i2 <- modifina(findf(q1i,2))
  f2i2 <- modifina(findf(q2i,2))
  makefinite <- function(x) x[is.finite(x)]
  fall1a <- fall1 <- sort(unique(c(y01,y11,y21,y31,f0i1,f1i1,f2i1)))
  fall2a <- fall2 <- sort(unique(c(y02,y12,y22,y32,f0i2,f1i2,f2i2)))
  fall1 <- makefinite(fall1)
  fall2 <- makefinite(fall2)
  getDs <- function(u,p,w){
    p0 <- 0*w
#    print(c(um,uM))
#    print(c(w[1],(rev(w))[1]))
    iu <- 1
    iw <- 1
    lu <- length(u)
    lw <- length(w)
    wa <- ua <- 0
    while(iu < lu){
       while(iw< lw &&w[iw] <= u[iu] ){
          p0[iw] <- if(w[iw]<=u[iu] && wa >= ua && w[iw]>ua )
                       (w[iw]-wa)/(u[iu]-ua)*p else 0
          wa <- w[iw]
          iw <- iw + 1
       }
       ua <- u[iu]
       iu <- iu + 1
    }
    return(p0)
  }
  i01 <- any(!is.finite(fall1a) & fall1a <0 )
  i11 <- any(!is.finite(fall1a) & fall1a >0 )
  i02 <- any(!is.finite(fall2a) & fall2a <0 )
  i12 <- any(!is.finite(fall2a) & fall2a >0 )
  p0i1 <- getDs(makefinite(f0i1),P1/gridN,fall1)
  print(c(sum(p0i1),P1))
  p1i1 <- getDs(makefinite(f1i1),P2/gridN,fall1)
  print(c(sum(p1i1),P2,sum(p0i1)+sum(p1i1)))
  p2i1 <- getDs(makefinite(f2i1),P3/gridN,fall1)
  print(c(sum(p2i1),P3,sum(p0i1)+sum(p1i1)+sum(p2i1)))
  p0i2 <- getDs(makefinite(f0i2),P1/gridN,fall2)
  print(c(sum(p0i2),P1))
  p1i2 <- getDs(makefinite(f1i2),P2/gridN,fall2)
  print(c(sum(p1i2),P2,sum(p0i2)+sum(p1i2)))
  p2i2 <- getDs(makefinite(f2i2),P3/gridN,fall2)
  print(c(sum(p2i2),P3,sum(p0i2)+sum(p1i2)+sum(p2i2)))
  grid1 <- cbind(fall1,cumsum(p0i1+p1i1+p2i1))
  grid2 <- cbind(fall2,cumsum(p0i2+p1i2+p2i2))
  grid1[,2]<- grid1[,2]/rev(grid1[,2])[1]*(1-eps)
  grid2[,2]<- grid2[,2]/rev(grid2[,2])[1]*(1-eps)
  grid1a <- grid1
  grid2a <- grid2
  if(i01) grid1a <- rbind(c(-Inf,0),grid1a)
  if(i11) grid1a <- rbind(grid1a, c(Inf,1))
  if(i02) grid2a <- rbind(c(-Inf,0),grid2a)
  if(i12) grid2a <- rbind(grid2a, c(Inf,1))
  pfun10  <- approxfun(grid1[,1],grid1[,2])
  dfun10  <- splinefun(grid1[,1],grid1[,2])
  qfun10  <- approxfun(grid1[,2],grid1[,1])
  pfun20  <- approxfun(grid2[,1],grid2[,2])
  dfun20  <- splinefun(grid2[,1],grid2[,2])
  qfun20  <- approxfun(grid2[,2],grid2[,1])
  pfun1 <- function(q, lower.tail= TRUE, log.p=FALSE){
     if(i01){i0 <- q<m1-eps} else {i0 <- q>m1-eps & q<m1+eps}
     if(i11){i1 <- q>M1+eps} else {i1 <- q>M1-eps & q<M1+eps}
     ina <- rep(FALSE,length.out=length(q))
     if(!i01){ina <- q<m1}
     if(!i11){ina <- ina | (q>M1)}
     ii <- q>=m1+eps & q <=M1-eps
     p <- q*0
     p[ina] <- NA
     p[i0] <- 0
     p[i1] <- 1
     p[ii] <- pfun10(q[ii])
     if(!lower.tail) p <- 1-p
     if(log.p) p <- log(p)
     return(p)
  }
  pfun2 <- function(q, lower.tail= TRUE, log.p=FALSE){
     if(i02){i0 <- q<m2-eps} else {i0 <- q>m2-eps & q<m2+eps}
     if(i12){i1 <- q>M2+eps} else {i1 <- q>M2-eps & q<M2+eps}
     ina <- rep(FALSE,length.out=length(q))
     if(!i02){ina <- q<m2}
     if(!i12){ina <- ina | (q>M2)}
     ina <- q<m2
     i0 <- q>=m2 & q<m2+eps
     i1 <- q>M2-eps
     ii <- q>=m2+eps & q <=M2-eps
     p <- q*0
     p[ina] <- NA
     p[i0] <- 0
     p[i1] <- 1
     p[ii] <- pfun20(q[ii])
     if(!lower.tail) p <- 1-p
     if(log.p) p <- log(p)
     return(p)
  }
  qfun1 <- function(p, lower.tail= TRUE, log.p=FALSE){
     if(log.p) p <- exp(p)
     if(!lower.tail) p <- 1-p
     ina <- p<0 | p>1
     i0 <- p>=0 & p<=eps
     i1 <- p>1-eps & p <= 1
     ii <- p>=eps & p <= 1-eps
     q <- p*0
     q[ina] <- NA
     q[i0] <- m1
     q[i1] <- M1
     q[ii] <- qfun10(p[ii])
     return(q)
  }
  qfun2 <- function(p, lower.tail= TRUE, log.p=FALSE){
     if(log.p) p <- exp(p)
     if(!lower.tail) p <- 1-p
     ina <- p<0 | p>1
     i0 <- p>=0 & p<=eps
     i1 <- p>1-eps & p <= 1
     ii <- p>=eps & p <= 1-eps
     q <- p*0
     q[ina] <- NA
     q[i0] <- m2
     q[i1] <- M2
     q[ii] <- qfun20(p[ii])
     return(q)
  }
  dfun1 <- function(x, log = FALSE){
     d <- dfun10(x,deriv=1)
     d <- pmax(d,0)
     if(log) d <- log(d)
     return(d)
  }
  dfun2 <- function(x, log = FALSE){
     d <- dfun20(x,deriv=1)
     d <- pmax(d,0)
     if(log) d <- log(d)
     return(d)
  }

  par(mfrow=c(2,2))
  plot(fall1,grid1[,2],type="l", xlab=expression(Lambda[1]), ylab="cdf", log= if(log) "y" else "")
  plot(fall2,grid2[,2],type="l", xlab=expression(Lambda[2]), ylab="cdf", log= if(log) "y" else "")
  plot(fall1,dfun1(fall1),type="l", xlab=expression(Lambda[1]), ylab="density", log= if(log) "y" else "")
  plot(fall2,dfun2(fall2),type="l", xlab=expression(Lambda[2]), ylab="density", log= if(log) "y" else "")
  par(mfrow=c(1,1))

  return(list(c(x1,x2,P1,P2,P3),grid1=grid1a,grid2=grid2a,minmax1=c(m1,M1),
              minmax2=c(m2,M2), y1=c(y01,y11,y21,y31), y2=c(y02,y12,y22,y32),
              d1=dfun1,p1=pfun1,q1=qfun1,d2=dfun2,p2=pfun2,q2=qfun2))
}
res <- getL2DistrGEV(xi=.3, p1=0.1, p2=0.3, gridN=40000,eps=1e-10)
res <- getL2DistrGEV(xi=3, p1=0.1, p2=0.3, gridN=1000)
