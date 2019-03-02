### NOT YET USED
if(FALSE){
#.my.do.pts(data=data, fct=absInfoEval(x,absInfo.f),
#           which.lbs=which.lbs, which.Order=which.Order)

# .my.do.points does
#     (1) it selects points by which.lbs and which.order
#     (2) it rescales their coordinates according to a transformed coordinate system
#     (3) it determines their individual size according to cexfun
#     (4) it plots them
#     (5) it labels them

.my.do.pts <- function(data, fct, fctC, which.lbs, which.Order, distr,
                       jit.fac, jit.tol,
                       col.sh, col.ns, pch.sh, pch.ns,
                       cex.sh, cex.ns, cexfun.sh=NULL, cexfun.ns=NULL,
                       with.lab, lab.pts, lab.font, lab.pos,
                       rfct = NULL, rfctC = NULL,
                       scaleX, scaleX.fct, scaleX.inv,
                       scaleY, scaleY.fct, xlim, ylim, dots, dots.points,
                       fct0=NULL, fct0C=NULL){
      
      if(is.null(fct0)) fct0 <- function(x) absInfoEval(x, fct)
      if(is.null(fct0C)) fct0C <- function(x) absInfoEval(x, fctC)
      if(is.null(rfct)) rfct <- fct0
      if(is.null(rfctC)) rfctC <- fct0C
      
      res <- .my.inner.do.pts(data,fct0,rfct, which.lbs, which.Order, distr,
                             jit.fac[1:2], jit.tol[1:2], scaleX, scaleX.fct, scaleX.inv,
                             scaleY, scaleY.fct, xlim, ylim, dots)
      
      resC <- .my.inner.do.pts(data,fct0C,rfctC, which.lbs, which.Order, distr,
                             jit.fac[3:4], jit.tol[3:4], scaleX, scaleX.fct, scaleX.inv,
                             scaleY, scaleY.fct, xlim, ylim, dots)
      id <- res[["sh"]]$i
      idC <- resC[["sh"]]$i
      nid <- res[["ns"]]$i
      nidC <- resC[["ns"]]$i

      inScl.sh <- .my.inner.scale(res,resC,cex.sh,cexfun=cexfun.sh, what="sh")
      inScl.ns <- .my.inner.scale(res,resC,cex.ns,cexfun=cexfun.ns, what="ns")

      .do.pts(res[["sh"]]$resc$x, res[["sh"]]$resc$y,
             inScl.sh$f1, col.sh[[1]], pch.sh[id],dots.points)
      .do.pts(resC[["sh"]]$resc$x, resC[["sh"]]$resc$y,
             inScl.sh$f1c, col.sh[[2]], pch.sh[idC],dots.points)
      .do.pts(res[["ns"]]$resc$x, res[["ns"]]$resc$y,
             inScl.ns$f1, col.ns[[1]], pch.ns[nid],dots.points)
      .do.pts(resC[["ns"]]$resc$x, resC[["ns"]]$resc$y,
             inScl.ns$f1c, col.ns[[2]], pch.ns[nidC],dots.points)
      if(with.lab0){
         lab.pts <- if(is.null(lab.pts)) list(id=id, idC=idC)
                       else list(id=lab.pts[id],idC=lab.pts[idC])
         .tx(res[["sh"]]$resc$x, res[["sh"]]$resc$y, lab.pts[["id"]], inScl.sh$f1/2, col.sh[1], lab.font[1], lab.pos[1])
         .tx(resC[["sh"]]$resc$x, resC[["sh"]]$resc$y, lab.pts[["idC"]], inScl.sh$f1c/2, col.sh[2], lab.font[2], lab.pos[2])
      }
}

.my.inner.scale <- function(res,resC,cex,cexfun=NULL, what){
      c1fun <- if(is.null(cexfun)) NULL else cexfun[[1]]
      c2fun <- if(is.null(cexfun)) NULL else cexfun[[2]]
      f1 <-  .cexscale(res[["what"]]$y,rescC[["what"]]$y,cex=cex[1], fun = c1fun)
      f1c <- .cexscale(resC[["what"]]$y,resc[["what"]]$y,cex=cex[2], fun = c2fun)
      return(list(f1=f1,f1c=f1c))
}

.my.inner.do.pts <- function(data, srtfct, rescfct, which.lbs, which.Order, distr,
                             jit.fac, jit.tol, scaleX, scaleX.fct, scaleX.inv,
                             scaleY, scaleY.fct, xlim, ylim, dots){

       sel <- .SelectOrderData(data, srtfct, which.lbs, which.Order)

       i.d <- sel$ind
       i0.d <- sel$ind1
       y.d <- sel$y
       x.d <- sel$data

       nd <- length(data)
       n <- length(i.d)
       n.ns <- nd-n

       y.d.ns <- sel$y.ns
       x.d.ns <- sel$data.ns
       i.d.ns <- sel$ind.ns
       resc.pts <- .my.resc.pts(x.d,y.d,distr,jit.fac[1], jit.tol[1],
                                rescfct, scaleX, scaleX.fct, #scaleX.inv,
                                scaleY, scaleY.fct, xlim, ylim, dots)
       resc.pts.ns <- .my.resc.pts(x.d.ns, y.d.ns, distr, jit.fac[2],
                                   jit.tol[2], rescfct, scaleX,
                                   scaleX.fct, # scaleX.inv, 
                                   scaleY, scaleY.fct,
                                   xlim, ylim, dots)
       return(list(sh=list(i=i.d,   i0=i0.d,x=x.d,   y=y.d,   n=n,   resc=resc.pts),
                   ns=list(i=i.d.ns,        x=x.d.ns,y=y.d.ns,n=n.ns,resc=resc.pts.ns)))
}

.my.resc.pts <- function(x,y,distr, jit.fac, jit.tol, resc.fct, scaleX,
                         scaleX.fct, # scaleX.inv, 
                         scaleY, scaleY.fct,
                         xlim, ylim, dots){
       resc.dat <-.rescalefct(x, resc.fct,
                              scaleX, scaleX.fct, # scaleX.inv,
                              scaleY, scaleY.fct, xlim, ylim, dots)
       x.r <- resc.dat$X
       y.r <- resc.dat$Y

       y0.r <- y.r
       if(is(distr, "DiscreteDistribution")){
          y0.r <- jitter(y0.r, factor = jit.fac)
       }else{
          iR <- .isReplicated(y0.r, jit.tol)
          if(any(iR)&&jit.fac>0)
             y0.r[iR] <- jitter(y0.r[iR], factor = jit.fac)
       }
       return(list(x=x.r,y=y.r, yj = y0.r))
}

.do.pts <- function(x,y,cxa,ca,pa, dots)
        do.call(points,args=c(list(x,y,cex=cxa,col=ca,pch=pa),dots))
.tx <- function(xa,ya,lb,cx,ca,ft,pos)
       text(x=xa,y=ya,labels=lb,cex=cx, col=ca, font=ft, pos=pos)

}
