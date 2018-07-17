#------------------------------------
#### utilities copied from package RobAStBase v.1.0  svn-rev 767
#------------------------------------

.evalListRec <- function(list0){ ## a list
    len <- length(list0)
    if(len==0L) return(list0)
    for(i in 1:len) {
        if(is.list(list0[[i]])){ list0[[i]] <- .evalListRec(list0[[i]])
           }else list0[[i]] <- eval(list0[[i]])
    }
    return(list0)
}


.makedotsLowLevel <- function(dots){
       dots$sub <- dots$xlab <- dots$ylab <- dots$main <- dots$type <- NULL
       dots$xlim <- dots$ylim <- dots$yaxt <- dots$axes <- dots$xaxt <- NULL
       dots$panel.last <- dots$panel.first <- dots$frame.plot <- dots$ann <-NULL
       dots$log <- dots$asp <- NULL
       return(dots)
}

.deleteDotsABLINE <- function(dots){
    dots$reg <- dots$a <- dots$b <- NULL
    dots$untf <- dots$h <- dots$v <- NULL
    dots
}

.deleteDotsTEXT <- function(dots){
   dots$labels <- dots$offset <- dots$vfont <- dots$pos <- dots$font <- NULL
   dots
}

.makedotsP <- function(dots){
    dots <- .makedotsLowLevel(dots)
    dots$lwd <- NULL
    dots$attr.pre <- NULL
    .deleteDotsABLINE(.deleteDotsTEXT(dots))
}

.SelectOrderData <- function(data, fct, which.lbs, which.Order, which.nonlbs = NULL){
   ## for data to be plotted in; performs two selections:
   ## on unordered (original) data (acc. to which.lbs)
   ## on data ordered acc. to fct a selection acc. to which.Order is done
   ## return value: list with elements
   #      data, the selected/thinned out data,
   #      y = fct(data)
   #      ind the indices of the selected data in the original data
   #      ind1 the indices of the data selected by which.lbs in the original data
     dimL <- !is.null(dim(data))
     d1  <- if(dimL) dim(data) else 1
     n   <- if(dimL) nrow(data) else length(data)
     ind <- 1:n

     ### function evaluation
     y <- if(dimL) apply(data, 1, fct) else sapply(data,fct)

#------------------------------------------------------------------------------
     ## selected data : data.t
#------------------------------------------------------------------------------

     ### first selection
     if(is.null(which.lbs)) which.lbs <- 1:n
     ## which.lbs0 is a logical of length the original data set, selecting
     ##     the remaining obs after first selection
     which.lbs0 <- ind %in% which.lbs
     # the remaining nb of obs after first selection
     n.s <- sum(which.lbs0)
     ## produce index for shown data after first selection
     ind.s <- ind[which.lbs0]
     ## function values after first selection
     y.s <- y[ind.s]

     ### ordering
     oN.s <- order(y.s)
     ## indices remaining after first selection ordered
     ##         from largest function value to smallest
     ind1.s <- rev(ind[oN.s])

     ### second selection
     ## selection of ordered
     if(is.null(which.Order))
          which.Order <- 1:n.s ## if no 2nd selection performed use all remaining obs.

     ## from ranks in remaining selection pick out those in which.order
     in.t <- (n.s+1)-which.Order
     in.t <- in.t[in.t>0]
     oN.t <-  oN.s[in.t] ## use largest ones in this order
     oN.t <- oN.t[!is.na(oN.t)]

     ## remaining number of observations after 2nd selection
     n.t <- length(oN.t)
     ## observations indices after 2nd selection
     ind.t <- ind.s[oN.t]
     ind.t <- ind.t[!is.na(ind.t)]
     ## function values after 2nd selection
     y.t <- y[ind.t]
     ## data after both selections
#     data.t <- if(dimL) data[ind.t,] else data[ind.t]
#     # if needed recast it to matrix/array
#     if(dimL) dim(data.t) <- c(n.t,d1[-1])
     data.t <- .SelectIndex(data,1,ind.t)

#------------------------------------------------------------------------------
     ## data not labelled: data.ns
#------------------------------------------------------------------------------
     if(is.null(which.nonlbs)) which.nonlbs <- 1:n
     #### non selected obs' indices after 1st selection
     ind.ns0 <- ind[!which.lbs0]
     #### non selected obs' indices in 2nd selection
     ind.nt <- if(length(oN.t)) ind.s[-oN.t] else numeric(0)
     #### non selected obs' in total is the union of both non-selected ones
     ind.ns1 <- unique(sort(c(ind.ns0, ind.nt)))
     ind.ns <- ind.ns1[ind.ns1 %in% which.nonlbs]
     ## number of non-selected obs'
     n.ns <- length(ind.ns)

#     which.lbns0 <-ind %in% ind.ns
#     which.lbnx <- rep(which.lbns0, length.out=length(data))

     ## non selected data
     data.ns <- .SelectIndex(data,1,ind.ns)
#     data.ns <- data[which.lbnx]
     # if needed recast it to matrix
#     if(dimL) dim(data.ns) <- c(n.ns,d1[-1])

     y.ns <- y[ind.ns]

     return(list(data=data.t, y=y.t, ind=ind.t, ind1=ind1.s, data.ns=data.ns, y.ns=y.ns, ind.ns = ind.ns))
}

.SelectIndex <- function(data,index,selection){
  dims <- dim(data)
  if(is.null(dims)) return(data[selection])
  datav <- data
  dimv <- dims
  if(index!=1){
     n <- length(dims)
     dims1 <- dims[-index]
     ind0 <- 1:n
     ind1 <- if(index<n) c((1:(index-1))+1,1,((index+1):n)) else c((1:(index-1))+1,1)
     ind2 <- c(index,ind0[-index])
     datav <- aperm(data,ind2)
     dimv <- dims[ind2]
  }
  len0 <- dimv[1]
  len1 <- prod(dimv[-1])
  lens <- length(selection)
  sel <- numeric(lens*len1)
  dimss <- dimv
  dimss[1] <- lens
  for(j in 1:len1)
     sel[1:lens+(j-1)*lens] <- selection+(j-1)*len0
  datas <- datav[sel]
  dim(datas) <- dimss
  if(index!=1){
     datas <- aperm(datas,ind1)
  }
  return(datas)
}

.plotRescaledAxis <- function(scaleX,scaleX.fct, scaleX.inv,
                              scaleY,scaleY.fct, scaleY.inv,
                              xlim, ylim, X, ypts = 400, n = 11,
                              finiteEndpoints = rep(FALSE,4),
                              x.ticks = NULL, y.ticks = NULL, withbox = TRUE){
# plots rescaled axes acc. to logicals scaleX, scaleY
# to this end uses trafos scaleX.fct with inverse scale.inv
# resp. scaleY.fct; it respects xlim and  ylim (given in orig. scale)
# return value: none
        if(scaleX){
           if(is.null(x.ticks)){
               x <- pretty(scaleX.inv(X))
               if(!is.null(xlim)) x <- pmax(x, xlim[1])
               if(!is.null(xlim)) x <- pmin(x, xlim[2])
               X <- .DistrCollapse(scaleX.fct(x),0*x)$supp
               x <- scaleX.inv(X)
               x <- x[is.finite(x)]
               x <- pretty(x,n=n)
               X <- .DistrCollapse(scaleX.fct(x),0*x)$supp
               x <- scaleX.inv(X)
               x <- x[is.finite(x)]
               x <- pretty(x,n=length(x))
               x[.isEqual01(x)&x<0.4] <- 0
               X <- scaleX.fct(x)
               xf <- prettyNum(x)
               i01 <- !.isEqual01(X)
               xf <- xf[i01]
               Xi <- X
               X <- X[i01]
               i0 <- any(!i01&Xi<0.5)
               i1 <- any(!i01&Xi>0.5)
               if(i0){ xf <- c(NA,xf); X <- c(0, X)}
               if(i1){ xf <- c(xf,NA); X <- c(X, 1)}
               axis(1,at=X,labels=xf)
               if(!finiteEndpoints[1]&i0) axis(1,at=0,labels=expression(-infinity))
               if(!finiteEndpoints[2]&i1) axis(1,at=1,labels=expression(infinity))
            }else{
               if(is.null(xlim)){ xlim <- c(-Inf,Inf)}else{
                  if(is.na(xlim[1])) xlim[1] <- -Inf
                  if(is.na(xlim[2])) xlim[2] <- Inf }
               x.ticks <- sort(unique(x.ticks[!is.na(x.ticks)]))
               xf <- pmin(pmax(x.ticks[is.finite(x.ticks)],xlim[1]),xlim[2])
               Xf <- scaleX.fct(xf)
               axis(1,at=Xf,labels=xf)
               if(-Inf %in% x.ticks) axis(1,at=0,labels=expression(-infinity))
               if(Inf %in% x.ticks)  axis(1,at=1,labels=expression(infinity))
            }
            if(withbox) box()
        }else{
            if(!is.null(x.ticks)){
               if(is.null(xlim)){ xlim <- c(-Inf,Inf)}else{
                  if(is.na(xlim[1])) xlim[1] <- -Inf
                  if(is.na(xlim[2])) xlim[2] <- Inf }
               x.ticks <- sort(unique(x.ticks[!is.na(x.ticks)]))
               xf <- pmin(pmax(x.ticks[is.finite(x.ticks)],xlim[1]),xlim[2])
               axis(1,at=xf,labels=xf)
               if(-Inf %in% x.ticks) axis(1,at=0,labels=expression(-infinity))
               if(Inf %in% x.ticks)  axis(1,at=1,labels=expression(infinity))
               if(withbox) box()
            }
        }
        if(scaleY){
           if(is.null(y.ticks)){
               Y0 <- if(!is.null(ylim)) max(0, scaleY.fct(ylim[1])) else 0
               Y1 <- if(!is.null(ylim)) min(1, scaleY.fct(ylim[2])) else 1
               Y <- seq(Y0,Y1, length=ypts)
               y <- pretty(scaleY.inv(Y),n=n)
               Y <- .DistrCollapse(scaleY.fct(y),0*y)$supp
               y <- scaleY.inv(Y)
               y <- y[is.finite(y)]
               y <- pretty(y,n=length(y))
               y[.isEqual01(y)&y<0.4] <- 0
               Y <- scaleX.fct(y)
               yf <- prettyNum(y)
               Y <- scaleY.fct(y)
               i01 <- !.isEqual01(Y)
               yf <- yf[i01]
               Yi <- Y
               Y <- Y[i01]
               i0 <- any(!i01&Yi<0.5)
               i1 <- any(!i01&Yi>0.5)
               if(i0){ yf <- c(NA,yf); Y <- c(0, Y)}
               if(i1){ yf <- c(yf,NA); Y <- c(Y, 1)}
               axis(2,at=Y,labels=yf)
               if(!finiteEndpoints[3]&i0) axis(2,at=0,labels=expression(-infinity))
               if(!finiteEndpoints[4]&i1) axis(2,at=1,labels=expression(infinity))
            }else{
               if(is.null(ylim)){ ylim <- c(-Inf,Inf)}else{
                  if(is.na(ylim[1])) ylim[1] <- -Inf
                  if(is.na(ylim[2])) ylim[2] <- Inf }
               y.ticks <- sort(unique(y.ticks[!is.na(y.ticks)]))
               yf <- pmin(pmax(y.ticks[is.finite(y.ticks)],ylim[1]),ylim[2])
               Yf <- scaleY.fct(yf)
               axis(2,at=Yf,labels=yf)
               if(-Inf %in% y.ticks) axis(2,at=0,labels=expression(-infinity))
               if(Inf %in% y.ticks)  axis(2,at=1,labels=expression(infinity))
            }
            if(withbox) box()
        }else{
            if(!is.null(y.ticks)){
               if(is.null(ylim)){ ylim <- c(-Inf,Inf)}else{
                  if(is.na(ylim[1])) ylim[1] <- -Inf
                  if(is.na(ylim[2])) ylim[2] <- Inf }
               y.ticks <- sort(unique(y.ticks[!is.na(y.ticks)]))
               yf <- pmin(pmax(y.ticks[is.finite(y.ticks)],ylim[1]),ylim[2])
               axis(2,at=yf,labels=yf)
               if(-Inf %in% y.ticks) axis(2,at=0,labels=expression(-infinity))
               if(Inf %in% y.ticks)  axis(2,at=1,labels=expression(infinity))
               if(withbox) box()
           }
        }
   return(invisible(NULL))
}

.rescalefct <- function(x, fct,
         scaleX = FALSE, scaleX.fct, scaleX.inv,
         scaleY = FALSE, scaleY.fct = pnorm,
         xlim, ylim, dots){

# if scaleX rescales x, if scaleY rescales fct(x);
# to this end uses trafos scaleX.fct with inverse scale.inv
# resp. scaleY.fct; it respects xlim and  ylim (given in orig. scale)
# thins out the scaled values if necessary and accordingly modifies
# slots xaxt, yaxt, axes of dots to indicate the new axes have to be drawn
#    paradigm small letters = orig. scale, capital letters = transformed scale
# return value: list with (thinned out) x and y, X and Y and modified dots
         if(length(x)==0) return(list(x=NULL,y=NULL,X=NULL,Y=NULL,scy=NA,dots=dots))
         if(!is.null(dots$log)){
             scaleX <- scaleX & !grepl("x", dots$log)
             scaleY <- scaleY & !grepl("y", dots$log)
         }
         X <- x
         wI <- 1:length(x)
         if(scaleX){
            if(!is.null(xlim)){
                   dots$xlim <- scaleX.fct(xlim)
                   x <- x[x>=xlim[1] & x<=xlim[2]]
            }
            Xo <- X <- scaleX.fct(x)
            X <- .DistrCollapse(X, 0*X)$supp
            wI <- sapply(X, function(uu){ w<- which(uu==Xo); if(length(w)>0) w[1] else NA})
            wI <- wI[!is.na(wI)]
            x <- scaleX.inv(X)
            dots$axes <- NULL
            dots$xaxt <- "n"
         }
         Y <- y <- if(is.function(fct)) fct(x) else fct[wI,1]
         scy <- if(is.function(fct)) NA else fct[wI,2]
         if(scaleY){
            Y <- scaleY.fct(y)
            if(!is.null(ylim)) dots$ylim <- scaleY.fct(ylim)
            dots$axes <- NULL
            dots$yaxt <- "n"
            }
         return(list(x=x,y=y,X=X,Y=Y,scy=scy,dots=dots))
}
