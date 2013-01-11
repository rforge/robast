.rescalefct <- function(x, fct, scaleX = FALSE, scaleX.fct, scaleX.inv,
         scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
         xlim, ylim, dots){
         X <- x
         if(scaleX){
            if(!is.null(xlim)){
                   dots$xlim <- scaleX.fct(xlim)
                   x <- x[x>=xlim[1] & x<=xlim[2]]
            }
            X <- scaleX.fct(x)
            X <- distr:::.DistrCollapse(X, 0*X)
            x <- scaleX.inv(X)
            dots$axes <- NULL
            dots$xaxt <- "n"
         }
         Y <- y <- sapply(fct,x)
         if(scaleY){
            Y <- scaleY.fct(y)
            if(!is.null(ylim)) dots$ylim <- scaleY.fct(ylim[,i])
            dots$axes <- NULL
            dots$yaxt <- "n"
            }
         return(list(x=x,y=y,X=X,Y=Y,dots=dots))
}

.plotRescaledAxis <- function(scaleX,scaleX.fct, scaleX.inv,
                              scaleY,scaleY.fct, scaleY.inv,
                              xlim, ylim, X, ypts = 400){
        if(scaleX){
               x <- pretty(scaleX.inv(X))
               if(!is.null(xlim)) x <- pmax(x, scaleY.fct(xlim[1]))
               if(!is.null(xlim)) x <- pmin(x, scaleY.fct(xlim[2]))
               X <- distr:::.DistrCollapse(scaleX.fct(x),0*x)
               x <- scaleX.inv(X)
               i01 <- !distr:::.isEqual01(X)
               x <- x[i01]
               Xi <- X
               X <- X[i01]
               i0 <- any(!i01&Xi<0.5)
               i1 <- any(!i01&Xi>0.5)
               if(i0){ x <- c(NA,x); X <- c(0, X)}
               if(i1){ x <- c(x,NA); X <- c(X, 1)}
               axis(1,at=X,labels=x)
               if(i0) axis(1,at=0,labels=expression(-infinity))
               if(i1) axis(1,at=1,labels=expression(infinity))
            }
        if(scaleY){
               Y0 <- if(!is.null(ylim)) max(0, scaleY.fct(ylim[1])) else 0
               Y1 <- if(!is.null(ylim)) min(1, scaleY.fct(ylim[2])) else 1
               Y <- seq(Y0,Y1, length=ypts)
               y <- pretty(scaleY.inv(Y))
               Y <- distr:::.DistrCollapse(scaleY.fct(y),0*y)
               y <- scaleY.inv(Y)
               i01 <- !distr:::.isEqual01(Y)
               y <- y[i01]
               Yi <- Y
               Y <- Y[i01]
               i0 <- any(!i01&Yi<0.5)
               i1 <- any(!i01&Yi>0.5)
               if(i0){ y <- c(NA,y); Y <- c(0, Y)}
               if(i1){ y <- c(y,NA); Y <- c(Y, 1)}
               axis(2,at=Y,labels=y)
               if(i0) axis(2,at=0,labels=expression(-infinity))
               if(i1) axis(2,at=1,labels=expression(infinity))
        }
   return(invisible(NULL))
}

.legendCoord <- function(x, scX, scX.fct, scY, scY.fct){
                if (is.character(x)) return(x)
                x1 <- if(scX) scX.fct(x[1]) else x[1]
                x2 <- if(scY) scY.fct(x[2]) else x[2]
                return(c(x1,x2))
            }
