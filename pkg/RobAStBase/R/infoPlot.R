setMethod("infoPlot", "IC",
    function(object, data = NULL,
             ..., withSweave = getdistrOption("withSweave"),
             col = par("col"), lwd = par("lwd"), lty, 
             colI = grey(0.5), lwdI = 0.7*par("lwd"), ltyI = "dotted",
             main = FALSE, inner = TRUE, sub = FALSE, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3], 
             with.automatic.grid = TRUE,
             with.legend = TRUE, legend = NULL, legend.bg = "white",
             legend.location = "bottomright", legend.cex = 0.8,
             x.vec = NULL, scaleX = FALSE, scaleX.fct, scaleX.inv,
             scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
             scaleN = 9, x.ticks = NULL, y.ticks = NULL,
             mfColRow = TRUE, to.draw.arg = NULL,
             cex.pts = 1, cex.pts.fun = NULL, col.pts = par("col"),
             pch.pts = 19,
             cex.npts = 1, cex.npts.fun = NULL, col.npts = grey(.5),
             pch.npts = 20,
             jitter.fac = 1, with.lab = FALSE, cex.lbs = 1, adj.lbs = c(0,0),
             col.lbs = col.pts, lab.pts = NULL,
             lab.font = NULL, alpha.trsp = NA,
             which.lbs = NULL, which.Order  = NULL, which.nonlbs = NULL,
             attr.pre = FALSE, return.Order = FALSE,
             ylab.abs = "absolute information", 
             ylab.rel= "relative information",
             withSubst = TRUE){

        args0 <- list(object = object, data = data, withSweave = withSweave,
             col = col, lwd = lwd,
             lty = if(!missing(lty)) lty else NULL,
             colI = colI, lwdI = lwdI,
             ltyI = ltyI, main = main, inner = inner, sub = sub,
             col.inner = col.inner, cex.inner = cex.inner,
             bmar = bmar, tmar = tmar, with.automatic.grid = with.automatic.grid,
             with.legend = with.legend, legend = legend, legend.bg = legend.bg,
             legend.location = legend.location, legend.cex = legend.cex,
             x.vec = x.vec, scaleX = scaleX,
             scaleX.fct = if(!missing(scaleX.fct)) scaleX.fct else NULL,
             scaleX.inv = if(!missing(scaleX.inv)) scaleX.inv else NULL,
             scaleY = scaleY, scaleY.fct = scaleY.fct,
             scaleY.inv = scaleY.inv, scaleN = scaleN, x.ticks = x.ticks,
             y.ticks = y.ticks, mfColRow = mfColRow, to.draw.arg = to.draw.arg,
             cex.pts = cex.pts, cex.pts.fun = cex.pts.fun, col.pts = col.pts,
             pch.pts = pch.pts, cex.npts = cex.npts, cex.npts.fun = cex.npts.fun,
             col.npts = col.npts, pch.npts = pch.npts,
             jitter.fac = jitter.fac, with.lab = with.lab,
             cex.lbs = cex.lbs, adj.lbs = adj.lbs,
             col.lbs = if(!missing(col.lbs)) col.lbs else if(!missing(col.pts)) col.pts else par("col"),
             lab.pts = lab.pts, lab.font = lab.font, alpha.trsp = alpha.trsp,
             which.lbs = which.lbs, which.Order  = which.Order,
             which.nonlbs = which.nonlbs, attr.pre = attr.pre,
             return.Order = return.Order, ylab.abs = ylab.abs, ylab.rel= ylab.rel,
             withSubst = withSubst)

        mc <- match.call(call = sys.call(sys.parent(1)))
        objectc <- mc$object
        dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
        plotInfo <- list(call = mc, dots=dots, args=args0)

        L2Fam <- eval(object@CallL2Fam)

        if(missing(scaleX.fct)){
           scaleX.fct <- p(L2Fam)
           scaleX.inv <- q.l(L2Fam)
        }


        dotsP <- dots
        dotsP$type <- dotsP$lty <- dotsP$col <- dotsP$lwd <- NULL
        dotsP$xlab <- dotsP$ylab <- NULL
        dotsP$axes <- NULL

        pF.0 <- expression({})
        if(!is.null(dots[["panel.first"]])){
            pF.0 <- .panel.mingle(dots,"panel.first")
        }
        pL.0 <- expression({})
        if(!is.null(dots[["panel.last"]])){
            pL.0 <- .panel.mingle(dots,"panel.last")
        }
        dotsP$panel.first <- NULL
        dotsP$panel.last <- NULL

        dotsLeg <- dotsT <- .makedotsLowLevel(dots)
        dotsT <- dotsL <- .makedotsLowLevel(dotsP)
        dotsT["main"] <- dotsT["cex.main"] <- dotsT["col.main"] <- NULL
        dotsT["line"] <- NULL

        dots.points <-   .makedotsPt(dots)

        withbox <- TRUE
        if(!is.null(dots[["withbox"]])) withbox <- dots[["withbox"]]
        dots["withbox"] <- NULL
        dots["type"] <- NULL
        xlab <- dots$xlab; if(is.null(xlab)) xlab <- "x"
        dots$xlab <- dots$ylab <- NULL

        trafO <- trafo(L2Fam@param)
        dimsA <- dims <- nrow(trafO)
        dimm <- ncol(trafO)
        
        if(missing(data)) data <- NULL

        to.draw <- .getToDraw(dims, trafO, L2Fam, to.draw.arg, "Abs")
        to.draw1 <- to.draw[to.draw>1]
        dims0 <- length(to.draw1)
        dims1 <- length(to.draw)
        nrows <- trunc(sqrt(dims0))
        ncols <- ceiling(dims0/nrows)
        in1to.draw <- (1%in%to.draw)

        yaxt0 <- xaxt0 <- rep("s",dims1)
        if(!is.null(dots$xaxt)){ xaxt1 <- eval(dots$xaxt); xaxt0 <- rep(xaxt1, length.out=dims1)}
        if(!is.null(dots$yaxt)){ yaxt1 <- eval(dots$yaxt); yaxt0 <- rep(yaxt1, length.out=dims1)}

        logArg <- NULL
        if(!is.null(dots[["log"]]))
            logArg <- rep(dots[["log"]], length.out=dims1)
        dotsP$log <- dots$log <- NULL

        dotsP0 <- vector("list",dims1)
        if(!is.null(dotsP)) for(i in 1:dims1) dotsP0[[i]] <- dotsP
        dotsP <- dotsP0

        for(i in 1:dims1){dotsP[[i]]$xaxt <- xaxt0[i];dotsP[[i]]$yaxt <- yaxt0[i]}

        if(!is.null(logArg))
            for(i in 1:dims1) dotsP[[i]]$log <- logArg[i]

        if(!is.null(x.ticks)){
           x.ticks <- .fillList(x.ticks, dims1)
           for(i in 1:dims1){
               if(!is.null(x.ticks[[i]]))
                   if(!is.null(logArg)) if(!grepl("x",logArg[i])) dotsP[[i]]$xaxt <- "n"
           }
        }
        if(!is.null(y.ticks)){
           y.ticks <- .fillList(y.ticks, dims1)
           for(i in 1:dims1){
               if(!is.null(y.ticks[[i]]))
                   if(!is.null(logArg)) if(!grepl("y",logArg[i])) dotsP[[i]]$yaxt <- "n"
           }
        }

        scaleX <- rep(scaleX, length.out=dims1)
        scaleY <- rep(scaleY, length.out=dims1)
        scaleX <- scaleX & (xaxt0!="n")
        scaleY <- scaleY & (yaxt0!="n")

        scaleX.fct <- .fillList(scaleX.fct, dims1)
        scaleX.inv <- .fillList(scaleX.inv, dims1)

        scaleY.fct <- .fillList(scaleY.fct, dims1)
        scaleY.inv <- .fillList(scaleY.inv, dims1)


        distr <- L2Fam@distribution
        if(!is(distr, "UnivariateDistribution") | is(distr, "CondDistribution"))
            stop("not yet implemented")

        xlim <- eval(dots$xlim)
        ylim <- eval(dots$ylim)
        .xylim <- .getXlimYlim(dots,dotsP, dims1, xlim, ylim)
           dots <- .xylim$dots; dotsP <- .xylim$dotsP
           xlim <- .xylim$xlim; ylim <- .xylim$ylim

        if(missing(x.vec)) x.vec <- NULL
        x.v.ret <- .getX.vec(distr, dims1, eval(dots$lty), x.vec, scaleX, scaleX.fct, scaleX.inv, .xylim$xm, .xylim$xM)
              lty <- x.v.ret$lty; plty <- x.v.ret$plty; x.vec <- x.v.ret$x.vec


        if(with.legend){
          if(missing(legend.location)){
             legend.location <- .fillList("topright", dims1   )
             if (in1to.draw) legend.location[[1]] <-  "bottomright"
          }else{
             legend.location <- as.list(legend.location)
             legend.location <- .fillList(legend.location, dims1  )
          }
          if(is.null(legend)){
             legend <- vector("list",dims1)
             legend <- .fillList(list(as.list(c("class. opt. IC", objectc))),
                                                 dims1)
          }
        }



         trafo <- trafo(L2Fam@param)
            
      .pT <- .prepareTitles(withSubst,
              presubArg2 = c("%C", "%D", "%A"),
              presubArg3 = c(as.character(class(object)[1]),
                             as.character(date()),
                             as.character(deparse(objectc))),
              dots,
              mainText = gettextf("Information Plot for IC %%A"),
              L2Fam, inner, dims0, dims, to.draw, trafO, L2Fam, type = "info", bmar, tmar)

        dots <- .pT$dots; main <- .pT$main; mainL <- .pT$mainL; lineT <- .pT$lineT
        sub <- .pT$sub; subL <- .pT$subL; bmar <- .pT$bmar; tmar <- .pT$tmar;
        innerT <- .pT$innerT; innerL <- .pT$innerL; .mpresubs <- .pT$.mpresubs


            QFc <- diag(dimsA)
            if(is(object,"ContIC") & dimsA>1 )
               {if (is(normtype(object),"QFNorm")) QFc <- QuadForm(normtype(object))
                QFc0 <- solve( trafo %*% solve(L2Fam@FisherInfo) %*% t(trafo ))
                if (is(normtype(object),"SelfNorm")|is(normtype(object),"InfoNorm")) 
                    QFc <- QFc0
               }

            absInfoEval <- function(x,y, withNorm = FALSE){
                       if(length(x)){
                       aI <- .msapply(x, y@Map[[1]])
                       if(withNorm) aI <- aI / max(aI)
                       }else aI <- NULL
                       return(aI)
            }


            QFc.5 <- sqrt(PosSemDefSymmMatrix(QFc))

            classIC <- as(trafo %*% solve(L2Fam@FisherInfo) %*% L2Fam@L2deriv, "EuclRandVariable")
            absInfoClass.f <- t(classIC) %*% QFc %*% classIC
#            absInfoClass <- absInfoEval(x.vec, absInfoClass.f)

            QF <- diag(dimsA)
            if(is(object,"ContIC") & dimsA>1 )
               {if (is(normtype(object),"QFNorm")) QF <- QuadForm(normtype(object))}
            QF.5 <- sqrt(PosSemDefSymmMatrix(QF))

            IC1 <- as(diag(dimsA) %*% object@Curve, "EuclRandVariable")
            absInfo.f <- t(IC1) %*% QF %*% IC1
#            absInfo <- absInfoEval(x.vec, absInfo.f)

            plotInfo$absInfoEval <- absInfoEval
            plotInfo$absInfoClass.f <- function(x) absInfoEval(x,absInfoClass.f)
            plotInfo$absInfo.f <- function(x) absInfoEval(x,absInfo.f)

            w0 <- getOption("warn")
            options(warn = -1)
            on.exit(options(warn = w0))
#            opar$cin <- opar$cra <- opar$csi <- opar$cxy <-  opar$din <- NULL
            opar <- par(no.readonly = TRUE)
            on.exit(par(opar))
            omar <- par("mar")


            wmar <- FALSE
            if(!missing(bmar)||!missing(tmar)){
                 lpA <- max(dims1,1)
                 parArgsL <- vector("list",lpA)
                 wmar <- TRUE
                 if(missing(bmar)) bmar <- omar[1]
                 if(missing(tmar)) bmar <- omar[3]
                 bmar <- rep(bmar, length.out=lpA)
                 tmar <- rep(tmar, length.out=lpA)
                 for( i in 1:lpA)
                      parArgsL[[i]] <- list(mar = c(bmar[i],omar[2],tmar[i],omar[4]))
                 plotInfo$parArgsL <- parArgsL
            }

           .pFL <- .preparePanelFirstLast(with.automatic.grid , dims1, pF.0, pL.0,
                             logArg, scaleX, scaleY, x.ticks, y.ticks,
                             scaleX.fct, scaleY.fct)

            pF <- .pFL$pF
            pF.abs <- if(in1to.draw) pF[[1]] else NULL
            pF.rel <- if(in1to.draw) pF[-1] else pF
            if(is.list(pL.0)){
               pL.abs <- if(in1to.draw) pL.0 else NULL
               pL.rel <- if(in1to.draw) pL.0 else pL
            }else{pL.rel <- pL.abs <- pL <- pL.0 }

            plotInfo$to.draw <- to.draw
            plotInfo$panelFirst <- pF
            plotInfo$panelLast <- pL
            plotInfo$gridS <- .pFL$gridS


            wmar <- FALSE
            if(!missing(bmar)||!missing(tmar)){
               wmar <- TRUE
               bmar <-
               nmar <- c(bmar[i],omar[2],tmar[i],omar[4])
            }
            trEnv <- new.env()

            if(!is.null(data)){

               n <- if(!is.null(dim(data))) nrow(data) else length(data)

               lab.pts <- if(is.null(lab.pts)) paste(1:n) else rep(lab.pts,length.out=n)

               if(!is.null(cex.pts.fun)){
                   cex.pts.fun <- .fillList(cex.pts.fun, (dims1)*2)
               }
               if(!is.null(cex.npts.fun)){
                   cex.npts.fun <- .fillList(cex.npts.fun, (dims1)*2)
               }

               if(missing(adj.lbs)) adj.lbs <- c(0,0)
               if(!is.array(adj.lbs) ||
                 (is.array(adj.lbs)&&!all.equal(dim(adj.lbs),c(2,2,dims1)))){
                  adj.lbs <- array(rep(adj.lbs, length.out= 2*dims1*2),
                                     dim=c(2,2,dims1))
               }
               adjC.lbs <- matrix(adj.lbs[,2,],nrow=2,ncol=dims1)
               adj.lbs <- matrix(adj.lbs[,1,],nrow=2,ncol=dims1)


               if(attr.pre){
                 if(missing(pch.pts)) pch.pts <- 1
                 if(!is.matrix(pch.pts))
                     pch.pts <- t(matrix(rep(pch.pts, length.out= 2*n),2,n))

                 if(missing(col.pts)) col.pts <- c(col, colI)
                 if(!is.matrix(col.pts))
                    col.pts <- t(matrix(rep(col.pts, length.out= 2*n),2,n))

                 if(missing(cex.pts)) cex.pts <- 1
                 if(!is.matrix(cex.pts))
                    cex.pts <- matrix(rep(cex.pts, length.out= 2*n),n,2)

                 if(missing(cex.lbs)) cex.lbs <- 1
                 if(!is.array(cex.lbs) ||
                   (is.array(cex.lbs)&&!all.equal(dim(cex.lbs),c(n,2,dims1)))){
                    cex.lbs <- array(rep(cex.lbs, length.out= n*dims1*2),
                                     dim=c(n,2,dims1))
                   }

                 if(missing(col.lbs)) col.lbs <- col.pts
                 if(!is.matrix(col.lbs))
                    col.lbs <- t(matrix(rep(col.lbs, length.out= 2*n),2,n))

                 }

               sel <- .SelectOrderData(data, function(x)absInfoEval(x,absInfo.f),
                                       which.lbs, which.Order, which.nonlbs)
               sel.C <- .SelectOrderData(data, function(x)absInfoEval(x,absInfoClass.f),
                                       which.lbs, which.Order, which.nonlbs)

               plotInfo$sel <- sel
               plotInfo$sel.C <- sel.C

               i.d <- sel$ind
               i.dC <- sel.C$ind
               i0.d <- sel$ind1
               i0.dC <- sel.C$ind1
               y.d <- sel$y
               y.dC <- sel.C$y
               x.d <- sel$data
               x.dC <- sel.C$data
               n.s <- length(i.d)

               i.d.ns <- sel$ind.ns
               i.dC.ns <- sel.C$ind.ns
               y.d.ns <- sel$y.ns
               y.dC.ns <- sel.C$y.ns
               x.d.ns <- sel$data.ns
               x.dC.ns <- sel.C$data.ns
               n.ns <- length(i.d.ns)

            selAlly <- c(sel$y,sel.C$y)
            selAlly.n <- c(sel$y.ns,sel.C$y.ns)


            plotInfo$IC <- i0.d
            plotInfo$IC.class <- i0.dC

            labC.pts <-  lab.pts[sel.C$ind]
            lab.pts <-  lab.pts[sel$ind]

            if(attr.pre){
               col0.pts <- col.pts[sel$ind,1]
               colC.pts <- col.pts[sel.C$ind,2]
               col0.npts <- col.npts[sel$ind.ns,1]
               colC.npts <- col.npts[sel.C$ind.ns,2]
               col.pts <- col0.pts; col.npts <- col0.npts

               pch0.pts <- pch.pts[sel$ind,1]
               pchC.pts <- pch.pts[sel.C$ind,2]
               pch0.npts <- pch.npts[sel$ind.ns,1]
               pchC.npts <- pch.npts[sel.C$ind.ns,2]
               pch.pts <- pch0.pts; pch.npts <- pch0.npts

               cex0.pts <- cex.pts[sel$ind,1]
               cexC.pts <- cex.pts[sel.C$ind,1]
               cex0.npts <- cex.npts[sel$ind.ns,1]
               cexC.npts <- cex.npts[sel.C$ind.ns,2]
               cex.pts <- cex0.pts; cex.npts <- cex0.npts

               cex0.lbs <- matrix(cex.lbs[sel$ind,1,],nrow=n.s,ncol=dims1)
               cexC.lbs <- matrix(cex.lbs[sel.C$ind,2,],nrow=n.s,ncol=dims1)
               cex.lbs <- cex0.lbs

               col0.lbs <- col.lbs[sel$ind,1]
               colC.lbs <- col.lbs[sel$ind,2]
               col.lbs <- col0.lbs
          }else{
               if(missing(pch.pts)) pch.pts <- 1
               if(!is.matrix(pch.pts))
                   pch.pts <- t(matrix(rep(pch.pts, length.out= 2*n.s),2,n.s))
               pchC.pts <- pch.pts[,2]
               pch.pts <- pch.pts[,1]

               if(missing(pch.npts)) pch.npts <- 2
               if(!is.matrix(pch.npts))
                   pch.npts <- t(matrix(rep(pch.npts, length.out= 2*n.ns),2,n.ns))
               pchC.npts <- pch.npts[,2]
               pch.npts <- pch.npts[,1]

               if(missing(col.pts)) col.pts <- c(col, colI)
               if(!is.matrix(col.pts))
                  col.pts <- t(matrix(rep(col.pts, length.out= 2*n.s),2,n.s))
               colC.pts <- col.pts[,2]
               col.pts <- col.pts[,1]

               if(missing(col.npts)) col.pts <- c(col, colI)
               if(!is.matrix(col.npts))
                  col.npts <- t(matrix(rep(col.npts, length.out= 2*n.ns),2,n.ns))
               colC.npts <- col.npts[,2]
               col.npts <- col.npts[,1]

               if(missing(cex.pts)) cex.pts <- 1
               if(!is.matrix(cex.pts))
                  cex.pts <- matrix(rep(cex.pts, length.out= 2*n.s),n.s,2)
               cexC.pts <- cex.pts[,2]
               cex.pts <- cex.pts[,1]

               if(missing(cex.npts)) cex.npts <- 1
               if(!is.matrix(cex.npts))
                  cex.npts <- matrix(rep(cex.npts, length.out= 2*n.ns),n.ns,2)
               cexC.npts <- cex.npts[,2]
               cex.npts <- cex.npts[,1]

               if(missing(cex.lbs)) cex.lbs <- 1
               if(!is.array(cex.lbs) ||
                   (is.array(cex.lbs)&&all.equal(dim(cex.lbs),c(n.s,2,dims1)))){
                    cex.lbs <- array(rep(cex.lbs, length.out= n.s*dims1*2),
                                     dim=c(n.s,2,dims1))
                   }
               cexC.lbs <- matrix(cex.lbs[,2,],nrow=n.s,ncol=dims1)
               cex.lbs <- matrix(cex.lbs[,1,],nrow=n.s,ncol=dims1)

               if(missing(col.lbs)) col.lbs <- col.pts
               if(!is.matrix(col.lbs))
                    col.lbs <- t(matrix(rep(col.lbs, length.out= 2*n.s),2,n.s))
               colC.lbs <- col.lbs[,2]
               col.lbs <- col.lbs[,1]
            }

               jitter.fac <- rep(jitter.fac, length.out=2)
               lab.font <- rep(lab.font, length.out=2)


               resc.dat <-.rescalefct(x.d, function(x) absInfoEval(x,absInfo.f),
                              scaleX[1], scaleX.fct[[1]], scaleX.inv[[1]],
                              scaleY[1], scaleY.fct[[1]], xlim[,1], ylim[,1], dotsP[[1]])
               resc.datC <-.rescalefct(x.dC, function(x) absInfoEval(x,absInfoClass.f),
                              scaleX[1], scaleX.fct[[1]], scaleX.inv[[1]],
                              scaleY[1], scaleY.fct[[1]], xlim[,1], ylim[,1], dotsP[[1]])
               resc.dat.ns <-.rescalefct(x.d.ns, function(x) absInfoEval(x,absInfo.f),
                              scaleX[1], scaleX.fct[[1]], scaleX.inv[[1]],
                              scaleY[1], scaleY.fct[[1]], xlim[,1], ylim[,1], dotsP[[1]])
               resc.datC.ns <-.rescalefct(x.dC.ns, function(x) absInfoEval(x,absInfoClass.f),
                              scaleX[1], scaleX.fct[[1]], scaleX.inv[[1]],
                              scaleY[1], scaleY.fct[[1]], xlim[,1], ylim[,1], dotsP[[1]])

               plotInfo$resc.dat.abs <- resc.dat
               plotInfo$resc.dat.abs.ns <- resc.dat.ns
               plotInfo$resc.datC.abs <- resc.datC
               plotInfo$resc.datC.abs.ns <- resc.datC.ns

               x.dr <- resc.dat$X
               x.dCr <- resc.datC$X
               y.dr <- resc.dat$Y
               y.dCr <- resc.datC$Y

               x.dr.ns <- resc.dat.ns$X
               x.dCr.ns <- resc.datC.ns$X
               y.dr.ns <- resc.dat.ns$Y
               y.dCr.ns <- resc.datC.ns$Y

#               lab.pts <- if(is.null(lab.pts))
#                               cbind(i.d, i.dC)
#                          else cbind(lab.pts[i.d],lab.pts[i.dC])


               dots.points <-   .makedotsPt(dots)

               do.pts <- function(x,y,cxa,ca,pa){
                 if(length(x)>0)
                    do.call(points,args=c(list(x,y,cex=cxa,col=ca,pch=pa),
                            dots.points))}
               tx <- function(xa,ya,lb,cx,ca,ad){
                 if(length(xa)>0)
                    if(!is.null(lb)) text(x=xa,y=ya,labels=lb,cex=cx, col=ca, adj=ad)
               }
               alp.v <- rep(alpha.trsp, length.out = dims1)


               pL.abs <- substitute({

                   pI <- get("plotInfo", envir = trEnv0)

                   if(length(ICy0)){
                      ICy0r1 <- ICy0r
                      ICy0cr1 <- ICy0cr
                      if(is(distr, "DiscreteDistribution")){
                        ICy0r1 <- jitter(ICy0r1, factor = jitter.fac0[1])
                        ICy0cr1 <- jitter(ICy0cr1, factor = jitter.fac0[2])
                      }
                      c1fun <- if(is.null(cexfun)) NULL else cexfun[[1]]
                      c2fun <- if(is.null(cexfun)) NULL else cexfun[[2]]
                      f1 <- .cexscale(ICy0,ICy0c,cex=cex0, fun = c1fun)
                      f1c <- .cexscale(ICy0c,ICy0,cex=cex0C, fun = c2fun)
                      col.pts <- if(!is.na(al0)) .msapply(col0,
                                 addAlphTrsp2col, alpha=al0) else col0
                      colC.pts <- if(!is.na(al0)) .msapply(col0C,
                                 addAlphTrsp2col, alpha=al0) else col0C

                      pI$doPtsAbs <- list(x = x0, y = ICy0r1, cex = f1,
                                             col = col.pts, pch = pch0)
                      pI$doPtsAbsC <- list(x = x0c, y = ICy0cr1, cex = f1c,
                                             col = colC.pts, pch = pch0C)
                      do.pts(x0, ICy0r1, f1,col.pts,pch0)
                      do.pts(x0c, ICy0cr1, f1c,colC.pts,pch0C)

                      if(with.lab0){
                         tx(x0, ICy0r1, lab.pts0, cex.lbs0, col.lbs0, adj.lbs0)
                         tx(x0c, ICy0cr1, labC.pts0, cexC.lbs0, colC.lbs0, adjC.lbs0)
                         pI$doLabsAbs <- list(x = x0, y = ICy0r1, adj = adj.lbs0,
                                         lab = lab.pts0, cex = cex.lbs0, col= col.lbs0)
                         pI$doLabsCAbs <- list(x = x0c, y = ICy0cr1, adj = adjC.lbs0,
                                         lab = labC.pts0, cex = cexC.lbs0, col= colC.lbs0)
                      }
                   }
                   if(length(ICy0.ns)){
                      ICy0r1.ns <- ICy0r.ns
                      ICy0cr1.ns <- ICy0cr.ns
                      if(is(distr, "DiscreteDistribution")){
                         ICy0r1.ns <- jitter(ICy0r1.ns, factor = jitter.fac0[1])
                         ICy0cr1.ns <- jitter(ICy0cr1.ns, factor = jitter.fac0[2])
                      }
                      c1fun.ns <- if(is.null(cexfun.ns)) NULL else cexfun.ns[[1]]
                      c2fun.ns <- if(is.null(cexfun.ns)) NULL else cexfun.ns[[2]]
                      f1.ns <- .cexscale(ICy0.ns,ICy0c.ns,cex=cex0.ns, fun = c1fun.ns)
                      f1c.ns <- .cexscale(ICy0c.ns,ICy0.ns,cex=cex0C.ns, fun = c2fun.ns)
                      col.npts <- if(!is.na(al0)) .msapply(col0.ns,
                                  addAlphTrsp2col, alpha=al0) else col0.ns
                      colC.npts <- if(!is.na(al0)) .msapply(col0C.ns,
                                   addAlphTrsp2col, alpha=al0) else col0C.ns

                      pI$doPtsAbs.ns <- list(x = x0.ns, y = ICy0r1.ns, cex = f1.ns,
                                             col = col.npts, pch = pch0.ns)
                      pI$doPtsAbsC.ns <- list(x = x0c.ns, y = ICy0cr1.ns, cex = f1c.ns,
                                             col = colC.npts, pch = pch0C.ns)
                      do.pts(x0.ns, ICy0r1.ns, f1.ns,col.npts,pch0.ns)
                      do.pts(x0c.ns, ICy0cr1.ns, f1c.ns,colC.npts,pch0C.ns)
                   }

                   assign("plotInfo", pI, envir = trEnv0)
                   pL0
                   }, list(ICy0 = y.d,  ICy0c = y.dC,
                           ICy0r = y.dr, ICy0cr = y.dCr,
                           ICy0c.ns = y.dC.ns, ICy0.ns = y.d.ns,
                           ICy0r.ns = y.dr.ns, ICy0cr.ns = y.dCr.ns,
                           pL0 = pL.abs,
                           x0 = x.dr, x0c = x.dCr,
                           x0.ns = x.dr.ns, x0c.ns = x.dCr.ns,
                           al0 = alp.v[1],
                           cex0 = cex.pts,
                           pch0 = pch.pts,
                           col0 = col.pts,
                           cex0.ns = cex.npts,
                           pch0.ns = pch.npts,
                           col0.ns = col.npts,
                           cex0C = cexC.pts,
                           pch0C = pchC.pts,
                           col0C = colC.pts,
                           cex0C.ns = cexC.npts,
                           pch0C.ns = pchC.npts,
                           col0C.ns = colC.npts,
                           lab.pts0 = lab.pts,
                           labC.pts0 = labC.pts,
                           with.lab0 = with.lab, n0 = n,
                           jitter.fac0 = jitter.fac, cexfun = cex.pts.fun,
                           cexfun.ns = cex.npts.fun,
                           cex.lbs0 = cex.lbs[,1],
                           cexC.lbs0 = cexC.lbs[,1],
                           adj.lbs0 = adj.lbs[,1],
                           adjC.lbs0 = adjC.lbs[,1],
                           col.lbs0 = col.lbs,
                           colC.lbs0 = colC.lbs,
                           trEnv0 = trEnv)
                           )

               pL.rel <- substitute({

                   pI <- get("plotInfo", envir = trEnv0)

                   dotsP0 <- dotsP[[i1]]

                   if(length(x0)){

                      y0.vec <- .msapply(x0,  IC1.i.5@Map[[indi]])
                      if(!is.null(y0.vec)) y0.vec <- y0.vec^2/ICy0
                      y0c.vec <- .msapply(x0c, classIC.i.5@Map[[indi]])
                      if(!is.null(y0c.vec)) y0c.vec <- y0c.vec^2/ICy0c

                      if(is(distr, "DiscreteDistribution")){
                         if(length(y0.vec)) y0.vec <- jitter(y0.vec, factor = jitter.fac0[1])
                         if(length(y0c.vec)) y0c.vec <- jitter(y0c.vec, factor = jitter.fac0[2])
                      }
                      resc.rel <- .rescalefct(x0, cbind(y0.vec,ICy0),
                              scaleX[i1], scaleX.fct[[i1]], scaleX.inv[[i1]],
                              scaleY[i1], scaleY.fct[[i1]], xlim[,i1], ylim[,i1], dotsP0)
                      resc.rel.c <- .rescalefct(x0c, cbind(y0c.vec,ICy0c),
                              scaleX[i1], scaleX.fct[[i1]], scaleX.inv[[i1]],
                              scaleY[i1], scaleY.fct[[i1]], xlim[,i1], ylim[,i1], dotsP0)

                      pI$resc.dat.rel[[i]] <- resc.rel
                      pI$resc.datC.rel[[i]] <- resc.rel.c

                      c1fun <- if(is.null(cexfun)) NULL else cexfun[[(i1-1)*2+1]]
                      c2fun <- if(is.null(cexfun)) NULL else cexfun[[(i1-1)*2+2]]
                      f1 <- .cexscale(resc.rel$scy,resc.rel.c$scy,cex=cex0, fun=c1fun)
                      f1c <- .cexscale(resc.rel.c$scy,resc.rel$scy,cex=cex0C, fun=c2fun)

                      col.pts <- if(!is.na(al0[i1])) .msapply(col0,
                              addAlphTrsp2col, alpha=al0[i1]) else col0
                      colC.pts <- if(!is.na(al0[i1])) .msapply(col0C,
                              addAlphTrsp2col, alpha=al0[i1]) else col0C

                      pI$doPtsRel[[i]] <- list(x = resc.rel$X, y = resc.rel$Y, cex = f1,
                                             col = col.pts, pch = pch0)
                      pI$doPtsRelC[[i]] <- list(x = resc.rel.c$X, y = resc.rel.c$Y, cex = f1c,
                                             col = colC.pts, pch = pch0C)
                      do.pts(resc.rel$X, resc.rel$Y, f1,col.pts,pch0)
                      do.pts(resc.rel.c$X, resc.rel.c$Y, f1c,colC.pts,pch0C)

                      if(with.lab0){
                        cexl <- cex.lbs0[,i1]; cexlC <- cexC.lbs0[,i1]
                        adjl <- adj.lbs0[,i1]; adjlC <- adjC.lbs0[,i1]
                        tx(resc.rel$X, resc.rel$Y, lab.pts0, cexl, col.lbs0, adjl)
                        tx(resc.rel.c$X, resc.rel.c$Y, labC.pts0, cexlC, colC.lbs0, adjlC)
                        pI$doLabsRel[[i]] <- list(x = resc.rel$X, y = resc.rel$Y,
                                         lab = lab.pts0, cex = cexl, col= col.lbs0, adj=adjl)
                        pI$doLabsCRel[[i]] <- list(x = resc.rel.c$X, y = resc.rel.c$Y,
                                         lab = labC.pts0, cex = cexlC, col= colC.lbs0, adj=adjl)
                      }
                   }
                   if(length(x0.ns)){

                      y0.vec.ns <- .msapply(x0.ns,  IC1.i.5@Map[[indi]])
                      if(!is.null(y0.vec.ns)) y0.vec.ns <- y0.vec.ns^2/ICy0.ns
                      y0c.vec.ns <- .msapply(x0c.ns, classIC.i.5@Map[[indi]])
                      if(!is.null(y0c.vec.ns)) y0c.vec.ns <- y0c.vec.ns^2/ICy0c.ns

                      if(is(distr, "DiscreteDistribution")){
                        if(length(y0.vec.ns)) y0.vec.ns <- jitter(y0.vec.ns, factor = jitter.fac0[1])
                        if(length(y0c.vec.ns)) y0c.vec.ns <- jitter(y0c.vec.ns, factor = jitter.fac0[2])
                      }

                      resc.rel.ns <- .rescalefct(x0.ns, cbind(y0.vec.ns,ICy0.ns),
                              scaleX[i1], scaleX.fct[[i1]], scaleX.inv[[i1]],
                              FALSE, scaleY.fct[[i1]], dots$xlim, dots$ylim, dotsP0)
                      resc.rel.c.ns <- .rescalefct(x0c.ns, cbind(y0c.vec.ns,ICy0c.ns),
                              scaleX[i1], scaleX.fct[[i1]], scaleX.inv[[i1]],
                              FALSE, scaleY.fct[[i1]], dots$xlim, dots$ylim, dotsP0)

                      pI$resc.dat.rel.ns[[i]] <- resc.rel.ns
                      pI$resc.datC.rel.ns[[i]] <- resc.rel.c.ns

                      c1fun.ns <- if(is.null(cexfun.ns)) NULL else cexfun.ns[[(i1-1)*2+1]]
                      c2fun.ns <- if(is.null(cexfun.ns)) NULL else cexfun.ns[[(i1-1)*2+2]]
                      f1.ns <- .cexscale(resc.rel.ns$scy,resc.rel.c.ns$scy,cex=cex0.ns, fun = c1fun.ns)
                      f1c.ns <- .cexscale(resc.rel.c.ns$scy,resc.rel.ns$scy,cex=cex0C.ns, fun = c2fun.ns)

                      col.npts <- if(!is.na(al0[i1])) .msapply(col0.ns,
                              addAlphTrsp2col, alpha=al0[i1]) else col0.ns
                      colC.npts <- if(!is.na(al0[i1])) .msapply(col0C.ns,
                              addAlphTrsp2col, alpha=al0[i1]) else col0C.ns

                      pI$doPtsRel.ns[[i]] <- list(x = resc.rel.ns$X, y = resc.rel.ns$Y, cex = f1.ns,
                                             col = col.npts, pch = pch0.ns)
                      pI$doPtsRelC.ns[[i]] <- list(x = resc.rel.c.ns$X, y = resc.rel.c.ns$Y, cex = f1c.ns,
                                             col = colC.npts, pch = pch0C.ns)
                      do.pts(resc.rel.ns$X, resc.rel.ns$Y, f1.ns,col.npts,pch0.ns)
                      do.pts(resc.rel.c.ns$X, resc.rel.c.ns$Y, f1c.ns,colC.npts,pch0C.ns)
                   }

                   assign("plotInfo", pI, envir = trEnv0)
                   pL0
                   }, list(ICy0 = y.d, ICy0c = y.dC,
                           ICy0c.ns = y.dC.ns, ICy0.ns = y.d.ns,
                           pL0 = pL.rel,
                           x0 = x.d, x0c = x.dC,
                           x0.ns = x.d.ns, x0c.ns = x.dC.ns,
                           cex0 = cex.pts,
                           pch0 = pch.pts,
                           col0 = col.pts,
                           cex0.ns = cex.npts,
                           pch0.ns = pch.npts,
                           col0.ns = col.npts,
                           cex0C = cexC.pts,
                           pch0C = pchC.pts,
                           col0C = colC.pts,
                           cex0C.ns = cexC.npts,
                           pch0C.ns = pchC.npts,
                           col0C.ns = colC.npts,
                           lab.pts0 = lab.pts,
                           labC.pts0 = labC.pts,
                           with.lab0 = with.lab, n0 = n, al0 = alp.v,
                           jitter.fac0 = jitter.fac, cexfun = cex.pts.fun,
                           cexfun.ns = cex.npts.fun,
                           cex.lbs0 = cex.lbs,
                           cexC.lbs0 = cexC.lbs,
                           adj.lbs0 = adj.lbs,
                           adjC.lbs0 = adjC.lbs,
                           col.lbs0 = col.lbs,
                           colC.lbs0 = colC.lbs,
                           trEnv0 = trEnv)
                  )

            }
            fac.leg <- if(dims0>1) 3/4 else .75/.8

            if(1 %in% to.draw){
               indi <- 1
               resc <-.rescalefct(x.vec[[1]], function(x) absInfoEval(x,absInfo.f),
                              scaleX[1], scaleX.fct[[1]], scaleX.inv[[1]],
                              scaleY[1], scaleY.fct[[1]], xlim[,1], ylim[,1], dotsP[[1]])
               resc.C <-.rescalefct(x.vec[[1]], function(x) absInfoEval(x,absInfoClass.f),
                              scaleX[1], scaleX.fct[[1]], scaleX.inv[[1]],
                              scaleY[1], scaleY.fct[[1]], xlim[,1], ylim[,1], dotsP[[1]])


               plotInfo$resc.abs <- resc
               plotInfo$resc.C.abs <- resc.C

               dotsP[[1]] <- resc$dots

               if(wmar) do.call(par, args = parArgsL[[1]])

               finiteEndpoints <- rep(FALSE,4)
               if(scaleX[1]){
                  finiteEndpoints[1] <- is.finite(scaleX.inv[[1]](min(resc.C$X, xlim[1,1],na.rm=TRUE)))
                  finiteEndpoints[2] <- is.finite(scaleX.inv[[1]](max(resc.C$X, xlim[2,1],na.rm=TRUE)))
               }
               if(scaleY[1]){
                  finiteEndpoints[3] <- is.finite(scaleY.inv[[1]](min(resc.C$Y, ylim[1,1],na.rm=TRUE)))
                  finiteEndpoints[4] <- is.finite(scaleY.inv[[1]](max(resc.C$Y, ylim[2,1],na.rm=TRUE)))
               }


               plotInfo$absPlotArgs <- c(list(resc.C$X, resc.C$Y, type = plty,
                   lty = ltyI, col = colI, lwd = lwdI,
                   xlab = .mpresubs(xlab), ylab = .mpresubs(ylab.abs), panel.last = pL.abs,
                   panel.first = pF.abs),
                   dotsP[[1]])
               assign("plotInfo", plotInfo, envir = trEnv)
               do.call(plot, args=c(list(resc.C$X, resc.C$Y, type = plty,
                   lty = ltyI, col = colI, lwd = lwdI,
                   xlab = .mpresubs(xlab), ylab = .mpresubs(ylab.abs), panel.last = pL.abs,
                   panel.first = pF.abs),
                   dotsP[[1]]))
               plotInfo <- get("plotInfo", envir = trEnv)
               plotInfo$absPlotUsr <- par("usr")

               do.call(lines, args=c(list(resc$X, resc$Y, type = plty,
                       lty = lty, lwd = lwd, col = col), dotsL))
               plotInfo$absPlotCArgs <- c(list(resc$X, resc$Y, type = plty,
                       lty = lty, lwd = lwd, col = col), dotsL)

               x.ticks0 <- if(xaxt0[1]!="n") x.ticks[[1]] else NULL
               y.ticks0 <- if(yaxt0[1]!="n") y.ticks[[1]] else NULL


               .plotRescaledAxis(scaleX[1], scaleX.fct[[1]], scaleX.inv[[1]],
                              scaleY[1],scaleY.fct[[1]], scaleY.inv[[1]],
                              xlim[,1], ylim[,1], resc$X, ypts = 400,
                              n = scaleN, x.ticks = x.ticks0,
                              y.ticks = y.ticks0, withbox = withbox)
               plotInfo$absAxis <- list(scaleX[1], scaleX.fct[[1]], scaleX.inv[[1]],
                              scaleY[1],scaleY.fct[[1]], scaleY.inv[[1]],
                              xlim[,1], ylim[,1], resc$X, ypts = 400,
                              n = scaleN, x.ticks = x.ticks0,
                              y.ticks = y.ticks0, withbox = withbox)

               if(with.legend){
                 legend(.legendCoord(legend.location[[1]], scaleX[1], scaleX.fct[[1]],
                        scaleY[1], scaleY.fct[[1]]), legend = legend[[1]], bg = legend.bg,
                     lty = c(ltyI, lty), col = c(colI, col), 
                     lwd = c(lwdI, lwd), cex = legend.cex*fac.leg)
                 plotInfo$absLegend <- list(.legendCoord(legend.location[[1]],
                        scaleX[1], scaleX.fct[[1]], scaleY[1], scaleY.fct[[1]]),
                        legend = legend[[1]], bg = legend.bg,
                     lty = c(ltyI, lty), col = c(colI, col),
                     lwd = c(lwdI, lwd), cex = legend.cex*fac.leg)
               }

               if(innerL){
                  do.call(title, args=c(list(main = innerT[[1]]),  dotsT,
                          line = lineT, cex.main = cex.inner, col.main = col.inner))
                  plotInfo$absTitle <- c(list(main = innerT[[1]]),  dotsT,
                          line = lineT, cex.main = cex.inner, col.main = col.inner)
               }
            }

            if(dims > 1 && length(to.draw[to.draw!=1])>0){
                nrows <- trunc(sqrt(dims0))
                ncols <- ceiling(dims0/nrows)
#                if (!withSweave||!mfColRow)
#                     dN <- substitute({devNew()}) else substitute({})

                IC1.i.5 <- QF.5%*%IC1
                classIC.i.5 <- QFc.5%*%classIC

                plotInfo$relInfo.f <- function(x,i){
                                         den <- sapply(x,IC1.i.5@Map[[i]])
                                         nom <- absInfoEval(x,absInfo.f)
                                         den^2/nom}
                plotInfo$relInfoClass.f <- function(x,i){
                                         den <- sapply(x,classIC.i.5@Map[[i]])
                                         nom <- absInfoEval(x,absInfoClass.f)
                                         den^2/nom}

                if(!is.null(data)){
                    plotInfo$resc.dat.rel <- plotInfo$resc.datC.rel <- vector("list", dims0)
                    plotInfo$resc.dat.rel.ns <- plotInfo$resc.datC.rel.ns <- vector("list", dims0)
                }
                plotInfo$relPlotUsr <- plotInfo$par.rel <- vector("list", dims0)
                plotInfo$relPlotArgs <- plotInfo$relPlotCArgs <- vector("list", dims0)
                plotInfo$relY <- plotInfo$relYc <-plotInfo$relAxis <- vector("list", dims0)
                plotInfo$relLegend <- plotInfo$relTitle <- vector("list", dims0)
                plotInfo$doLabsRel <- plotInfo$doLabsCRel <- vector("list", dims0)

                if(mfColRow){
                  if(!withSweave&&in1to.draw && length(dev.list())>0) devNew()
                  par(mfrow = c(nrows, ncols))
                  plotInfo$rel.mfrow <- c(nrows, ncols)
                }

                for(i in 1:dims0){
                    indi <- to.draw1[i]-1
                    i1 <- i + in1to.draw


                    y.vecn <- absInfoEval(x.vec[[i1]],absInfo.f)
                    y.vec1 <- .msapply(x.vec[[i1]], IC1.i.5@Map[[indi]])
                    if(!is.null(y.vec1))
                       y.vec1 <- y.vec1^2/y.vecn

                    y.vecnC <- absInfoEval(x.vec[[i1]],absInfoClass.f)
                    y.vec1C <- .msapply(x.vec[[i1]], classIC.i.5@Map[[indi]])
                    if(!is.null(y.vec1C))
                       y.vec1C <- y.vec1C^2/y.vecnC


                    resc <-.rescalefct(x.vec[[i1]], cbind(y.vec1,y.vecn),
                              scaleX[i1], scaleX.fct[[i1]], scaleX.inv[[i1]],
                              scaleY[i1], scaleY.fct[[i1]], xlim[,i1], ylim[,i1], dotsP[[i1]])
                    resc.C <-.rescalefct(x.vec[[i1]], cbind(y.vec1C,y.vecnC),
                              scaleX[i1], scaleX.fct[[i1]], scaleX.inv[[i1]],
                              scaleY[i1], scaleY.fct[[i1]], xlim[,i1], ylim[,i1], dotsP[[i1]])

                    plotInfo$relY[[i]] <- resc$Y
                    plotInfo$relYc[[i]] <- resc.C$Y

                    if(wmar) do.call(par, args = parArgsL[[i+in1to.draw]])

                    finiteEndpoints <- rep(FALSE,4)
                    if(scaleX[i1]){
                      finiteEndpoints[1] <- is.finite(scaleX.inv[[i1]](min(resc$X, xlim[1,i1],na.rm=TRUE)))
                      finiteEndpoints[2] <- is.finite(scaleX.inv[[i1]](max(resc$X, xlim[2,i1],na.rm=TRUE)))
                    }
                    if(scaleY[i1]){
                       finiteEndpoints[3] <- is.finite(scaleY.inv[[i1]](min(resc$Y, resc.C$Y, ylim[1,i1],na.rm=TRUE)))
                       finiteEndpoints[4] <- is.finite(scaleY.inv[[i1]](max(resc$Y, resc.C$Y, ylim[2,i1],na.rm=TRUE)))
                    }


                    plotInfo$relPlotArgs[[i]] <- c(list(resc$X, resc$Y, type = plty,
                                  lty = lty, xlab = .mpresubs(xlab), ylab = .mpresubs(ylab.rel),
                                  col = col, lwd = lwd, panel.last = pL.rel,
                                  panel.first = pF.rel[[i]]),  dotsP[[i1]])

                    assign("plotInfo", plotInfo, envir = trEnv)
                    do.call(plot, args=c(list(resc$X, resc$Y, type = plty,
                                  lty = lty, xlab = .mpresubs(xlab), ylab = .mpresubs(ylab.rel),
                                  col = col, lwd = lwd, panel.last = pL.rel,
                                  panel.first = pF.rel[[i]]),  dotsP[[i1]]))
                    plotInfo <- get("plotInfo", envir = trEnv)
                    plotInfo$relPlotUsr[[i]] <- par("usr")

                    plotInfo$relPlotCArgs[[i]] <- c(list(resc.C$X, resc.C$Y, type = plty,
                            lty = ltyI, col = colI, lwd = lwdI), dotsL)
                    do.call(lines, args = c(list(resc.C$X, resc.C$Y, type = plty,
                            lty = ltyI, col = colI, lwd = lwdI), dotsL))

                    x.ticks0 <- if(xaxt0[i1]!="n") x.ticks[[i1]] else NULL
                    y.ticks0 <- if(yaxt0[i1]!="n") y.ticks[[i1]] else NULL


                    .plotRescaledAxis(scaleX[i1], scaleX.fct[[i1]], scaleX.inv[[i1]],
                              scaleY[i1],scaleY.fct[[i1]],
                              scaleY.inv[[i1]], dots$xlim,
                              dots$ylim, resc$X, ypts = 400, n = scaleN,
                              finiteEndpoints = finiteEndpoints,
                              x.ticks = x.ticks0,
                              y.ticks = y.ticks0, withbox = withbox)
                    plotInfo$relAxis[[i]] <- list(scaleX[i1], scaleX.fct[[i1]], scaleX.inv[[i1]],
                              scaleY[i1],scaleY.fct[[i1]],
                              scaleY.inv[[i1]], dots$xlim,
                              dots$ylim, resc$X, ypts = 400, n = scaleN,
                              finiteEndpoints = finiteEndpoints,
                              x.ticks = x.ticks0,
                              y.ticks = y.ticks0, withbox = withbox)

                    if(with.legend){
                      legend(.legendCoord(legend.location[[i1]],
                                 scaleX[i1], scaleX.fct[[i1]], scaleY[i1], scaleY.fct[[i1]]),
                           bg = legend.bg, legend = legend[[i1]],
                           col = c(colI, col), lwd = c(lwdI, lwd),
                           lty = c(ltyI, lty), cex = legend.cex*fac.leg)
                      plotInfo$relLegend[[i]] <- list(.legendCoord(legend.location[[i1]],
                                 scaleX[i1], scaleX.fct[[i1]], scaleY[i1], scaleY.fct[[i1]]),
                           bg = legend.bg, legend = legend[[i1]],
                           col = c(colI, col), lwd = c(lwdI, lwd),
                           lty = c(ltyI, lty), cex = legend.cex*fac.leg)
                    }
                    if(innerL){
                       do.call(title, args = c(list(main = innerT[[i1]]),
                               dotsT, line = lineT, cex.main = cex.inner, 
                               col.main = col.inner))
                       plotInfo$relTitle[[i]] <- c(list(main = innerT[[i1]]),
                               dotsT, line = lineT, cex.main = cex.inner,
                               col.main = col.inner)
                    }
                }
            }
        cex.main <- if(!hasArg(cex.main)) par("cex.main") else dots$"cex.main"
        col.main <- if(!hasArg(col.main)) par("col.main") else dots$"col.main"
        if (mainL){
            mtext(text = main, side = 3, cex = cex.main, adj = .5,
                  outer = TRUE, padj = 1.4, col = col.main)
            plotInfo$mainL <- list(text = main, side = 3, cex = cex.main, adj = .5,
               outer = TRUE, padj = 1.4, col = col.main)
        }
        cex.sub <- if(!hasArg(cex.sub)) par("cex.sub") else dots$"cex.sub"
        col.sub <- if(!hasArg(col.sub)) par("col.sub") else dots$"col.sub"
        if (subL){
            mtext(text = sub, side = 1, cex = cex.sub, adj = .5,
                  outer = TRUE, line = -1.6, col = col.sub)
            plotInfo$subL <- list(text = sub, side = 1, cex = cex.sub, adj = .5,
               outer = TRUE, line = -1.6, col = col.sub)
        }
        class(plotInfo) <- c("plotInfo","DiagnInfo")
        if(return.Order){ whichRet <- names(plotInfo) %in% c("sel","sel.C")
                            return(plotInfo[whichRet])}
        invisible(plotInfo)
        }
    )
 
