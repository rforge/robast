setMethod("comparePlot", signature("IC","IC"),
    function(obj1,obj2, obj3 = NULL, obj4 = NULL, data = NULL,
             ..., withSweave = getdistrOption("withSweave"),
             forceSameModel = FALSE,
             main = FALSE, inner = TRUE, sub = FALSE,
             col = par("col"), lwd = par("lwd"), lty,
             col.inner = par("col.main"), cex.inner = 0.8,
             bmar = par("mar")[1], tmar = par("mar")[3],
             with.automatic.grid = TRUE,
             with.legend = FALSE, legend = NULL, legend.bg = "white",
             legend.location = "bottomright", legend.cex = 0.8,
             withMBR = FALSE, MBRB = NA, MBR.fac = 2, col.MBR = par("col"),
             lty.MBR = "dashed", lwd.MBR = 0.8,
             x.vec = NULL, scaleX = FALSE, scaleX.fct, scaleX.inv,
             scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
             scaleN = 9, x.ticks = NULL, y.ticks = NULL,
             mfColRow = TRUE, to.draw.arg = NULL,
             cex.pts = 1, cex.pts.fun = NULL, col.pts = par("col"),
             pch.pts = 19,
             cex.npts = 1, cex.npts.fun = NULL, col.npts = par("col"),
             pch.npts = 20,
             jitter.fac = 1, with.lab = FALSE, cex.lbs = 1, adj.lbs = c(0,0),
             col.lbs = col.pts, lab.pts = NULL,
             lab.font = NULL, alpha.trsp = NA,
             which.lbs = NULL, which.Order  = NULL, which.nonlbs = NULL,
             attr.pre = FALSE, return.Order = FALSE,
             withSubst = TRUE){

        args0 <- list(obj1 = obj1, obj2 = obj2, obj3 = obj3, obj4 = obj4,
             data = data, withSweave = withSweave, forceSameModel = forceSameModel,
             main = main, inner = inner, sub = sub, col = col, lwd = lwd,
             lty = if(!missing(lty)) lty else NULL,
             col.inner = col.inner, cex.inner = cex.inner,
             bmar = bmar, tmar = tmar, with.automatic.grid = with.automatic.grid,
             with.legend = with.legend, legend = legend, legend.bg = legend.bg,
             legend.location = legend.location, legend.cex = legend.cex,
             withMBR = withMBR, MBRB = MBRB, MBR.fac = MBR.fac, col.MBR = col.MBR,
             lty.MBR = lty.MBR, lwd.MBR = lwd.MBR,
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
             return.Order = return.Order, withSubst = withSubst)

        .mc <- match.call(call = sys.call(sys.parent(1)))
        dots <- match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)$"..."
        plotInfo <- list(call = .mc, dots=dots, args=args0)
        .xc<- function(obj) as.character(deparse(.mc[[obj]]))
        xc <- c(.xc("obj1"), .xc("obj2"))
        if(!is.null(obj3)) xc <- c(xc, .xc("obj3"))
        if(!is.null(obj4)) xc <- c(xc, .xc("obj4"))

        dotsP <- dots
        dotsP$type <- dotsP$lty <- dotsP$col <- dotsP$lwd <- NULL
        dotsP$xlab <- dotsP$ylab <- NULL

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

        dotsLeg <- dotsT <- dotsL <- .makedotsLowLevel(dots)
        dotsLeg$lty <- dotsLeg$lwd <- dotsLeg$col <- NULL
        dots.points <-   .makedotsPt(dots)
        
        ncomp <- 2+ (!missing(obj3)|!is.null(obj3)) +
                    (!missing(obj4)|!is.null(obj4))

        if(missing(col)) col <- 1:ncomp
           else col <- rep(col, length.out = ncomp)
        if(missing(lwd))  lwd <- rep(1,ncomp)
           else lwd <- rep(lwd, length.out = ncomp)
        if(!missing(lty)) rep(lty, length.out = ncomp)

        dots["type"] <- NULL
        xlab <- dots$xlab; if(is.null(xlab)) xlab <- "x"
        ylab <- dots$ylab; if(is.null(ylab)) ylab <- "(partial) IC"
        dots$xlab <- dots$ylab <- NULL

        L2Fam <- eval(obj1@CallL2Fam)
        if(forceSameModel)
        if(!identical(CallL2Fam(obj1),CallL2Fam(obj2)))
            stop("ICs need to be defined for the same model")

        if(missing(scaleX.fct)){
           scaleX.fct <- p(L2Fam)
           scaleX.inv <- q.l(L2Fam)
        }

        trafO <- trafo(L2Fam@param)
        dims  <- nrow(trafO)
        dimm <- ncol(trafO)

        to.draw <- .getToDraw(dims, trafO, L2Fam, to.draw.arg)
        dims0 <- length(to.draw)
        nrows <- trunc(sqrt(dims0))
        ncols <- ceiling(dims0/nrows)

        yaxt0 <- xaxt0 <- rep("s",dims0)
        if(!is.null(dots$xaxt)){ xaxt1 <- eval(dots$xaxt); xaxt0 <- rep(xaxt1, length.out=dims0)}
        if(!is.null(dots$yaxt)){ yaxt1 <- eval(dots$yaxt); yaxt0 <- rep(yaxt1, length.out=dims0)}

        logArg <- NULL
        if(!is.null(dots[["log"]]))
            logArg <- rep(dots[["log"]], length.out=dims0)
        dotsP$log <- dots$log <- NULL

        dotsP0 <- vector("list",dims0)
        if(!is.null(dotsP)) for(i in 1:dims0) dotsP0[[i]] <- dotsP
        dotsP <- dotsP0
        for(i in 1:dims0){dotsP[[i]]$xaxt <- xaxt0[i];dotsP[[i]]$yaxt <- yaxt0[i]}

        if(!is.null(logArg))
            for(i in 1:dims0) dotsP[[i]]$log <- logArg[i]

        if(!is.null(x.ticks)){
           x.ticks <- .fillList(x.ticks, dims0)
           for(i in 1:dims0){
               if(!is.null(x.ticks[[i]]))
                   if(!is.null(logArg)) if(!grepl("x",logArg[i])) dotsP[[i]]$xaxt <- "n"
           }
        }
        if(!is.null(y.ticks)){
           y.ticks <- .fillList(y.ticks, dims0)
           for(i in 1:dims0){
               if(!is.null(y.ticks[[i]]))
                   if(!is.null(logArg)) if(!grepl("y",logArg[i])) dotsP[[i]]$yaxt <- "n"
           }
        }

        scaleX <- rep(scaleX, length.out=dims0)
        scaleY <- rep(scaleY, length.out=dims0)
        scaleX <- scaleX & (xaxt0!="n")
        scaleY <- scaleY & (yaxt0!="n")

        scaleX.fct <- .fillList(scaleX.fct, dims0)
        scaleX.inv <- .fillList(scaleX.inv, dims0)

        scaleY.fct <- .fillList(scaleY.fct, dims0)
        scaleY.inv <- .fillList(scaleY.inv, dims0)

        MBRB <- matrix(rep(t(MBRB), length.out=dims0*2),ncol=2, byrow=T)
        MBRB <- MBRB * MBR.fac

        distr <- L2Fam@distribution
        if(!is(distr, "UnivariateDistribution")) stop("not yet implemented")


        xlim <- eval(dots$xlim)
        ylim <- eval(dots$ylim)
        .xylim <- .getXlimYlim(dots,dotsP, dims0, xlim, ylim)
          dots <- .xylim$dots; dotsP <- .xylim$dotsP
          xlim <- .xylim$xlim; ylim <- .xylim$ylim;

        if(missing(x.vec)) x.vec <- NULL
        x.v.ret <- .getX.vec(distr, dims0, eval(dots$lty), x.vec, scaleX, scaleX.fct, scaleX.inv, .xylim$xm, .xylim$xM)
              lty <- x.v.ret$lty; plty <- x.v.ret$plty; x.vec <- x.v.ret$x.vec

        dims <- nrow(trafo(L2Fam@param)); ID <- diag(dims)
        IC1 <- as(ID %*% obj1@Curve, "EuclRandVariable")
        IC2 <- as(ID %*% obj2@Curve, "EuclRandVariable")

        if(is(obj3, "IC")){
          if(forceSameModel)
           if(!identical(CallL2Fam(obj1),CallL2Fam(obj3)))
               stop("ICs need to be defined for the same model")
           IC3 <- as(ID %*% obj3@Curve, "EuclRandVariable")
        }

        if(is(obj4, "IC")){
          if(forceSameModel)
           if(!identical(CallL2Fam(obj1),CallL2Fam(obj4)))
               stop("ICs need to be defined for the same model")
           IC4 <- as(ID %*% obj4@Curve, "EuclRandVariable")
        }


        .pT <- .prepareTitles(withSubst,
              presubArg2 = c(paste("%C",1:ncomp,sep=""),"%D",
                             paste("%A",1:ncomp,sep="")),
              presubArg3 = c(as.character(class(obj1)[1]),
                    as.character(class(obj2)[1]),
                    if(is.null(obj3))NULL else as.character(class(obj3)[1]),
                    if(is.null(obj4))NULL else as.character(class(obj4)[1]),
                    as.character(date()),
                    xc), dots,
              mainText = paste(gettextf("Plot for ICs"),
                                paste("%A", 1:ncomp, sep="", collapse=", "),
                                sep=" "),
              L2Fam, inner, dims0, dims, to.draw, trafO, L2Fam, type = "compare", bmar, tmar)
        dots <- .pT$dots; main <- .pT$main; mainL <- .pT$mainL; lineT <- .pT$lineT
        sub <- .pT$sub; subL <- .pT$subL; bmar <- .pT$bmar; tmar <- .pT$tmar;
        innerT <- .pT$innerT; innerL <- .pT$innerL; .mpresubs <- .pT$.mpresubs

        w0 <- getOption("warn"); options(warn = -1); on.exit(options(warn = w0))

        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        omar <- par("mar")

        if(mfColRow) par(mfrow = c(nrows, ncols)) else{
          if(!withSweave && length(dev.list())>0) devNew()
        }

        wmar <- FALSE
        if(!missing(bmar)||!missing(tmar)){
             lpA <- max(dims0,1)
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

        if(is(distr, "DiscreteDistribution")){
            x.vecD <- vector("list", dims0)
            for(i in 1:dims0)
                x.vecD[[i]] <- seq(from = min(x.vec[[i]]), to = max(x.vec[[i]]), length = 1000)
        }
        dotsT$main <- dotsT$cex.main <- dotsT$col.main <- dotsT$line <- NULL

        .pFL <- .preparePanelFirstLast(with.automatic.grid , dims0, pF.0, pL.0,
                             logArg, scaleX, scaleY, x.ticks, y.ticks,
                             scaleX.fct, scaleY.fct)
        pF <- .pFL$pF; pL <- .pFL$pL; gridS <- .pFL$gridS

        plotInfo$to.draw <- to.draw
        plotInfo$panelFirst <- pF
        plotInfo$panelLast <- pL
        plotInfo$gridS <- gridS


        plotInfo$IC1.f <- IC1.f <- function(x,i) .msapply(x, IC1@Map[[i]])
        plotInfo$IC2.f <- IC2.f <- function(x,i) .msapply(x, IC2@Map[[i]])

        if(is(obj3, "IC")){
           plotInfo$IC3.f <- IC3.f <- function(x,i) .msapply(x, IC3@Map[[i]])
        }
        if(is(obj4, "IC")){
           plotInfo$IC4.f <- IC4.f <- function(x,i) .msapply(x, IC4@Map[[i]])
         }

        trEnv <- new.env()

        if(!is.null(data)){

            sel1 <- sel2 <- sel3 <- sel4 <- NULL

            plotInfo$resc.dat <- vector("list", ncomp*dims0)
            plotInfo$resc.dat.ns <- vector("list", ncomp*dims0)
            plotInfo$doPts <- plotInfo$doPts.ns <- vector("list", ncomp*dims0)
            plotInfo$doLabs <- vector("list", ncomp*dims0)

            n <- if(!is.null(dim(data))) nrow(data) else length(data)

            lab.pts <- if(is.null(lab.pts)) paste(1:n) else rep(lab.pts,length.out=n)

            if(!is.null(cex.pts.fun)){
                  cex.pts.fun <- .fillList(cex.pts.fun, dims0*ncomp)}
            if(!is.null(cex.npts.fun)){
                  cex.npts.fun <- .fillList(cex.npts.fun, dims0*ncomp)}

            if(missing(adj.lbs)) cex.lbs <- c(0,0)
            if(!is.array(adj.lbs) ||
                (is.array(adj.lbs)&&!all.equal(dim(adj.lbs),c(2,ncomp,dims0)))){
                 adj.lbs <- array(rep(adj.lbs, length.out= 2*dims0*ncomp),
                                      dim=c(2,ncomp,dims0))
            }

            if(attr.pre){
               if(missing(pch.pts)) pch.pts <- 1
               if(!is.matrix(pch.pts))
                   pch.pts <- t(matrix(rep(pch.pts, length.out= ncomp*n),ncomp,n))

               if(missing(col.pts)) col.pts <- 1:ncomp
               if(!is.matrix(col.pts))
                  col.pts <- t(matrix(rep(col.pts, length.out= ncomp*n),ncomp,n))

               if(missing(cex.pts)) cex.pts <- 1
               if(!is.matrix(cex.pts))
                  cex.pts <- matrix(rep(cex.pts, length.out= ncomp*n),n,ncomp)

               if(missing(cex.lbs)) cex.lbs <- 1
               if(!is.array(cex.lbs) ||
                   (is.array(cex.lbs)&&!all.equal(dim(cex.lbs),c(n,ncomp,dims0)))){
                    cex.lbs <- array(rep(cex.lbs, length.out= n*dims0*ncomp),
                                     dim=c(n,ncomp,dims0))
                   }

               if(missing(col.lbs)) col.lbs <- col.pts
               if(!is.matrix(col.lbs))
                  col.lbs <- t(matrix(rep(col.lbs, length.out= ncomp*n),ncomp,n))
               }


            absInfoEval <- function(x,IC){
                  QF <- ID
                  if(is(IC,"ContIC") & dims>1 ){
                     if (is(normtype(IC),"QFNorm"))
                          QF <- QuadForm(normtype(IC))
                  }
                  absInfo.f <- t(IC) %*% QF %*% IC
                  return(.msapply(x, absInfo.f@Map[[1]]))
            }

            plotInfo$ICabs.f <- absInfoEval
            plotInfo$IC1abs.f <- function(x) absInfoEval(x,IC1)
            plotInfo$IC2abs.f <- function(x) absInfoEval(x,IC2)


            def.sel <- function(IC){
                 fct.aI <- function(x) absInfoEval(x,IC)
                 return(.SelectOrderData(data, fct.aI, which.lbs, which.Order,
                                         which.nonlbs))}
                 
            sel1 <- def.sel(IC1); sel2 <- def.sel(IC2)
            plotInfo$sel1 <- sel1
            plotInfo$sel2 <- sel2
            plotInfo$obj1 <- sel1$ind
            plotInfo$obj2 <- sel2$ind
            selAlly.s <- c(sel1$y,sel2$y)
            selAlly.ns <- c(sel1$y.ns,sel2$y.ns)

            n.s <- length(sel1$ind)
            n.ns <- length(sel1$ind.ns)

            lab0.pts <- matrix(NA, n.s, ncomp)
            lab0.pts[,1] <- lab.pts[sel1$ind]
            lab0.pts[,2] <- lab.pts[sel2$ind]

            if(attr.pre){
               col0.pts     <- col.pts[sel1$ind,]
               col0.pts[,2] <- col.pts[sel2$ind,2]

               pch0.pts     <- pch.pts[sel1$ind,]
               pch0.pts[,2] <- pch.pts[sel2$ind,2]

               cex0.pts     <- cex.pts[sel1$ind,]
               cex0.pts[,2] <- cex.pts[sel2$ind,2]


               cex0.lbs      <- cex.lbs[sel1$ind,,,drop=FALSE]
               cex0.lbs[,2,] <- cex.lbs[sel2$ind,2,]

               col0.lbs     <- col.lbs[sel1$ind,]
               col0.lbs[,2] <- col.lbs[sel2$ind,2]

               col.npts     <- col.pts[sel1$ind.ns,]
               col.npts[,2] <- col.pts[sel2$ind.ns,2]

               pch.npts     <- pch.pts[sel1$ind.ns,]
               pch.npts[,2] <- pch.pts[sel2$ind.ns,2]

               cex.npts     <- cex.pts[sel1$ind.ns,]
               cex.npts[,2] <- cex.pts[sel2$ind.ns,2]
            }


            if(is(obj3, "IC")){ sel3 <- def.sel(IC3)
                                plotInfo$sel3 <- sel3
                                plotInfo$obj3 <- sel3$ind
                                selAlly.s <- c(selAlly.s,sel3$y)
                                selAlly.ns <- c(selAlly.ns,sel3$y.ns)
                                plotInfo$IC3abs.f <- function(x) absInfoEval(x,IC3)
                                lab0.pts[,3] <- lab.pts[sel3$ind]
                                if(attr.pre){
                                   col0.pts[,3] <- col.pts[sel3$ind,3]
                                   pch0.pts[,3] <- pch.pts[sel3$ind,3]
                                   cex0.pts[,3] <- cex.pts[sel3$ind,3]
                                   cex0.lbs[,3,] <- cex.lbs[sel3$ind,3,]
                                   col0.lbs[,3] <- col.lbs[sel3$ind,3]
                                   col.npts[,3] <- col.pts[sel3$ind.ns,3]
                                   pch.npts[,3] <- pch.pts[sel3$ind.ns,3]
                                   cex.npts[,3] <- cex.pts[sel3$ind.ns,3]
                                }
                              }
            if(is(obj4, "IC")){ sel4 <- def.sel(IC4)
                                plotInfo$sel4 <- sel4
                                plotInfo$obj4 <- sel4$ind
                                selAlly.s <- c(selAlly.s,sel4$y)
                                selAlly.ns <- c(selAlly.ns,sel4$y.ns)
                                plotInfo$IC4abs.f <- function(x) absInfoEval(x,IC4)
                                lab0.pts[,4] <- lab.pts[sel4$ind]
                                if(attr.pre){
                                   col0.pts[,4] <- col.pts[sel4$ind,4]
                                   pch0.pts[,4] <- pch.pts[sel4$ind,4]
                                   cex0.pts[,4] <- cex.pts[sel4$ind,4]
                                   cex0.lbs[,4,] <- cex.lbs[sel4$ind,4,]
                                   col0.lbs[,4] <- col.lbs[sel4$ind,4]
                                   col.npts[,4] <- col.pts[sel4$ind.ns,4]
                                   pch.npts[,4] <- pch.pts[sel4$ind.ns,4]
                                   cex.npts[,4] <- cex.pts[sel4$ind.ns,4]
                                }
                              }

            lab.pts <- lab0.pts

            if(attr.pre){
               col.pts <- col0.pts
               pch.pts <- pch0.pts
               cex.pts <- cex0.pts
               cex.lbs <- cex0.lbs
               col.lbs <- col0.lbs
            }else{
               if(missing(pch.pts)) pch.pts <- 1
               if(!is.matrix(pch.pts))
                   pch.pts <- t(matrix(rep(pch.pts, length.out= ncomp*n.s),ncomp,n.s))
               if(missing(pch.npts)) pch.npts <- 2
               if(!is.matrix(pch.npts))
                   pch.npts <- t(matrix(rep(pch.npts, length.out= ncomp*n.ns),ncomp,n.ns))

               if(missing(col.pts)) col.pts <- 1:ncomp
               if(!is.matrix(col.pts))
                  col.pts <- t(matrix(rep(col.pts, length.out= ncomp*n.s),ncomp,n.s))
               if(missing(col.npts)) col.npts <- 1:ncomp
               if(!is.matrix(col.npts))
                  col.npts <- t(matrix(rep(col.npts, length.out= ncomp*n.ns),ncomp,n.ns))

               if(missing(cex.pts)) cex.pts <- 1
               if(!is.matrix(cex.pts))
                  cex.pts <- matrix(rep(cex.pts, length.out= ncomp*n.s),n.s,ncomp)
               if(missing(cex.npts)) cex.npts <- 1
               if(!is.matrix(cex.npts))
                  cex.npts <- matrix(rep(cex.npts, length.out= ncomp*n.ns),n.ns,ncomp)

               if(missing(cex.lbs)) cex.lbs <- 1
               if(!is.array(cex.lbs) ||
                   (is.array(cex.lbs)&&all.equal(dim(cex.lbs),c(n.s,ncomp,dims0)))){
                    cex.lbs <- array(rep(cex.lbs, length.out= n.s*dims0*dims0),
                                     dim=c(n.s,ncomp,dims0))
                   }

               if(missing(col.lbs)) col.lbs <- col.pts
               if(!is.matrix(col.lbs))
                  col.lbs <- t(matrix(rep(col.lbs, length.out= ncomp*n.s),ncomp,n.s))
            }


            dots.points <- .makedotsLowLevel(dots)
            dots.points$col <- dots.points$cex <- dots.points$pch <- NULL
            alp.v <- rep(alpha.trsp,length.out = ncomp)

            plotInfo$resc.D <- plotInfo$resc <- vector("list", ncomp*dims0)
            plotInfo$resc.dat <- plotInfo$resc.dat.ns <- vector("list", ncomp*dims0)
            plotInfo$doPts <- plotInfo$doPts.ns <- plotInfo$doLabs <- vector("list", ncomp*dims0)
            plotInfo$PlotLines <- plotInfo$PlotPoints <- plotInfo$PlotUsr <- vector("list", dims0)
            plotInfo$PlotLinesD <- plotInfo$PlotArgs <- plotInfo$Axis <- vector("list", dims0)
            plotInfo$MBR <- plotInfo$Legend <- plotInfo$innerTitle <- vector("list", dims0)


            pL.o <- pL
            pL <- substitute({

                 doIt <- function(sel.l,fct.l,j.l, trEnv1){

                     pI <- get("plotInfo", envir = trEnv1)

                     rescd <- .rescalefct(sel.l$data, fct.l, scaleX[i], scaleX.fct[[i]],
                                   scaleX.inv[[i]], scaleY[i], scaleY.fct[[i]], xlim[,i],
                                   ylim[,i], dotsP[[i]])
                     if(is(distr, "DiscreteDistribution")){
                        if(length(rescd$Y))
                           rescd$Y <- jitter(rescd$Y, factor = jitter.fac0[j.l])
                     }
                     i.l <- sel.l$ind
                     n.l <- length(i.l)

                     i.l.ns <- sel.l$ind.ns
                     n.l.ns <- length(i.l.ns)



                     if(n.l){
                        pI$resc.dat[[(i-1)*ncomp+j.l]] <- rescd[i.l]
                        lab.pts.l <- if(is.null(lab0)) paste(i.l) else lab0[,j.l]
                        col.l <- if(is.na(al0[j.l])) col0[,j.l] else
                                 addAlphTrsp2col(col0[,j.l], al0[j.l])
                        pch.l <- pch0[,j.l]
                        cfun <- if(is.null(cexfun)) NULL else cexfun[[(i-1)*ncomp+j.l]]
                        cex.l    <- .cexscale(sel.l$y,selAlly.s,cex=cex0[,j.l], fun = cfun)   ##.cexscale in infoPlot.R
                        if(length(rescd$X[i.l])){
                           pI$doPts[[(i-1)*ncomp+j.l]] <- c(list(rescd$X[i.l], rescd$Y[i.l], cex = cex.l,
                              col = col.l, pch = pch.l), dwo0)
                           do.call(points, args=c(list(rescd$X[i.l], rescd$Y[i.l], cex = cex.l,
                              col = col.l, pch = pch.l), dwo0))
                           if(with.lab0){
                              text(rescd$X[i.l], rescd$Y[i.l], labels = lab.pts.l,
                                   cex = cexl0[,j.l,i], col = coll0[,j.l], adj=adjl0[,j.l,i])
                              pI$doLabs[[(i-1)*ncomp+j.l]] <- list(rescd$X[i.l], rescd$Y[i.l], labels = lab.pts.l,
                                   cex = cexl0[,j.l,i], col = coll0[,j.l],adj=adjl0[,j.l,i])
                           }
                        }
                     }
                     if(n.l.ns){
                        pI$resc.dat.ns[[(i-1)*ncomp+j.l]] <- rescd[i.l.ns]
                        col.l.ns <- if(is.na(al0[j.l])) coln0[,j.l] else
                                 addAlphTrsp2col(coln0[,j.l], al0[j.l])
                        pch.l.ns <- pchn0[,j.l]
                        cfun.ns <- if(is.null(cexnfun)) NULL else cexnfun[[(i-1)*ncomp+j.l]]
                        cex.l.ns <- .cexscale(sel.l$y.ns,selAlly.ns, cex=cexn0[,j.l], fun = cfun.ns)   ##.cexscale in infoPlot.R
                        if(length(rescd$X[i.l.ns])){
                           pI$doPts.ns[[(i-1)*ncomp+j.l]] <- c(list(rescd$X[i.l.ns], rescd$Y[i.l.ns], cex = cex.l.ns,
                              col = col.l.ns, pch = pch.l.ns), dwo0)
                           do.call(points, args=c(list(rescd$X[i.l.ns], rescd$Y[i.l.ns], cex = cex.l.ns,
                              col = col.l.ns, pch = pch.l.ns), dwo0))
                        }
                     }
                     assign("plotInfo", pI, envir = trEnv1)
                 }
                 doIt(sel1,fct1,1, trEnv0);  doIt(sel2,fct2,2, trEnv0)
                 if(is(obj3, "IC")) doIt(sel3,fct3,3, trEnv0)
                 if(is(obj4, "IC")) doIt(sel4,fct4,4, trEnv0)
                 pL0
              }, list(pL0 = pL.o, cex0 = cex.pts, pch0 = pch.pts, col0 = col.pts,
                      jitter.fac0 = jitter.fac, dwo0 = dots.points, al0 = alp.v,
                      with.lab0 = with.lab, lab0 = lab.pts, cexfun=cex.pts.fun,
                      cexn0 = cex.npts, pchn0 = pch.npts, coln0 = col.npts,
                      cexnfun=cex.npts.fun, trEnv0 = trEnv, cexl0 = cex.lbs,
                      adjl0 = adj.lbs, coll0 = col.lbs)
                      #,scaleX = scaleX, scaleX.fct = scaleX.fct,
                      #scaleX.inv = scaleX.inv, scaleY = scaleY,
                      #scaleY.fct = scaleY.fct, scaleY.inv = scaleY.inv)
            )
        }


        plotInfo$resc.D <- plotInfo$resc <- vector("list", ncomp*dims0)
        plotInfo$PlotArgs <- plotInfo$PlotPoints <- vector("list", dims0)
        plotInfo$PlotLinesD <- plotInfo$MBR <- plotInfo$PlotUsr <- vector("list", dims0)
        plotInfo$Axis <- plotInfo$PlotLines <- vector("list", dims0)
        plotInfo$Legend <- plotInfo$innerTitle <- vector("list", dims0)


        for(i in 1:dims0){
            indi <- to.draw[i]
            if(!is.null(ylim)) dotsP[[i]]$ylim <- ylim[,i]

            fct1 <- function(x) IC1.f(x,indi)

            resc.args <- c(list(x.vec[[i]], "fc"=fct1, scaleX[i], scaleX.fct[[i]],
                                scaleX.inv[[i]], scaleY[i], scaleY.fct[[i]], xlim[,i],
                                ylim[,i], dotsP[[i]]))

            resc1 <- do.call(.rescalefct, resc.args)
            resc.args$fc <- fct2 <- function(x) IC2.f(x,indi)
            resc2 <- do.call(.rescalefct, resc.args)

            plotInfo$resc[[(i-1)*ncomp+1]] <- resc1
            plotInfo$resc[[(i-1)*ncomp+2]] <- resc2

            dotsP[[i]] <- resc1$dots
            matp  <- cbind(resc1$Y, resc2$Y)

            if(is(obj3, "IC")){
                resc.args$fc <- fct3 <- function(x) IC3.f(x,indi)
                resc3 <- do.call(.rescalefct, resc.args)
                plotInfo$resc[[(i-1)*ncomp+3]] <- resc3
                matp  <- cbind(matp,resc3$Y)
            }
            if(is(obj4, "IC")){
                resc.args$fc <- fct4 <- function(x) IC4.f(x,indi)
                resc4 <- do.call(.rescalefct, resc.args)
                plotInfo$resc[[(i-1)*ncomp+4]] <- resc4
                matp  <- cbind(matp,resc4$Y)
            }

            ym <- min(matp,na.rm=T)
            yM <- max(matp,na.rm=T)
            y0 <- matp[,1]
            y0[1:2] <- c(ym,yM)

            finiteEndpoints <- rep(FALSE,4)
            if(scaleX[i]){
               finiteEndpoints[1] <- is.finite(scaleX.inv[[i]](min(x.vec[[i]], xlim[1,i],na.rm=TRUE)))
               finiteEndpoints[2] <- is.finite(scaleX.inv[[i]](max(x.vec[[i]], xlim[2,i],na.rm=TRUE)))
            }
            if(scaleY[i]){
               finiteEndpoints[3] <- is.finite(scaleY.inv[[i]](min(ym, ylim[1,i],na.rm=TRUE)))
               finiteEndpoints[4] <- is.finite(scaleY.inv[[i]](max(yM, ylim[2,i],na.rm=TRUE)))
            }

            if(wmar) do.call(par,args=parArgsL[[i]])

            assign("plotInfo", plotInfo, envir = trEnv)
            do.call(plot, args=c(list(x = resc1$X, y = y0,
                 type = "n", xlab = .mpresubs(xlab), ylab = .mpresubs(ylab),
                 lty = lty[1], col = addAlphTrsp2col(col[1],0),
                 lwd = lwd[1]), dotsP[[i]], list(panel.last = pL, panel.first=pF[[i]])))
            plotInfo <- get("plotInfo", envir = trEnv)

            plotInfo$PlotUsr[[i]] <- par("usr")

            if(plty=="p"){
               do.call(matpoints, args = c(list( x = resc1$X, y = matp,
                    col = col), dots.points))
               plotInfo$PlotPoints[[i]] <- c(list( x = resc1$X, y = matp,
                    col = col), dots.points)
            }

            do.call(matlines, args = c(list( x = resc1$X, y = matp,
                    lty = lty, col = col, lwd = lwd), dotsL))
            plotInfo$PlotLines[[i]] <- c(list( x = resc1$X, y = matp,
                    lty = lty, col = col, lwd = lwd), dotsL)

            x.ticks0 <- if(xaxt0[i]!="n") x.ticks[[i]] else NULL
            y.ticks0 <- if(yaxt0[i]!="n") y.ticks[[i]] else NULL

            .plotRescaledAxis(scaleX[i], scaleX.fct[[i]], scaleX.inv[[i]],
                              scaleY[i],scaleY.fct[[i]], scaleY.inv[[i]], xlim[,i],
                              ylim[,i], resc1$X, ypts = 400, n = scaleN,
                              finiteEndpoints = finiteEndpoints,
                              x.ticks = x.ticks0, y.ticks = y.ticks0)
            plotInfo$Axis[[i]] <- list( scaleX[i], scaleX.fct[[i]], scaleX.inv[[i]],
                              scaleY[i],scaleY.fct[[i]], scaleY.inv[[i]], xlim[,i],
                              ylim[,i], resc1$X, ypts = 400, n = scaleN,
                              finiteEndpoints = finiteEndpoints,
                              x.ticks = x.ticks0, y.ticks = y.ticks0)
            if(withMBR){
                MBR.i <- MBRB[i,]
                if(scaleY[i]) MBR.i <- scaleY.fct[[i]](MBR.i)
                abline(h=MBR.i, col=col.MBR, lty=lty.MBR, lwd = lwd.MBR)
                plotInfo$MBR[[i]] <- list(h=MBR.i, col=col.MBR, lty=lty.MBR, lwd = lwd.MBR)
            }

            if(is(distr, "DiscreteDistribution")){
                 rescD.args <- c(list(x.vecD, "fc"=fct1, scaleX[i], scaleX.fct[[i]],
                                scaleX.inv[[i]], scaleY[i], scaleY.fct[[i]], xlim[,i],
                                ylim[,i], dotsP[[i]]))
                 resc1D <- do.call(.rescalefct, rescD.args)
                 rescD.args$fc <- fct2
                 resc2D <- do.call(.rescalefct, rescD.args)
                 matpD  <- cbind(resc1D$Y, resc2D$Y)
                 plotInfo$resc.D[[(i-1)*ncomp+1]] <- resc1D
                 plotInfo$resc.D[[(i-1)*ncomp+2]] <- resc2D

                 if(is(obj3, "IC")){
                    rescD.args$fc <- fct3
                    resc3D <- do.call(.rescalefct, rescD.args)
                    plotInfo$resc.D[[(i-1)*ncomp+3]] <- resc3D
                    matpD  <- cbind(matpD, resc3D$Y)
                 }
                 if(is(obj4, "IC")){
                    rescD.args$fc <- fct4
                    resc4D <- do.call(.rescalefct, rescD.args)
                    plotInfo$resc.D[[(i-1)*ncomp+4]] <- resc4D
                    matpD  <- cbind(matpD, resc4D$Y)
                 }
                 do.call(matlines, c(list(resc1D$X, matpD, lty = lty,
                         col = col, lwd = lwd), dotsL))
                 plotInfo$PlotLinesD[[i]] <- c(list(resc1D$X, matpD, lty = lty,
                         col = col, lwd = lwd), dotsL)
            }

           if(innerL){
              do.call(title, args=c(list(main = innerT[[indi]]),  dotsT,
                      line = lineT, cex.main = cex.inner, col.main = col.inner))
              plotInfo$innerTitle[[i]] <- c(list(main = innerT[[indi]]),  dotsT,
                      line = lineT, cex.main = cex.inner, col.main = col.inner)
           }
        }

        if(with.legend){
           if(is.null(legend)) legend <- xc
           legend(.legendCoord(legend.location, scaleX, scaleX.fct,
                        scaleY, scaleY.fct[[i]]), col = col, bg = legend.bg,
                      legend = legend, dotsLeg, cex = legend.cex, lwd = lwd, lty = lty)
           plotInfo$Legend[[i]] <- list(.legendCoord(legend.location, scaleX, scaleX.fct,
                        scaleY, scaleY.fct[[i]]), col = col, bg = legend.bg,
                      legend = legend, dotsLeg, cex = legend.cex, lwd = lwd, lty = lty)
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
        if(return.Order){ whichRet <- names(plotInfo) %in% c("obj1","obj2","obj3","obj4")
                            return(plotInfo[whichRet])}
        return(invisible(plotInfo))
    })

.comparePlot.bkp <- getMethod("comparePlot", signature("IC","IC"))
