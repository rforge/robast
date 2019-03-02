setMethod("plot", signature(x = "IC", y = "missing"),
    function(x, ...,withSweave = getdistrOption("withSweave"),
             main = FALSE, inner = TRUE, sub = FALSE, 
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
             mfColRow = TRUE, to.draw.arg = NULL, withSubst = TRUE){

        args0 <- list(x = x, withSweave = withSweave,
             main = main, inner = inner, sub = sub,
             col.inner = col.inner, cex.inner = cex.inner,
             bmar = bmar, tmar = tmar, with.automatic.grid = with.automatic.grid,
             with.legend = with.legend, legend = legend, legend.bg = legend.bg,
             legend.location = legend.location, legend.cex = legend.cex,
             withMBR = withMBR, MBRB = MBRB, MBR.fac = MBR.fac, col.MBR = col.MBR,
             lty.MBR = lty.MBR, lwd.MBR = lwd.MBR,
             x.vec = x.vec, scaleX = scaleX,
             scaleX.fct = if(!missing(scaleX.fct)) scaleX.fct else NULL,
             scaleX.inv = if(!missing(scaleX.inv)) scaleX.inv else NULL,
             scaleY = scaleY,
             scaleY.fct = scaleY.fct,
             scaleY.inv = scaleY.inv, scaleN = scaleN, x.ticks = x.ticks,
             y.ticks = y.ticks, mfColRow = mfColRow, to.draw.arg = to.draw.arg,
             withSubst = withSubst)
        mc <- match.call(call = sys.call(sys.parent(1)))
        dots <- match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)$"..."
        plotInfo <- list(call = mc, dots=dots, args=args0)

        xc <- mc$x
        xcc <- as.character(deparse(xc))
        dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
        dotsLeg <- dotsT <- dotsL <- .makedotsLowLevel(dots)
        dotsLeg$lty <- dotsLeg$lwd <- dotsLeg$col <- NULL

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

        L2Fam <- eval(x@CallL2Fam)
        if(missing(scaleX.fct)){
           scaleX.fct <- p(L2Fam)
           scaleX.inv <- q.l(L2Fam)
        }

        trafO <- trafo(L2Fam@param)
        dims  <- nrow(trafO)
        
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

        distr <- L2Fam@distribution
        if(!is(distr, "UnivariateDistribution")) stop("not yet implemented")


        MBRB <- matrix(rep(t(MBRB), length.out=dims0*2),ncol=2, byrow=T)
        MBRB <- MBRB * MBR.fac


        xlim <- eval(dots$xlim)
        ylim <- eval(dots$ylim)
        .xylim <- .getXlimYlim(dots,dotsP, dims0, xlim, ylim)
           dots <- .xylim$dots; dotsP <- .xylim$dotsP
           xlim <- .xylim$xlim; ylim <- .xylim$ylim

        if(missing(x.vec)) x.vec <- NULL
        x.v.ret <- .getX.vec(distr, dims0, eval(dots$lty), x.vec, scaleX, scaleX.fct, scaleX.inv, .xylim$xm, .xylim$xM)
              lty <- x.v.ret$lty; plty <- x.v.ret$plty; x.vec <- x.v.ret$x.vec
        lwd <- if(is.null(dots$lwd))  1 else eval(dots$lwd)

        .pFL <- .preparePanelFirstLast(with.automatic.grid , dims0, pF.0, pL.0,
                             logArg, scaleX, scaleY, x.ticks, y.ticks,
                             scaleX.fct, scaleY.fct)
           pF <- .pFL$pF; pL <- .pFL$pL; gridS <- .pFL$gridS


        plotInfo$to.draw <- to.draw
        plotInfo$panelFirst <- pF
        plotInfo$panelLast <- pL
        plotInfo$gridS <- gridS

        
        if(!is.null(dots[["lty"]]))  dots["lty"] <- NULL
        if(!is.null(dots[["type"]])) dots["type"] <- NULL
        xlab <- dots$xlab; if(is.null(xlab)) xlab <- "x"
        ylab <- dots$ylab; if(is.null(ylab)) ylab <- "(partial) IC"
        dots$xlab <- dots$ylab <- NULL

        IC1 <- as(diag(dims) %*% x@Curve, "EuclRandVariable")

        .pT <- .prepareTitles(withSubst,
                  presubArg2 = c("%C", "%D", "%A"),
                  presubArg3 = c(as.character(class(x)[1]),
                                 as.character(date()),
                                 xcc),
                  dots,
                  mainText =  gettextf("Plot for IC %%A"), ###
                  L2Fam, inner, dims0, dims, to.draw, trafO, x, type = "all", bmar, tmar)

           dots <- .pT$dots; main <- .pT$main; mainL <- .pT$mainL; lineT <- .pT$lineT
           sub <- .pT$sub; subL <- .pT$subL; bmar <- .pT$bmar; tmar <- .pT$tmar;
           innerT <- .pT$innerT; innerL <- .pT$innerL; .mpresubs <- .pT$.mpresubs

        if(with.legend){
          fac.leg <- if(dims0>1) 3/4 else .75/.8
          if(missing(legend.location)){
             legend.location <- .fillList("bottomright", dims0)
          }else{
             legend.location <- as.list(legend.location)
             legend.location <- .fillList(legend.location, dims0)
          }
          if(is.null(legend)){
             legend <- vector("list",dims0)
#             legend <- .fillList(as.list(xc),dims0)
             legend <- .fillList(as.list(xc),dims0)
          }else{
             if(!is.list(legend)) legend <- .fillList(legend,dims0)
          }
        }


        w0 <- getOption("warn")
        options(warn = -1)
        on.exit(options(warn = w0))

        opar <- par(no.readonly = TRUE)
        omar <- par("mar")
        on.exit(par(opar))

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

        dotsT$main <- dotsT$cex.main <- dotsT$col.main <- dotsT$line <- NULL

        dotsT["pch"] <- dotsT["cex"] <- NULL
        dotsT["col"] <- dotsT["lwd"] <- NULL
        dotsL["cex"] <- dotsLeg["bg"] <- dotsLeg["cex"] <- NULL
        dots$ylim <- NULL

        plotInfo$resc.D <- plotInfo$resc <- vector("list", dims0)
        plotInfo$PlotLinesD <- plotInfo$PlotUsr <- vector("list", dims0)
        plotInfo$PlotArgs <- plotInfo$Axis <- vector("list", dims0)
        plotInfo$MBR <- plotInfo$Legend <- plotInfo$innerTitle <- vector("list", dims0)

        IC.f <- function(x,i) .msapply(x, IC1@Map[[i]])

        plotInfo$IC.f <- IC.f

        for(i in 1:dims0){

            indi <- to.draw[i]
            if(!is.null(ylim)) dots$ylim <- ylim[,i]       

            IC.f.i <- function(x) IC.f(x,indi)

            resc <-.rescalefct(x.vec[[i]], IC.f.i, scaleX[i], scaleX.fct[[i]],
                              scaleX.inv[[i]], scaleY[i], scaleY.fct[[i]], xlim[,i],
                              ylim[,i], dots)

            plotInfo$resc[[i]] <- resc
            dots <- resc$dots
            dots$xlim <- xlim[,i]
            dots$ylim <- ylim[,i]
            x.vec1 <- resc$X
            y.vec1 <- resc$Y

            finiteEndpoints <- rep(FALSE,4)
            if(scaleX[i]){
               finiteEndpoints[1] <- is.finite(scaleX.inv[[i]](min(x.vec1, xlim[1,i])))
               finiteEndpoints[2] <- is.finite(scaleX.inv[[i]](max(x.vec1, xlim[2,i])))
            }
            if(scaleY[i]){
               finiteEndpoints[3] <- is.finite(scaleY.inv[[i]](min(y.vec1, ylim[1,i])))
               finiteEndpoints[4] <- is.finite(scaleY.inv[[i]](max(y.vec1, ylim[2,i])))
            }


            if(wmar) do.call(par,args=parArgsL[[i]])

            plotInfo$PlotArgs[[i]] <- c(list(x=x.vec1, y=y.vec1, type = plty,
                                      lty = lty, lwd = lwd,
                                      xlab = .mpresubs(xlab), ylab = .mpresubs(ylab),
                                      panel.first = pF[[i]],
                                      panel.last = pL), dotsP[[i]])

            do.call(plot, args=c(list(x=x.vec1, y=y.vec1, type = plty,
                                      lty = lty, lwd = lwd,
                                      xlab = .mpresubs(xlab), ylab = .mpresubs(ylab),
                                      panel.first = pF[[i]],
                                      panel.last = pL), dotsP[[i]]))

            x.ticks0 <- if(xaxt0[i]!="n") x.ticks[[i]] else NULL
            y.ticks0 <- if(yaxt0[i]!="n") y.ticks[[i]] else NULL


            plotInfo$PlotUsr[[i]] <- par("usr")
            .plotRescaledAxis(scaleX[i], scaleX.fct[[i]], scaleX.inv[[i]],
                              scaleY[i],scaleY.fct[[i]], scaleY.inv[[i]],
                              xlim[,i], ylim[,i], x.vec1, ypts = 400, n = scaleN,
                              finiteEndpoints = finiteEndpoints,
                              x.ticks = x.ticks[[i]], y.ticks = y.ticks[[i]])
            plotInfo$Axis[[i]] <- list(scaleX[i], scaleX.fct[[i]], scaleX.inv[[i]],
                              scaleY[i],scaleY.fct[[i]], scaleY.inv[[i]],
                              xlim[,i], ylim[,i], x.vec1, ypts = 400, n = scaleN,
                              finiteEndpoints = finiteEndpoints,
                              x.ticks = x.ticks[[i]], y.ticks = y.ticks[[i]])
            if(withMBR){
                MBR.i <- MBRB[i,]
                if(scaleY[i]) MBR.i <- scaleY.fct[[i]](MBR.i)
                abline(h=MBR.i, col=col.MBR, lty=lty.MBR, lwd = lwd.MBR)
                plotInfo$MBR[[i]] <- list(h=MBR.i, col=col.MBR, lty=lty.MBR, lwd = lwd.MBR)
            }
            if(is(distr, "DiscreteDistribution")){
                x.vec1D <- seq(from = min(x.vec[[i]]), to = max(x.vec[[i]]), length = 1000)
                rescD <-.rescalefct(x.vec1D, IC.f.i, scaleX[i], scaleX.fct[[i]],
                                scaleX.inv[[i]], scaleY[i], scaleY.fct[[i]], xlim[,i],
                                ylim[,i], dotsP[[i]])
                plotInfo$resc.D[[i]] <- rescD
                x.vecD <- rescD$X
                y.vecD <- rescD$Y

                dotsL$lty <- NULL

                if(is.null(dotsL$lwd)) dotsL$lwd <- lwd
                do.call(lines,args=c(list(x.vecD, y.vecD,
                                          lty = "dotted"), dotsL))
                plotInfo$PlotLinesD[[i]] <- c(list(x.vecD, y.vecD,
                                          lty = "dotted"), dotsL)
            }
            do.call(title,args=c(list(main = innerT[i]), dotsT, line = lineT,
                    cex.main = cex.inner, col.main = col.inner))
            plotInfo$innerTitle[[i]] <- c(list(main = innerT[i]), dotsT, line = lineT,
                    cex.main = cex.inner, col.main = col.inner)

            if(with.legend){
               legend(.legendCoord(legend.location[[i]], scaleX[i], scaleX.fct[[i]],
                        scaleY[i], scaleY.fct[[i]]), bg = legend.bg,
                      legend = legend[[i]], dotsLeg, cex = legend.cex*fac.leg,
                      lwd = lwd, lty = lty)
               plotInfo$Legend[[i]] <- list(.legendCoord(legend.location[[i]],
                      scaleX[i], scaleX.fct[[i]], scaleY[i], scaleY.fct[[i]]),
                      bg = legend.bg, legend = legend[[i]], dotsLeg,
                      cex = legend.cex*fac.leg, lwd = lwd, lty = lty)
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
        return(invisible(plotInfo))
    })


setMethod("plot", signature(x = "IC",y = "numeric"),
          function(x, y, ...,
          cex.pts = 1, cex.pts.fun = NULL, col.pts = par("col"),
          pch.pts = 19,
          cex.npts = 1, cex.npts.fun = NULL, col.npts = par("col"),
          pch.npts = 20,
          jitter.fac = 1, with.lab = FALSE, cex.lbs = 1, adj.lbs = c(0,0),
          col.lbs = col.pts,
          lab.pts = NULL, lab.font = NULL, alpha.trsp = NA,
          which.lbs = NULL, which.Order  = NULL, which.nonlbs = NULL,
          attr.pre = FALSE, return.Order = FALSE){

        args0 <- list(x = x, y = y, cex.pts = cex.pts, cex.pts.fun = cex.pts.fun,
             col.pts = col.pts, pch.pts = pch.pts, cex.npts = cex.npts,
             cex.npts.fun = cex.npts.fun, col.npts = col.npts, pch.npts = pch.npts,
             jitter.fac = jitter.fac, with.lab = with.lab,
             cex.lbs = cex.lbs, adj.lbs = adj.lbs,
             col.lbs = if(!missing(col.lbs)) col.lbs else if(!missing(col.pts)) col.pts else par("col"),
             lab.pts = lab.pts,
             lab.font = lab.font, alpha.trsp = alpha.trsp,
             which.lbs = which.lbs, which.Order  = which.Order,
             which.nonlbs = which.nonlbs, attr.pre = attr.pre,
             return.Order = return.Order)
        mc <- match.call(call = sys.call(sys.parent(1)))
        dots <- match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)$"..."
        plotInfo <- list(call = mc, dots=dots, args=args0)

    n <- if(!is.null(dim(y))) nrow(y) else length(y)

    L2Fam <- eval(x@CallL2Fam)
    trafO <- trafo(L2Fam@param)
    dims0 <- length(.getToDraw(nrow(trafO), trafO, L2Fam, eval(dots$to.draw.arg)))

    if(missing(adj.lbs)) adj.lbs <- c(0,0)
    if(!is.matrix(adj.lbs) ||
          (is.matrix(adj.lbs)&&!all.equal(dim(adj.lbs),c(2,dims0)))){
          adj.lbs <- matrix(rep(adj.lbs, length.out= dims0*2),nrow=2,ncol=dims0)
    }

    lab.pts <- if(is.null(lab.pts)) paste(1:n) else rep(lab.pts,length.out=n)

    if(attr.pre){
       if(missing(pch.pts)) pch.pts <- 1
       if(!length(pch.pts)==n)
          pch.pts <- rep(pch.pts, length.out= n)
       if(missing(col.pts)) col.pts <- par("col")
       if(!length(col.pts)==n)
          col.pts <- rep(col.pts, length.out= n)
       if(missing(cex.pts)) cex.pts <- 1
       if(!length(cex.pts)==n)
          cex.pts <- rep(cex.pts, length.out= n)
       if(missing(cex.lbs)) cex.lbs <- 1
       if(!is.matrix(cex.lbs) ||
          (is.matrix(cex.lbs)&&!all.equal(dim(cex.lbs),c(n,dims0)))){
          cex.lbs <- matrix(rep(cex.lbs, length.out= n*dims0),nrow=n,ncol=dims0)
       }
       if(missing(col.lbs)) col.lbs <- col.pts
       if(!length(col.lbs)==n)
          col.lbs <- rep(col.lbs, length.out= n)
    }


    L2Fam <- eval(x@CallL2Fam)
    trafO <- trafo(L2Fam@param)
    dims <- nrow(trafO)
    dimm <- length(L2Fam@param)
    QF <- diag(dims)

    if(is(x,"ContIC") & dims>1 )
      {if (is(normtype(x),"QFNorm")) QF <- QuadForm(normtype(x))}

    IC1 <- as(diag(dims) %*% x@Curve, "EuclRandVariable")
    absInfo <- t(IC1) %*% QF %*% IC1
    ICMap <- IC1@Map

    ICabs.f <- function(x) .msapply(x, absInfo@Map[[1]])
    plotInfo$ICabs.f <- ICabs.f

    IC.f <- function(x,i) .msapply(x, IC1@Map[[i]])
    plotInfo$IC.f <- IC.f

    sel <- .SelectOrderData(y, ICabs.f, which.lbs, which.Order, which.nonlbs)
    plotInfo$sel <- sel
    plotInfo$obj <- sel$ind

    i.d <- sel$ind
    i0.d <- sel$ind1
    n.s <- length(i.d)

    i.d.ns <- sel$ind.ns
    n.ns <- length(i.d.ns)

    lab.pts <- lab.pts[sel$ind]
    if(attr.pre){
       col.pts <- col.pts[sel$ind]
       col.npts <- col.pts[sel$ind.ns]
       pch.npts <- pch.pts[sel$ind.ns]
       pch.pts <- pch.pts[sel$ind]
       cex.npts <- cex.pts[sel$ind.ns]
       cex.pts <- cex.pts[sel$ind]
       cex.lbs <-  cex.lbs[sel$ind,]
       col.lbs <-  col.lbs[sel$ind]
    }else{
       if(missing(pch.pts)) pch.pts <- 1
       if(!length(pch.pts)==n.s)
          pch.pts <- rep(pch.pts, length.out= n.s)
       if(missing(col.pts)) col.pts <- par("col")
       if(!length(col.pts)==n.s)
          col.pts <- rep(col.pts, length.out= n.s)
       if(missing(cex.pts)) cex.pts <- 1
       if(!length(cex.pts)==n.s)
          cex.pts <- rep(cex.pts, length.out= n.s)

       if(missing(pch.npts)) pch.npts <- 1
       if(!length(pch.npts)==n.ns)
          pch.npts <- rep(pch.npts, length.out= n.ns)
       if(missing(col.npts)) col.npts <- par("col")
       if(!length(col.npts)==n.ns)
          col.npts <- rep(col.npts, length.out= n.ns)
       if(missing(cex.npts)) cex.npts <- 1
       if(!length(cex.npts)==n.ns)
          cex.npts <- rep(cex.npts, length.out= n.ns)
       if(!is.matrix(cex.lbs) ||
          (is.matrix(cex.lbs)&&!all.equal(dim(cex.lbs),c(n.s,dims0)))){
          cex.lbs <- matrix(rep(cex.lbs, length.out= n.s*dims0),nrow=n.s,ncol=dims0)
       }
       if(missing(col.lbs)) col.lbs <- col.pts
       if(!length(col.lbs)==n.s)
          col.lbs <- rep(col.lbs, length.out= n.s)
    }


    dots.without <- dots
    dots.without$col <- dots.without$cex <- dots.without$pch <- NULL

    dims0 <- .getDimsTD(L2Fam,dots[["to.draw.arg"]])

    if(!is.null(cex.pts.fun)){
                  cex.pts.fun <- .fillList(cex.pts.fun)}
    if(!is.null(cex.npts.fun)){
                  cex.npts.fun <- .fillList(cex.npts.fun)}

    pL <- expression({})
    if(!is.null(dots$panel.last))
        pL <- .panel.mingle(dots,"panel.last")
    if(is.list(pL)){
       pL <- .fillList(pL, dims0)

       if(dims0) for(i in 1:dims0){
          if(is.null(pL[[i]])) pL[[i]] <- expression({})
       }
       pL <- substitute({pL1 <- pL0
                         pL1[[i]]},
                         list(pL0=pL))
    }

    dots$panel.last <- NULL

    plotInfo$resc.dat <- plotInfo$resc.dat.ns <- vector("list", dims0)
    plotInfo$doPts <- plotInfo$doPts.ns <- plotInfo$doLabs <- vector("list", dims0)

    trEnv <- new.env()

    pL <- substitute({
        pI <- get("plotInfo", envir = trEnv0)

        IC.f.i <- function(x) IC.f.0(x,indi)

        if(length(y0s)){
            resc.dat <-.rescalefct(y0s, IC.f.i,
                              scaleX[i], scaleX.fct[[i]], scaleX.inv[[i]],
                              scaleY[i], scaleY.fct[[i]], xlim[,i], ylim[,i],
                              dwo0)
            pI$resc.dat[[i]] <- resc.dat
            y1 <- resc.dat$X
            ICy <- resc.dat$Y
            if(is(distr, "DiscreteDistribution")){
               if(length(ICy)) ICy <- jitter(ICy, factor = jitter.fac0) }
            col.pts <- if(!is.na(al0)) .msapply(col0, addAlphTrsp2col,alpha=al0) else col0
            cfun <- if(is.null(cexfun)) NULL else cexfun[[i]]
            cex.l    <- .cexscale(absy0,absy0,cex=cex0, fun = cfun)   ##.cexscale in infoPlot.R

            pI$doPts[[i]] <- c(list(y1, ICy, cex = cex.l,
                        col = col.pts, pch = pch0), dwo0)
            do.call(points, args=c(list(y1, ICy, cex = cex.l,
                        col = col.pts, pch = pch0), dwo0))

            if(with.lab0){
               text(x = y0s, y = ICy, adj=adj.lb0[,i], labels = lab.pts0,
                    cex = cex.lb0[,i], col = col.lb0)
               pI$doLabs[[i]] <- list(x = y0s, y = ICy, adj=adj.lb0[,i],
                    labels = lab.pts0, cex = cex.lb0[,i], col = col.lb0)
            }
        }

        if(length(y0s.ns)){
            resc.dat.ns <-.rescalefct(y0s.ns, IC.f.i,
                              scaleX[i], scaleX.fct[[i]], scaleX.inv[[i]],
                              scaleY[i], scaleY.fct[[i]], xlim[,i], ylim[,i],
                              dwo0)
            pI$resc.dat.ns[[i]] <- resc.dat.ns
            y1.ns <- resc.dat.ns$X
            ICy.ns <- resc.dat.ns$Y
            if(is(distr, "DiscreteDistribution"))
               {if(length(ICy.ns)) ICy.ns <- jitter(ICy.ns, factor = jitter.fac0) }

           col.npts <- if(!is.na(al0)) .msapply(col0.ns, addAlphTrsp2col,alpha=al0) else col0.ns
           cfun.ns <- if(is.null(cexnfun)) NULL else cexnfun[[i]]
           cex.l.ns <- .cexscale(absy0.ns,absy0.ns, cex=cex0.ns, fun = cfun.ns)   ##.cexscale in infoPlot.R

           pI$doPts.ns[[i]] <- c(list(y1.ns, ICy.ns, cex = cex.l.ns,
                        col = col.npts, pch = pch0.ns), dwo0)
           do.call(points, args=c(list(y1.ns, ICy.ns, cex = cex.l.ns,
                        col = col.npts, pch = pch0.ns), dwo0))
        }

        assign("plotInfo", pI, envir = trEnv0)
        pL0
        }, list(pL0 = pL, IC.f.0 = IC.f,
                y0s = sel$data, absy0 = sel$y,
                y0s.ns = sel$data.ns, absy0.ns = sel$y.ns,
                dwo0 = dots.without,
                cex0 = cex.pts, pch0 = pch.pts, col0 = col.pts,
                cex0.ns = cex.npts, pch0.ns = pch.npts, col0.ns = col.npts,
                with.lab0 = with.lab, lab.pts0 = lab.pts,
                al0 = alpha.trsp, jitter.fac0 = jitter.fac,
                cexfun=cex.pts.fun, cexnfun=cex.npts.fun,
                trEnv0 = trEnv, cex.lb0 = cex.lbs, adj.lb0 = adj.lbs,
                col.lb0=col.lbs
                ))

  assign("plotInfo", plotInfo, envir = trEnv)
  ret <- do.call("plot", args = c(list(x = x, panel.last = pL), dots))
  plotInfo <- get("plotInfo", envir = trEnv)
  ret$dots <- ret$args <- ret$call <- NULL
  plotInfo <- c(plotInfo, ret)
  class(plotInfo) <- c("plotInfo","DiagnInfo")
  if(return.Order){ whichRet <- names(plotInfo) == "obj"
                    return(plotInfo[whichRet])}
  return(invisible(plotInfo))
})


