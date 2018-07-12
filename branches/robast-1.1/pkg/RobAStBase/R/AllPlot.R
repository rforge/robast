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

       .mpresubs <- if(withSubst){
                     function(inx) 
                      .presubs(inx, c("%C", "%A", "%D" ),
                          c(as.character(class(x)[1]), 
                            as.character(date()), 
                            xcc))
                     }else function(inx)inx

        if(!is.logical(inner)){
          if(!is.list(inner))
              inner <- as.list(inner)
            #stop("Argument 'inner' must either be 'logical' or a 'list'")
           inner <- .fillList(inner,4)
           innerD <- inner[1:3]
           innerL <- inner[4] 
        }else{innerD <- innerL <- inner}


        L2Fam <- eval(x@CallL2Fam)
        if(missing(scaleX.fct)){
           scaleX.fct <- p(L2Fam)
           scaleX.inv <- q.l(L2Fam)
        }

        trafO <- trafo(L2Fam@param)
        dims  <- nrow(trafO)
        
        to.draw <- 1:dims
        dimnms  <- c(rownames(trafO))
        if(is.null(dimnms))
           dimnms <- paste("dim",1:dims,sep="")
        if(! is.null(to.draw.arg)){
            if(is.character(to.draw.arg)) 
                 to.draw <- pmatch(to.draw.arg, dimnms)
            else if(is.numeric(to.draw.arg)) 
                 to.draw <- to.draw.arg
        }
        dims0 <- length(to.draw)
        nrows <- trunc(sqrt(dims0))
        ncols <- ceiling(dims0/nrows)

        if(!is.null(x.ticks)) dots$xaxt <- "n"
        if(!is.null(y.ticks)){
           y.ticks <- .fillList(y.ticks, dims0)
           dots$yaxt <- "n"
        }

        scaleY.fct <- .fillList(scaleY.fct, dims0)
        scaleY.inv <- .fillList(scaleY.inv, dims0)

        pF <- expression({})
        if(!is.null(dots[["panel.first"]])){
            pF <- .panel.mingle(dots,"panel.first")
        }
        ..panelFirst <- .fillList(pF,dims0)
        if(with.automatic.grid)
            ..panelFirst <- .producePanelFirstS(
                  ..panelFirst,x, to.draw.arg, FALSE,
                  x.ticks = x.ticks, scaleX = scaleX, scaleX.fct = scaleX.fct,
                  y.ticks = y.ticks, scaleY = scaleY, scaleY.fct = scaleY.fct)
        gridS <- if(with.automatic.grid)
                 substitute({grid <- function(...){}}) else expression({})
        pF <- vector("list",dims0)
        if(dims0>0)
           for(i in 1:dims0){
               pF[[i]] <- substitute({ gridS0
                                        pF0},
                          list(pF0=..panelFirst[[i]], gridS0=gridS))
           }

        pL <- expression({})
        if(!is.null(dots[["panel.last"]])){
            pL <- .panel.mingle(dots,"panel.last")
        }
        ..panelLast <- .fillList(pL,dims0)
        pL <- vector("list",dims0)
        if(dims0>0)
           for(i in 1:dims0)
               pL[[i]] <- if(is.null(..panelLast[[i]])) expression({}) else ..panelLast[[i]]

        dots$panel.last <- dots$panel.first <- NULL

        plotInfo$to.draw <- to.draw
        plotInfo$panelFirst <- pF
        plotInfo$panelLast <- pL
        plotInfo$gridS <- gridS


        MBRB <- matrix(rep(t(MBRB), length.out=dims0*2),ncol=2, byrow=T)
        MBRB <- MBRB * MBR.fac

        e1 <- L2Fam@distribution
        if(!is(e1, "UnivariateDistribution")) stop("not yet implemented")

        if(is(e1, "UnivariateDistribution")){
           xlim <- eval(dots$xlim)
           if(!is.null(xlim)){ 
               xm <- min(xlim)
               xM <- max(xlim)
               if(!length(xlim) %in% c(2,2*dims0))
                  stop("Wrong length of Argument xlim");
               xlim <- matrix(xlim, 2,dims0)
            }
            if(is(e1, "AbscontDistribution")){
                lower0 <- getLow(e1, eps = getdistrOption("TruncQuantile")*2)
                upper0 <- getUp(e1, eps = getdistrOption("TruncQuantile")*2)
                me <- median(e1); s <- IQR(e1)
                lower1 <- me - 6 * s
                upper1 <- me + 6 * s
                lower <- max(lower0, lower1)
                upper <- min(upper0, upper1)
                if(!is.null(xlim)){ 
                  lower <- min(lower,xm)
                  upper <- max(upper,xM)
                }
                h <- upper - lower
                if(is.null(x.vec)){
                   if(scaleX){
                      xpl <- scaleX.fct(lower - 0.1*h)
                      xpu <- scaleX.fct(upper + 0.1*h)
                      xp.vec <- seq(from = xpl, to = xpu, length = 1000)
                      x.vec <- scaleX.inv(xp.vec)
                   }else{
                      x.vec <- seq(from = lower - 0.1*h, to = upper + 0.1*h, length = 1000)
                   }
                }
                plty <- "l"
                lty <- "solid"
            }else{
                if(!is.null(x.vec)){
                   if(is(e1, "DiscreteDistribution"))
                      x.vec <- intersect(x.vec,support(e1))
                }else{
                   if(is(e1, "DiscreteDistribution")) x.vec <- support(e1)
                   else{
                      x.vec <- r(e1)(1000)
                      x.vec <- sort(unique(x.vec))
                   }
                }
                plty <- "p"
                lty <- "dotted"
                if(!is.null(dots$xlim)) x.vec <- x.vec[(x.vec>=xm) & (x.vec<=xM)]

            }
         }
         ylim <- eval(dots$ylim)
         if(!is.null(ylim)){ 
               if(!length(ylim) %in% c(2,2*dims0)) 
                  stop("Wrong length of Argument ylim"); 
               ylim <- matrix(ylim, 2,dims0)
         }

        
        if(!is.null(dots[["lty"]]))  dots["lty"] <- NULL
        if(!is.null(dots[["type"]])) dots["type"] <- NULL
        xlab <- dots$xlab; if(is.null(xlab)) xlab <- "x"
        ylab <- dots$ylab; if(is.null(ylab)) ylab <- "(partial) IC"
        dots$xlab <- dots$ylab <- NULL

        IC1 <- as(diag(dims) %*% x@Curve, "EuclRandVariable")

        mainL <- FALSE
        subL <- FALSE
        lineT <- NA


     if (hasArg(main)){
         mainL <- TRUE
         if (is.logical(main)){
             if (!main) mainL <-  FALSE
             else
                  main <- gettextf("Plot for IC %%A") ###
                          ### double  %% as % is special for gettextf
             }
         main <- .mpresubs(main)
         if (mainL) {
             if(missing(tmar))
                tmar <- 5
             if(missing(cex.inner))
                cex.inner <- .65
             lineT <- 0.6
             }
     }
     if (hasArg(sub)){
         subL <- TRUE
         if (is.logical(sub)){
             if (!sub) subL <-  FALSE
             else       sub <- gettextf("generated %%D")
                          ### double  %% as % is special for gettextf
         }
         sub <- .mpresubs(sub)
         if (subL)
             if (missing(bmar)) bmar <- 6
     }

     if(is.logical(innerL)){
        tnm  <- c(rownames(trafO))
        tnms <- if(is.null(tnm)) paste(1:dims) else 
                                 paste("'", tnm, "'", sep = "") 
        mnm <- names(L2Fam@param@main)
        mnms <- if(is.null(mnm)) NULL else paste("'", mnm, "' = ", sep = "") 
        mss  <- paste(mnms, round(L2Fam@param@main, 3), collapse=", ",sep="")
        innerT <- paste(gettextf("Component "),  tnms, 
                        gettextf("\nof"), #gettextf(" of L_2 derivative\nof"),
                        name(x)[1],
                        gettextf("\nwith main parameter ("), mss,")")
        if(!is.null(L2Fam@param@nuisance)){
            nnm <- names(L2Fam@param@nuisance)
            nnms <- if(is.null(nnm)) NULL else paste("'", nnm, "' = ", sep = "") 
            innerT <- paste(innerT,
                        gettextf("\nand nuisance parameter ("),
                        paste(nnms,round(L2Fam@param@nuisance, 3), collapse = ", "),
                        ")",
                        sep=""  )
        }
        if(!is.null(L2Fam@param@fixed)){
            fnm <- names(L2Fam@param@fixed)
            fnms <- if(is.null(fnm)) NULL else paste("'", fnm, "' = ", sep = "") 
            innerT <- paste(innerT,
                        gettextf("\nand fixed known parameter ("),
                        paste(fnms, round(L2Fam@param@fixed, 3), collapse = ", "),
                        ")",
                        sep=""  )
        }
     }else{
        innerT <- lapply(inner, .mpresubs)
        innerT <- .fillList(innerT,dims)
        if(dims0<dims){
           innerT0 <- innerT
           for(i in 1:dims0) innerT[to.draw[i]] <- innerT0[i]          
        }
     }

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
             legend <- .fillList(as.list(xc),dims0)
          }
        }


        w0 <- getOption("warn")
        options(warn = -1)
        on.exit(options(warn = w0))
        opar <- par(no.readonly = TRUE)
#        opar$cin <- opar$cra <- opar$csi <- opar$cxy <-  opar$din <- NULL
        on.exit(par(opar))
        if (!withSweave)
             devNew()
        
        parArgs <- NULL
        if(mfColRow)
           parArgs <- list(mfrow = c(nrows, ncols))

        omar <- par("mar")
        parArgs <- c(parArgs,list(mar = c(bmar,omar[2],tmar,omar[4])))

        do.call(par,args=parArgs)


        dotsT["pch"] <- dotsT["cex"] <- NULL
        dotsT["col"] <- dotsT["lwd"] <- NULL
        dotsL["cex"] <- dotsLeg["bg"] <- dotsLeg["cex"] <- NULL
        dots$ylim <- NULL

        plotInfo$resc.D <- plotInfo$resc <- vector("list", dims0)
        plotInfo$PlotLinesD <- plotInfo$PlotUsr <- vector("list", dims0)
        plotInfo$PlotArgs <- plotInfo$Axis <- vector("list", dims0)
        plotInfo$MBR <- plotInfo$Legend <- plotInfo$innerTitle <- vector("list", dims0)


        for(i in 1:dims0){
            indi <- to.draw[i]
            if(!is.null(ylim)) dots$ylim <- ylim[,i]       
            fct <- function(x) .msapply(x, IC1@Map[[indi]])
            print(xlim[,i])
            resc <-.rescalefct(x.vec, fct, scaleX, scaleX.fct,
                              scaleX.inv, scaleY, scaleY.fct[[i]], xlim[,i],
                              ylim[,i], dots)

            plotInfo$resc[[i]] <- resc
            dots <- resc$dots
            dots$xlim <- xlim[,i]
            dots$ylim <- ylim[,i]
            x.vec1 <- resc$X
            y.vec1 <- resc$Y

            finiteEndpoints <- rep(FALSE,4)
            if(scaleX){
               finiteEndpoints[1] <- is.finite(scaleX.inv(min(x.vec1, xlim[1,i])))
               finiteEndpoints[2] <- is.finite(scaleX.inv(max(x.vec1, xlim[2,i])))
            }
            if(scaleY){
               finiteEndpoints[3] <- is.finite(scaleY.inv[[i]](min(y.vec1, ylim[1,i])))
               finiteEndpoints[4] <- is.finite(scaleY.inv[[i]](max(y.vec1, ylim[2,i])))
            }


            plotInfo$PlotArgs[[i]] <- c(list(x=x.vec1, y=y.vec1, type = plty, lty = lty,
                                      xlab = .mpresubs(xlab), ylab = .mpresubs(ylab),
                                      panel.first = pF[[i]],
                                      panel.last = pL[[i]]), dots)
            do.call(plot, args=c(list(x=x.vec1, y=y.vec1, type = plty, lty = lty,
                                      xlab = .mpresubs(xlab), ylab = .mpresubs(ylab),
                                      panel.first = pF[[i]],
                                      panel.last = pL[[i]]), dots))

            plotInfo$PlotUsr[[i]] <- par("usr")
            .plotRescaledAxis(scaleX, scaleX.fct, scaleX.inv,
                              scaleY,scaleY.fct[[i]], scaleY.inv[[i]],
                              xlim[,i], ylim[,i], x.vec1, ypts = 400, n = scaleN,
                              finiteEndpoints = finiteEndpoints,
                              x.ticks = x.ticks, y.ticks = y.ticks[[i]])
            plotInfo$Axis[[i]] <- list(scaleX, scaleX.fct, scaleX.inv,
                              scaleY,scaleY.fct[[i]], scaleY.inv[[i]],
                              xlim[,i], ylim[,i], x.vec1, ypts = 400, n = scaleN,
                              finiteEndpoints = finiteEndpoints,
                              x.ticks = x.ticks, y.ticks = y.ticks[[i]])
            if(withMBR){
                MBR.i <- MBRB[i,]
                if(scaleY) MBR.i <- scaleY.fct[[i]](MBR.i)
                abline(h=MBR.i, col=col.MBR, lty=lty.MBR, lwd = lwd.MBR)
                plotInfo$MBR[[i]] <- list(h=MBR.i, col=col.MBR, lty=lty.MBR, lwd = lwd.MBR)
            }
            if(is(e1, "DiscreteDistribution")){
                x.vec1D <- seq(from = min(x.vec), to = max(x.vec), length = 1000)
                rescD <-.rescalefct(x.vec1D, fct, scaleX, scaleX.fct,
                                scaleX.inv, scaleY, scaleY.fct[[i]], xlim[,i],
                                ylim[,i], dots)
                plotInfo$resc.D[[i]] <- rescD
                x.vecD <- rescD$X
                y.vecD <- rescD$Y

                dotsL$lty <- NULL
                do.call(lines,args=c(list(x.vecD, y.vecD,
                                          lty = "dotted"), dotsL))
                plotInfo$PlotLinesD[[i]] <- c(list(x.vecD, y.vecD,
                                          lty = "dotted"), dotsL)
            }
            do.call(title,args=c(list(main = innerT[indi]), dotsT, line = lineT,
                    cex.main = cex.inner, col.main = col.inner))
            plotInfo$innerTitle[[i]] <- c(list(main = innerT[indi]), dotsT, line = lineT,
                    cex.main = cex.inner, col.main = col.inner)

            if(with.legend){
               legend(.legendCoord(legend.location[[i]], scaleX, scaleX.fct,
                        scaleY, scaleY.fct[[i]]), bg = legend.bg,
                      legend = legend[[i]], dotsLeg, cex = legend.cex*fac.leg)
               plotInfo$Legend[[i]] <- list(.legendCoord(legend.location[[i]],
                      scaleX, scaleX.fct, scaleY, scaleY.fct[[i]]), bg = legend.bg,
                      legend = legend[[i]], dotsLeg, cex = legend.cex*fac.leg)
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
        invisible(plotInfo)
    })


setMethod("plot", signature(x = "IC",y = "numeric"),
          function(x, y, ...,
          cex.pts = 1, cex.pts.fun = NULL, col.pts = par("col"),
          pch.pts = 1,
          cex.npts = 1, cex.npts.fun = NULL, col.npts = par("col"),
          pch.npts = 2,
          jitter.fac = 1, with.lab = FALSE,
          lab.pts = NULL, lab.font = NULL, alpha.trsp = NA,
          which.lbs = NULL, which.Order  = NULL, which.nonlbs = NULL,
          attr.pre = FALSE, return.Order = FALSE){

        args0 <- list(x = x, y = y, cex.pts = cex.pts, cex.pts.fun = cex.pts.fun,
             col.pts = col.pts, pch.pts = pch.pts, cex.npts = cex.npts,
             cex.npts.fun = cex.npts.fun, col.npts = col.npts, pch.npts = pch.npts,
             jitter.fac = jitter.fac, with.lab = with.lab, lab.pts = lab.pts,
             lab.font = lab.font, alpha.trsp = alpha.trsp,
             which.lbs = which.lbs, which.Order  = which.Order,
             which.nonlbs = which.nonlbs, attr.pre = attr.pre,
             return.Order = return.Order)
        mc <- match.call(call = sys.call(sys.parent(1)))
        dots <- match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)$"..."
        plotInfo <- list(call = mc, dots=dots, args=args0)

    n <- if(!is.null(dim(y))) nrow(y) else length(y)
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
       lab.pts <- if(is.null(lab.pts)) paste(1:n) else rep(lab.pts,length.out=n)
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

    sel <- .SelectOrderData(y, function(x).msapply(x, absInfo@Map[[1]]),
                            which.lbs, which.Order, which.nonlbs)
    plotInfo$sel <- sel
    plotInfo$obj <- sel$ind1

    i.d <- sel$ind
    i0.d <- sel$ind1
    n <- length(i.d)

    i.d.ns <- sel$ind.ns
    n.ns <- length(i.d.ns)

    if(attr.pre){
       col.pts <- col.pts[sel$ind]
       col.npts <- col.pts[sel$ind.ns]
       pch.npts <- pch.pts[sel$ind.ns]
       pch.pts <- pch.pts[sel$ind]
       cex.npts <- cex.pts[sel$ind.ns]
       cex.pts <- cex.pts[sel$ind]
       lab.pts <- lab.pts[sel$ind]
    }else{
       if(missing(pch.pts)) pch.pts <- 1
       if(!length(pch.pts)==n)
          pch.pts <- rep(pch.pts, length.out= n)
       if(missing(col.pts)) col.pts <- par("col")
       if(!length(col.pts)==n)
          col.pts <- rep(col.pts, length.out= n)
       if(missing(cex.pts)) cex.pts <- 1
       if(!length(cex.pts)==n)
          cex.pts <- rep(cex.pts, length.out= n)
       lab.pts <- if(is.null(lab.pts)) paste(1:n) else rep(lab.pts,length.out=n)

       if(missing(pch.npts)) pch.npts <- 1
       if(!length(pch.npts)==n.ns)
          pch.npts <- rep(pch.npts, length.out= n.ns)
       if(missing(col.npts)) col.npts <- par("col")
       if(!length(col.npts)==n.ns)
          col.npts <- rep(col.npts, length.out= n.ns)
       if(missing(cex.npts)) cex.npts <- 1
       if(!length(cex.npts)==n.ns)
          cex.npts <- rep(cex.npts, length.out= n.ns)
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
    pL <- .fillList(pL, dims0)
    if(dims0) for(i in 1:dims0){
       if(is.null(pL[[i]])) pL[[i]] <- expression({})
    }
    dots$panel.last <- NULL

    plotInfo$resc.dat <- plotInfo$resc.dat.ns <- vector("list", dims0)
    plotInfo$doPts <- plotInfo$doPts.ns <- plotInfo$doLabs <- vector("list", dims0)

    trEnv <- new.env()

    pL <- substitute({
        pI <- get("plotInfo", envir = trEnv0)

        y1 <- y0s
        ICy <- .msapply(y0s,ICMap0[[indi]])
        resc.dat <-.rescalefct(y0s, function(x) .msapply(x,ICMap0[[indi]]),
                              scaleX, scaleX.fct, scaleX.inv,
                              scaleY, scaleY.fct[[i]], xlim[,i], ylim[,i],
                              dwo0)
        pI$resc.dat[[i]] <- resc.dat
        y1 <- resc.dat$X
        ICy <- resc.dat$Y
        if(is(e1, "DiscreteDistribution")){
           if(length(ICy)) ICy <- jitter(ICy, factor = jitter.fac0) }

        y1.ns <- y0s.ns
        ICy.ns <- .msapply(y0s.ns,ICMap0[[indi]])
        resc.dat.ns <-.rescalefct(y0s.ns, function(x) .msapply(x,ICMap0[[indi]]),
                              scaleX, scaleX.fct, scaleX.inv,
                              scaleY, scaleY.fct[[i]], xlim[,i], ylim[,i],
                              dwo0)
        pI$resc.dat.ns[[i]] <- resc.dat.ns
        y1.ns <- resc.dat.ns$X
        ICy.ns <- resc.dat.ns$Y
        if(is(e1, "DiscreteDistribution"))
           {if(length(ICy.ns)) ICy.ns <- jitter(ICy.ns, factor = jitter.fac0) }

        col.pts <- if(!is.na(al0)) .msapply(col0, addAlphTrsp2col,alpha=al0) else col0
        col.npts <- if(!is.na(al0)) .msapply(col0.ns, addAlphTrsp2col,alpha=al0) else col0.ns

        cfun <- if(is.null(cexfun)) NULL else cexfun[[i]]
        cfun.ns <- if(is.null(cexnfun)) NULL else cexnfun[[i]]

        cex.l    <- .cexscale(absy0,absy0,cex=cex0, fun = cfun)   ##.cexscale in infoPlot.R
        cex.l.ns <- .cexscale(absy0.ns,absy0.ns, cex=cex0.ns, fun = cfun.ns)   ##.cexscale in infoPlot.R

        if(length(y1)){
        pI$doPts[[i]] <- c(list(y1, ICy, cex = cex.l,
                        col = col.pts, pch = pch0), dwo0)
        do.call(points, args=c(list(y1, ICy, cex = cex.l,
                        col = col.pts, pch = pch0), dwo0))
        }
        if(length(y1.ns)){
        pI$doPts.ns[[i]] <- c(list(y1.ns, ICy.ns, cex = cex.l.ns,
                        col = col.npts, pch = pch0.ns), dwo0)
        do.call(points, args=c(list(y1.ns, ICy.ns, cex = cex.l.ns,
                        col = col.npts, pch = pch0.ns), dwo0))
        }
        if(with.lab0 && length(y0s)){
           text(x = y0s, y = ICy, labels = lab.pts0,
                cex = cex.l/2, col = col0)
           pI$doLabs[[i]] <- list(x = y0s, y = ICy, labels = lab.pts0,
                cex = cex.l/2, col = col0)
        }
        assign("plotInfo", pI, envir = trEnv0)
        pL0
        }, list(pL0 = pL, ICMap0 = ICMap,
                y0s = sel$data, absy0 = sel$y,
                y0s.ns = sel$data.ns, absy0.ns = sel$y.ns,
                dwo0 = dots.without,
                cex0 = cex.pts, pch0 = pch.pts, col0 = col.pts,
                cex0.ns = cex.npts, pch0.ns = pch.npts, col0.ns = col.npts,
                with.lab0 = with.lab, lab.pts0 = lab.pts,
                al0 = alpha.trsp, jitter.fac0 = jitter.fac,
                cexfun=cex.pts.fun, cexnfun=cex.npts.fun,
                trEnv0 = trEnv
                ))
  assign("plotInfo", plotInfo, envir = trEnv)
  ret <- do.call("plot", args = c(list(x = x, panel.last = pL), dots))
  plotInfo <- get("plotInfo", envir = trEnv)
  ret$dots <- ret$args <- ret$call <- NULL
  plotInfo <- c(plotInfo, ret)
  class(plotInfo) <- c("plotInfo","DiagnInfo")
  if(return.Order) return(plotInfo)
  return(invisible(plotInfo))
})

