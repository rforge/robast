.oldcomparePlot <- getMethod("comparePlot", signature("IC","IC"))

setMethod("comparePlot", signature("IC","IC"),
    function(obj1,obj2, obj3 = NULL, obj4 = NULL, data = NULL,
             ..., withSweave = getdistrOption("withSweave"),
             forceSameModel =  FALSE,
             main = FALSE, inner = TRUE, sub = FALSE,
             col = par("col"), lwd = par("lwd"), lty,
             col.inner = par("col.main"), cex.inner = 0.8,
             bmar = par("mar")[1], tmar = par("mar")[3],
             with.automatic.grid = TRUE,
             with.legend = FALSE, legend = NULL, legend.bg = "white",
             legend.location = "bottomright", legend.cex = 0.8,
             withMBR = FALSE, MBRB = NA, MBR.fac = 2, col.MBR = par("col"),
             lty.MBR = "dashed", lwd.MBR = 0.8,  n.MBR = 10000,
             x.vec = NULL, scaleX = FALSE, scaleX.fct, scaleX.inv,
             scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
             scaleN = 9, x.ticks = NULL, y.ticks = NULL,
             mfColRow = TRUE, to.draw.arg = NULL,
             cex.pts = 1, cex.pts.fun = NULL,
             col.pts = par("col"),  pch.pts = 1,
             cex.npts = 1, cex.npts.fun = NULL,
             col.npts = par("col"),  pch.npts = 2,
             jitter.fac = 1, with.lab = FALSE,
             lab.pts = NULL, lab.font = NULL, alpha.trsp = NA,
             which.lbs = NULL, which.Order  = NULL, which.nonlbs = NULL,
             attr.pre = FALSE, return.Order = FALSE, withSubst = TRUE){

        args0 <- list(obj1 = obj1, obj2 = obj2, obj3 = obj3, obj4 = obj4,
             data = data, withSweave = withSweave, forceSameModel = forceSameModel,
             main = main, inner = inner, sub = sub, col = col, lwd = lwd,
             lty = if(!missing(lty)) lty else NULL,
             col.inner = col.inner, cex.inner = cex.inner,
             bmar = bmar, tmar = tmar, with.automatic.grid = with.automatic.grid,
             with.legend = with.legend, legend = legend, legend.bg = legend.bg,
             legend.location = legend.location, legend.cex = legend.cex,
             withMBR = withMBR, MBRB = MBRB, MBR.fac = MBR.fac, col.MBR = col.MBR,
             lty.MBR = lty.MBR, lwd.MBR = lwd.MBR,  n.MBR = n.MBR,
             x.vec = x.vec, scaleX = scaleX,
             scaleX.fct = if(!missing(scaleX.fct)) scaleX.fct else NULL,
             scaleX.inv = if(!missing(scaleX.inv)) scaleX.inv else NULL,
             scaleY = scaleY, scaleY.fct = scaleY.fct,
             scaleY.inv = scaleY.inv, scaleN = scaleN, x.ticks = x.ticks,
             y.ticks = y.ticks, mfColRow = mfColRow, to.draw.arg = to.draw.arg,
             cex.pts = cex.pts, cex.pts.fun = cex.pts.fun, col.pts = col.pts,
             pch.pts = pch.pts, cex.npts = cex.npts, cex.npts.fun = cex.npts.fun,
             col.npts = col.npts, pch.npts = pch.npts,
             jitter.fac = jitter.fac, with.lab = with.lab, lab.pts = lab.pts,
             lab.font = lab.font, alpha.trsp = alpha.trsp,
             which.lbs = which.lbs, which.Order  = which.Order,
             which.nonlbs = which.nonlbs, attr.pre = attr.pre,
             return.Order = return.Order, withSubst = withSubst)

        .mc <- match.call(call = sys.call(sys.parent(1)))
        dots <- match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)$"..."
        plotInfo <- list(call = .mc, dots=dots, args=args0)

        mcl <- match.call(call = sys.call(sys.parent(1)), expand.dots = TRUE)

        L2Fam <- eval(obj1@CallL2Fam); trafO <- trafo(L2Fam@param)
        dims  <- nrow(trafO); to.draw <- 1:dims
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

        MBRB <- matrix(rep(t(MBRB), length.out=dims0*2),ncol=2, byrow=T)
        if(withMBR && all(is.na(MBRB))){
           ICmbr <- try(getStartIC(model = L2Fam, risk = MBRRisk()), silent=TRUE)
           if(is(ICmbr,"try-error")){
              robModel <- InfRobModel(center = L2Fam, neighbor =
                             ContNeighborhood(radius = 0.5))
              ICmbr <- try(optIC(model = robModel, risk = asBias()), silent=TRUE)
           }
           if(!is(ICmbr,"try-error"))
              MBRB <- .getExtremeCoordIC(ICmbr, distribution(L2Fam), to.draw,
                                         n=n.MBR)
           else withMBR <- FALSE
        }
        mcl$MBRB <- MBRB
        mcl$withMBR <- withMBR
        ret <- do.call(.oldcomparePlot, as.list(mcl[-1]),
                        envir=parent.frame(2))
        ret$dots <- ret$args <- ret$call <- NULL
        plotInfo <- c(plotInfo, ret)
        class(plotInfo) <- c("plotInfo","DiagnInfo")
        if(return.Order){ whichRet <- names(plotInfo) %in% c("obj1","obj2","obj3","obj4")
                            return(plotInfo[whichRet])}
        return(invisible(plotInfo))
      })
