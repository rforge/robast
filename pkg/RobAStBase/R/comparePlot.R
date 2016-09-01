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
             pch.pts = 1, jitter.fac = 1, with.lab = FALSE,
             lab.pts = NULL, lab.font = NULL, alpha.trsp = NA,
             which.lbs = NULL, which.Order  = NULL, return.Order = FALSE,
             withSubst = TRUE){

        .mc <- match.call(call = sys.call(sys.parent(1)))
        .xc<- function(obj) as.character(deparse(.mc[[obj]]))
        xc <- c(.xc("obj1"), .xc("obj2"))
        if(!is.null(obj3)) xc <- c(xc, .xc("obj3"))
        if(!is.null(obj4)) xc <- c(xc, .xc("obj4"))
        dots <- match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)$"..."
        dotsP <- dots
        dotsLeg <- dotsT <- dotsL <- .makedotsLowLevel(dots)
        dots.points <-   .makedotsPt(dots)
        
        ncomp <- 2+ (!missing(obj3)|!is.null(obj3)) +
                    (!missing(obj4)|!is.null(obj4))

        if(missing(cex.pts)) cex.pts <- 1
        cex.pts <- rep(cex.pts, length.out= ncomp)

        if(missing(col)) col <- 1:ncomp
           else col <- rep(col, length.out = ncomp)
        if(missing(lwd))  lwd <- rep(1,ncomp)
           else lwd <- rep(lwd, length.out = ncomp)
        if(!missing(lty)) rep(lty, length.out = ncomp)
        if(missing(col.pts)) col.pts <- 1:ncomp

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
           scaleX.inv <- q(L2Fam)
        }

        trafO <- trafo(L2Fam@param)
        dims  <- nrow(trafO)
        dimm <- ncol(trafO)

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

        if(!is.null(x.ticks)) dotsP$xaxt <- "n"
        if(!is.null(y.ticks)){
           y.ticks <- .fillList(y.ticks, dims0)
           dotsP$yaxt <- "n"
        }

        if(!is.null(cex.pts.fun)){
           cex.pts.fun <- .fillList(cex.pts.fun, dims0*ncomp)
        }


        scaleY.fct <- .fillList(scaleY.fct, dims0)
        scaleY.inv <- .fillList(scaleY.inv, dims0)

        MBRB <- matrix(rep(t(MBRB), length.out=dims0*2),ncol=2, byrow=T)
        MBRB <- MBRB * MBR.fac

        distr <- L2Fam@distribution
        if(!is(distr, "UnivariateDistribution")) stop("not yet implemented")

        xlim <- dotsP$xlim <- eval(dots$xlim)
        if(!is.null(xlim)){
               xm <- min(xlim)
               xM <- max(xlim)
               xlim <- matrix(xlim, 2,dims0)
            }
        if(is(distr, "AbscontDistribution")){
            lower0 <- getLow(distr, eps = getdistrOption("TruncQuantile")*2)
            upper0 <- getUp(distr, eps = getdistrOption("TruncQuantile")*2)
            me <- median(distr); s <- IQR(distr)
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
            if(missing(lty)) lty <- "solid"
        }else{
            if(!is.null(x.vec)){
               if(is(distr, "DiscreteDistribution"))
                   x.vec <- intersect(x.vec,support(distr))
            }else{
               if(is(distr, "DiscreteDistribution")) x.vec <- support(distr) else{
                   x.vec <- r(distr)(1000)
                   x.vec <- sort(unique(x.vec))
               }
            }
            plty <- "p"
            if(missing(lty)) lty <- "dotted"
            if(!is.null(xlim)) x.vec <- x.vec[(x.vec>=xm) & (x.vec<=xM)]
        }
        ylim <- eval(dots$ylim)
        if(!is.null(ylim)){
               if(! length(ylim) %in% c(2,2*dims0))
                  stop("Wrong length of Argument ylim");
               ylim <- matrix(ylim, 2,dims0)
        }
        dots$ylim <- dots$xlim <- NULL

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

      lineT <- NA

      
      .mpresubs <- if(withSubst){function(inx)
            .presubs(inx, c(paste("%C",1:ncomp,sep=""),
                                     "%D",
                                    paste("%A",1:ncomp,sep="")),
                  c(as.character(class(obj1)[1]),
                    as.character(class(obj2)[1]),
                    if(is.null(obj3))NULL else as.character(class(obj3)[1]),
                    if(is.null(obj4))NULL else as.character(class(obj4)[1]),
                    as.character(date()),
                    xc))} else function(inx)inx

        mainL <- FALSE
        if (hasArg(main)){
            mainL <- TRUE
            if (is.logical(main)){
                if (!main) mainL <- FALSE else
                     main <- paste(gettextf("Plot for ICs"),
                                paste("%A", 1:ncomp, sep="", collapse=", "),
                                sep=" ")
            }
            main <- .mpresubs(main)
            if (mainL) {
                if(missing(tmar)) tmar <- 5
                if(missing(cex.inner)) cex.inner <- .65
                lineT <- 0.6
            }
        }
        subL <- FALSE
        if (hasArg(sub)){
            subL <- TRUE
            if (is.logical(sub)){
                if (!sub) subL <-  FALSE  else sub <- gettextf("generated %%D")
            }
            sub <- .mpresubs(sub)
            if (subL)  if (missing(bmar)) bmar <- 6
        }

        .mknam <- function(val, iP = "", txt){
            nm <- names(val)
            nms <- if(is.null(nm)) NULL else paste("'", nm, "' = ", sep = "")
            iP <- paste(iP, "\n", gettext(txt), " (",
                        paste(nms, round(val,3), collapse = ", "), ")", sep = "")
        }

        innerParam <- .mknam(L2Fam@param@main, "", "with main parameter")
        if(!is.null(L2Fam@param@nuisance))
            innerParam <- .mknam(L2Fam@param@nuisance, innerParam,
                                 "and nuisance parameter")
        if(!is.null(L2Fam@param@fixed))
            innerParam <- .mknam(L2Fam@param@fixed, innerParam,
                             "and fixed known parameter")

        if(!is.logical(inner)){
            if(!is.list(inner))
                inner <- as.list(inner)
            innerT <- .fillList(inner,dims)
            if(dims0<dims){
               innerT0 <- innerT
               for(i in 1:dims0) innerT[to.draw[i]] <- innerT0[i]
            }
            innerL <- TRUE
        }else{if(any(is.na(inner))||any(!inner)) {
             innerT <- as.list(rep("",dims)); innerL <- FALSE
            }else{innerL <- TRUE
                  tnm  <- c(rownames(trafO))
                  tnms <- if(is.null(tnm)) paste(1:dims) else
                                           paste("'", tnm, "'", sep = "")
                  innerT <- as.list(paste(paste(gettext("Component "),  tnms,
                                   gettext(" of (partial) IC\nfor "),
                                   name(L2Fam)[1], sep =""), innerParam))
               }
          }


        w0 <- getOption("warn"); options(warn = -1); on.exit(options(warn = w0))

        opar <- par(no.readonly = TRUE)
        if(mfColRow){ on.exit(par(opar)); par(mfrow = c(nrows, ncols)) }

        if(is(distr, "DiscreteDistribution"))
                x.vecD <- seq(from = min(x.vec), to = max(x.vec), length = 1000)

        dotsT$main <- dotsT$cex.main <- dotsT$col.main <- dotsT$line <- NULL

        pF <- expression({})
        if(!is.null(dots[["panel.first"]])){
            pF <- .panel.mingle(dots,"panel.first")
        }
        ..panelFirst <- .fillList(pF,dims0)
        if(with.automatic.grid)
           ..panelFirst <- .producePanelFirstS(
                ..panelFirst,obj1 , to.draw.arg, FALSE,
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
        dots$panel.first <- NULL
        pL <- expression({})
        if(!is.null(dots[["panel.last"]])){
            pL <- .panel.mingle(dots,"panel.last")
        }
        pL <- .fillList(pL, dims0)
        dotsP$panel.last <- NULL

        sel1 <- sel2 <- sel3 <- sel4 <- NULL
        if(!is.null(data)){
            n <- if(!is.null(dim(data))) nrow(data) else length(data)
            lab.pts <- rep(lab.pts, length.out=n)

            absInfoEval <- function(x,IC){
                  QF <- ID
                  if(is(IC,"ContIC") & dims>1 ){
                     if (is(normtype(IC),"QFNorm"))
                          QF <- QuadForm(normtype(IC))
                  }
                  absInfo.f <- t(IC) %*% QF %*% IC
                  return(sapply(x, absInfo.f@Map[[1]]))
            }
            def.sel <- function(IC){
                 fct.aI <- function(x) absInfoEval(x,IC)
                 return(.SelectOrderData(data, fct.aI, which.lbs, which.Order))}
                 
            sel1 <- def.sel(IC1); sel2 <- def.sel(IC2)
            selAlly <- c(sel1$y,sel2$y)

            if(is(obj3, "IC")){ sel3 <- def.sel(IC3)
                                selAlly <- c(selAlly,sel3$y)
                              }
            if(is(obj4, "IC")){ sel4 <- def.sel(IC4)
                                selAlly <- c(selAlly,sel4$y)
                              }

            dots.points <- .makedotsLowLevel(dots)
            dots.points$col <- dots.points$cex <- dots.points$pch <- NULL
            alp.v <- rep(alpha.trsp,length.out = ncomp)

            pL <- substitute({
                 doIt <- function(sel.l,fct.l,j.l){
                     rescd <- .rescalefct(sel.l$data, fct.l, scaleX, scaleX.fct,
                                   scaleX.inv, scaleY, scaleY.fct[[i]], xlim[,i],
                                   ylim[,i], dotsP)
                     if(is(distr, "DiscreteDistribution"))
                        rescd$Y <- jitter(rescd$Y, factor = jitter.fac0[j.l])
                     i.l <- sel.l$ind
                     n.l <- length(i.l)
                     pch.pts.l <- rep(pch0, length.out=n.l)
                     lab.pts.l <- if(is.null(lab0)) paste(i.l) else lab0[i.l]

                     col.l <- if(is.na(al0[j.l])) col0[j.l] else
                                 addAlphTrsp2col(col0[j.l], al0[j.l])

                     cfun <- if(is.null(cexfun)) NULL else cexfun[[(i-1)*ncomp+j.l]]

                     cex.l <- .cexscale(sel.l$y,selAlly,cex=cex0[j.l], fun = cfun)   ##.cexscale in infoPlot.R
                     do.call(points, args=c(list(rescd$X, rescd$Y, cex = cex.l,
                             col = col.l, pch = pch.pts.l), dwo0))
                     if(with.lab0)
                        text(rescd$X, rescd$Y, labels = lab.pts.l,
                             cex = cex.l/2, col = col.l)
                 }
                 doIt(sel1,fct1,1);  doIt(sel2,fct2,2)
                 if(is(obj3, "IC")) doIt(sel3,fct3,3)
                 if(is(obj4, "IC")) doIt(sel4,fct4,4)
                 pL0
              }, list(pL0 = pL, cex0 = cex.pts, pch0 = pch.pts, col0 = col.pts,
                      jitter.fac0 = jitter.fac, dwo0 = dots.points, al0 = alp.v,
                      with.lab0 = with.lab, lab0 = lab.pts, cexfun=cex.pts.fun)
            )
        }

        dotsP$type <- dotsP$lty <- dotsP$col <- dotsP$lwd <- NULL
        dotsP$xlab <- dotsP$ylab <- NULL

        for(i in 1:dims0){
            indi <- to.draw[i]
            if(!is.null(ylim)) dotsP$ylim <- ylim[,i]

            fct1 <- function(x) sapply(x, IC1@Map[[indi]])

            resc.args <- c(list(x.vec, "fc"=fct1, scaleX, scaleX.fct,
                                scaleX.inv, scaleY, scaleY.fct[[i]], xlim[,i],
                                ylim[,i], dotsP))
            resc1 <- do.call(.rescalefct, resc.args)
            resc.args$fc <- fct2 <- function(x) sapply(x, IC2@Map[[indi]])
            resc2 <- do.call(.rescalefct, resc.args)

            dotsP <- resc1$dots
            matp  <- cbind(resc1$Y, resc2$Y)

            if(is(obj3, "IC")){
                resc.args$fc <- fct3 <- function(x) sapply(x, IC3@Map[[indi]])
                resc3 <- do.call(.rescalefct, resc.args)
                matp  <- cbind(matp,resc3$Y)
            }
            if(is(obj4, "IC")){
                resc.args$fc <- fct4 <- function(x) sapply(x, IC4@Map[[indi]])
                resc4 <- do.call(.rescalefct, resc.args)
                matp  <- cbind(matp,resc4$Y)
            }

            ym <- min(matp,na.rm=T)
            yM <- max(matp,na.rm=T)
            y0 <- matp[,1]
            y0[1:2] <- c(ym,yM)

            finiteEndpoints <- rep(FALSE,4)
            if(scaleX){
               finiteEndpoints[1] <- is.finite(scaleX.inv(min(x.vec, xlim[1],na.rm=TRUE)))
               finiteEndpoints[2] <- is.finite(scaleX.inv(max(x.vec, xlim[2],na.rm=TRUE)))
            }
            if(scaleY){
               finiteEndpoints[3] <- is.finite(scaleY.inv[[i]](min(ym, ylim[1,i],na.rm=TRUE)))
               finiteEndpoints[4] <- is.finite(scaleY.inv[[i]](max(yM, ylim[2,i],na.rm=TRUE)))
            }

            do.call(plot, args=c(list(x = resc1$X, y = y0,
                 type = "n", xlab = .mpresubs(xlab), ylab = .mpresubs(ylab),
                 lty = lty[1], col = addAlphTrsp2col(col[1],0),
                 lwd = lwd[1]), dotsP, list(panel.last = pL[[i]], panel.first=pF[[i]])))
            if(plty=="p")
               do.call(matpoints, args = c(list( x = resc1$X, y = matp,
                    col = col), dots.points))

            do.call(matlines, args = c(list( x = resc1$X, y = matp,
                    lty = lty, col = col, lwd = lwd), dotsL))


            .plotRescaledAxis(scaleX, scaleX.fct, scaleX.inv,
                              scaleY,scaleY.fct[[i]], scaleY.inv[[i]], xlim[,i],
                              ylim[,i], resc1$X, ypts = 400, n = scaleN,
                              finiteEndpoints = finiteEndpoints,
                              x.ticks = x.ticks, y.ticks = y.ticks[[i]])
            if(withMBR){
                MBR.i <- MBRB[i,]
                if(scaleY) MBR.i <- scaleY.fct[[i]](MBR.i)
                abline(h=MBR.i, col=col.MBR, lty=lty.MBR, lwd = lwd.MBR)
            }

            if(is(distr, "DiscreteDistribution")){
                 rescD.args <- c(list(x.vecD, "fc"=fct1, scaleX, scaleX.fct,
                                scaleX.inv, scaleY, scaleY.fct[[i]], xlim[,i],
                                ylim[,i], dotsP))
                 resc1D <- do.call(.rescalefct, rescD.args)
                 rescD.args$fc <- fct2
                 resc2D <- do.call(.rescalefct, rescD.args)
                 matpD  <- cbind(resc1D$Y, resc2D$Y)
                 if(is(obj3, "IC")){
                    rescD.args$fc <- fct3
                    resc3D <- do.call(.rescalefct, rescD.args)
                    matpD  <- cbind(matpD, resc3D$Y)
                 }
                 if(is(obj4, "IC")){
                    rescD.args$fc <- fct4
                    resc4D <- do.call(.rescalefct, rescD.args)
                    matpD  <- cbind(matpD, resc4D$Y)
                 }
                 do.call(matlines, c(list(resc1D$X, matpD, lty = lty,
                         col = col, lwd = lwd), dotsL))
            }

           if(innerL)
              do.call(title, args=c(list(main = innerT[[indi]]),  dotsT,
                      line = lineT, cex.main = cex.inner, col.main = col.inner))
        }

        if(with.legend){
           if(is.null(legend)) legend <- xc
           legend(.legendCoord(legend.location, scaleX, scaleX.fct,
                        scaleY, scaleY.fct[[i]]), col = col, bg = legend.bg,
                      legend = legend, dotsLeg, cex = legend.cex)
        }

        cex.main <- if(!hasArg(cex.main)) par("cex.main") else dots$"cex.main"
        col.main <- if(!hasArg(col.main)) par("col.main") else dots$"col.main"
        if (mainL)
            mtext(text = main, side = 3, cex = cex.main, adj = .5,
                  outer = TRUE, padj = 1.4, col = col.main)

        cex.sub <- if(!hasArg(cex.sub)) par("cex.sub") else dots$"cex.sub"
        col.sub <- if(!hasArg(col.sub)) par("col.sub") else dots$"col.sub"
        if (subL)
            mtext(text = sub, side = 1, cex = cex.sub, adj = .5,
                  outer = TRUE, line = -1.6, col = col.sub)

        if(return.Order) return(list(obj1=sel1$ind1, obj2=sel2$ind1,
                                     obj3=sel3$ind1, obj4=sel4$ind1))
        invisible()
    })

