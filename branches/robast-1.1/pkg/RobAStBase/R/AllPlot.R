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

################################################################################
## 1. preparation: fingle around with arguments:
################################################################################
#  1.1 read out dots, object, L2Fam, scaling
################################################################################
        mc <- match.call(call = sys.call(sys.parent(1)))
        xc <- mc$x
        xcc <- as.character(deparse(xc))
        dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."


        L2Fam <- eval(x@CallL2Fam)
        if(missing(scaleX.fct)){
           scaleX.fct <- p(L2Fam)
           scaleX.inv <- q(L2Fam)
        }
        if(missing(scaleY.fct)){
           scaleY.fct <- pnorm
           scaleY.inv <- qnorm
        }

################################################################################
#  1.2  clean up dots arguments
################################################################################
        dotsLeg <- dotsT <- dotsL <- .makedotsLowLevel(dots)

################################################################################
#  1.3 parameter trafo and dimensions of the panels
################################################################################
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

################################################################################
#  1.4 preparation of cex, scaling  per panel, legend
################################################################################
        if(!is.null(x.ticks)) dots$xaxt <- "n"
        if(!is.null(y.ticks)){
           y.ticks <- .fillList(y.ticks, dims0)
           dots$yaxt <- "n"
        }

        scaleY.fct <- .fillList(scaleY.fct, dims0)
        scaleY.inv <- .fillList(scaleY.inv, dims0)

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

################################################################################
#  1.5  prepare titles
################################################################################
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


################################################################################
#  2. pre- and posthooks per panel (panel last -> pL, panel first -> pF)
################################################################################
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
        if(dims0>0)
           pL <- vector("list",dims0)
           for(i in 1:dims0)
               pL[[i]] <- if(is.null(..panelLast[[i]])) expression({}) else ..panelLast[[i]]


        dots$panel.last <- dots$panel.first <- NULL


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


################################################################################
#  2.2. preparation: what is to be done "on exit"
################################################################################

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

################################################################################
#  3. creating the panel plots
################################################################################
        icpInfo <- vector("list",0)
        icpInfo$panels <- vector("list",dims0)
        for(i in 1:dims0){
            indi <- to.draw[i]
            if(!is.null(ylim)) dots$ylim <- ylim[,i]       
            fct <- function(x) sapply(x, IC1@Map[[indi]])
            print(xlim[,i])
            resc <-.rescalefct(x.vec, fct, scaleX, scaleX.fct,
                              #scaleX.inv, 
                              scaleY, scaleY.fct[[i]], xlim[,i],
                              ylim[,i], dots)
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


            plot.args <- c(list(x=x.vec1, y=y.vec1, type = plty, lty = lty,
                                      xlab = .mpresubs(xlab), ylab = .mpresubs(ylab),
                                      panel.first = pF[[i]],
                                      panel.last = pL[[i]]), dots)
            do.call(plot, args=plot.args)
            icpInfo$panels[[i]]$plot.args <- plot.args
            rm(plot.args)
            
            .plotRescaledAxis.args <- list(scaleX, scaleX.fct, scaleX.inv,
                              scaleY,scaleY.fct[[i]], scaleY.inv[[i]],
                              xlim[,i], ylim[,i], x.vec1, ypts = 400, n = scaleN,
                              finiteEndpoints = finiteEndpoints,
                              x.ticks = x.ticks, y.ticks = y.ticks[[i]])
            do.call(.plotRescaledAxis, args=.plotRescaledAxis.args)
            icpInfo$panels[[i]]$.plotRescaledAxis.args <- .plotRescaledAxis.args
            rm(.plotRescaledAxis.args)                  

            if(withMBR){
                MBR.i <- MBRB[i,]
                if(scaleY) MBR.i <- scaleY.fct[[i]](MBR.i)
                MBR.args <- list(h=MBR.i, col=col.MBR, lty=lty.MBR, lwd = lwd.MBR)
                do.call(abline, args=MBR.args)
                icpInfo$panels[[i]]$MBR.args <- MBR.args
                rm(MBR.args)
                
            }
            if(is(e1, "DiscreteDistribution")){
                x.vec1D <- seq(from = min(x.vec), to = max(x.vec), length = 1000)
                rescD <-.rescalefct(x.vec1D, fct, scaleX, scaleX.fct,
                                #scaleX.inv, 
                                scaleY, scaleY.fct[[i]], xlim[,i],
                                ylim[,i], dots)
                x.vecD <- rescD$X
                y.vecD <- rescD$Y

                dotsL$lty <- NULL
                lines.args <- c(list(x.vecD, y.vecD, lty = "dotted"), dotsL)
                do.call(lines, args = lines.args)
                icpInfo$panels[[i]]$lines.args <- lines.args
                rm(lines.args)
            }
            
            title.args <- c(list(main = innerT[indi]), dotsT, line = lineT,
                            cex.main = cex.inner, col.main = col.inner)
            do.call(title, args=title.args)
            icpInfo$panels[[i]]$title.args <- title.args
            rm(title.args)        

            if(with.legend){
               legend.args <- c(list(.legendCoord(legend.location[[i]], scaleX, 
                        scaleX.fct, scaleY, scaleY.fct), bg = legend.bg,
                        legend = legend[[i]]), dotsLeg, cex = legend.cex*fac.leg)
               do.call(graphics::legend, args=legend.args)
               icpInfo$panels[[i]]$legend.args <- legend.args
               rm(legend.args)                      
            }          

        }
################################################################################
#  4. outer titles
################################################################################
        cex.main <- if(!hasArg(cex.main)) par("cex.main") else dots$"cex.main"
        col.main <- if(!hasArg(col.main)) par("col.main") else dots$"col.main"
        if (mainL){
            main.args <- list(text = main, side = 3, cex = cex.main, adj = .5,
                  outer = TRUE, padj = 1.4, col = col.main)
            do.call(mtext, args=main.args)
            icpInfo$main.args <- main.args 
            rm(main.args)                
        }


        cex.sub <- if(!hasArg(cex.sub)) par("cex.sub") else dots$"cex.sub"
        col.sub <- if(!hasArg(col.sub)) par("col.sub") else dots$"col.sub"
        if (subL){
            sub.args <- list(text = sub, side = 1, cex = cex.sub, adj = .5,
                  outer = TRUE, line = -1.6, col = col.sub)
            do.call(mtext, args=sub.args)
            ipInfo$sub.args <- sub.args 
            rm(sub.args)
        }

  class(icpInfo) <- c("ICPlotInfo","DiagnInfo")
  retv <- list(call=mc, ICPlotInfo = icpInfo)       
  invisible(return(retv))
})


setMethod("plot", signature(x = "IC",y = "numeric"),
          function(x, y, ..., 
####             
             cex.pts = 1, 
             cex.pts.fun = NULL, 
             col.pts = par("col"),
             pch.pts = 1,              
             jit.fac = 1, 
             jit.tol = .Machine$double.eps, 
             with.lab = FALSE,       
             lab.pts = NULL, lab.col = par("col"), lab.font = NULL, lab.adj = NULL,            
             alpha.trsp = NA,          
             which.lbs = NULL, which.Order = NULL, which.nonlbs = NULL, return.Order = FALSE,             
             draw.nonlbl = TRUE,  ## should non-labelled observations also be drawn?             
             cex.nonlbl = 0.3,    ## character expansion(s) for non-labelled observations
             cex.nonlbl.fun = NULL, ## like cex.pts.fun for non-labelled observations
             col.nonlbl = par("col"),   
             pch.nonlbl = "."    ## plotting symbol(s) for non-labelled observations
          ){

    mc <- match.call(call = sys.call(sys.parent(1)))
    dots <- match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)$"..."

    n <- if(!is.null(dim(y))) nrow(y) else length(y)
    lab.pts <- if(is.null(lab.pts)) paste(1:n) else rep(lab.pts,n)

    L2Fam <- eval(x@CallL2Fam)
    trafO <- trafo(L2Fam@param)
    dims <- nrow(trafO)
    dimm <- length(L2Fam@param)
    QF <- diag(dims)


################################################################################
#  2.1. preparation: norm, function to evaluate it for both robust and classic
################################################################################
    if(is(x,"ContIC") & dims>1 )
      {if (is(normtype(x),"QFNorm")) QF <- QuadForm(normtype(x))}

    IC1 <- as(diag(dims) %*% x@Curve, "EuclRandVariable")
    absInfo <- t(IC1) %*% QF %*% IC1
    ICMap <- IC1@Map

    sel <- .SelectOrderData(y, function(x)sapply(x, absInfo@Map[[1]]),
                            which.lbs, which.Order, which.nonlbs)
    i.d <- sel[["ind"]]
    i0.d <- sel[["ind1"]]
    n.s <- length(i.d)

    i.d.ns <- sel[["ind.ns"]]
    n.ns <- length(i.d.ns)

    if(length(col.pts)==n){
       col.pts0 <- col.pts
       col.pts <- col.pts0[i.d]
       col.nonlbl <- if(draw.nonlbl && n.ns > 0 ) col.pts0[i.d.ns] else NULL      
    }else{
       col.pts <- rep(col.pts,length.out=n.s)
       col.nonlbl <- if(draw.nonlbl && n.ns > 0 ) 
                        rep(col.nonlbl,length.out=n.ns) else NULL
    }
    if(length(pch.pts)==n){
       pch.pts0 <- pch.pts
       pch.pts <- pch.pts0[i.d]
       pch.nonlbl <- if(draw.nonlbl && n.ns > 0 ) pch.pts0[i.d.ns] else NULL     
    }else{
       pch.pts <- rep(pch.pts,length.out=n.s)
       pch.nonlbl <- if(draw.nonlbl && n.ns > 0 ) 
                        rep(pch.nonlbl,length.out=n.ns) else NULL      
    }
    if(length(cex.pts)==n){
       cex.pts0 <- cex.pts
       cex.pts <- cex.pts0[i.d]
       cex.nonlbl <- if(draw.nonlbl && n.ns > 0 ) cex.pts0[i.d.ns] else NULL     
    }else{
       cex.pts <- rep(cex.pts,length.out=n.s)
       cex.nonlbl <- if(draw.nonlbl && n.ns > 0 ) 
                        rep(cex.nonlbl,length.out=n.ns) else NULL     
    }
    if(length(lab.col)==n){
       lab.col <- lab.col[i.d]
    }else{
       lab.col <- rep(lab.col,length.out=n.s)
    }

    dots.without <- dots
    dots.without$col <- dots.without$cex <- dots.without$pch <- NULL

    dims0 <- .getDimsTD(L2Fam,dots[["to.draw.arg"]])
    alp.v <- matrix(rep(alpha.trsp, length.out = (n.s+n.ns)*dims0),
                    (n.s+n.ns),dims0)
    alp.v.s <- alp.v[i.d,,drop=FALSE]
    alp.v.ns <- if(draw.nonlbl && n.ns > 0 ) alp.v[i.d.ns,,drop=FALSE] else NULL

    if(!is.null(cex.pts.fun)){
       cex.pts.fun <- .fillList(cex.pts.fun, dims0)
    }
    if(!is.null(cex.nonlbl.fun)&& draw.nonlbl && n.ns > 0 ) {
       cex.nonlbl.fun <- .fillList(cex.nonlbl.fun, dims0)
    }


    lab.adj <- if(is.null(lab.adj)){ matrix(0.5,n,dims0)
               }else{
                  if(length(lab.adj)%in%c(1,2))
                     lab.adj <- rep(lab.adj, length.out=2*dims0) 
                  if(length(lab.adj)==2*dims0){
                     lab.adj <- matrix( rep(lab.adj,
                                 times=rep(n,times=2*dims0)), n,2*dims0) 
                  }else{
                      if(!is.matrix(lab.adj))
                          lab.adj <- matrix(rep(lab.adj, length.out=n),n,1)
                      if(ncol(lab.adj)==1)
                          lab.adj <- cbind(lab.adj,lab.adj)  
                      if(ncol(lab.adj)==2) 
                         lab.adj <- lab.adj[,rep(1:2,dims0)] 
                      if(ncol(lab.adj)!=2*dims0) 
                         stop("Wrong number of columns in arg 'lab.adj'.")
                  }             
               }
    
################################################################################
#  2.5 plotting in data : preparation
################################################################################
    pL <- expression({})
    if(!is.null(dots$panel.last))
        pL <- .panel.mingle(dots,"panel.last")
    pL <- .fillList(pL, dims0)
    if(dims0) for(i in 1:dims0){
       if(is.null(pL[[i]])) pL[[i]] <- expression({})
    }
    dots$panel.last <- NULL


################################################################################
#  2.6 inserting the code to plot in data into panel last
################################################################################
    pL <- substitute({
        y1 <- y0s
        ICy <- sapply(y0s,ICMap0[[indi]])
        #print(xlim[,i])
        resc.dat <-.rescalefct(y0s, function(x) sapply(x,ICMap0[[indi]]),
                              scaleX, scaleX.fct, #scaleX.inv,
                              scaleY, scaleY.fct[[i]], xlim[,i], ylim[,i],
                              dwo0)
        y1 <- resc.dat$X
        ICy <- resc.dat$Y
        sel <- resc.dat$idx

        if(is(e1, "DiscreteDistribution")){
           ICy <- jitter(ICy, factor = jit.fac0)
        }else{if(any(.isReplicated(ICy, jit.tol0))&&jit.fac0>0)
                 ICy <- jitter(ICy, factor = jit.fac0)
        }

        al0.si <- al0.s[sel,i]
        col.s <- .alphTrspWithNA(col0[sel],al0.si)

        cfun <- if(is.null(cexfun)) NULL else cexfun[[i]]
        cex.s <- .cexscale(resc.dat$scy,resc.dat$scy,cex=cex0, fun=cfun)

        do.call(points, args=c(list(y1, ICy, cex = cex.s,
                        col = col.s, pch = pch0[sel]), dwo0))
        if(with.lab0){
           for(kk in 1:length(y0s))
               text(x = y0s[kk], y = ICy[kk], labels = lab.pts0[kk],
                cex = log(absy0[kk]+1)*1.5*cex0, col = lab.col0[sel],
                font= lab.ft, adj=lab.ad0[kk,(i-1)*2+(1:2)])
        }

        if(dononlb){
            resc.dat.ns <-.rescalefct(y0s.ns, function(x) sapply(x,ICMap0[[indi]]),
                                  scaleX, scaleX.fct,# scaleX.inv,
                                  scaleY, scaleY.fct[[i]], xlim[,i], ylim[,i],
                                  dwo0)
            y1.ns <- resc.dat.ns$X
            ICy.ns <- resc.dat.ns$Y
            sel.ns <- resc.dat.ns$idx
    
            if(is(e1, "DiscreteDistribution")){
               ICy.ns <- jitter(ICy.ns, factor = jit.fac0)
            }else{if(any(.isReplicated(ICy.ns, jit.tol0))&&jit.fac0>0)
                     ICy.ns <- jitter(ICy.ns, factor = jit.fac0)
            }
    
            al0.nsi <- al0.ns[sel.ns,i]
            col.ns <- .alphTrspWithNA(col0.ns[sel.ns],al0.nsi)
    
            cfun.ns <- if(is.null(cexfun.ns)) NULL else cexfun.ns[[i]]
            cex.ns <- .cexscale(resc.dat.ns$scy,resc.dat.ns$scy,cex=cex0.ns, fun=cfun.ns)
    
    
            do.call(points, args=c(list(y1.ns, ICy.ns, cex = cex.ns,
                            col = col.ns, pch = pch0.ns[sel.ns]), dwo0))
        }
        pL0
        }, list(pL0 = pL, ICMap0 = ICMap, y0s = sel[["data"]], absy0 = sel$y,
                y0s.ns = sel[["data.ns"]], dwo0 = dots.without, 
                cex0 = cex.pts, pch0 = pch.pts, col0 = col.pts, 
                cex0.ns = cex.nonlbl, pch0.ns = pch.nonlbl, col0.ns = col.nonlbl, 
                cexfun = cex.pts.fun ,cexfun.ns=cex.nonlbl.fun,
                with.lab0 = with.lab, lab.pts0 = lab.pts[i.d],
                al0.s = alp.v.s, al0.ns = alp.v.ns, lab.ft=lab.font,
                jit.fac0 = jit.fac, jit.tol0=jit.tol, lab.ad0=lab.adj,
                dononlb = draw.nonlbl&(n.ns>0) 
                ))

  plotArgs <- c(list(x = x, panel.last = pL), dots)
  retvPlot <- do.call("plot", args = plotArgs)
  retvPlot$call <- NULL

  class(plotArgs) <- c("ICPlotInfo","DiagnInfo")
  retv <- list(call=mc, ICPlotInfo = c(retvPlot, dataArgs=plotArgs))       

  if(return.Order){ 
     retOrder <- i0.d
     retv$retOrder <- retOrder         
  }
  invisible(return(retv))
})

