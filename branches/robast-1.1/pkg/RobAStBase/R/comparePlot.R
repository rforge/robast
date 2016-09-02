setMethod("comparePlot", signature("IC","IC"),
    function(obj1,obj2, obj3 = NULL, obj4 = NULL, data = NULL,
             ..., withSweave = getdistrOption("withSweave"),
             forceSameModel = FALSE,
             main = FALSE, inner = TRUE, sub = FALSE,
             col = par("col"), lwd = par("lwd"), lty,
             col.inner = par("col.main"), cex.inner = 0.8,
             bmar = par("mar")[1], tmar = par("mar")[3],

             with.automatic.grid = TRUE, ##new
             
             with.legend = FALSE, 
             legend = NULL,       ##new
             legend.bg = "white",
             legend.location = "bottomright", 
             legend.cex = 0.8,
             
             withMBR = FALSE, MBRB = NA, MBR.fac = 2, col.MBR = par("col"),
             lty.MBR = "dashed", lwd.MBR = 0.8,

             #new: scaling 
             x.vec = NULL, scaleX = FALSE, scaleX.fct, scaleX.inv,
             scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
             scaleN = 9, x.ticks = NULL, y.ticks = NULL,
             
             mfColRow = TRUE, to.draw.arg = NULL,
             
             alpha.trsp = NA,
             
             cex.pts = 1, 
             cex.pts.fun = NULL, 
             col.pts = par("col"),
             pch.pts = 1, 
             
             jit.fac = 1, 
             jit.tol = .Machine$double.eps, 
             
             with.lab = FALSE,
             
             lab.pts = NULL, lab.col = par("col"), lab.font = NULL, lab.adj = NULL,
             
             which.lbs = NULL, which.Order = NULL, which.nonlbs = NULL, return.Order = FALSE,
             
             draw.nonlbl = TRUE,  ## should non-labelled observations also be drawn?
             
             cex.nonlbl = 0.3,    ## character expansion(s) for non-labelled observations
             cex.nonlbl.fun = NULL, ## like cex.pts.fun for non-labelled observations
             col.nonlbl = par("col"),   
             pch.nonlbl = ".",    ## plotting symbol(s) for non-labelled observations
             
             withSubst = TRUE){

################################################################################
## 1. preparation: fingle around with arguments:
################################################################################
#  1.1 read out dots, object, L2Fam, scaling
################################################################################
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

        if(missing(cex.nonlbl)) cex.nonlbl <- 1
        cex.nonlbl <- rep(cex.nonlbl,length.out=ncomp)

        if(missing(col)) col <- 1:ncomp
           else col <- rep(col, length.out = ncomp)
        if(missing(lwd))  lwd <- rep(1,ncomp)
           else lwd <- rep(lwd, length.out = ncomp)
        if(!missing(lty)) rep(lty, length.out = ncomp)

        L2Fam <- eval(obj1@CallL2Fam)
        if(forceSameModel)
        if(!identical(CallL2Fam(obj1),CallL2Fam(obj2)))
            stop("ICs need to be defined for the same model")

        if(missing(scaleX.fct)){
           scaleX.fct <- p(L2Fam)
           scaleX.inv <- q(L2Fam)
        }
        if(missing(scaleY.fct)){
           scaleY.fct <- pnorm
           scaleY.inv <- qnorm
        }

################################################################################
#  1.2 axis type: withbox xlab, ylab
################################################################################
        dots["type"] <- NULL
        xlab <- dots$xlab; if(is.null(xlab)) xlab <- "x"
        ylab <- dots$ylab; if(is.null(ylab)) ylab <- "(partial) IC"
        dots$xlab <- dots$ylab <- NULL


################################################################################
#  1.3 parameter trafo and dimensions of the panels
################################################################################
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

################################################################################
#  1.4 preparation of cex, scaling  per panel, legend
################################################################################

        if(!is.null(x.ticks)) dotsP$xaxt <- "n"
        if(!is.null(y.ticks)){
           y.ticks <- .fillList(y.ticks, dims0)
           dotsP$yaxt <- "n"
        }

        if(!is.null(cex.pts.fun)){
           cex.pts.fun <- .fillList(cex.pts.fun, dims0*ncomp)
        }
        if(!is.null(cex.nonlbl.fun)){
           cex.nonlbl.fun <- .fillList(cex.nonlbl.fun, dims0*ncomp)
        }


        scaleY.fct <- .fillList(scaleY.fct, dims0)
        scaleY.inv <- .fillList(scaleY.inv, dims0)

        MBRB <- matrix(rep(t(MBRB), length.out=dims0*2),ncol=2, byrow=T)
        MBRB <- MBRB * MBR.fac


################################################################################
#  1.5  determine x and y grid
################################################################################
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


################################################################################
#  1.6  clean up dots arguments
################################################################################
      dots$ylim <- dots$xlim <- NULL
      lineT <- NA
      dotsT$main <- dotsT$cex.main <- dotsT$col.main <- dotsT$line <- NULL

      
################################################################################
#  1.7  prepare titles
################################################################################
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


################################################################################
#  2. computation of x,y-points to be plotted
################################################################################
#  2.1. preparation: ICs to evaluate 
################################################################################
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

################################################################################
#  2.2. preparation: what is to be done "on exit"
################################################################################
        w0 <- getOption("warn"); options(warn = -1); on.exit(options(warn = w0))

################################################################################
#  2.3. panel axis : mar, xaxt, yaxt, ...
################################################################################
        opar <- par(no.readonly = TRUE)
        if(mfColRow){ on.exit(par(opar)); par(mfrow = c(nrows, ncols)) }

        if(is(distr, "DiscreteDistribution"))
                x.vecD <- seq(from = min(x.vec), to = max(x.vec), length = 1000)


################################################################################
#  2.4. pre- and posthooks per panel (panel last -> pL, panel first -> pF)
################################################################################
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
        dotsP$panel.last <- NULL

        ..panelLast <- .fillList(pL, dims0)
        if(dims0 > 0){
           pL <- vector("list",dims0)
           for(i in 1:dims0)
               pL[[i]] <- if(is.null(..panelLast[[i]])) expression({}) else ..panelLast[[i]]
        }

################################################################################
#  2.5 plotting in data : preparation
################################################################################
        sel1 <- sel2 <- sel3 <- sel4 <- NULL
        if(!is.null(data)){
            n <- if(!is.null(dim(data))) nrow(data) else length(data)

            lab.pts <- rep(lab.pts, length.out=n)
            jit.fac <- rep(jit.fac, length.out=ncomp)
            jit.tol <- rep(jit.tol, length.out=ncomp)

            with.lab <- rep(with.lab, length.out=ncomp)
            lab.font <- rep(lab.font, length.out=ncomp)

            lab.adj <- if(is.null(lab.adj)){ matrix(0.5,n,2*ncomp)
                          }else{
                             if(length(lab.adj)%in%c(1,2))
                                lab.adj <- rep(lab.adj, length.out=2*ncomp) 
                             if(length(lab.adj)==2*ncomp){
                                lab.adj <- matrix( rep(lab.adj,
                                            times=rep(n,times=2*ncomp)), n,2*ncomp) 
                             }else{
                                 if(!is.matrix(lab.adj))
                                     lab.adj <- matrix(rep(lab.adj, length.out=n),n,1)
                                 if(ncol(lab.adj)==1) 
                                     lab.adj <- cbind(lab.adj,lab.adj)  
                                 if(ncol(lab.adj)==2) 
                                    lab.adj <- lab.adj[,rep(1:2,ncomp)] 
                                 if(ncol(lab.adj)!=2*ncomp) 
                                    stop("Wrong number of columns in arg 'lab.adj'.")
                             }             
                          }
                          
            absInfoEval <- function(x,IC,obj0=obj1){
                  QF <- ID
                  if(is(IC,"ContIC") & dims>1 ){
                     if (is(normtype(obj0),"QFNorm"))
                          QF <- QuadForm(normtype(obj0))
                  }
                  absInfo.f <- t(IC) %*% QF %*% IC
                  return(sapply(x, absInfo.f@Map[[1]]))
            }
            def.sel <- function(IC, obj00){
                 fct.aI <- function(x) absInfoEval(x,IC, obj0=obj00)
                 return(.SelectOrderData(data, fct.aI, which.lbs, 
                                         which.Order, which.nonlbs))}
                 
            sel1 <- def.sel(IC1,obj1); sel2 <- def.sel(IC2,obj2)
            selAlly <- c(sel1[["y"]],sel2[["y"]],sel1[["y.ns"]],sel2[["y.ns"]])

            if(is(obj3, "IC")){ sel3 <- def.sel(IC3,obj3)
                                selAlly <- c(selAlly,sel3[["y"]],sel3[["y.ns"]])
                              }
            if(is(obj4, "IC")){ sel4 <- def.sel(IC4,obj4)
                                selAlly <- c(selAlly,sel4[["y"]],sel4[["y.ns"]])
                              }

            n.s <- length(sel1[["ind"]]) #nrow(selAlly)
            n.ns<- length(sel1[["ind.ns"]])

            if(is.matrix(col.pts)){
               col.pts0 <- matrix(rep(col.pts,length.out=n*ncomp),n,ncomp)
               col.pts <- .cbindSel4(sel1,sel2,if(is(obj3, "IC")) sel3 else NULL, 
                                     if(is(obj4, "IC")) sel4 else NULL, 
                                     matr=col.pts0,sel=TRUE)
               col.nonlbl <- if(n.ns>0 && draw.nonlbl)
                  .cbindSel4(sel1,sel2,if(is(obj3, "IC")) sel3 else NULL, 
                                        if(is(obj4, "IC")) sel4 else NULL, 
                                        matr=col.pts0,sel=FALSE) else NULL
            }else{ 
               col.pts <- rep(col.pts, length.out=ncomp)
               col.nonlbl <- if(n.ns>0 && draw.nonlbl) 
                                rep(col.nonlbl, length.out=ncomp) else NULL
            }
            
            if(is.matrix(pch.pts)){
               pch.pts0 <- matrix(rep(pch.pts,length.out=n*ncomp),n,ncomp)
               pch.pts <- .cbindSel4(sel1,sel2,if(is(obj3, "IC")) sel3 else NULL, 
                                     if(is(obj4, "IC")) sel4 else NULL, 
                                     matr=pch.pts0,sel=TRUE)
               pch.nonlbl <- if(n.ns>0 && draw.nonlbl)
                  .cbindSel4(sel1,sel2,if(is(obj3, "IC")) sel3 else NULL, 
                                       if(is(obj4, "IC")) sel4 else NULL, 
                             matr=pch.pts0,sel=FALSE)  else NULL
            }else{ 
               pch.pts <- matrix(rep(pch.pts, length.out=ncomp*n.s),n.s,ncomp)
               pch.nonlbl <- if(n.ns>0 && draw.nonlbl)
                  matrix(rep(pch.nonlbl, length.out=ncomp*n.ns),
                         n.ns,ncomp) else NULL
            }
            
            if(is.matrix(cex.pts)){
               cex.pts0 <- matrix(rep(cex.pts,length.out=n*ncomp),n,ncomp)
               cex.pts <- .cbindSel4(sel1,sel2,if(is(obj3, "IC")) sel3 else NULL, 
                                     if(is(obj4, "IC")) sel4 else NULL, 
                                     matr=cex.pts0,sel=TRUE)
               cex.nonlbl <- if(n.ns>0 && draw.nonlbl)
                       .cbindSel4(sel1,sel2,if(is(obj3, "IC")) sel3 else NULL, 
                                            if(is(obj4, "IC")) sel4 else NULL, 
                                  matr=cex.pts0,sel=FALSE)  else NULL
            }else{ 
               cex.pts <- rep(cex.pts,length.out=ncomp)
               cex.nonlbl <- if(n.ns>0 && draw.nonlbl) 
                          rep(cex.nonlbl,length.out=ncomp) else NULL
            }
     
               ## label to be plotted at points


            if(is.matrix(lab.col)){
               lab.col0 <- matrix(rep(lab.col,length.out=n*ncomp),n,ncomp)
               lab.col  <- .cbindSel4(sel1,sel2,if(is(obj3, "IC")) sel3 else NULL, 
                                     if(is(obj4, "IC")) sel4 else NULL, 
                                     matr=lab.col0,sel=TRUE)
            }else{ 
               lab.col <- rep(lab.col, length.out=ncomp)
            }

            
            dots.points <- .makedotsLowLevel(dots)
            dots.points$col <- dots.points$cex <- dots.points$pch <- NULL

            alp.v <- array(rep(alpha.trsp, length.out = ncomp*n*dims0),
                               dim=c(n,dims0,ncomp))

################################################################################
#  2.6 inserting the code to plot in data into panel last
################################################################################
            pL <- substitute({
                 doIt <- function(sel.l,fct.l,j.l){

                     rescd <- .rescalefct(sel.l[["data"]], fct.l, scaleX, scaleX.fct,
                                   #scaleX.inv, 
                                   scaleY, scaleY.fct[[i]], xlim[,i],
                                   ylim[,i], dotsP)
                     if(is(distr, "DiscreteDistribution")){
                        rescd$Y <- jitter(rescd$Y, factor = jit.fac0[j.l])
                     }else{ if(any(.isReplicated(rescd$Y, jit.tol0[j.l]))&&jit.fac0[j.l]>0)
                             rescd$Y <- jitter(rescd$Y, factor = jit.fac0[j.l])                    
                     }   

                     sel <- rescd$idx
                     i.l <- (sel.l[["ind"]])[sel]
                     n.l <- length(i.l)
                     pch.s.l <- pch0.s[sel,j.l]
                     lab.pts.l <- if(is.null(lab0)) paste(i.l) else lab0[i.l]
                     
                     ##local alpha transp
                     al0.s.l <- al0[i.l,indi,j.l]
                     col.s.l <- .alphTrspWithNA(col0.s[sel,j.l],al0.s.l)

                     cfun.s <- if(is.null(cexfun.s)) NULL else cexfun.s[[(i-1)*ncomp+j.l]]
                     ##.cexscale in plotUtils.R
                     cex.s.l <- .cexscale(sel.l[["y"]][sel], selAlly, 
                                          cex=cex0.s[sel,j.l], fun = cfun.s[[j.l]])   

                     do.call(points, args=c(list(rescd$X, rescd$Y, cex = cex.s.l,
                             col = col.s.l, pch = pch.s.l), dwo0))
                     if(with.lab0){
                        for(kk in 1:n.l)
                            text(rescd$X[kk], rescd$Y[kk], labels = lab.pts.l[kk],
                               cex = cex.s.l/2, col = lab.col0[sel,j.l], 
                               font = lab.ft0[j.l], 
                               adj = lab.ad0[kk,(j.l-1)*2+(1:2)])
                     }

                     if(dononlb){
                         rescd.ns <- .rescalefct(sel.l[["data.ns"]], fct.l, 
                                       scaleX, scaleX.fct, #scaleX.inv, 
                                       scaleY, scaleY.fct[[i]], xlim[,i],
                                       ylim[,i], dotsP)

                         if(is(distr, "DiscreteDistribution")){
                            rescd.ns$Y <- jitter(rescd.ns$Y, 
                                                 factor = jit.fac0[j.l])
                         }else{ if(any(.isReplicated(rescd.ns$Y, 
                                        jit.tol0[j.l])) && jit.fac0[j.l]>0)
                                 rescd.ns$Y <- jitter(rescd.ns$Y, 
                                                      factor = jit.fac0[j.l])                    
                         }   
                  
                         sel.ns <- rescd.ns$idx
                         i.nl <- (sel.l[["ind.ns"]])[rescd.ns$idx] 
                         pch.ns.l <- pch0.ns[sel.ns,j.l]
                  
                  
                         ##local alpha transp
                         al0.ns.l <- al0[i.nl,indi,j.l]
                         col.ns.l <- .alphTrspWithNA(col0.ns[sel.ns,j.l],al0.ns.l)
                  
                         cfun.ns <- if(is.null(cexfun.ns)) NULL else 
                                       cexfun.ns[[(i-1)*ncomp+j.l]]
                  
                         ##.cexscale in plotUtils.R
                         cex.ns.l <- .cexscale(sel.l[["y.ns"]][sel.ns], selAlly, 
                                               cex=cex0.ns[sel.ns,j.l], 
                                               fun = cfun.ns[[j.l]])   
                  
                  
                         do.call(points, args=c(list(rescd.ns$X, rescd.ns$Y, 
                                 cex = cex.ns.l, col = col.ns.l, pch = pch.ns.l), 
                                 dwo0))
                     }        
                 }
                 doIt(sel1,fct1,1);  doIt(sel2,fct2,2)
                 if(is(obj3, "IC")) doIt(sel3,fct3,3)
                 if(is(obj4, "IC")) doIt(sel4,fct4,4)
                 pL0
              }, list(pL0 = pL, 
                      cex0.s = cex.pts, cex0.ns = cex.nonlbl,
                      pch0.s = pch.pts, pch0.ns = pch.nonlbl,
                      col0.s = col.pts, col0.ns = col.nonlbl,
                      cfun.s=cex.pts.fun, cfun.ns=cex.nonlbl.fun,
                      lab0 = lab.pts, lab.ad0 = lab.adj,
                      lab.col0 = lab.col, lab.ft0 = lab.font, 
                      jit.fac0 = jit.fac, jit.tol0=jit.tol,
                      dwo0 = dots.points, al0 = alp.v,
                      with.lab0 = with.lab, lab0 = lab.pts,
                      dononlb = draw.nonlbl&(n.ns>0) 
                      )
            )
        }

        dotsP$type <- dotsP$lty <- dotsP$col <- dotsP$lwd <- NULL
        dotsP$xlab <- dotsP$ylab <- NULL

################################################################################
#  3. creating the panel plots in a for loop
################################################################################

        cpInfo <- vector("list",0)
        cpInfo$panels <- vector("list",dims0)
        for(i in 1:dims0){
            indi <- to.draw[i]
            if(!is.null(ylim)) dotsP$ylim <- ylim[,i]

            fct1 <- function(x) sapply(x, IC1@Map[[indi]])

            resc.args <- c(list(x.vec, "fc"=fct1, scaleX, scaleX.fct, # scaleX.inv, 
                                scaleY, scaleY.fct[[i]], xlim[,i],
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

            plot.args <- c(list(x = resc1$X, y = y0,
                 type = "n", xlab = .mpresubs(xlab), ylab = .mpresubs(ylab),
                 lty = lty[1], col = addAlphTrsp2col(col[1],0),
                 lwd = lwd[1]), dotsP, list(panel.last = pL[[i]], panel.first=pF[[i]]))
            do.call(plot, args=plot.args)
            cpInfo$panels[[i]]$plot.args <- plot.args
            rm(plot.args)
            
            if(plty=="p"){
               matpoints.args <- c(list( x = resc1$X, y = matp,
                    col = col), dots.points)
               do.call(matpoints, args = matpoints.args)
               cpInfo$panels[[i]]$matpoints.args <- matpoints.args 
               rm(matpoints.args)
            }
            
            matlines.args <- c(list( x = resc1$X, y = matp,
                    lty = lty, col = col, lwd = lwd), dotsL)
            do.call(matlines, args = matlines.args)
            cpInfo$panels[[i]]$matlines.args <- matlines.args
            rm(matlines.args)

            .plotRescaledAxis.args <- list(scaleX, scaleX.fct, scaleX.inv,
                              scaleY,scaleY.fct[[i]], scaleY.inv[[i]], xlim[,i],
                              ylim[,i], resc1$X, ypts = 400, n = scaleN,
                              finiteEndpoints = finiteEndpoints,
                              x.ticks = x.ticks, y.ticks = y.ticks[[i]])
            do.call(.plotRescaledAxis, args=.plotRescaledAxis.args)
            cpInfo$panels[[i]]$.plotRescaledAxis.args <- .plotRescaledAxis.args
            rm(.plotRescaledAxis.args)                  

            if(withMBR){
                MBR.i <- MBRB[i,]
                if(scaleY) MBR.i <- scaleY.fct[[i]](MBR.i)
                MBR.args <- list(h=MBR.i, col=col.MBR, lty=lty.MBR, lwd = lwd.MBR)
                do.call(abline, args=MBR.args)
                cpInfo$panels[[i]]$MBR.args <- MBR.args
                rm(MBR.args)
            }

            if(is(distr, "DiscreteDistribution")){
                 rescD.args <- c(list(x.vecD, "fc"=fct1, scaleX, scaleX.fct,
                                #scaleX.inv, 
                                scaleY, scaleY.fct[[i]], xlim[,i],
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
                 matlinesD.args <- c(list(resc1D$X, matpD, lty = lty,
                         col = col, lwd = lwd), dotsL)
                 do.call(matlines, matlinesD.args)
                 cpInfo$panels[[i]]$matlinesD.args <- matlinesD.args 
                 rm(matlinesD.args)        
            }

           if(innerL){
              title.args <- c(list(main = innerT[[indi]]),  dotsT,
                      line = lineT, cex.main = cex.inner, col.main = col.inner)
              do.call(title, args=title.args)
              cpInfo$panels[[i]]$title.args <- title.args
              rm(title.args)        
           }           
        }

################################################################################
#  4. outer titles
################################################################################
        if(with.legend){
           if(is.null(legend)) legend <- xc
           legend.args <- c(list(.legendCoord(legend.location, scaleX, scaleX.fct,
                                scaleY, scaleY.fct[[i]]), col = col, 
                                bg = legend.bg, legend = legend), dotsLeg, 
                                cex = legend.cex)
           do.call(graphics::legend, args=legend.args)
           cpInfo$legend.args <- legend.args
           rm(legend.args)                      
        }

        cex.main <- if(!hasArg(cex.main)) par("cex.main") else dots$"cex.main"
        col.main <- if(!hasArg(col.main)) par("col.main") else dots$"col.main"
        if (mainL){
            main.args <- list(text = main, side = 3, cex = cex.main, adj = .5,
                  outer = TRUE, padj = 1.4, col = col.main)
            do.call(mtext, args = main.args)
            cpInfo$main.args <- main.args 
            rm(main.args)
        }
        cex.sub <- if(!hasArg(cex.sub)) par("cex.sub") else dots$"cex.sub"
        col.sub <- if(!hasArg(col.sub)) par("col.sub") else dots$"col.sub"
        if (subL){
            sub.args <- list(text = sub, side = 1, cex = cex.sub, adj = .5,
                  outer = TRUE, line = -1.6, col = col.sub)
            do.call(mtext, args = sub.args)
            cpInfo$sub.args <- sub.args 
            rm(sub.args)
        }

        class(cpInfo) <- c("comparePlotInfo","DiagnInfo")
        retv <- list(call=.mc, comparePlotInfo = cpInfo)       
        if(return.Order){ 
           retOrder <- list(obj1=sel1$ind1, obj2=sel2$ind1,
                                     obj3=sel3$ind1, obj4=sel4$ind1)
           retv$retOrder <- retOrder         
        }
        invisible(return(retv))
    })

