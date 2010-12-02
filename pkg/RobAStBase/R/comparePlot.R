setMethod("comparePlot", signature("IC","IC"),
    function(obj1,obj2, obj3 = NULL, obj4 = NULL, data = NULL,
             ..., withSweave = getdistrOption("withSweave"), 
             main = FALSE, inner = TRUE, sub = FALSE, 
             col = par("col"), lwd = par("lwd"), lty, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3], 
             with.legend = TRUE, legend.bg = "white",
             legend.location = "bottomright", legend.cex = 0.8,
             mfColRow = TRUE, to.draw.arg = NULL,
             cex.pts = 1, col.pts = par("col"),
             pch.pts = 1, jitter.fac = 1, with.lab = FALSE,
             lab.pts = NULL, lab.font = NULL,
             which.lbs = NULL, which.Order  = NULL, return.Order = FALSE){

        xc1 <- as.character(deparse(match.call(call = sys.call(sys.parent(1)))$obj1))
        xc2 <- as.character(deparse(match.call(call = sys.call(sys.parent(1)))$obj2))
        xc <- c(xc1,xc2)
        if(!is.null(obj3))
            xc <- c(xc,as.character(deparse(match.call(call = sys.call(sys.parent(1)))$obj3)))
        if(!is.null(obj4))
            xc <- c(xc,as.character(deparse(match.call(call = sys.call(sys.parent(1)))$obj4)))
        
        dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."

        ncomp <- 2+ (!missing(obj3)|!is.null(obj3)) +  
                    (!missing(obj4)|!is.null(obj4))
         
        if(missing(col)) col <- 1:ncomp
           else col <- rep(col, length.out = ncomp)
        if(missing(lwd))  lwd <- rep(1,ncomp)
           else lwd <- rep(lwd, length.out = ncomp)
        if(!missing(lty)) rep(lty, length.out = ncomp)
        if(missing(col.pts)) col.pts <- 1:ncomp

        
        if(!is.null(dots[["type"]])) dots["type"] <- NULL
        if(!is.null(dots[["xlab"]])) dots["xlab"] <- NULL
        if(!is.null(dots[["ylab"]])) dots["ylab"] <- NULL
        
        dotsP <- dotsL <- dotsT <- dots

        L2Fam <- eval(obj1@CallL2Fam)
        L2Fam1c <- obj1@CallL2Fam
        L2Fam2c <- obj2@CallL2Fam
        if(!identical(L2Fam1c,L2Fam2c))
            stop("ICs need to be defined for the same model")

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

        e1 <- L2Fam@distribution
        if(!is(e1, "UnivariateDistribution")) stop("not yet implemented")

        xlim <- eval(dots$xlim)
        if(!is.null(xlim)){ 
               xm <- min(xlim)
               xM <- max(xlim)
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
            x.vec <- seq(from = lower - 0.1*h, to = upper + 0.1*h, length = 1000)
            plty <- "l"
            if(missing(lty)) lty <- "solid"
        }else{
            if(is(e1, "DiscreteDistribution")) x.vec <- support(e1)
            else{
                x.vec <- r(e1)(1000)
                x.vec <- sort(unique(x.vec))
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
        dots$ylim <- NULL
        dotsP$xlim <- xlim
        dots$xlim <- NULL

        dims <- nrow(trafo(L2Fam@param))
        IC1 <- as(diag(dims) %*% obj1@Curve, "EuclRandVariable")
        IC2 <- as(diag(dims) %*% obj2@Curve, "EuclRandVariable")


        obj <- obj3
        if(is(obj, "IC"))
           {
           if(!identical(L2Fam1c,obj@CallL2Fam))
               stop("ICs need to be defined for the same model")
           IC3 <- as(diag(dims) %*% obj3@Curve, "EuclRandVariable")
           }

        obj <- obj4
        if(is(obj, "IC"))
           {
           if(!identical(L2Fam1c,obj@CallL2Fam))
               stop("ICs need to be defined for the same model")
           IC4 <- as(diag(dims) %*% obj4@Curve, "EuclRandVariable")
           }

      lineT <- NA

      .mpresubs <- function(inx)
                    distr:::.presubs(inx, c(paste("%C",1:ncomp,sep=""),
                                             "%D", 
                                            paste("%A",1:ncomp,sep="")),
                          c(as.character(class(obj1)[1]),
                            as.character(class(obj2)[1]),
                            if(is.null(obj3))NULL else as.character(class(obj3)[1]),
                            if(is.null(obj4))NULL else as.character(class(obj4)[1]),
                            as.character(date()),
                            xc))
            
        mainL <- FALSE
        if (hasArg(main)){
                 mainL <- TRUE
                 if (is.logical(main)){
                     if (!main) mainL <-  FALSE
                     else
                          main <- paste(gettextf("Plot for ICs"), 
                                        paste("%A", 1:ncomp, sep="", collapse=", "),
                                        sep=" ") ###
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
        subL <- FALSE
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
        mnm <- names(L2Fam@param@main)
        mnms <- if(is.null(mnm)) NULL else paste("'", mnm, "' = ", sep = "") 
        innerParam <-  paste(gettext("\nwith main parameter ("), 
                                    paste(mnms, round(L2Fam@param@main, 3), 
                                          collapse = ", "),
                                 ")", sep = "")
        if(!is.null(L2Fam@param@nuisance)){            
            nnm <- names(L2Fam@param@nuisance)
            nnms <- if(is.null(nnm)) NULL else paste("'", nnm, "' = ", sep = "") 
            innerParam <- paste(innerParam,
                                gettext("\nand nuisance parameter ("), 
                                    paste(nnms, round(L2Fam@param@nuisance, 3), 
                                           collapse = ", "),
                                ")", sep ="")
        }
        if(!is.null(L2Fam@param@fixed)){
            fnm <- names(L2Fam@param@fixed)
            fnms <- if(is.null(fnm)) NULL else paste("'", fnm, "' = ", sep = "") 
            innerParam <- paste(innerParam,
                                gettext("\nand fixed known parameter ("), 
                                    paste(fnms, round(L2Fam@param@fixed, 3), 
                                           collapse = ", "),
                                ")", sep ="")
        }    
            if(!is.logical(inner)){
#                if(!is.character(inner))
#                    stop("Argument 'inner' must either be 'logical' or a character vector")
                if(!is.list(inner))
                    inner <- as.list(inner)                
                innerT <- distr:::.fillList(inner,dims)
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


        w0 <- getOption("warn")
        options(warn = -1)
        on.exit(options(warn = w0))
        opar <- par(no.readonly = TRUE)
#        opar$cin <- opar$cra <- opar$csi <- opar$cxy <-  opar$din <- NULL
        if(mfColRow) on.exit(par(opar))
        
        if(mfColRow)
             par(mfrow = c(nrows, ncols))

        if(is(e1, "DiscreteDistribution"))
                x.vec1 <- seq(from = min(x.vec), to = max(x.vec), length = 1000)
        
            dotsT["main"] <- NULL
            dotsT["cex.main"] <- NULL
            dotsT["col.main"] <- NULL
            dotsT["line"] <- NULL

        pL <- expression({})
        if(!is.null(dotsP$panel.last))
            pL <- dotsP$panel.last
        dotsP$panel.last <- NULL

        if(!is.null(data)){
               n <- if(!is.null(dim(data))) nrow(data) else length(data)
               oN01 <- oN02 <- oN03 <- oN04 <- NULL
               if(is.null(which.lbs))
                  which.lbs <- 1:n
               which.lbs0 <- (1:n) %in% which.lbs
               which.lbx <- rep(which.lbs0, length.out=length(data))
               data0 <- data[which.lbx]
               n <- if(!is.null(dim(data0))) nrow(data0) else length(data0)
               oN <- (1:n)[which.lbs0]

               cex.pts <- rep(cex.pts, length.out=ncomp)
               col.pts <- rep(col.pts, length.out=ncomp)
               pch.pts <- matrix(rep(pch.pts, length.out=ncomp*n),n,ncomp)
               jitter.fac <- rep(jitter.fac, length.out=ncomp)
               with.lab <- rep(with.lab, length.out=ncomp)
               lab.font <- rep(lab.font, length.out=ncomp)


               absInfoEval <- function(x,object){
                        QF <- diag(dims)
                        if(is(object,"ContIC") & dims>1 )
                           {if (is(normtype(object),"QFNorm"))
                                QF <- QuadForm(normtype(object))}

                        IC1 <- as(diag(dims) %*% object@Curve, "EuclRandVariable")
                        absInfo.f <- t(IC1) %*% QF %*% IC1
                        return(sapply(x, absInfo.f@Map[[1]]))}

               lD <- length(data)+1
               aI1 <- absInfoEval(x=data,object=obj1)
               aI2 <- absInfoEval(x=data,object=obj2)
               aI3 <- if(is.null(obj3)) NULL else absInfoEval(x=data,object=obj3)
               aI4 <- if(is.null(obj4)) NULL else absInfoEval(x=data,object=obj4)


               aI10 <- aI1[which.lbs]
               aI20 <- aI2[which.lbs]
               aI30 <- if(is.null(obj3)) NULL else aI3[which.lbs]
               aI40 <- if(is.null(obj4)) NULL else aI4[which.lbs]

               data01 <- data02 <- data03 <- data04 <- data0

               if (n==length(data0)) {
                   oN1 <-  order(aI10)
                   oN2 <-  order(aI20)
                   oN3 <-  if(is.null(obj3)) NULL else order(aI30)
                   oN4 <-  if(is.null(obj4)) NULL else order(aI40)

                   oN01 <- order(aI1)
                   oN01 <- oN01[oN01 %in% which.lbs]
                   data01 <- data0[oN1]
                   aI10 <- aI1[oN1]

                   oN02 <- order(aI2)
                   oN02 <- oN02[oN02 %in% which.lbs]
                   data02 <- data0[oN2]
                   aI20 <- aI2[oN2]

                   if(!is.null(obj3)){
                     oN03 <- order(aI3)
                     oN03 <- oN03[oN03 %in% which.lbs]
                     data03 <- data0[oN3]
                     aI30 <- aI3[oN3]
                   }
                   if(!is.null(obj4)){
                     oN04 <- order(aI4)
                     oN04 <- oN04[oN04 %in% which.lbs]
                     data04 <- data0[oN4]
                     aI40 <- aI4[oN4]
                   }
                   if(!is.null(which.Order)){
                       lD <- length(oN01)
                       oN1 <- oN01[lD-which.Order]
                       data01 <- data[oN1]
                       aI10 <- aI1[oN1]

                       oN2 <- oN02[lD-which.Order]
                       data02 <- data[oN2]
                       aI20 <- aI2[oN2]

                       if(!is.null(obj3)){
                           oN3 <- oN03[lD-which.Order]
                           data03 <- data[oN3]
                           aI30 <- aI3[oN3]
                       }

                       if(!is.null(obj4)){
                           oN4 <- oN04[lD-which.Order]
                           data04 <- data[oN4]
                           aI40 <- aI4[oN4]
                       }
                       n <- length(oN1)
                   }
               }
               if(is.null(lab.pts)){
                    lab.pts <- matrix(paste(c(oN1,oN2)),n,2)
                    if(is.null(obj3)) lab.pts <- cbind(lab.pts, paste(oN3))
                    if(is.null(obj4)) lab.pts <- cbind(lab.pts, paste(oN4))
               }else
                  lab.pts <- matrix(rep(lab.pts, length.out = ncomp*n), n, ncomp)

               dots.points <- dots
               dots.points$col <- dots.points$cex <- dots.points$pch <- NULL


               pL <- substitute({
                   ICy1 <- sapply(y01,IC1@Map[[indi]])
                   ICy2 <- sapply(y02,IC2@Map[[indi]])
                   if(!is.null(obj30))
                      ICy3 <- sapply(y03,IC3@Map[[indi]])
                   if(!is.null(obj40))
                      ICy3 <- sapply(y04,IC4@Map[[indi]])

                   if(is(e1, "DiscreteDistribution")){
                      ICy1 <- jitter(ICy1, factor = jitter.fac0[1])
                      ICy2 <- jitter(ICy2, factor = jitter.fac0[2])
                      if(!is.null(obj30))
                          ICy3 <- jitter(ICy3, factor = jitter.fac0[3])
                      if(!is.null(obj40))
                          ICy4 <- jitter(ICy4, factor = jitter.fac0[4])
                   }
                   do.call(points, args=c(list(y01, ICy1, cex = log(aI10+1)*3*cex0[1],
                                   col = col0[1], pch = pch0[,1]), dwo0))
                   do.call(points, args=c(list(y02, ICy2, cex = log(aI20+1)*3*cex0[2],
                                   col = col0[2], pch = pch0[,2]), dwo0))
                   if(!is.null(obj30))
                       do.call(points, args=c(list(y03, ICy3, cex = log(aI30+1)*3*cex0[3],
                                   col = col0[3], pch = pch0[,3]), dwo0))
                   if(!is.null(obj40))
                       do.call(points, args=c(list(y04, ICy4, cex = log(aI40+1)*3*cex0[4],
                                   col = col0[4], pch = pch0[,4]), dwo0))
                   if(with.lab0){
                      text(x = y01, y = ICy1, labels = lab.pts0[,1],
                           cex = log(aI10s+1)*1.5*cex0[1], col = col0[1])
                      text(x = y02, y = ICy2, labels = lab.pts0[,2],
                           cex = log(aI20s+1)*1.5*cex0[2], col = col0[2])
                      if(!is.null(obj30))
                          text(x = y03, y = ICy3, labels = lab.pts0[,3],
                               cex = log(aI30s+1)*1.5*cex0[3], col = col0[3])
                      if(!is.null(obj40))
                          text(x = y04, y = ICy4, labels = lab.pts0[,4],
                               cex = log(aI40s+1)*1.5*cex0[4], col = col0[4])
                   }
                   pL0
                   }, list(pL0 = pL,
                           y01 = data01, y02 = data02, y03 = data03, y04 = data04,
                           aI10s = aI10, aI20s = aI20, aI30s = aI30, aI40s = aI40,
                           obj10 = obj1, obj20 = obj2, obj30 = obj3, obj40 = obj4,
                           dwo0 = dots.points, cex0 = cex.pts, pch0 = pch.pts,
                           col0 = col.pts, with.lab0 = with.lab,
                           lab.pts0 = lab.pts, n0 = n,
                           jitter.fac0 = jitter.fac
                           ))

            }


        for(i in 1:dims0){
            indi <- to.draw[i]
            if(!is.null(ylim)) dotsP$ylim <- ylim[,i]       
            matp  <- cbind(sapply(x.vec, IC1@Map[[indi]]),
                           sapply(x.vec, IC2@Map[[indi]]))

            if(is(obj3, "IC"))
                matp  <- cbind(matp,sapply(x.vec, IC3@Map[[indi]]))
            if(is(obj4, "IC"))
                matp  <- cbind(matp,sapply(x.vec, IC4@Map[[indi]]))

            do.call(plot, args=c(list( x = x.vec, y = matp[,1],
                 type = plty, lty = lty, col = col[1], lwd = lwd,
                 xlab = "x", ylab = "(partial) IC"), dotsP, list(panel.last = pL)))
            do.call(matlines, args = c(list( x = x.vec, y = matp[,-1],
                    lty = lty, col = col[-1], lwd = lwd), dotsL))

            if(is(e1, "DiscreteDistribution")){
                 matp1 <- cbind(sapply(x.vec1, IC1@Map[[indi]]),
                                sapply(x.vec1, IC2@Map[[indi]]))
                 if(is(obj3, "IC"))
                    matp1  <- cbind(matp1,sapply(x.vec1, IC3@Map[[indi]]))
                 if(is(obj4, "IC"))
                    matp1  <- cbind(matp1,sapply(x.vec1, IC4@Map[[indi]]))
                 do.call(matlines, c(list(x.vec1, matp1, lty = lty, 
                         col = col, lwd = lwd), dotsL))
                 }

           if(innerL)
              do.call(title, args=c(list(main = innerT[[indi]]),  dotsT,
                      line = lineT, cex.main = cex.inner, col.main = col.inner))
        }
        
        if(with.legend)
           legend(legend.location, legend = xc, col = col, bg = legend.bg,
               lwd = lwd*1.5, lty = lty, cex = legend.cex)

        if(!hasArg(cex.main)) cex.main <- par("cex.main") else cex.main <- dots$"cex.main"
        if(!hasArg(col.main)) col.main <- par("col.main") else col.main <- dots$"col.main"
        if (mainL)
            mtext(text = main, side = 3, cex = cex.main, adj = .5,
                  outer = TRUE, padj = 1.4, col = col.main)

        if(!hasArg(cex.sub)) cex.sub <- par("cex.sub") else cex.sub <- dots$"cex.sub"
        if(!hasArg(col.sub)) col.sub <- par("col.sub") else col.sub <- dots$"col.sub"
        if (subL)
            mtext(text = sub, side = 1, cex = cex.sub, adj = .5,
                  outer = TRUE, line = -1.6, col = col.sub)

        if(return.Order) return(list(obj1=oN01,obj2=oN02,obj3=oN03,obj4=oN04))
        invisible()
    })
