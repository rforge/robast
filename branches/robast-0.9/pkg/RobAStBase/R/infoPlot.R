setMethod("infoPlot", "IC",
    function(object, data = NULL,
             ..., withSweave = getdistrOption("withSweave"),
             col = par("col"), lwd = par("lwd"), lty, 
             colI = grey(0.5), lwdI = 0.7*par("lwd"), ltyI = "dotted",
             main = FALSE, inner = TRUE, sub = FALSE, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3], 
             with.legend = TRUE, legend = NULL, legend.bg = "white",
             legend.location = "bottomright", legend.cex = 0.8,
             scaleX = FALSE, scaleX.fct, scaleX.inv,
             scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
             mfColRow = TRUE, to.draw.arg = NULL,
             cex.pts = 1, col.pts = par("col"),
             pch.pts = 1, jitter.fac = 1, with.lab = FALSE,
             lab.pts = NULL, lab.font = NULL,
             which.lbs = NULL, which.Order  = NULL, return.Order = FALSE,
             ylab.abs = "absolute information", 
             ylab.rel= "relative information"){

        objectc <- match.call(call = sys.call(sys.parent(1)))$object
        dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
                   
        if(DEBUG.INFO) print("DATA")
        if(DEBUG.INFO) print(str(data))
        L2Fam <- eval(object@CallL2Fam)

        if(missing(scaleX.fct)){
           scaleX.fct <- p(L2Fam)
           scaleX.inv <- q(L2Fam)
        }


        dots["type"] <- NULL
        if(!is.null(dots[["xlab"]])) xlab0 <- dots[["xlab"]]
        dots["ylab"] <- NULL
        
        trafO <- trafo(L2Fam@param)
        dims <- nrow(trafO)
        dimm <- ncol(trafO)
        
        to.draw <- 1:(dims+1)
        dimnms <- rownames(trafO)
        if(is.null(dimnms))
           dimnms <- names(main(L2Fam@param))# paste("dim",1:dims,sep="")
        pdimnms <- c("Abs",dimnms)
        if(! is.null(to.draw.arg)){
            if(is.character(to.draw.arg)) 
                 to.draw <- pmatch(to.draw.arg, pdimnms)
            else if(is.numeric(to.draw.arg)) 
                 to.draw <- to.draw.arg
        }
        
        to.draw1 <- to.draw[to.draw>1]
        dims0 <- length(to.draw1)
        nrows <- trunc(sqrt(dims0))
        ncols <- ceiling(dims0/nrows)
        in1to.draw <- (1%in%to.draw)

        if(with.legend){
          if(missing(legend.location)){
             legend.location <- distr:::.fillList(list("topright"), dims0+in1to.draw   )
             if (in1to.draw) legend.location[[1]] <-  "bottomright"
          }else{
             legend.location <- as.list(legend.location)
             legend.location <- distr:::.fillList(legend.location, dims0+in1to.draw   )
          }
          if(is.null(legend)){
             legend <- vector("list",dims0+in1to.draw)
             legend <- distr:::.fillList(as.list(c("class. opt. IC", objectc)),
                                                 dims0+in1to.draw)
          }
        }
        e1 <- L2Fam@distribution
        if(!is(e1, "UnivariateDistribution") | is(e1, "CondDistribution"))
            stop("not yet implemented")

        if(is(e1, "UnivariateDistribution")){
           xlim <- eval(dots$xlim)
           if(!is.null(xlim)){ 
               xm <- min(xlim)
               xM <- max(xlim)
               dots$xlim <- NULL
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
         }
         ylim <- eval(dots$ylim)
         if(!is.null(ylim)){ 
               if(!length(ylim) %in% c(2,2*(dims0+in1to.draw))) 
                  stop("Wrong length of Argument ylim"); 
               ylim <- matrix(ylim, nrow=2,ncol=dims0+in1to.draw)
               dots$ylim <- NULL
         }

         dotsL <- .makeLowLevel(dots)
         dotsL$lwd <- dotsL$col <- dotsL$lty <- NULL
         dotsT <- dotsP <- dotsL
         dotsT["main"] <- dotsT["cex.main"] <- dotsT["col.main"] <- NULL
         dotsT["line"] <- NULL
         dotsP$xlim <- xlim
         
         trafo <- trafo(L2Fam@param)
            
            
            mainL <- FALSE
            subL <- FALSE
            lineT <- NA
       
           .mpresubs <- function(inx)
                    distr:::.presubs(inx, c("%C", "%D", "%A"),
                          c(as.character(class(object)[1]),
                            as.character(date()),
                            as.character(deparse(objectc))))
           
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
                #if(!is.character(inner))
                #stop("Argument 'inner' must either be 'logical' or a 'list'")
                if(!is.list(inner))
                    inner <- as.list(inner)                
                innerT <- distr:::.fillList(inner,1+dims)
                if(dims0<dims){
                   innerT0 <- innerT
                   for(i in 1:dims0) innerT[1+to.draw[i]] <- innerT0[1+i]          
                }
                innerL <- TRUE
            }else{if(any(is.na(inner))||any(!inner)) {
                     innerT <- as.list(rep("",1+dims)); innerL <- FALSE
                }else{innerL <- TRUE
                      tnm  <- rownames(trafO)
                      tnms <- if(is.null(tnm)) paste(1:dims) else 
                                               paste("'", tnm, "'", sep = "") 
                      innerT <- as.list(paste(c( paste(gettext("Absolute information of (partial) IC for "), 
                                       name(L2Fam)[1], sep =""),
                                   paste(gettext("Relative information of \ncomponent "),
                                       tnms, 
                                       gettext(" of (partial) IC\nfor "), 
                                       name(L2Fam)[1], sep ="")), innerParam))
                   }
              }


            QFc <- diag(dims)
            if(is(object,"ContIC") & dims>1 )
               {if (is(normtype(object),"QFNorm")) QFc <- QuadForm(normtype(object))
                QFc0 <- solve( trafo %*% solve(L2Fam@FisherInfo) %*% t(trafo ))
                if (is(normtype(object),"SelfNorm")|is(normtype(object),"InfoNorm")) 
                    QFc <- QFc0
               }

            absInfoEval <- function(x,y, withNorm = FALSE){
                       aI <- sapply(x, y@Map[[1]])
                       if(withNorm) aI <- aI / max(aI)
                       return(aI)
            }


            QFc.5 <- sqrt(PosSemDefSymmMatrix(QFc))

            classIC <- as(trafo %*% solve(L2Fam@FisherInfo) %*% L2Fam@L2deriv, "EuclRandVariable")
            absInfoClass.f <- t(classIC) %*% QFc %*% classIC
            absInfoClass <- absInfoEval(x.vec, absInfoClass.f)

            QF <- diag(dims)
            if(is(object,"ContIC") & dims>1 )
               {if (is(normtype(object),"QFNorm")) QF <- QuadForm(normtype(object))}
            QF.5 <- sqrt(PosSemDefSymmMatrix(QF))

            IC1 <- as(diag(dims) %*% object@Curve, "EuclRandVariable")
            absInfo.f <- t(IC1) %*% QF %*% IC1
            absInfo <- absInfoEval(x.vec, absInfo.f)


            w0 <- getOption("warn")
            options(warn = -1)
            on.exit(options(warn = w0))
            opar <- par(no.readonly = TRUE)
#            opar$cin <- opar$cra <- opar$csi <- opar$cxy <-  opar$din <- NULL
            if(mfColRow) on.exit(par(opar))
#            if (!withSweave)
#               devNew()

            omar <- par("mar")
            parArgs <- list(mar = c(bmar,omar[2],tmar,omar[4]))
            do.call(par,args=parArgs)

            
            pL.rel <- pL.abs <- pL <- expression({})
            if(!is.null(dots$panel.last))
               {pL.rel <- pL.abs <- pL <- dots$panel.last}

            if(!is.null(data)){

               n <- if(!is.null(dim(data))) nrow(data) else length(data)
               if(!is.null(lab.pts))
                    lab.pts <-  matrix(rep(lab.pts, length.out=2*n),n,2)

               sel <- .SelectOrderData(data, function(x)absInfoEval(x,absInfo.f),
                                       which.lbs, which.Order)
               sel.C <- .SelectOrderData(data, function(x)absInfoEval(x,absInfoClass.f),
                                       which.lbs, which.Order)
               i.d <- sel$ind
               i.dC <- sel.C$ind
               i0.d <- sel$ind1
               i0.dC <- sel.C$ind1
               y.d <- sel$y
               y.dC <- sel.C$y
               x.d <- sel$data
               x.dC <- sel.C$data
               n <- length(i.d)
               
               if(missing(col.pts)) col.pts <- c(col, colI)
               col.pts <- rep(col.pts, length.out=2)
               pch.pts <- matrix(rep(pch.pts, length.out=2*n),n,2)
               jitter.fac <- rep(jitter.fac, length.out=2)
               with.lab <- rep(with.lab, length.out=2)
               lab.font <- rep(lab.font, length.out=2)


               resc.dat <-.rescalefct(x.d, function(x) absInfoEval(x,absInfo.f),
                              scaleX, scaleX.fct, scaleX.inv,
                              scaleY, scaleY.fct, scaleY.inv,
                              dots$xlim, dots$ylim, dots)
               resc.datC <-.rescalefct(x.d, function(x) absInfoEval(x,absInfoClass.f),
                              scaleX, scaleX.fct, scaleX.inv,
                              scaleY, scaleY.fct, scaleY.inv,
                              dots$xlim, dots$ylim, dots)

               x.d <- resc.dat$X
               x.dC <- resc.datC$X
               y.d <- resc.dat$Y
               y.dC <- resc.datC$Y

               lab.pts <- if(is.null(lab.pts))
                               cbind(i.d, i.dC)
                          else cbind(lab.pts[i.d],lab.pts[i.dC])


               dots.points <-   .makedotsPt(dots)

               do.pts <- function(x,y,cxa,ca,pa)
                    do.call(points,args=c(list(x,y,cex=cxa,col=ca,pch=pa),
                            dots.points))
               tx <- function(xa,ya,lb,cx,ca)
                     text(x=xa,y=ya,labels=lb,cex=cx, col=ca)
               pL.abs <- substitute({
                   if(is(e1, "DiscreteDistribution")){
                      ICy0 <- jitter(ICy0, factor = jitter.fac0[1])
                      ICy0c <- jitter(ICy0c, factor = jitter.fac0[2])
                   }
                   f1 <- log(ICy0+1)*3*cex0[1]
                   f1c <- log(ICy0c+1)*3*cex0[2]
                   do.pts(y0, ICy0, f1,col0[1],pch0[,1])
                   do.pts(y0c, ICy0c, f1c,col0[2],pch0[,2])
                   if(with.lab0){
                      tx(y0, ICy0, lab.pts0, f1/2, col0[1])
                      tx(y0c, ICy0c, lab.pts0C, f1c/2, col0[2])
                   }
                   pL0
                   }, list(ICy0 = y.d, ICy0c = y.dC,
                           pL0 = pL, y0 = x.d, y0c = x.dC,
                           cex0 = cex.pts, pch0 = pch.pts,
                           col0 = col.pts, with.lab0 = with.lab, n0 = n,
                           lab.pts0 = lab.pts[i.d], lab.pts0C = lab.pts[i.dC],
                           jitter.fac0 = jitter.fac)
                           )

               pL.rel <- substitute({
                     y0.vec <- sapply(y0,  IC1.i.5@Map[[indi]])^2/ICy0
                     y0c.vec <- sapply(y0c, classIC.i.5@Map[[indi]])^2/ICy0c
                   if(is(e1, "DiscreteDistribution")){
                      y0.vec <- jitter(y0.vec, factor = jitter.fac0[1])
                      y0c.vec <- jitter(y0c.vec, factor = jitter.fac0[2])
                   }
                   f1 <- log(ICy0+1)*3*cex0[1]
                   f1c <- log(ICy0c+1)*3*cex0[2]
                   do.pts(y0, y0.vec, f1,col0[1],pch0[,1])
                   do.pts(y0c, y0c.vec, f1c,col0[2],pch0[,2])
                   if(with.lab0){
                      text(y0, y0.vec, lab.pts0, f1/2, col0[1])
                      text(y0c, y0c.vec, lab.pts0C, f1c/2, col0[2])
                   }
                   pL0
                   }, list(ICy0c = y.dC, ICy0 = y.d,
                           pL0 = pL, y0 = x.d, y0c = x.dC,
                           cex0 = cex.pts, pch0 = pch.pts,
                           col0 = col.pts, with.lab0 = with.lab,n0 = n,
                           lab.pts0 = lab.pts[i.d], lab.pts0C = lab.pts[i.dC],
                           jitter.fac0 = jitter.fac
                           ))
            }

            if(!is.null(ylim))
                dotsP$ylim <- ylim[,1]       
            
            fac.leg <- if(dims0>1) 3/4 else .75/.8 

            
            dotsP$axes <- NULL
            if(1 %in% to.draw){
               resc <-.rescalefct(x.vec, function(x) absInfoEval(x,absInfo.f),
                              scaleX, scaleX.fct, scaleX.inv,
                              scaleY, scaleY.fct, scaleY.inv,
                              dots$xlim, dots$ylim, dotsP)
               resc.C <-.rescalefct(x.vec, function(x) absInfoEval(x,absInfoClass.f),
                              scaleX, scaleX.fct, scaleX.inv,
                              scaleY, scaleY.fct, scaleY.inv,
                              dots$xlim, dots$ylim, dotsP)
               dotsP1 <- dotsP <- resc$dots
               dotsP$yaxt <- dots$yaxt
               x.vec0 <- resc$x
               x.vec0C <- resc.C$x
               x.vec1 <- resc$X
               x.vec1C <- resc.C$X
               y.vec1 <- resc$Y
               y.vec1C <- resc.C$Y

               do.call(plot, args=c(list(x.vec1, y.vec1, type = plty,
                   lty = ltyI, col = colI, lwd = lwdI,
                   xlab = xlab0, ylab = ylab.abs, panel.last = pL.abs),
                   dotsP1))
               do.call(lines, args=c(list(x.vec1C, y.vec1C, type = plty,
                       lty = lty, lwd = lwd, col = col), dotsL))
               .plotRescaledAxis(scaleX, scaleX.fct, scaleX.inv,
                              scaleY,scaleY.fct, scaleY.inv,
                              dots$xlim, dots$xlim, x.vec1, ypts = 400)
               if(with.legend)
                 legend(.legendCoord(legend.location[[1]], scaleX, scaleX.fct,
                        scaleY, scaleY.fct), legend = legend[[1]], bg = legend.bg,
                     lty = c(ltyI, lty), col = c(colI, col), 
                     lwd = c(lwdI, lwd), cex = legend.cex*fac.leg)


               if(innerL)
                  do.call(title, args=c(list(main = innerT[[1]]),  dotsT,
                          line = lineT, cex.main = cex.inner, col.main = col.inner))
            }
            
            if(dims > 1 && length(to.draw[to.draw!=1])>0){
                nrows <- trunc(sqrt(dims))
                ncols <- ceiling(dims/nrows)
                if (!withSweave||!mfColRow)
                     devNew()
                if(mfColRow)
                   parArgs <- c(parArgs,list(mfrow = c(nrows, ncols)))

                do.call(par,args=parArgs)

                IC1.i.5 <- QF.5%*%IC1
                classIC.i.5 <- QFc.5%*%classIC
                for(i in 1:dims0){
                    indi <- to.draw1[i]-1
                    if(!is.null(ylim)) 
                         dotsP$ylim <- ylim[,in1to.draw+i]       
                    else dotsP$ylim <- c(0,1)

                    y.vec1 <- sapply(x.vec0, IC1.i.5@Map[[indi]])^2/
                              absInfoEval(x,absInfo.f)
                    y.vec1C <- sapply(x.vec0C, classIC.i.5@Map[[indi]])^2/
                              absInfoEval(x,absInfoClass.f)

                    do.call(plot, args=c(list(x.vec1, y.vec1, type = plty,
                                  lty = lty, xlab = xlab0, ylab = ylab.rel,
                                  col = col, lwd = lwd, panel.last = pL.rel),
                                  dotsP))

                    do.call(lines, args = c(list(x.vec1C, y.vec1C, type = plty,
                            lty = ltyI, col = colI, lwd = lwdI), dotsL))
                    .plotRescaledAxis(scaleX, scaleX.fct, scaleX.inv,
                              FALSE,scaleY.fct, scaleY.inv,
                              dots$xlim, dots$xlim, x.vec1, ypts = 400)
                    if(with.legend)
                      legend(.legendCoord(legend.location[[i+in1to.draw]],
                                 scaleX, scaleX.fct, scaleY, scaleY.fct),
                           bg = legend.bg, legend = legend[[i+in1to.draw]],
                           col = c(colI, col), lwd = c(lwdI, lwd),
                           lty = c(ltyI, lty), cex = legend.cex*fac.leg)
                    if(innerL)
                       do.call(title, args = c(list(main = innerT[[1+indi]]),  
                               dotsT, line = lineT, cex.main = cex.inner, 
                               col.main = col.inner))
                }
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

        if(return.Order) return(list(IC=i0.d,IC.class=i0.dC))

        invisible()
        }
    )
 