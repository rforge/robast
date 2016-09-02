setMethod("infoPlot", "IC",
    function(object, data = NULL,
             ..., withSweave = getdistrOption("withSweave"),
             col = par("col"), lwd = par("lwd"), lty, 
             colI = grey(0.5), lwdI = 0.7*par("lwd"), ltyI = "dotted",
             main = FALSE, inner = TRUE, sub = FALSE, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3], 
             
             with.automatic.grid = TRUE, ##new
             
             with.legend = TRUE, 
             legend = NULL,       ##new
             legend.bg = "white",
             legend.location = "bottomright", 
             legend.cex = 0.8,
             
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
             
             which.lbs = NULL, which.Order  = NULL, which.nonlbs = NULL, return.Order = FALSE,
             
             draw.nonlbl = TRUE,  ## should non-labelled observations also be drawn?
             
             cex.nonlbl = 0.3,    ## character expansion(s) for non-labelled observations
             cex.nonlbl.fun = NULL, ## like cex.pts.fun for non-labelled observations
             col.nonlbl = par("col"),   
             pch.nonlbl = ".",    ## plotting symbol(s) for non-labelled observations
             
             ylab.abs = "absolute information", 
             ylab.rel= "relative information",
             
             withSubst = TRUE){

################################################################################
## 1. preparation: fingle around with arguments:
################################################################################
#  1.1 read out dots, object, L2Fam, scaling
################################################################################
        
        mc <- match.call(call = sys.call(sys.parent(1)))
        ## for later purposes in axis annotation save name of argument object 
        objectc <- mc$object
        
        ## for later purposes store the unevaluated objects in ...
        dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
        ## create respective L2 Family           
        L2Fam <- eval(object@CallL2Fam)

        ## for rescaling use PIT, i.e. the cdf/quantile in the L2Fam
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
        ## read out, if present, argument withbox and, if present, 
        #    delete it from dots
        withbox <- TRUE
        if(!is.null(dots[["withbox"]])) withbox <- dots[["withbox"]]
        dots["withbox"] <- NULL
        
        # delete argument "type"
        dots["type"] <- NULL
        ## read out, if present xlab, delete xlab, ylab from dots 
        xlab <- dots$xlab; if(is.null(xlab)) xlab <- "x"
        dots$xlab <- dots$ylab <- NULL

################################################################################
#  1.3 parameter trafo and dimensions of the panels
################################################################################
        ## read out trafo of the parameter and dimensions
        trafO <- trafo(L2Fam@param)
        dimsA <- dims <- nrow(trafO)
        dimm <- ncol(trafO)

        ## read out dimension names
        dimnms <- rownames(trafO)
        if(is.null(dimnms))
           dimnms <- names(main(L2Fam@param))
        pdimnms <- c("Abs",dimnms)
        
        ## find out how many panels are to be plotted
        to.draw <- 1:(dims+1)
        if(! is.null(to.draw.arg)){
            if(is.character(to.draw.arg)) 
                 to.draw <- pmatch(to.draw.arg, pdimnms)
            else if(is.numeric(to.draw.arg)) 
                 to.draw <- to.draw.arg
        }
        
        ## write to to.draw1 the coordinates to be plotted
        to.draw1 <- to.draw[to.draw>1]
        dims0 <- length(to.draw1)
        ## find out a default partition into panels
        nrows <- trunc(sqrt(dims0))
        ncols <- ceiling(dims0/nrows)

        ## is absInfo to be plotted?
        in1to.draw <- (1%in%to.draw)

################################################################################
#  1.4 preparation of cex, scaling  per panel, legend
################################################################################
        ## recycle, if necessary, the functions to rescale the point sizes
        ##  (in x and y coordinates, separately)
        if(!is.null(cex.pts.fun)){
           cex.pts.fun <- .fillList(cex.pts.fun, (dims0+in1to.draw)*2)
        }
        if(!is.null(cex.nonlbl.fun)){
           cex.nonlbl.fun <- .fillList(cex.nonlbl.fun, (dims0+in1to.draw)*2)
        }

        ## recycle, if necessary, the rescaling functions for y axis
        scaleY.fct <- .fillList(scaleY.fct, length(to.draw))
        scaleY.inv <- .fillList(scaleY.inv, length(to.draw))

        ## are x.ticks given? if not produce them automatically...
        if(!is.null(x.ticks)) dots$xaxt <- "n"

        ## are y.ticks given? if not produce them automatically...
        if(!is.null(y.ticks)){
           y.ticks <- .fillList(y.ticks, dims0+in1to.draw)
           dots$yaxt <- "n"
        }

        ## should a legend be drawn?
        if(with.legend){
          ## if so find a location where to place it
          if(missing(legend.location)){
             legend.location <- .fillList("topright", dims0+in1to.draw   )
             if (in1to.draw) legend.location[[1]] <-  "bottomright"
          }else{
             legend.location <- as.list(legend.location)
             legend.location <- .fillList(legend.location, dims0+in1to.draw  )
          }
          ## if so find what to writed into legend
          if(is.null(legend)){
             legend <- vector("list",dims0+in1to.draw)
             legend <- .fillList(list(as.list(c(gettext("class. opt. IC"), objectc))),
                                                 dims0+in1to.draw)
          }
        }
        
################################################################################
#  1.5  determine x and y grid
################################################################################
        ## find out the distribution of the observations
        distr <- L2Fam@distribution
        if(!is(distr, "UnivariateDistribution") | is(distr, "CondDistribution"))
            stop("not yet implemented")

        if(is(distr, "UnivariateDistribution")){

           xlim <- eval(dots$xlim)
           if(!is.null(xlim)){ 
               xm <- min(xlim)
               xM <- max(xlim)
               dots$xlim <- NULL
            }
            

            ## get x-coordinates to be plotted                
            if(is(distr, "AbscontDistribution")){

            ## case of continuous distribution

                if(is.null(x.vec)){
                   x.vec <- .getDefaultXgrid(distr = distr, xlim = xlim, 
                                scaleX=scaleX, scaleX.fct = scaleX.fct,
                                scaleX.inv=scaleX.inv, nXpts = 1000) 
                }
                plty <- "l"
                if(missing(lty)) lty <- "solid"
            }else{

            ## case of discrete distribution

                if(!is.null(x.vec)){
                   if(is(distr, "DiscreteDistribution"))
                      x.vec <- intersect(x.vec,support(distr))
                }else{
                   if(is(distr, "DiscreteDistribution")) x.vec <- support(distr)
                   else{
                      x.vec <- r(distr)(1000)
                      x.vec <- sort(unique(x.vec))
                   }
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

################################################################################
#  1.6  clean up dots arguments
################################################################################
         ## dotsP: elements from dots for points
         ## dotsL: elements from dots for lines
         ## dotsT: elements from dots for text

         dotsP <- dots
         dotsP$type <- dotsP$lty <- dotsP$col <- dotsP$lwd <- NULL
         dotsP$xlab <- dotsP$ylab <- NULL

         dotsL <- .makedotsLowLevel(dotsP)
         dotsT <- dotsL
         dotsT["main"] <- dotsT["cex.main"] <- dotsT["col.main"] <- NULL
         dotsT["line"] <- NULL
         dotsP$xlim <- xlim
         
         ## get trafo
         trafo <- trafO
            
         
         #-----------------------

################################################################################
#  1.7  prepare titles
################################################################################
         ## titles:
            mainL <- FALSE
            subL <- FALSE
            lineT <- NA
       
           .mpresubs <- if(withSubst){function(inx)
                    .presubs(inx, c("%C", "%D", "%A"),
                          c(as.character(class(object)[1]),
                            as.character(date()),
                            as.character(deparse(objectc))))
                    } else function(inx)inx
           
            if (hasArg("main")){
                 mainL <- TRUE
                 if (is.logical(main)){
                     if (!main) mainL <-  FALSE
                     else
                          main <- gettextf("Information Plot for IC %%A") ###
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
             if (hasArg("sub")){
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
                innerT <- .fillList(inner,1+dims)
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
         #-----------------------

################################################################################
#  2. computation of x,y-points to be plotted
################################################################################
#  2.1. preparation: norm, function to evaluate it for both robust and classic
################################################################################

            QFc <- diag(dimsA)
            if(is(object,"ContIC") & dimsA>1 )
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

            QF <- diag(dimsA)
            if(is(object,"ContIC") & dimsA>1 )
               {if (is(normtype(object),"QFNorm")) QF <- QuadForm(normtype(object))}
            QF.5 <- sqrt(PosSemDefSymmMatrix(QF))

            IC1 <- as(diag(dimsA) %*% object@Curve, "EuclRandVariable")
            absInfo.f <- t(IC1) %*% QF %*% IC1
            absInfo <- absInfoEval(x.vec, absInfo.f)


################################################################################
#  2.2. preparation: what is to be done "on exit"
################################################################################
            w0 <- getOption("warn")
            options(warn = -1)
            on.exit(options(warn = w0))
            opar <- par(no.readonly = TRUE)
#            opar$cin <- opar$cra <- opar$csi <- opar$cxy <-  opar$din <- NULL
            if(mfColRow) on.exit(par(opar))
#            if (!withSweave)
#               devNew()

################################################################################
#  2.3. panel axis : mar, xaxt, yaxt, ...
################################################################################
            omar <- par("mar")
            lpA <- max(length(to.draw),1)
            parArgsL <- vector("list",lpA)
            bmar <- rep(bmar, length.out=lpA)
            tmar <- rep(tmar, length.out=lpA)
            xaxt0 <- if(is.null(dots$xaxt)) {
                      if(is.null(dots$axes)||eval(dots$axes))
                         rep(par("xaxt"),lpA) else rep("n",lpA)
                      }else rep(eval(dots$xaxt),lpA)
            yaxt0 <- if(is.null(dots$yaxt)) {
                      if(is.null(dots$axes)||eval(dots$axes))
                         rep(par("yaxt"),lpA) else rep("n",lpA)
                      }else rep(eval(dots$yaxt),lpA)

            for( i in 1:lpA){
                 parArgsL[[i]] <- list(mar = c(bmar[i],omar[2],tmar[i],omar[4])
                                      ,xaxt=xaxt0[i], yaxt= yaxt0[i]
                                      )
            }

            
################################################################################
#  2.4. pre- and posthooks per panel (panel last -> pL, panel first -> pF)
################################################################################
            pL <- expression({})
            if(!is.null(dots[["panel.last"]])){
                pL <- .panel.mingle(dots,"panel.last")
            }
            dotsP$panel.last <- NULL

            ..panelLast <- .fillList(pL, length(to.draw))
            pL.abs <- if(is.null(..panelLast[[1]])) expression({}) else pL[[1]]}
            pL.rel <- pL
            
            if(in1to.draw)
               pL.rel <- if(length(pL)>1) pL[-1] else NULL
            if(length(pL.rel)>0){
               for(i in 1:length(pL.rel))
                   pL.rel[[i]] <- if(is.null(..panelLast[[i]])) expression({}) else ..panelLast[[i]]
            }
            

            pF <- expression({})
            if(!is.null(dots[["panel.first"]])){
               pF <- .panel.mingle(dots,"panel.first")
            }
            dotsP$panel.first <- NULL
            ..panelFirst <- .fillList(pF, length(to.draw))
            if(with.automatic.grid)
                ..panelFirst <- .producePanelFirstS(
                    ..panelFirst,object, to.draw.arg, TRUE,
                    x.ticks = x.ticks, scaleX = scaleX, scaleX.fct = scaleX.fct,
                    y.ticks = y.ticks, scaleY = scaleY, scaleY.fct = scaleY.fct)
            gridS <- if(with.automatic.grid)
                  substitute({grid <- function(...){}}) else expression({})
            pF.rel <- NULL
            if(in1to.draw){
               if(length(to.draw)>1){
                  for(j in 1:(length(to.draw)-1))
                      pF.rel[[j]] <- substitute({ gridS0
                                      .absInd <- FALSE
                                       pF0 }, list(pF0=..panelFirst[[j+1]], gridS0=gridS))
                  } else pF.rel <- NULL        
               pF.abs <- substitute({ gridS0
                                      .absInd <- TRUE
                                      pF0
                                      }, list(pF0=..panelFirst[[1]], gridS0=gridS))
            }else{
               pF.abs <- NULL
               if(length(to.draw>0)){
                  for(j in 1:(length(to.draw)))
                      pF.rel[[j]] <- substitute({ gridS0
                                      .absInd <- FALSE
                                      pF0 
                                      }, list(pF0=..panelFirst[[j]], gridS0=gridS))
                  } else pF.rel <- NULL        
            }
            dotsP$panel.last <- dotsP$panel.first <- NULL
            
################################################################################
#  2.5 plotting in data : preparation
################################################################################
            if(!is.null(data)){

               n <- if(!is.null(dim(data))) nrow(data) else length(data)
               if(!is.null(lab.pts))
                    lab.pts <-  matrix(rep(lab.pts, length.out=2*n),n,2)

               sel <- .SelectOrderData(data, function(x)absInfoEval(x,absInfo.f),
                                       which.lbs, which.Order, which.nonlbs)
               sel.C <- .SelectOrderData(data, function(x)absInfoEval(x,absInfoClass.f),
                                       which.lbs, which.Order, which.nonlbs)
               i.d <- sel[["ind"]]
               i.dC <- sel.C[["ind"]]
               i0.d <- sel[["ind1"]]
               i0.dC <- sel.C[["ind1"]]
               y.d <- sel[["y"]]
               y.dC <- sel.C[["y"]]
               x.d <- sel[["data"]]
               x.dC <- sel.C[["data"]]
               n0 <- length(data)
               n.s <- length(i.d)
               n.ns <- 0
               
               if(draw.nonlbl&&n.ns>0){
                  i.d.ns <- sel[["ind.ns"]]
                  i.dC.ns <- sel.C[["ind.ns"]]
                  x.d.ns <- sel[["data.ns"]]
                  x.dC.ns <- sel.C[["data.ns"]]
                  y.d.ns <- sel[["y.ns"]]
                  y.dC.ns <- sel.C[["y.ns"]]
                  n.ns <- length(i.d.ns)
               }else{
                     i.d.ns <- i.dC.ns <- x.d.ns <- NULL
                     x.dC.ns <- y.d.ns <- y.dC.ns <- NULL               
               }
               
               if(missing(col.pts)) col.pts <- c(col, colI)
               if(is.matrix(col.pts)){
                  col.pts0 <- matrix(rep(col.pts,length.out=n0*2),n0,2)
                  col.pts <- cbind(col.pts[i.d,1],col.pts[i.dC,2])
                  col.nonlbl <- if(draw.nonlbl&&n.ns>0)
                     cbind(col.pts[i.d.ns,1],col.pts[i.dC.ns,2]) else NULL
               }else{ 
                  col.pts <- t(matrix(rep(col.pts, length.out=2),2,n.s))
                  col.nonlbl <- if(draw.nonlbl&&n.ns>0)
                     t(matrix(rep(col.nonlbl, length.out=2),2,n.ns)) else NULL
               }
               
               if(is.matrix(pch.pts)){
                  pch.pts0 <- matrix(rep(pch.pts,length.out=n0*2),n0,2)
                  pch.pts <- cbind(pch.pts[i.d,1],pch.pts[i.dC,2])
                  pch.nonlbl <- if(draw.nonlbl&&n.ns>0)
                     cbind(pch.pts[i.d.ns,1],pch.pts[i.dC.ns,2]) else NULL
               }else{ 
                  pch.pts <- matrix(rep(pch.pts, length.out=2*n),n.s,2)
                  pch.nonlbl <- if(draw.nonlbl&&n.ns>0)
                     matrix(rep(pch.nonlbl, length.out=2*n.ns),n.ns,2) else NULL
               }
               
               if(is.matrix(cex.pts)){
                  cex.pts0 <- matrix(rep(cex.pts,length.out=n0*2),n0,2)
                  cex.pts <- cbind(cex.pts[i.d,1],cex.pts[i.dC,2])
                  cex.nonlbl <- if(draw.nonlbl&&n.ns>0)
                     cbind(cex.pts[i.d.ns,1],cex.pts[i.dC.ns,2]) else NULL
               }else{ 
                  cex.pts <- t(matrix(rep(cex.pts,length.out=2),2,n.s))
                  cex.nonlbl <- if(draw.nonlbl&&n.ns>0)
                     t(matrix(rep(cex.nonlbl,length.out=2),2,n.ns)) else NULL
               }

               ## label to be plotted at points

               jit.fac <- rep(jit.fac, length.out=2)
               jit.tol <- rep(jit.tol, length.out=2)

               with.lab <- rep(with.lab, length.out=2)
               lab.font <- rep(lab.font, length.out=2)

               if(is.matrix(lab.col)){
                  lab.col0 <- matrix(rep(lab.col,length.out=n0*2),n0,2)
                  lab.col  <- cbind(lab.col[i.d,1],lab.col[i.dC,2])
               }else{ 
                  lab.col <- t(matrix(rep(lab.col, length.out=2),2,n.s))
               }

               lab.pts <- if(is.null(lab.pts))
                               cbind(i.d,i.dC)
                          else cbind(lab.pts[i.d],lab.pts[i.dC])

               lab.adj <- if(is.null(lab.adj)){ matrix(0.5,n0,4)
                          }else{ 
                             if(length(lab.adj)==3)
                                stop("Wrong length of arg 'lab.adj'.")
                             if(length(lab.adj)<5){
                                matrix( rep( rep(lab.adj, length.out=4),
                                              times=rep(n0,times=4)), n0,4) 
                             }else{
                                 if(!is.matrix(lab.adj))
                                     lab.adj <- matrix(rep(lab.adj, length.out=n0),n0,1)
                                 if(ncol(lab.adj)==1)
                                     lab.adj <- cbind(lab.adj,lab.adj)  
                                 if(ncol(lab.adj)==2)
                                     lab.adj <- cbind(lab.adj,lab.adj)
                                 if(ncol(lab.adj)!=4) 
                                    stop("Wrong number of columns in arg 'lab.adj'.")
                             }             
                          }

               resc.dat <-.rescalefct(x.d, function(x) absInfoEval(x,absInfo.f),
                              scaleX, scaleX.fct, #scaleX.inv,
                              scaleY, scaleY.fct[[1]], dots$xlim, dots$ylim, dots)
               resc.datC <-.rescalefct(x.d, function(x) absInfoEval(x,absInfoClass.f),
                              scaleX, scaleX.fct, #scaleX.inv,
                              scaleY, scaleY.fct[[1]], dots$xlim, dots$ylim, dots)

               x.dr <- resc.dat$X
               x.dCr <- resc.datC$X
               y.dr <- resc.dat$Y
               y.dCr <- resc.datC$Y
               n.s <- length(resc.dat$X)
              
               pch.ptsC <- pch.pts[resc.datC$idx,2]
               pch.pts <- pch.pts[resc.dat$idx,1]
               col.ptsC <- col.pts[resc.datC$idx,2]
               col.pts <- col.pts[resc.dat$idx,1]
               cex.ptsC <- cex.pts[resc.datC$idx,2]
               cex.pts <- cex.pts[resc.dat$idx,1]
               lab.colC <- lab.col[resc.datC$idx,2]
               lab.col <- lab.col[resc.dat$idx,1]
               lab.adjC <-(lab.adj[i.dC,,drop=FALSE])[resc.datC$idx,3:4]
               lab.adj <- (lab.adj[i.d,,drop=FALSE])[resc.datC$idx,1:2]
               lab.ptsC <- lab.pts[resc.datC$idx,2]
               lab.pts  <- lab.pts[resc.datC$idx,1]

               if(draw.nonlbl&&n.ns>0){
                  resc.dat.ns <-.rescalefct(x.d.ns, function(x) absInfoEval(x,absInfo.f),
                              scaleX, scaleX.fct, #scaleX.inv,
                              scaleY, scaleY.fct[[1]], dots$xlim, dots$ylim, dots)
                  resc.datC.ns <-.rescalefct(x.dC.ns, function(x) absInfoEval(x,absInfoClass.f),
                              scaleX, scaleX.fct,# scaleX.inv,
                              scaleY, scaleY.fct[[1]], dots$xlim, dots$ylim, dots)

                  x.dr.ns <- resc.dat.ns$X
                  x.dCr.ns <- resc.datC.ns$X
                  y.dr.ns <- resc.dat.ns$Y
                  y.dCr.ns <- resc.datC.ns$Y
                  pch.nonlbl <- pch.nonlbl[resc.dat$idx,1]
                  col.nonlbl <- col.nonlbl[resc.dat$idx,1]
                  cex.nonlbl <- cex.nonlbl[resc.dat$idx,1]
                  pch.nonlblC <- pch.nonlbl[resc.datC$idx,2]
                  col.nonlblC <- col.nonlbl[resc.datC$idx,2]
                  cex.nonlblC <- cex.nonlbl[resc.datC$idx,2]
               }else{
                  x.dr.ns <- x.dCr.ns <- y.dr.ns <- y.dCr.ns <- NULL             
                  pch.nonlbl <- col.nonlbl <- cex.nonlbl <- NULL
                  pch.nonlblC <- col.nonlblC <- cex.nonlblC <- NULL              
               }
            
               dots.points <-   .makedotsPt(dots)

               do.pts <- function(x,y,cxa,ca,pa, dots)
                    do.call(points,args=c(list(x,y,cex=cxa,col=ca,pch=pa),
                            dots))
               tx <- function(xa,ya,lb,cx,ca,ft,ad){
                     nxa <- length(xa)
                     for(i in 1:nxa)
                         text(x=xa[i], y=ya[i], labels=lb[i], cex=cx[i], 
                              col=ca[i], font=ft, adj=ad[i,])
               }
               alp.v <- matrix(
                        rep(alpha.trsp, length.out = n0*(dims0+in1to.draw)),
                        n0,(dims0+in1to.draw))
               alp.v.s <- alp.v[i.d,,drop=FALSE]         
#               print(dim(alp.v.s))
               alp.vC.s <- alp.v[i.dC,,drop=FALSE]         
               if(draw.nonlbl&&n.ns>0){
                  alp.v.ns <- alp.v[i.d.ns,,drop=FALSE]         
                  alp.vC.ns <- alp.v[i.dC.ns,,drop=FALSE]         
               }else  alp.v.ns <- alp.vC.ns <- NULL
################################################################################
#  2.6 inserting the code to plot in data into panel last
################################################################################

               pL.abs <- substitute({
                   ICy0r1 <- ICy0r
                   ICy0cr1 <- ICy0cr
                   
                   if(is(distr, "DiscreteDistribution")){
                      ICy0r1 <- jitter(ICy0r1, factor = jit.fac0[1])
                      ICy0cr1 <- jitter(ICy0cr1, factor = jit.fac0[2])
                   }else{
                      if(any(.isReplicated(ICy0r1, jit.tol0[1]))&&jit.fac0[1]>0)
                         ICy0r1 <- jitter(ICy0r1, factor = jit.fac0[1])
                      if(any(.isReplicated(ICy0cr1, jit.tol0[2]))&&jit.fac0[2]>0)
                         ICy0cr1 <- jitter(ICy0cr1, factor = jit.fac0[2])
                   }

                   n.s <- length(ICy0r)

                   c1fun <- if(is.null(cexfun)) NULL else cexfun[[1]]
                   c2fun <- if(is.null(cexfun)) NULL else cexfun[[2]]
                   f1 <- .cexscale(ICy0,ICy0c,cex=cex0, fun = c1fun)
                   f1c <- .cexscale(ICy0c,ICy0,cex=cex0C, fun = c2fun)

                   
                   col1.pts <- c(.alphTrspWithNA(col0,al0))
                   col1.ptsC <- c(.alphTrspWithNA(col0C,al0C))

#                   print(length(y0))
#                   print(length(y0c))
#                   print(length(ICy0r1))
#                   print(length(ICy0cr1))
                   do.pts(y0,  ICy0r1,  f1, col1.pts, pch0,dots.points)
                   do.pts(y0c, ICy0cr1, f1c,col1.ptsC,pch0C,dots.points)
                   
                   if(with.lab0){
                      tx(y0,  ICy0r1,  lab.pts0,  f1/2,  lab.col0, lab.ft0[1], lab.ad0)
                      tx(y0c, ICy0cr1, lab.pts0C, f1c/2, lab.col0C, lab.ft0[2], lab.ad0C)
                   }
#-------------------------------------------
                   if(dononlb){
                       ICy0r1.ns <- ICy0r.ns
                       ICy0cr1.ns <- ICy0cr.ns
                       if(is(distr, "DiscreteDistribution")){
                          ICy0r1.ns <- jitter(ICy0r1.ns, factor = jit.fac0[1])
                          ICy0cr1.ns <- jitter(ICy0cr1.ns, factor = jit.fac0[2])
                       }else{
                          if(any(.isReplicated(ICy0r1.ns, jit.tol0[1]))&&jit.fac0[1]>0)
                             ICy0r1.ns <- jitter(ICy0r1.ns, factor = jit.fac0[1])
                          if(any(.isReplicated(ICy0cr1.ns, jit.tol0[2]))&&jit.fac0[2]>0)
                             ICy0cr1.ns <- jitter(ICy0cr1.ns, factor = jit.fac0[2])
                       }
    
                       n.ns <- length(ICy0r.ns)
    
                       c1fun.ns <- if(is.null(cexfun.ns)) NULL else cexfun.ns[[1]]
                       c2fun.ns <- if(is.null(cexfun.ns)) NULL else cexfun.ns[[2]]
                       f1.ns <-  .cexscale(ICy0.ns,ICy0c.ns,cex=cex0.ns, fun = c1fun.ns)
                       f1c.ns <- .cexscale(ICy0c.ns,ICy0.ns,cex=cex0.nsC, fun = c2fun.ns)
    
                       col1.ns  <- c(.alphTrspWithNA(col0.ns,al0.ns))
                       col1.nsC <- c(.alphTrspWithNA(col0.nsC,al0C.ns))
    
#                   print(length(y0.ns))
#                   print(length(y0c.ns))
#                   print(length(ICy0r1.ns))
#                   print(length(ICyc0r1.ns))
                       do.pts(y0.ns,  ICy0r1.ns,  f1.ns, col1.ns, pch0.ns)
                       do.pts(y0c.ns, ICy0cr1.ns, f1c.ns,col1.nsC,pch0.nsC)
                   }

                   pL0
                   }, list(ICy0c = y.dC, ICy0 = y.d,
                           ICy0r = y.dr, ICy0cr = y.dCr,
                           ICy0c.ns = y.dC.ns, ICy0.ns = y.d.ns,
                           ICy0r.ns = y.dr.ns, ICy0cr.ns = y.dCr.ns,
                           pL0 = pL, 
                           y0 = x.dr, y0c = x.dCr,
                           y0.ns = x.dr.ns, y0c.ns=x.dCr.ns,
                           col0 = col.pts, col0.ns = col.nonlbl,
                           pch0 = pch.pts, pch0.ns = pch.nonlbl,
                           cex0 = cex.pts, cex0.ns = cex.nonlbl,
                           col0C = col.ptsC, col0.nsC = col.nonlblC,
                           pch0C = pch.ptsC, pch0.nsC = pch.nonlblC,
                           cex0C = cex.ptsC, cex0.nsC = cex.nonlblC,
                           cexfun = cex.pts.fun,cexfun.ns = cex.nonlbl.fun, 
                           al0 = alp.v.s[,1], 
                           al0C = alp.vC.s[,1], 
                           al0.ns = alp.v.ns[,1], 
                           al0C.ns = alp.vC.ns[,1], 
                           n0 = n, with.lab0 = with.lab, 
                           lab.col0 = lab.col, 
                           lab.col0C = lab.colC,lab.ft0 = lab.font, 
                           lab.pts0 = lab.pts, lab.pts0C = lab.ptsC, 
                           lab.ad0 = lab.adj,lab.ad0C = lab.adjC,
                           jit.fac0 = jit.fac, jit.tol0=jit.tol,
                           dononlb = draw.nonlbl&(n.ns>0)
                           )
                   )
               pL.rel <- substitute({
                     y0.vec <- sapply(y0,  IC1.i.5@Map[[indi]])^2/ICy0r
                     y0c.vec <- sapply(y0c, classIC.i.5@Map[[indi]])^2/ICy0cr
                   if(is(distr, "DiscreteDistribution")){
                      y0.vec <- jitter(y0.vec, factor = jit.fac0[1])
                      y0c.vec <- jitter(y0c.vec, factor = jit.fac0[2])
                   }else{
                      if(any(.isReplicated(y0.vec, jit.tol0[1]))&&jit.fac0[1]>0)
                         y0.vec <- jitter(y0.vec, factor = jit.fac0[1])
                      if(any(.isReplicated(y0c.vec, jit.tol0[2]))&&jit.fac0[2]>0)
                         y0c.vec <- jitter(y0c.vec, factor = jit.fac0[2])
                   }


                   n.s <- length(y0.vec)

                   col1.pts  <- c(.alphTrspWithNA(col0,al0[,i]))
                   col1.ptsC <- c(.alphTrspWithNA(col0C,al0C[,i]))

                   dotsP0 <- dotsP
                   resc.rel <- .rescalefct(y0, cbind(y0.vec,ICy0),
                              scaleX, scaleX.fct, #scaleX.inv,
                              FALSE, scaleY.fct[[i]], dots$xlim, dots$ylim, dotsP0)
                   resc.rel.c <- .rescalefct(y0c, cbind(y0c.vec,ICy0c),
                              scaleX, scaleX.fct,# scaleX.inv,
                              FALSE, scaleY.fct[[i]], dots$xlim, dots$ylim, dotsP0)

                   c1fun <- if(is.null(cexfun)) NULL else cexfun[[(i1-1)*2+1]]
                   c2fun <- if(is.null(cexfun)) NULL else cexfun[[(i1-1)*2+2]]


                   f1 <- .cexscale(resc.rel$scy,resc.rel.c$scy,cex=cex0, fun=c1fun)
                   f1c <- .cexscale(resc.rel.c$scy,resc.rel$scy,cex=cex0C, fun=c2fun)

                   do.pts(resc.rel$X,   resc.rel$Y,   f1, col1.pts, pch0, dots.points)
                   do.pts(resc.rel.c$X, resc.rel.c$Y, f1c,col1.ptsC,pch0C,dots.points)


                   if(with.lab0){
                      tx(resc.rel$X,   resc.rel$Y,   lab.pts0,  f1/2,  lab.col0, lab.ft0[1],  lab.ad0)
                      tx(resc.rel.c$X, resc.rel.c$Y, lab.pts0C, f1c/2, lab.col0C, lab.ft0[2], lab.ad0C)
                   }

#--------------------------------------------------------------------------------
                   if(dononlb){
                       y0.vec.ns <- sapply(y0.ns,  IC1.i.5@Map[[indi]])^2/ICy0r.ns
                       y0c.vec.ns <- sapply(y0c.ns, classIC.i.5@Map[[indi]])^2/ICy0cr.ns
                       if(is(distr, "DiscreteDistribution")){
                          y0.vec.ns <- jitter(y0.vec.ns, factor = jit.fac0[1])
                          y0c.vec.ns <- jitter(y0c.vec.ns, factor = jit.fac0[2])
                       }else{
                          if(any(.isReplicated(y0.vec.ns, jit.tol0[1]))&&jit.fac0[1]>0)
                             y0.vec.ns <- jitter(y0.vec.ns, factor = jit.fac0[1])
                          if(any(.isReplicated(y0c.vec.ns, jit.tol0[2]))&&jit.fac0[2]>0)
                             y0c.vec.ns <- jitter(y0c.vec.ns, factor = jit.fac0[2])
                       }
    
                       n.ns <- length(y0.vec.ns)
    
                       col1.ns <- c(.alphTrspWithNA(col.ns,al0.ns[,i]))
                       col1.nsC <- c(.alphTrspWithNA(col.nsC,al0C.ns[,i]))
    
                       resc.rel.ns <- .rescalefct(y0.ns, cbind(y0.vec.ns,ICy0.ns),
                                  scaleX, scaleX.fct, #scaleX.inv,
                                  FALSE, scaleY.fct[[i]], dots$xlim, dots$ylim, dotsP0)
                       resc.rel.c.ns <- .rescalefct(y0c.ns, cbind(y0c.vec.ns,ICy0c.ns),
                                  scaleX, scaleX.fct,# scaleX.inv,
                                  FALSE, scaleY.fct[[i]], dots$xlim, dots$ylim, dotsP0)
    
                       c1fun.ns <- if(is.null(cexfun.ns)) NULL else cexfun.ns[[1]]
                       c2fun.ns <- if(is.null(cexfun.ns)) NULL else cexfun.ns[[2]]
      
                       f1.ns <- .cexscale(resc.rel.ns$scy,resc.rel.c.ns$scy,cex=cex0.ns, fun=c1fun.ns)
                       f1c.ns <- .cexscale(resc.rel.c.ns$scy,resc.rel.ns$scy,cex=cex0.nsC, fun=c2fun.ns)
    
                       do.pts(resc.rel.ns$X, resc.rel.ns$Y, f1.ns,col.ns,pch0.ns)
                       do.pts(resc.rel.c.ns$X, resc.rel.c.ns$Y, f1c.ns,col.nsC,pch0.nsC)
                   }

                   pL0
                   }, list(ICy0c = y.dC, ICy0 = y.d,
                           ICy0r = y.dr, ICy0cr = y.dCr,
                           ICy0c.ns = y.dC.ns, ICy0.ns = y.d.ns,
                           ICy0r.ns = y.dr.ns, ICy0cr.ns = y.dCr.ns,
                           pL0 = pL, y0 = x.dr, y0c = x.dCr, 
                           y0.ns = x.dr.ns, y0c.ns = x.dCr.ns,
                           cex0 = cex.pts, cex0.ns = cex.nonlbl,
                           col0 = col.pts, col0.ns = col.nonlbl,
                           pch0 = pch.pts, pch0.ns = pch.nonlbl, 
                           cex0C = cex.ptsC, cex0.nsC = cex.nonlblC,
                           col0C = col.ptsC, col0.nsC = col.nonlblC,
                           pch0C = pch.ptsC, pch0.nsC = pch.nonlblC, 
                           al0 = alp.v.s, al0C = alp.vC.s, 
                           al0.ns = alp.v.ns, al0C.ns = alp.vC.ns, 
                           n0 = n, with.lab0 = with.lab, 
                           lab.col0 = lab.col, lab.col0C = lab.colC,lab.ft0 = lab.font, 
                           lab.pts0 = lab.pts, lab.pts0C = lab.pts, 
                           lab.ad0 = lab.adj, lab.ad0C = lab.adjC,
                           jit.fac0 = jit.fac, jit.tol0=jit.tol, 
                           cexfun = cex.pts.fun, cexfun.ns = cex.nonlbl.fun, 
                           dononlb = draw.nonlbl&(n.ns>0)
                           )
                   )
            }

################################################################################
#  3. creating the panel plots
################################################################################
#  3.1. if present: the absInfo-Plot
################################################################################
            if(!is.null(ylim))
                dotsP$ylim <- ylim[,1]       
            
            fac.leg <- if(dims0>1) 3/4 else .75/.8 


            dotsP$axes <- NULL

            ipInfo <- vector("list",0)

            if(1 %in% to.draw){
               resc <-.rescalefct(x.vec, function(x) absInfoEval(x,absInfo.f),
                              scaleX, scaleX.fct, #scaleX.inv,
                              scaleY, scaleY.fct[[1]], dots$xlim, dots$ylim, dotsP)
               resc.C <-.rescalefct(x.vec, function(x) absInfoEval(x,absInfoClass.f),
                              scaleX, scaleX.fct, #scaleX.inv,
                              scaleY, scaleY.fct[[1]], dots$xlim, dots$ylim, dotsP)
               dotsP1 <- dotsP <- resc$dots
               dotsP$yaxt <- dots$yaxt

               ipInfo$Abs <- vector("list",0)
               ipInfo$Abs$parArgsL <- parArgsL[[1]]               
               do.call(par, args = parArgsL[[1]])

               plot.args <- c(list(x=resc.C$X, y=resc.C$Y, type = plty,
                   lty = ltyI, col = colI, lwd = lwdI,
                   xlab = .mpresubs(xlab), ylab = .mpresubs(ylab.abs), 
                   panel.last = pL.abs,
                   panel.first = pF.abs
                   ),
                   dotsP1)

#                print(summary(plot.args$x))
#                print(summary(plot.args$y))
               ipInfo$Abs$plot.args <- plot.args
               do.call(plot, args=plot.args)
               rm(plot.args)
               
               lines.args <- c(list(resc$X, resc$Y, type = plty,
                       lty = lty, lwd = lwd, col = col), dotsL)
               ipInfo$Abs$lines.args <- lines.args
               do.call(lines, args=lines.args)
               rm(lines.args)
               
               
               scaleX0 <- scaleX & (xaxt0[1]!="n")
               scaleY0 <- scaleY & (yaxt0[1]!="n")
               x.ticks0 <- if(xaxt0[1]!="n") x.ticks else NULL
               y.ticks0 <- if(yaxt0[1]!="n") y.ticks[[1]] else NULL

               finiteEndpoints <- rep(FALSE,4)
               if(scaleX){
                  finiteEndpoints[1] <- is.finite(scaleX.inv(min(resc.C$X, xlim[1],na.rm=TRUE)))
                  finiteEndpoints[2] <- is.finite(scaleX.inv(max(resc.C$X, xlim[2],na.rm=TRUE)))
               }
               if(scaleY){
                  finiteEndpoints[3] <- is.finite(scaleY.inv[[1]](min(resc.C$Y, ylim[1,1],na.rm=TRUE)))
                  finiteEndpoints[4] <- is.finite(scaleY.inv[[1]](max(resc.C$Y, ylim[2,1],na.rm=TRUE)))
               }

               .plotRescaledAxis.args <- list(scaleX0, scaleX.fct, scaleX.inv,
                              scaleY0,scaleY.fct, scaleY.inv,
                              dots$xlim, dots$ylim, resc$X, ypts = 400,
                              n = scaleN, x.ticks = x.ticks0,
                              y.ticks = y.ticks0, withbox = withbox)
               ipInfo$Abs$.plotRescaledAxis.args <- .plotRescaledAxis.args
               do.call(.plotRescaledAxis,args=.plotRescaledAxis.args)
               rm(.plotRescaledAxis.args)
               
               if(with.legend){
                  legend.args <- list(.legendCoord(legend.location[[1]], scaleX, 
                        scaleX.fct, scaleY, scaleY.fct[[1]]), legend = legend[[1]], 
                     bg = legend.bg, lty = c(ltyI, lty), col = c(colI, col), 
                     lwd = c(lwdI, lwd), cex = legend.cex*fac.leg)
                  ipInfo$Abs$legend.args <- legend.args
                  do.call(graphics::legend, args=legend.args)
                  rm(legend.args)
               }

               if(innerL){
                  title.args <- c(list(main = innerT[[1]]), dotsT, 
                                  list(line = lineT, 
                                  cex.main = cex.inner, col.main = col.inner))
                  ipInfo$Abs$title.args <- title.args
                  do.call(title, args=title.args)
                  rm(title.args)       
               }   
            }
            
################################################################################
#  3.2. if present: the relInfoplots in a for loop
################################################################################
            ipInfo$Rel <- vector("list",length(to.draw[to.draw!=1]))
            if(dims > 1 && length(to.draw[to.draw!=1])>0){
                nrows <- trunc(sqrt(dims0))
                ncols <- ceiling(dims0/nrows)
                if (!withSweave||!mfColRow)
                     dN <- substitute({devNew()}) else substitute({})

                IC1.i.5 <- QF.5%*%IC1
                classIC.i.5 <- QFc.5%*%classIC
                if(!in1to.draw){
                   resc <-.rescalefct(x.vec, function(x) absInfoEval(x,absInfo.f),
                                  scaleX, scaleX.fct, #scaleX.inv,
                                  scaleY, scaleY.fct[[1]], dots$xlim, dots$ylim, dotsP)
                   resc.C <-.rescalefct(x.vec, function(x) absInfoEval(x,absInfoClass.f),
                                  scaleX, scaleX.fct, #scaleX.inv,
                                  scaleY, scaleY.fct[[1]], dots$xlim, dots$ylim, dotsP)
                   dotsP1 <- dotsP <- resc$dots
                   dotsP$yaxt <- dots$yaxt
                    
                }
                for(i in 1:dims0){
                    indi <- to.draw1[i]-1
                    i1 <- i + in1to.draw
                    if(!is.null(ylim)) 
                         dotsP$ylim <- ylim[,in1to.draw+i]       
                    else dotsP$ylim <- c(0,1)

                    y.vec1 <- sapply(resc$x, IC1.i.5@Map[[indi]])^2/
                              absInfoEval(resc$x,absInfo.f)
                    y.vec1C <- sapply(resc.C$x, classIC.i.5@Map[[indi]])^2/
                              absInfoEval(resc.C$x,absInfoClass.f)

                    if(mfColRow){
                       parArgsL[[i+in1to.draw]] <- c(parArgsL[[i+in1to.draw]],
                                                  list(mfrow = c(nrows, ncols)))
                       eval(dN)
                       if(i==1) do.call(par,args=parArgsL[[i+in1to.draw]])
                    }else{do.call(par,args=parArgsL[[i+in1to.draw]])}
                    
                    ipInfo$Rel[[i]]$parArgsL <- parArgsL[[i+in1to.draw]]               

                    if(!is.null(dotsP$log)) if(dotsP$log=="y") dotsP$log <- NULL
                    if(!is.null(dotsP$log)) if(dotsP$log=="xy") dotsP$log <- "x"
                    
                    plot.args <- c(list(resc$X, y.vec1, type = plty,
                                  lty = lty, xlab = .mpresubs(xlab), ylab = .mpresubs(ylab.rel),
                                  col = col, lwd = lwd, panel.last = pL.rel,
                                  panel.first = pF.rel[[i]]),  dotsP)
                    do.call(plot, args=plot.args)
                    ipInfo$Rel[[i]]$plot.args <- plot.args 
                    rm(plot.args)

                    lines.args <- c(list(resc.C$X, y.vec1C, type = plty,
                            lty = ltyI, col = colI, lwd = lwdI), dotsL)
                    do.call(lines, args = lines.args)
                    ipInfo$Rel[[i]]$lines.args <- lines.args 
                    rm(lines.args)
                    
                    scaleX0 <- scaleX & (xaxt0[i+in1to.draw]!="n")
                    scaleY0 <- scaleY & (yaxt0[i+in1to.draw]!="n")
                    x.ticks0 <- if(xaxt0[i+in1to.draw]!="n") x.ticks else NULL
                    y.ticks0 <- if(yaxt0[i+in1to.draw]!="n") y.ticks[[i+in1to.draw]] else NULL

                    finiteEndpoints <- rep(FALSE,4)
                    if(scaleX){
                      finiteEndpoints[1] <- is.finite(scaleX.inv(min(resc$X, xlim[1],na.rm=TRUE)))
                      finiteEndpoints[2] <- is.finite(scaleX.inv(max(resc$X, xlim[2],na.rm=TRUE)))
                    }
                    if(scaleY){
                       finiteEndpoints[3] <- is.finite(scaleY.inv[[i]](min(y.vec1, ylim[1,i+in1to.draw],na.rm=TRUE)))
                       finiteEndpoints[4] <- is.finite(scaleY.inv[[i]](max(y.vec1, ylim[2,i+in1to.draw],na.rm=TRUE)))
                    }

                    .plotRescaledAxis.args <- list(scaleX0, scaleX.fct, scaleX.inv,
                              FALSE,scaleY.fct[[i]],
                              scaleY.inv[[i]], dots$xlim,
                              dots$ylim, resc$X, ypts = 400, n = scaleN,
                              finiteEndpoints = finiteEndpoints,
                              x.ticks = x.ticks0,
                              y.ticks = y.ticks0, withbox = withbox)                              
                    do.call(.plotRescaledAxis, args=.plotRescaledAxis.args)
                    ipInfo$Rel[[i]]$.plotRescaledAxis.args <- .plotRescaledAxis.args 
                    rm(.plotRescaledAxis.args)          

                    if(with.legend){
                      legend.args <- list(.legendCoord(legend.location[[i1]],
                                 scaleX, scaleX.fct, scaleY, scaleY.fct[[i]]),
                           bg = legend.bg, legend = legend[[i1]],
                           col = c(colI, col), lwd = c(lwdI, lwd),
                           lty = c(ltyI, lty), cex = legend.cex*fac.leg)
                      do.call(graphics::legend, args= legend.args)
                      ipInfo$Rel[[i]]$legend.args <-legend.args 
                      rm(legend.args)                           
                    }
                    if(innerL){
                       title.args  <- c(list(main = innerT[[1+indi]]),  
                               dotsT, line = lineT, cex.main = cex.inner, 
                               col.main = col.inner)
                       do.call(title, args = title.args)
                       ipInfo$Rel[[i]]$title.args <- title.args         
                       rm(title.args)        
                    }           
                }
            }
################################################################################
#  4. outer titles
################################################################################
        cex.main <- if(!hasArg("cex.main")) par("cex.main") else dots$"cex.main"
        col.main <- if(!hasArg("col.main")) par("col.main") else dots$"col.main"
        if (mainL){
            main.args <- list(text = main, side = 3, cex = cex.main, adj = .5,
                  outer = TRUE, padj = 1.4, col = col.main)
            do.call(mtext, args=main.args)
            ipInfo$main.args <- main.args 
            rm(main.args)
                  
        }
        cex.sub <- if(!hasArg("cex.sub")) par("cex.sub") else dots$"cex.sub"
        col.sub <- if(!hasArg("col.sub")) par("col.sub") else dots$"col.sub"
        if (subL){
            sub.args <- list(text = sub, side = 1, cex = cex.sub, adj = .5,
                  outer = TRUE, line = -1.6, col = col.sub)
            do.call(mtext, args=sub.args)
            ipInfo$sub.args <- sub.args 
            rm(sub.args)
        }
        class(ipInfo) <- c("infoPlotInfo","DiagnInfo")
        retv <- list(call=mc, infoPlotInfo = ipInfo)
        
        if(return.Order){ 
           retOrder <- list(IC=i0.d,IC.class=i0.dC)
           retv$retOrder <- retOrder         
        }
        invisible(return(retv))
        }
    )
 
