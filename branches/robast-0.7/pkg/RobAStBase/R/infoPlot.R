setMethod("infoPlot", "IC",
    function(object, ..., withSweave = getdistrOption("withSweave"), 
             colI = grey(0.5), lwdI = 0.7*par("lwd"),
             main = FALSE, inner = TRUE, sub = FALSE, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3], 
             mfColRow = TRUE, to.draw.arg = NULL){

        objectc <- match.call(call = sys.call(sys.parent(1)))$object
        dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
                   

        L2Fam <- eval(object@CallL2Fam)
        
        if(!hasArg(col)) col <- par("col") else col <- dots$col
        if(!hasArg(lwd)) lwd <- par("lwd") else lwd <- dots$lwd
        
        if(!is.null(dots[["lty"]]))  dots["lty"] <- NULL
        if(!is.null(dots[["type"]])) dots["type"] <- NULL
        if(!is.null(dots[["xlab"]])) dots["xlab"] <- NULL
        if(!is.null(dots[["ylab"]])) dots["ylab"] <- NULL
        
        trafO <- trafo(L2Fam@param)
        dims <- nrow(trafO)
        dimm <- length(L2Fam@param)
        
        to.draw <- 1:(dims+1)
        dimnms  <- c(rownames(trafO))
        if(is.null(dimnms))
           dimnms <- paste("dim",1:dims,sep="")
        pdimnms <- c("Abs",dimnms)
        if(!mfColRow && ! is.null(to.draw.arg)){
            if(is.character(to.draw.arg)) 
                 to.draw <- pmatch(to.draw.arg, pdimnms)
            else if(is.numeric(to.draw.arg)) 
                 to.draw <- to.draw.arg
        }
        
        to.draw1 <- to.draw[to.draw>1]
        dims0 <- length(to.draw1)
        nrows <- trunc(sqrt(dims0))
        ncols <- ceiling(dims0/nrows)

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
                lower <- if(is.finite(q(e1)(0)))
                     q(e1)(0) else q(e1)(getdistrOption("TruncQuantile"))
                upper <- if(is.finite(q(e1)(1))) 
                     q(e1)(1) else q(e1)(1 - getdistrOption("TruncQuantile"))
                if(!is.null(xlim)){ 
                  lower <- min(lower,xm)
                  upper <- max(upper,xM)
                }
                h <- upper - lower
                x.vec <- seq(from = lower - 0.1*h, to = upper + 0.1*h, length = 1000)
                plty <- "l"
                lty <- "solid"
            }else{
                if(is(e1, "DiscreteDistribution")) x.vec <- support(e1)
                else{
                   x.vec <- r(e1)(1000)
                   x.vec <- sort(unique(x.vec))
                }
                plty <- "p"
                lty <- "dotted"
                if(!is.null(xlim)) x.vec <- x.vec[(x.vec>=xm) & (x.vec<=xM)]
            }
         }
         ylim <- eval(dots$ylim)
         if(!is.null(ylim)){ 
               if(!length(ylim) %in% c(2,2*(dims0+(1%in%to.draw)))) 
                  stop("Wrong length of Argument ylim"); 
               ylim <- matrix(ylim, nrow=2,ncol=dims0+(1%in%to.draw))
               dots$ylim <- NULL
         }

         dotsP <- dotsL <- dotsT <- dots
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
                      tnm  <- c(rownames(trafO))
                      tnms <- if(is.null(tnm)) paste(1:dims) else 
                                               paste("'", tnm, "'", sep = "") 
                      innerT <- as.list(paste(c( paste(gettext("Absolute information of (partial) IC for"), 
                                       name(L2Fam)[1], sep =""),
                                   paste(gettext("Relative information of \ncomponent "),
                                       tnms, 
                                       gettext(" of (partial) IC\nfor "), 
                                       name(L2Fam)[1], sep ="")), innerParam))
                   }
              }


            QFc <- diag(dimm)
            if(is(object,"ContIC") & dims>1 )
               {if (is(normtype(object),"QFNorm")) QFc <- QuadForm(normtype(object))
                QFc0 <- solve( trafo %*% solve(L2Fam@FisherInfo) %*% t(trafo ))
                if (is(normtype(object),"SelfNorm")|is(normtype(object),"InfoNorm")) 
                    QFc <- QFc0
               }
            QFc.5 <- sqrt(PosSemDefSymmMatrix(QFc))

            classIC <- as(trafo %*% solve(L2Fam@FisherInfo) %*% L2Fam@L2deriv, "EuclRandVariable")
            absInfoClass <- t(classIC) %*% QFc %*% classIC
            absInfoClass <- sapply(x.vec, absInfoClass@Map[[1]])

            QF <- diag(dimm)
            if(is(object,"ContIC") & dims>1 )
               {if (is(normtype(object),"QFNorm")) QF <- QuadForm(normtype(object))}
            QF.5 <- sqrt(PosSemDefSymmMatrix(QF))

            IC1 <- as(diag(dimm) %*% object@Curve, "EuclRandVariable")
            absInfo <- t(IC1) %*% QF %*% IC1
            absInfo <- sapply(x.vec, absInfo@Map[[1]])

            
            w0 <- getOption("warn")
            options(warn = -1)
            on.exit(options(warn = w0))
            opar <- par()
            on.exit(par(opar))
#            if (!withSweave)
#               devNew()

            omar <- par("mar")
            parArgs <- list(mar = c(bmar,omar[2],tmar,omar[4]))
            do.call(par,args=parArgs)

            
            
            dotsP["col"] <- NULL
            dotsP["lwd"] <- NULL
            if(!is.null(ylim)) 
                dotsP$ylim <- ylim[,1]       
            if(1 %in% to.draw){
               do.call(plot, args=c(list(x.vec, absInfoClass, type = plty, 
                   lty = "dashed", col = colI, lwd = lwdI,
                   xlab = "x", ylab = "absolute information"), dotsP))
               do.call(lines, args=c(list(x.vec, absInfo, type = plty, lty = lty), 
                       dotsL))
               legend("top",
                     legend = c("class. opt. IC", objectc), 
                     lty = c(lty,"dashed"), col = c(colI, col), 
                     lwd=c(lwdI, lwd), cex = 0.75)

               dotsT["main"] <- NULL
               dotsT["cex.main"] <- NULL
               dotsT["col.main"] <- NULL
               dotsT["line"] <- NULL
               if(innerL)
                  do.call(title, args=c(list(main = innerT[[1]]),  dotsT,
                          line = lineT, cex.main = cex.inner, col.main = col.inner))
            }
            
            if(dims0 > 1){
                dotsP["ylim"] <- NULL
                dotsL["ylim"] <- NULL
                dotsT["ylim"] <- NULL
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
                         dotsP$ylim <- ylim[,(1%in%to.draw)+i]       
                    else dotsP$ylim <- c(0,1)
                    y.vec <- sapply(x.vec, IC1.i.5@Map[[indi]])^2/absInfo
                    do.call(plot, args=c(list(x.vec, y.vec, type = plty, 
                                  lty = lty, xlab = "x", 
                                  ylab = "relative information", 
                                  col = colI, lwd = lwdI), dotsP))

                    yc.vec <- sapply(x.vec, classIC.i.5@Map[[indi]])^2/absInfoClass
                    do.call(lines, args=c(list(x.vec, yc.vec, type = plty, 
                          lty = "dashed"), dotsL))
                    legend("topright",
                           legend = c("class. opt. IC", objectc), lty = c(lty,"dashed"), 
                               col = c(colI, col), lwd=c(lwdI, lwd),
                               cex = 0.6)
                    if(innerL)
                       do.call(title, args=c(list(main = innerT[[1+indi]]),  dotsT,
                               line = lineT, cex.main = cex.inner, col.main = col.inner))
                }
            }
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


        invisible()
        }
    )
 