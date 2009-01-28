setMethod("comparePlot", signature("IC","IC"),
    function(obj1,obj2, obj3 = NULL, obj4 = NULL, 
             ..., withSweave = getdistrOption("withSweave"), 
             main = FALSE, inner = TRUE, sub = FALSE, 
             col = par("col"), lwd = par("lwd"), lty, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3], 
             legend.location = "bottomright", 
             mfColRow = TRUE, to.draw.arg = NULL){

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
        dimm <- length(L2Fam@param)
        
        to.draw <- 1:dims
        dimnms  <- c(rownames(trafO))
        if(is.null(dimnms))
           dimnms <- paste("dim",1:dims,sep="")
        if(!mfColRow && ! is.null(to.draw.arg)){
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
        IC1 <- as(diag(dimm) %*% obj1@Curve, "EuclRandVariable")
        IC2 <- as(diag(dimm) %*% obj2@Curve, "EuclRandVariable")


        obj <- obj3
        if(is(obj, "IC"))
           {
           if(!identical(L2Fam1c,obj@CallL2Fam))
               stop("ICs need to be defined for the same model")
           IC3 <- as(diag(dimm) %*% obj3@Curve, "EuclRandVariable")
           }

        obj <- obj4
        if(is(obj, "IC"))
           {
           if(!identical(L2Fam1c,obj@CallL2Fam))
               stop("ICs need to be defined for the same model")
           IC4 <- as(diag(dimm) %*% obj4@Curve, "EuclRandVariable")
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
        opar <- par()
        on.exit(par(opar))
        
        if(mfColRow)
             par(mfrow = c(nrows, ncols))

        if(is(e1, "DiscreteDistribution"))
                x.vec1 <- seq(from = min(x.vec), to = max(x.vec), length = 1000)
        
            dotsT["main"] <- NULL
            dotsT["cex.main"] <- NULL
            dotsT["col.main"] <- NULL
            dotsT["line"] <- NULL

        for(i in 1:dims0){
            indi <- to.draw[i]
            if(!is.null(ylim)) dotsP$ylim <- ylim[,i]       
            matp  <- cbind(sapply(x.vec, IC1@Map[[indi]]),
                           sapply(x.vec, IC2@Map[[indi]]))
            if(is(obj3, "IC"))
                matp  <- cbind(matp,sapply(x.vec, IC3@Map[[indi]]))
            if(is(obj4, "IC"))
                matp  <- cbind(matp,sapply(x.vec, IC4@Map[[indi]]))

            do.call(matplot, args=c(list( x= x.vec, y=matp,
                 type = plty, lty = lty, col = col, lwd = lwd,
                 xlab = "x", ylab = "(partial) IC"), dotsP))

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
        
        legend(legend.location, legend = xc, col = col, 
               cex = 0.75, lwd = lwd*1.5, lty = lty)

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
    })
