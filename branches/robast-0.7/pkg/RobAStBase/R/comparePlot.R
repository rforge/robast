setMethod("comparePlot", signature("IC","IC"),
    function(obj1,obj2, obj3 = NULL, obj4 = NULL, 
             ..., withSweave = getdistrOption("withSweave"), 
             main = FALSE, inner = TRUE, sub = FALSE, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3], 
             mfColRow = TRUE){

        xc1 <- as.character(deparse(match.call(call = sys.call(sys.parent(1)))$obj1))
        xc2 <- as.character(deparse(match.call(call = sys.call(sys.parent(1)))$obj2))
        xc <- c(xc1,xc2)
        if(!is.null(obj3))
            xc <- c(xc,as.character(deparse(match.call(call = sys.call(sys.parent(1)))$obj3)))
        if(!is.null(obj4))
            xc <- c(xc,as.character(deparse(match.call(call = sys.call(sys.parent(1)))$obj4)))
        
        dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."

        ncomp <- 2+ !is.null(obj3) +  !is.null(obj4)
         
        if(is.null(dots[["col"]]))   dots$"col" <- 1:ncomp
        if(is.null(dots[["lwd"]]))   dots$"lwd" <- 1
        
        col <- dots[["col"]]
        lwd <- dots[["lwd"]]
        
        if(!is.null(dots[["lty"]]))  dots["lty"] <- NULL
        if(!is.null(dots[["type"]])) dots["type"] <- NULL
        if(!is.null(dots[["xlab"]])) dots["xlab"] <- NULL
        if(!is.null(dots[["ylab"]])) dots["ylab"] <- NULL
        
        dotsP <- dotsL <- dotsT <- dots


        L2Fam <- eval(obj1@CallL2Fam)
        L2Fam1c <- obj1@CallL2Fam
        L2Fam2c <- obj2@CallL2Fam
        if(!identical(L2Fam1c,L2Fam2c))
            stop("ICs need to be defined for the same model")

        e1 <- L2Fam@distribution
        if(!is(e1, "UnivariateDistribution")) stop("not yet implemented")

        if(is(e1, "AbscontDistribution")){
            lower <- ifelse(is.finite(q(e1)(0)), q(e1)(0), q(e1)(getdistrOption("TruncQuantile")))
            upper <- ifelse(is.finite(q(e1)(1)), q(e1)(1), q(e1)(1 - getdistrOption("TruncQuantile")))
            h <- upper - lower
            x.vec <- seq(from = lower - 0.1*h, to = upper + 0.1*h, length = 1000)
            plty <- "l"
            lty <- "solid"
        }else{
            if(is(e1, "DiscreteDistribution")){
                x.vec <- support(e1)
                plty <- "p"
                lty <- "dotted"
            }else{
                x.vec <- r(e1)(1000)
                x.vec <- sort(unique(x.vec))
                plty <- "p"
                lty <- "dotted"
            }
        }

        dims <- nrow(L2Fam@param@trafo)
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
            innerParam <-  paste(gettext("\nwith main parameter ("), 
                                    paste(round(L2Fam@param@main, 3), 
                                          collapse = ", "),
                                 ")", sep = "")
            if(!is.null(L2Fam@param@nuisance))
                innerParam <- paste(innerParam,
                                gettext("\nand nuisance parameter ("), 
                                    paste(round(L2Fam@param@nuisance, 3), 
                                           collapse = ", "),
                                ")", sep ="")
            if(!is.null(L2Fam@param@fixed))
                innerParam <- paste(innerParam,
                                gettext("\nand fixed known parameter ("), 
                                    paste(round(L2Fam@param@fixed, 3), 
                                           collapse = ", "),
                                ")", sep ="")
            
            if(!is.logical(inner)){
                if(!is.character(inner))
                    stop("Argument 'inner' must either be 'logical' or a character vector")
                innerT <- rep(inner,length.out=dims)
                innerL <- TRUE
            }else{if(any(is.na(inner))||any(!inner)) {
                 innerT <- ""; innerL <- FALSE
                }else{innerL <- TRUE
                      innerT <- paste(paste(gettext("Component "),  1:dims, 
                                       gettext(" of (partial) IC\nfor "), 
                                       name(L2Fam)[1], sep =""), innerParam)
                   }
              }


        w0 <- getOption("warn")
        options(warn = -1)
        on.exit(options(warn = w0))
        opar <- par()
        on.exit(par(opar))
        nrows <- trunc(sqrt(dims))
        ncols <- ceiling(dims/nrows)
        par(mfrow = c(nrows, ncols))

        if(is(e1, "DiscreteDistribution"))
                x.vec1 <- seq(from = min(x.vec), to = max(x.vec), length = 1000)

            dotsT["main"] <- NULL
            dotsT["cex.main"] <- NULL
            dotsT["col.main"] <- NULL
            dotsT["line"] <- NULL

        for(i in 1:dims){
            matp  <- cbind(sapply(x.vec, IC1@Map[[i]]),sapply(x.vec, IC2@Map[[i]]))
            if(is(obj3, "IC"))
                matp  <- cbind(matp,sapply(x.vec, IC3@Map[[i]]))
            if(is(obj4, "IC"))
                matp  <- cbind(matp,sapply(x.vec, IC4@Map[[i]]))

            do.call(matplot, args=c(list( x= x.vec, y=matp,
                 type = plty, lty = lty,
                 xlab = "x", ylab = "(partial) IC"), dotsP))

            if(is(e1, "DiscreteDistribution")){
                 matp1 <- cbind(sapply(x.vec1, IC1@Map[[i]]),sapply(x.vec1, IC2@Map[[i]]))
                 if(is(obj3, "IC"))
                    matp1  <- cbind(matp1,sapply(x.vec1, IC3@Map[[i]]))
                 if(is(obj4, "IC"))
                    matp1  <- cbind(matp1,sapply(x.vec1, IC4@Map[[i]]))
                 do.call(matlines, c(list(x.vec1, matp1, lty = "dotted"),dotsL))
                 }

           if(innerL)
              do.call(title, args=c(list(main = innerT[i]),  dotsT,
                      line = lineT, cex.main = cex.inner, col.main = col.inner))
        }
        
        legend("bottomright", 
               legend = xc, col = eval(dots[["col"]]), 
               cex=0.75, lwd=eval(dots[["lwd"]])*1.5)

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
