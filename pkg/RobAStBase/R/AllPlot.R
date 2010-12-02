setMethod("plot", signature(x = "IC", y = "missing"),
    function(x,...,withSweave = getdistrOption("withSweave"), 
             main = FALSE, inner = TRUE, sub = FALSE, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3],
             mfColRow = TRUE, to.draw.arg = NULL){

        xc <- match.call(call = sys.call(sys.parent(1)))$x
        dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."


        if(!is.logical(inner)){
          if(!is.list(inner))
              inner <- as.list(inner)
            #stop("Argument 'inner' must either be 'logical' or a 'list'")
           inner <- distr:::.fillList(inner,4)          
           innerD <- inner[1:3]
           innerL <- inner[4] 
        }else{innerD <- innerL <- inner}


        L2Fam <- eval(x@CallL2Fam)

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


        e1 <- L2Fam@distribution
        if(!is(e1, "UnivariateDistribution")) stop("not yet implemented")

        if(is(e1, "UnivariateDistribution")){
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
                lty <- "solid"
            }else{
                if(is(e1, "DiscreteDistribution")) x.vec <- support(e1)
                else{
                   x.vec <- r(e1)(1000)
                   x.vec <- sort(unique(x.vec))
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
        if(!is.null(dots[["xlab"]])) dots["xlab"] <- NULL
        if(!is.null(dots[["ylab"]])) dots["ylab"] <- NULL

        IC1 <- as(diag(dims) %*% x@Curve, "EuclRandVariable")

        mainL <- FALSE
        subL <- FALSE
        lineT <- NA

     .mpresubs <- function(inx)
                    distr:::.presubs(inx, c("%C", "%D", "%A"),
                          c(as.character(class(x)[1]),
                            as.character(date()),
                            as.character(deparse(xc))))

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
        innerT <- distr:::.fillList(innerT,dims)
        if(dims0<dims){
           innerT0 <- innerT
           for(i in 1:dims0) innerT[to.draw[i]] <- innerT0[i]          
        }
     }


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

        dotsT <- dots
        dotsT["main"] <- NULL
        dotsT["cex.main"] <- NULL
        dotsT["col.main"] <- NULL
        dotsT["line"] <- NULL


        dots$ylim <- NULL
        for(i in 1:dims0){
            indi <- to.draw[i]
            if(!is.null(ylim)) dots$ylim <- ylim[,i]       
            do.call(plot, args=c(list(x.vec, sapply(x.vec, IC1@Map[[indi]]), 
                                      type = plty, lty = lty,
                                      xlab = "x", ylab = "(partial) IC"),
                                 dots))     
            if(is(e1, "DiscreteDistribution")){
                x.vec1 <- seq(from = min(x.vec), to = max(x.vec), length = 1000)
                do.call(lines,args=c(list(x.vec1, sapply(x.vec1, IC1@Map[[indi]]), 
                                          lty = "dotted"), dots))
            }
            do.call(title,args=c(list(main = innerT[indi]), dotsT, line = lineT,
                    cex.main = cex.inner, col.main = col.inner))
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
    })


setMethod("plot", signature(x = "IC",y = "numeric"),
          function(x, y, ..., cex.pts = 1, col.pts = par("col"),
          pch.pts = 1, jitter.fac = 1, with.lab = FALSE,
          lab.pts = NULL, lab.font = NULL,
             which.lbs = NULL, which.Order  = NULL, return.Order = FALSE){
    dots <- match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)$"..."

    n <- if(!is.null(dim(y))) nrow(y) else length(y)
    oN0 <- NULL
    if(is.null(which.lbs))
       which.lbs <- 1:n
    which.lbs0 <- (1:n) %in% which.lbs
    which.lbx <- rep(which.lbs0, length.out=length(y))
    y0 <- y[which.lbx]
    n <- if(!is.null(dim(y0))) nrow(y0) else length(y0)
    oN <- (1:n)[which.lbs0]


    L2Fam <- eval(x@CallL2Fam)
    trafO <- trafo(L2Fam@param)
    dims <- nrow(trafO)
    dimm <- length(L2Fam@param)
    QF <- diag(dims)

    if(is(x,"ContIC") & dims>1 )
      {if (is(normtype(x),"QFNorm")) QF <- QuadForm(normtype(x))}

    IC1 <- as(diag(dims) %*% x@Curve, "EuclRandVariable")
    absInfo <- t(IC1) %*% QF %*% IC1
    ICMap <- IC1@Map

    absInfo <- sapply(y, absInfo@Map[[1]])
    absInfo0 <- absInfo[which.lbs]/max(absInfo)

    if (n==length(y0)) {
        oN <-  order(absInfo0)
        oN0 <- order(absInfo)
        oN0 <- oN0[oN0 %in% which.lbs]
        y0 <- y0[oN]
        if(!is.null(which.Order)){
            oN <- oN0[(n+1)-which.Order]
            y0 <- y[oN]
            absInfo0 <- absInfo[oN]/max(absInfo[oN])
        }
    }
    if(is.null(lab.pts)) lab.pts <- paste(oN)
    else {lab.pts <- rep(lab.pts, length.out=length(y))
          lab.pts <- lab.pts[oN]}

    dots.without <- dots
    dots.without$col <- dots.without$cex <- dots.without$pch <- NULL

    pL <- expression({})
    if(!is.null(dots$panel.last))
        pL <- dots$panel.last
    dots$panel.last <- NULL

    pL <- substitute({
        ICy <- sapply(y0s,ICMap0[[indi]])
        if(is(e1, "DiscreteDistribution"))
           ICy <- jitter(ICy, factor = jitter.fac0)
        do.call(points, args=c(list(y0s, ICy, cex = log(absy0+1)*3*cex0,
                        col = col0, pch = pch0), dwo0))
        if(with.lab0){
           text(x = y0s, y = ICy, labels = lab.pts0,
                cex = log(absy0+1)*1.5*cex0, col = col0)
        }
        pL0
        }, list(pL0 = pL, ICMap0 = ICMap, y0s = y0, absy0 = absInfo0,
                dwo0 = dots.without, cex0 = cex.pts, pch0 = pch.pts,
                col0 = col.pts, with.lab0 = with.lab, lab.pts0 = lab.pts,
                jitter.fac0 = jitter.fac
                ))

  do.call("plot", args = c(list(x = x, panel.last = pL), dots))
  if(return.Order) return(oN0)
  invisible()
})
