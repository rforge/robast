.isReplicated <- function(x, tol=.Machine$double.eps){
  tx <- table(x)
  rx <- as.numeric(names(tx[tx>1]))
  sapply(x, function(y) any(abs(y-rx)<tol))
}

.plotData <- function(
  ## helper function for cniper-type plots to plot in data
   data, # data to be plot in
   dots, # dots from the calling function
   fun, # function to determine risk difference
   L2Fam, # L2Family
   IC, # IC1 in cniperContPlot and eta in cniperPointPlot
   jit.fac,
   jit.tol
){
               dotsP <- .makedotsP(dots)

               al <- dotsP$alpha.trsp
               if(!is.null(al)) if(!is.na(al))
                   dotsP$col <- sapply(dotsP$col,
                                            addAlphTrsp2col, alpha=al)

               n <- if(!is.null(dim(data))) nrow(data) else length(data)
               if(!is.null(dots$lab.pts))
                    lab.pts <- rep(lab.pts, length.out=n)


               sel <- .SelectOrderData(data, function(x)sapply(x,fun),
                                       dots$which.lbs, dots$which.Order)
               i.d <- sel$ind
               i0.d <- sel$ind1
               y.d <- sel$y
               x.d <- sel$data
               n <- length(i.d)

               resc.dat <- .rescalefct(x.d, function(x) sapply(x,fun),
                              dots$scaleX, dots$scaleX.fct, dots$scaleX.inv,
                              dots$scaleY, dots$scaleY.fct,
                              dots$xlim, dots$ylim, dots)

               if(any(.isReplicated(resc.dat$X, jit.tol))&&jit.fac>0)
                       resc.dat$X <- jitter(resc.dat$X, factor = jit.fac)
               if(any(.isReplicated(resc.dat$Y, jit.tol))&&jit.fac>0)
                       resc.dat$Y <- jitter(resc.dat$Y, factor = jit.fac)

               dotsP$scaleX <- dotsP$scaleY <- dotsP$scaleN <-NULL
               dotsP$scaleX.fct <- dotsP$scaleY.fct <- NULL
               dotsP$scaleX.inv <- dotsP$scaleY.inv <- NULL
               dotsP$cex.pts <- dotsP$col.pts <- dotsP$lab.pts <- dotsP$pch.pts <- NULL
               dotsP$jit.fac <- dotsP$with.lab <- dotsP$alpha.trsp <- NULL
               dotsP$return.Order <- dotsP$cex.pts.fun <- NULL
               dotsP$x.ticks <- dotsP$y.ticks <- NULL
               dotsP$lab.font <- dotsP$which.lbs <- dotsP$which.lbs <- NULL

               dotsP$x <- resc.dat$X
               dotsP$y <- resc.dat$Y

               trafo <- trafo(L2Fam@param)
               dims <- nrow(trafo)
               QF <- diag(dims)
               if(is(IC,"ContIC") & dims>1 )
                      {if (is(normtype(IC),"QFNorm"))
                           QF <- QuadForm(normtype(IC))}

               absInfoEval <- function(x,y) sapply(x, y@Map[[1]])
               IC.rv <- as(diag(dims) %*% IC@Curve, "EuclRandVariable")
               absy.f <- t(IC.rv) %*% QF %*% IC.rv
               absy <- absInfoEval(x.d, absy.f)

               if(is.null(dots$cex.pts)) dots$cex.pts <- par("cex")

               dotsT <- dotsP
               dotsT$cex <- dotsP$cex/2
               dotsP$cex <- .cexscale(absy,absy,cex=dots$cex.pts, fun = dots$cex.pts.fun)
               dotsP$col <- dots$col.pts

               dotsT$pch <- NULL
               dotsT$labels <- if(is.null(dots$lab.pts)) i.d else dots$lab.pts[i.d]
               do.call(points,dotsP)
               if(!is.null(dots$with.lab))
                   if(dots$with.lab)  do.call(text,dotsT)
               if(!is.null(dots$return.Order))
                   if(dots$return.Order) return(i0.d)
        return(invisible(NULL))
        }


.getFunCnip <- function(IC1,IC2, risk, L2Fam, r, b20=NULL){

        bType <- biastype(risk)
        nType <- normtype(risk)
        fnorm <- fct(nType)
        riskfct <- getRiskFctBV(risk, bType)

       .getTrVar <- function(IC){
           R <- Risks(IC)[["trAsCov"]]
           if(is.null(R)) R <- getRiskIC(IC, risk = trAsCov(), L2Fam = L2Fam)
           if("trAsCov" %in% names(R)) R <- R[["trAsCov"]]
           if(length(R) > 1) R <- R$value
           return(c(R))
        }
        R1 <- .getTrVar (IC1)
        R2 <- .getTrVar (IC2)

        fun <- function(x){
            y1 <- sapply(x, function(x1)evalIC(IC1,as.matrix(x1,ncol=1)))
            b1 <- r*fnorm(y1)
            r1 <- riskfct(var=R1,bias=b1)
            if(!is.null(b20)){
               r2 <- riskfct(var=R2, bias=b20)
            }else{
               y2 <- sapply(x,function(x0) evalIC(IC2,x0))
               b2 <- r*fnorm(y2)
               r2 <- riskfct(var=R2,bias=b2)
            }
            r1 - r2
            return(r1-r2)
        }
        return(fun)
}

cniperCont <- function(IC1, IC2, data = NULL, ...,
                           neighbor, risk, lower=getdistrOption("DistrResolution"),
                           upper=1-getdistrOption("DistrResolution"), n = 101,
                           scaleX = FALSE, scaleX.fct, scaleX.inv,
                           scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
                           scaleN = 9, x.ticks = NULL, y.ticks = NULL,
                           cex.pts = 1, cex.pts.fun = NULL, col.pts = par("col"),
                           pch.pts = 19, jit.fac = 1, jit.tol = .Machine$double.eps, with.lab = FALSE,
                           lab.pts = NULL, lab.font = NULL, alpha.trsp = NA,
                           which.lbs = NULL, which.Order  = NULL,
                           return.Order = FALSE, 
             draw.nonlbl = TRUE,  ## should non-labelled observations also be drawn?
             cex.nonlbl = 0.3,    ## character expansion(s) for non-labelled observations
             pch.nonlbl = ".",    ## plotting symbol(s) for non-labelled observations
                           withSubst = TRUE){

        mcD <- match.call(expand.dots = FALSE)
        dots <- as.list(mcD$"...")
        mc <- match.call(#call = sys.call(sys.parent(1)),
                       expand.dots = TRUE)
        mcl <- as.list(mc[-1])
        IC1c <- as.character(deparse(IC1))
        IC2c <- as.character(deparse(IC2))

       .mpresubs <- if(withSubst){
                     function(inx) 
                      .presubs(inx, c("%C1", "%A1", "%C2", "%A2", "%D" ),
                          c(as.character(class(IC1)[1]), 
                            IC1c,
                            as.character(class(IC2)[1]), 
                            IC2c,
                            as.character(date())
                            ))
                     }else function(inx)inx

        if(!is.null(dots$main)) dots$main <- .mpresubs(dots$main)
        if(!is.null(dots$sub)) dots$sub <- .mpresubs(dots$sub)
        if(!is.null(dots$xlab)) dots$xlab <- .mpresubs(dots$xlab)
        if(!is.null(dots$ylab)) dots$ylab <- .mpresubs(dots$ylab)

        if(!is(IC1,"IC")) stop ("IC1 must be of class 'IC'")
        if(!is(IC2,"IC")) stop ("IC2 must be of class 'IC'")
        if(!identical(IC1@CallL2Fam, IC2@CallL2Fam))
        stop("IC1 and IC2 must be defined on the same model")

        L2Fam <- eval(IC1@CallL2Fam)

        b20 <- NULL
        fCpl <- eval(dots$fromCniperPlot)
        if(!is.null(fCpl)&&length(Risks(IC2)))
            if(fCpl) b20 <- neighbor@radius*Risks(IC2)$asBias$value
        dots$fromCniperPlot <- NULL
        
        fun <- .getFunCnip(IC1,IC2, risk, L2Fam, neighbor@radius, b20)

        if(missing(scaleX.fct)){
           scaleX.fct <- p(L2Fam)
           scaleX.inv <- q(L2Fam)
        }

        if("lower" %in% names(as.list(mc))) lower <- p(L2Fam)(lower)
        if("upper" %in% names(as.list(mc))) upper <- p(L2Fam)(upper)

        x <-  q(L2Fam)(seq(lower,upper,length=n))
        if(is(distribution(L2Fam), "DiscreteDistribution"))
           x <- seq(q(L2Fam)(lower),q(L2Fam)(upper),length=n)

        resc <- .rescalefct(x, fun, scaleX, scaleX.fct,
                     scaleX.inv, scaleY, scaleY.fct, dots$xlim, dots$ylim, dots)

        dotsPl <- dots
        dotsPl$x <- resc$X
        dotsPl$y <- resc$Y
        dotsPl$type <- "l"
        if(is.null(dotsPl$main)) dotsPl$main <- gettext("Cniper region plot")
        if(is.null(dotsPl$xlab)) dotsPl$xlab <- gettext("Dirac point")
        if(is.null(dotsPl$ylab))
           dotsPl$ylab <- gettext("Asymptotic Risk difference (IC1 - IC2)")

        colSet <- ltySet <- lwdSet <- FALSE
        if(!is.null(dotsPl$col)) {colSet <- TRUE; colo <- eval(dotsPl$col)}
        if(colSet) {
           colo <- rep(colo,length.out=2)
           dotsPl$col <- colo[1]
        }
        if(!is.null(dotsPl$lwd)) {lwdSet <- TRUE; lwdo <- eval(dotsPl$lwd)}
        if(lwdSet) {
           lwdo <- rep(lwdo,length.out=2)
           dotsPl$lwd <- lwdo[1]
        }
        if(!is.null(dotsPl$lty)) {ltySet <- TRUE; ltyo <- eval(dotsPl$lty)}
        if(ltySet && ((!is.numeric(ltyo) && length(ltyo)==1)||
                        is.numeric(ltyo))){
           ltyo <- list(ltyo,ltyo)
           dotsPl$lty <- ltyo[[1]]
        }else{ if (ltySet && !is.numeric(ltyo) && length(ltyo)==2){
                   dotsPl$lty <- ltyo[[1]]
            }
        }
        do.call(plot,dotsPl)

        dots$x <- dots$y <- NULL
        dotsl <- .makedotsLowLevel(dots)
        if(colSet) dotsl$col <- colo[2]
        if(lwdSet) dotsl$lwd <- lwdo[2]
        if(ltySet) dotsl$lty <- ltyo[[2]]

        dotsl$h <- if(scaleY) scaleY.fct(0) else 0
        do.call(abline, dotsl)

        .plotRescaledAxis(scaleX, scaleX.fct, scaleX.inv, scaleY,scaleY.fct,
                          scaleY.inv, dots$xlim, dots$ylim, resc$X, ypts = 400,
                          n = scaleN, x.ticks = x.ticks, y.ticks = y.ticks)
        if(!is.null(data)){
           dots$scaleX <- scaleX
           dots$scaleX.fct <-  scaleX.fct
           if(!is.null(mcl$scaleX.inv)) dots$scaleX.inv <-  scaleX.inv
           dots$scaleY <- scaleY
           dots$scaleY.fct <- scaleY.fct
           dots$scaleY.inv <- scaleY.inv
           dots$scaleN <- scaleN
           dots$x.ticks <- x.ticks
           dots$y.ticks <- y.ticks
           dots$cex.pts <- cex.pts
           dots$cex.pts.fun <- cex.pts.fun
           dots$col.pts <- col.pts
           dots$pch.pts <- pch.pts
           dots$jit.fac <- jit.fac
           dots$jit.tol <- jit.tol
           dots$with.lab <- with.lab
           dots$lab.pts <- lab.pts
           dots$lab.font <- lab.font
           dots$alpha.trsp <- alpha.trsp
           dots$which.lbs <- which.lbs
           dots$which.Order  <- which.Order
           dots$return.Order <- return.Order

           return(.plotData(data=data, dots=dots, fun=fun, L2Fam=L2Fam, IC=IC1, jit.fac=jit.fac, jit.tol=jit.tol))
        }
        invisible(NULL)
}

cniperPoint <- function(L2Fam, neighbor, risk= asMSE(),
                        lower=getdistrOption("DistrResolution"),
                        upper=1-getdistrOption("DistrResolution")){


        mc <- match.call(expand.dots = FALSE)

        if(is.null(as.list(mc)$lower)) lower <- q(L2Fam)(lower)
        if(is.null(as.list(mc)$upper)) upper <- q(L2Fam)(upper)
#        lower <- q(L2Fam)(lower)
#        upper <- q(L2Fam)(upper)

        robMod <- InfRobModel(center = L2Fam, neighbor = neighbor)

        psi <- optIC(model = L2Fam, risk = asCov())
        eta <- optIC(model = robMod, risk = risk)

        fun <- .getFunCnip(psi,eta, risk, L2Fam, neighbor@radius)
        res <- uniroot(fun, lower = lower, upper = upper)$root
        names(res) <- "cniper point"
        res
}

cniperPointPlot <- function(L2Fam, data=NULL, ..., neighbor, risk= asMSE(),
                        lower=getdistrOption("DistrResolution"),
                        upper=1-getdistrOption("DistrResolution"), n = 101,
                        withMaxRisk = TRUE,
                           scaleX = FALSE, scaleX.fct, scaleX.inv,
                           scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
                           scaleN = 9, x.ticks = NULL, y.ticks = NULL,
                           cex.pts = 1, cex.pts.fun = NULL, col.pts = par("col"),
                           pch.pts = 19, jit.fac = 1, jit.tol = .Machine$double.eps, 
                           with.lab = FALSE,
                           lab.pts = NULL, lab.font = NULL, alpha.trsp = NA,
                           which.lbs = NULL, which.Order  = NULL,
                           return.Order = FALSE, 
             draw.nonlbl = TRUE,  ## should non-labelled observations also be drawn?
             cex.nonlbl = 0.3,    ## character expansion(s) for non-labelled observations
             pch.nonlbl = ".",    ## plotting symbol(s) for non-labelled observations
                           withSubst = TRUE){

        mc0 <- match.call(#call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)
        mc <- match.call(#call = sys.call(sys.parent(1)),
                       expand.dots = TRUE)
        mcl <- as.list(mc[-1])
        dots <- as.list(mc0$"...")
        L2Famc <- as.character(deparse(L2Fam))

       .mpresubs <- if(withSubst){
                     function(inx) 
                      .presubs(inx, c("%C", "%A", "%D" ),
                          c(as.character(class(L2Fam)[1]), 
                            L2Famc,
                            as.character(date())
                            ))
                     }else function(inx)inx

        if(!is.null(dots$main)) dots$main <- .mpresubs(dots$main)
        if(!is.null(dots$sub)) dots$sub <- .mpresubs(dots$sub)
        if(!is.null(dots$xlab)) dots$xlab <- .mpresubs(dots$xlab)
        if(!is.null(dots$ylab)) dots$ylab <- .mpresubs(dots$ylab)
        if(is.null(mcl$risk)) mcl$risk <- asMSE()

        robMod <- InfRobModel(center = L2Fam, neighbor = neighbor)

        mcl$IC1 <- optIC(model = L2Fam, risk = asCov())
        mcl$IC2 <- if(is(risk,"interpolRisk")){
                     getStartIC(model=L2Fam, risk = risk)
                   }else optIC(model = robMod, risk = risk)
        mcl$L2Fam <- NULL
        if(is.null(dots$ylab))
           mcl$ylab <- gettext("Asymptotic Risk difference (classic - robust)")
        if(is.null(dots$main))
           mcl$main <- gettext("Cniper point plot")

        if(withMaxRisk) mcl$fromCniperPlot <- TRUE
        mcl$withMaxRisk <- NULL
        mcl$withSubst <- FALSE
        do.call(cniperCont, mcl)
}



 .cexscale <- function(y, y1=y, maxcex=4,mincex=0.05,cex, fun=NULL){
         if(is.null(fun)) fun <- function(x) log(1+abs(x))
         ly <- fun(y)
         ly1 <- fun(unique(c(y,y1)))
         my <- min(ly1,na.rm=TRUE)
         My <- max(ly1,na.rm=TRUE)
         ly0 <- (ly-my)/My
         ly1 <- ly0*(maxcex-mincex)+mincex
         return(cex*ly1)
 }
