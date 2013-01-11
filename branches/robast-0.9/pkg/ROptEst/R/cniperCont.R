PETER <- FALSE
setMethod("cniperCont", signature(IC1 = "IC",
                                  IC2 = "IC",
                                  L2Fam = "L2ParamFamily",
                                  neighbor = "ContNeighborhood",
                                  risk = "asMSE"),
    function(IC1, IC2, L2Fam, neighbor, risk, lower, upper, n = 101){
        R1 <- Risks(IC1)[["trAsCov"]]
        if(is.null(R1)) R1 <- getRiskIC(IC1, risk = trAsCov(), L2Fam = L2Fam)
        if(length(R1) > 1) R1 <- R1$value

        R2 <- Risks(IC2)[["trAsCov"]]
        if(is.null(R2)) R2 <- getRiskIC(IC2, risk = trAsCov(), L2Fam = L2Fam)
        if(length(R2) > 1) R2 <- R2$value

        r <- neighbor@radius

        fun <- function(x){
            y1 <- evalIC(IC1, x) 
            y2 <- evalIC(IC2, x)
            R1 - R2 + r^2*(as.vector(y1 %*% y1) - as.vector(y2 %*% y2))
        }
        x <- seq(from = lower, to = upper, length = n)
        y <- sapply(x, fun)
        plot(x, y, type = "l", main = "Cniper region plot", 
             xlab = "Dirac point", ylab = "Asymptotic MSE difference (IC1 - IC2)")
#        text(min(x), max(y)/2, "IC2", pos = 4)
#        text(min(x), min(y)/2, "IC1", pos = 4)
        abline(h = 0)
        invisible()
    })

setMethod("cniperPoint", signature(L2Fam = "L2ParamFamily",
                                   neighbor = "ContNeighborhood",
                                   risk = "asMSE"),
    function(L2Fam, neighbor, risk, lower, upper){
        D <- trafo(L2Fam@param)
        tr.invF <- sum(diag(D %*% solve(FisherInfo(L2Fam)) %*% t(D)))
        psi <- optIC(model = L2Fam, risk = asCov())
        robMod <- InfRobModel(center = L2Fam, neighbor = neighbor)
        eta <- optIC(model = robMod, risk = asMSE())
        maxMSE <- Risks(eta)$asMSE$value
        Delta <- sqrt(maxMSE - tr.invF)/neighbor@radius
        fun <- function(x){
            y <- evalIC(psi, x)
            sqrt(as.vector(y %*% y)) - Delta
        }
        res <- uniroot(fun, lower = lower, upper = upper)$root
        names(res) <- "cniper point"
        res
    })

setMethod("cniperPointPlot", signature(L2Fam = "L2ParamFamily",
                                   neighbor = "ContNeighborhood",
                                   risk = "asMSE"),
    function(L2Fam, neighbor, risk, lower, upper, n = 101, ...,
             scaleX = FALSE, scaleX.fct, scaleX.inv,
             scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
    ){
        dots <- as.list(match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)$"...")
        D <- trafo(L2Fam@param)
        tr.invF <- sum(diag(D %*% solve(FisherInfo(L2Fam)) %*% t(D)))
        psi <- optIC(model = L2Fam, risk = asCov())
        robMod <- InfRobModel(center = L2Fam, neighbor = neighbor)
        eta <- optIC(model = robMod, risk = asMSE())
        maxMSE <- Risks(eta)$asMSE$value
        fun <- function(x){
            y <- evalIC(psi, x)
            tr.invF + as.vector(y %*% y)*neighbor@radius^2 - maxMSE
        }
        dots$x <- x <- seq(from = lower, to = upper, length = n)
        dots$y <- sapply(x, fun)
        colSet <- ltySet <- lwdSet <- FALSE
        if(!is.null(dots$col)) {colSet <- TRUE; colo <- eval(dots$col)}
        if(colSet) {
           colo <- rep(colo,length.out=2)
           dots$col <- colo[1]
        }
        if(!is.null(dots$lwd)) {lwdSet <- TRUE; lwdo <- eval(dots$lwd)}
        if(lwdSet) {
           lwdo <- rep(lwdo,length.out=2)
           dots$lwd <- lwdo[1]
        }
        if(!is.null(dots$lty)) {ltySet <- TRUE; ltyo <- eval(dots$lty)}
        if(ltySet && ((!is.numeric(ltyo) && length(ltyo)==1)||
                        is.numeric(ltyo))){
           ltyo <- list(ltyo,ltyo)
           dots$lty <- ltyo[[1]]
        }else{ if (ltySet && !is.numeric(ltyo) && length(ltyo)==2){
                   dots$lty <- ltyo[[1]]
            }
        }
        if(is.null(dots$main)) dots$main <- gettext("Cniper point plot")
        if(is.null(dots$xlab)) dots$xlab <- gettext("Dirac point")
        if(is.null(dots$ylab))
            dots$ylab <- gettext("Asymptotic MSE difference (classic - robust)")
        dots$type <- "l"




        do.call(plot, dots)
#        text(min(x), max(y)/2, "Robust", pos = 4)
#        text(min(x), min(y)/2, "Classic", pos = 4)
        dots$x <- dots$y <- dots$xlab <- dots$ylab <- dots$main <- dots$type <- NULL
        dots$h <- 0
        if(colSet) dots$col <- colo[2]
        if(lwdSet) dots$lwd <- lwdo[2]
        if(ltySet) dots$lty <- ltyo[[2]]
        do.call(abline, dots)
        invisible()
    })

if(PETER){
.rescalefct <- RobAStBase:::.rescalefct
.makedotsP <- RobAStBase:::.makedotsP
.SelectOrderData <- RobAStBase:::.SelectOrderData

.plotData <- function(
  ## helper function for cniper-type plots to plot in data
   data, # data to be plot in
   dots, # dots from the calling function
   origCl, # call from the calling function
   fun, # function to determine risk difference
   L2Fam, # L2Family
   IC # IC1 in cniperContPlot and eta in cniperPointPlot
){
               dotsP <- .makedotsP(dots)
               dotsP$col <- rep(origCl$col.pts, length.out=n)
               dotsP$pch <- rep(origCl$pch.pts, length.out=n)

               n <- if(!is.null(dim(data))) nrow(data) else length(data)
               if(!is.null(lab.pts))
                    lab.pts <- rep(origCl$lab.pts, length.out=n)

               sel <- .SelectOrderData(data, function(x)sapply(x,fun),
                                       origCl$which.lbs, origCl$which.Order)
               i.d <- sel$ind
               i0.d <- sel$ind1
               y.d <- sel$y
               x.d <- sel$data
               x.d <- sel.C$data
               n <- length(i.d)

               resc.dat <- .rescalefct(x.d, function(x) sapply(x,fun),
                              origCl$scaleX, origCl$scaleX.fct, origCl$scaleX.inv,
                              origCl$scaleY, origCl$scaleY.fct,
                              dots$xlim, dots$ylim, dots)

               dotsP$x <- resc.dat$X
               dotsP$y <- resc.dat$Y

               trafo <- trafo(L2Fam@param)
               dims <- nrow(trafo)
               QF <- diag(dims)
               if(is(IC1,"ContIC") & dims>1 )
                      {if (is(normtype(object),"QFNorm"))
                           QF <- QuadForm(normtype(object))}

               absInfoEval <- function(x,y) sapply(x, y@Map[[1]])
               IC1.rv <- as(diag(dims) %*% IC@Curve, "EuclRandVariable")
               absy.f <- t(IC1.rv) %*% QF %*% IC1.rv
               absy <- absInfoEval(x.d, absy.f)

               dotsP$cex <-  log(absy+1)*3*rep(cex.pts, length.out=n)

               dotsT <- dotsP
               dotsT$pch <- NULL
               dotsT$cex <- dotsP$cex/2
               dotsT$labels <- if(is.null(lab.pts)) i.d else lab.pts[i.d]
               do.call(points,dotsP)
               if(with.lab)  do.call(text,dotsT)
               if(return.Order) return(i0.d)
        return(invisible(NULL))
        }

cniperContPlot <- function(IC1, IC2, data = NULL, ...,
                           neighbor, risk, lower=getdistrOption("DistrResolution"),
                           upper=1-getdistrOption("DistrResolution"), n = 101,
                           scaleX = FALSE, scaleX.fct, scaleX.inv,
                           scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
                           cex.pts = 1, col.pts = par("col"),
                           pch.pts = 1, jitter.fac = 1, with.lab = FALSE,
                           lab.pts = NULL, lab.font = NULL,
                           which.lbs = NULL, which.Order  = NULL,
                           return.Order = FALSE){

        mc <- match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)
        dots <- mc$"..."

        if(!is(IC1,"IC")) stop ("IC1 must be of class 'IC'")
        if(!is(IC2,"IC")) stop ("IC2 must be of class 'IC'")
        if(!identical(IC1@CallL2Fam, IC2@CallL2Fam))
        stop("IC1 and IC2 must be defined on the same model")

        L2Fam <- eval(IC1@CallL2Fam)
        riskfct <- getRiskBV(risk, biastype(risk))


        if(missing(scaleX.fct)){
           scaleX.fct <- p(L2Fam)
           scaleX.inv <- q(L2Fam)
        }

        R1 <- Risks(IC1)[["trAsCov"]]
        if(is.null(R1)) R1 <- getRiskIC(IC1, risk = trAsCov(), L2Fam = L2Fam)
        if(length(R1) > 1) R1 <- R1$value

        R2 <- Risks(IC2)[["trAsCov"]]
        if(is.null(R2)) R2 <- getRiskIC(IC2, risk = trAsCov(), L2Fam = L2Fam)
        if(length(R2) > 1) R2 <- R2$value

        r <- neighbor@radius

        fun <- function(x){
            y1 <- evalIC(IC1, x)
            y2 <- evalIC(IC2, x)
            riskfct(R1,r*fct(normtype(risk))(y1))-
            riskfct(R2,r*fct(normtype(risk))(y2))
        }
        x <-  q(L2Fam)(seq(lower,upper,length=n))
        resc <- .rescalefct(x, function(u) sapply(u,fun), scaleX, scaleX.fct,
                     scaleX.inv, scaleY, scaleY.fct, dots$xlim, dots$ylim, dots)
        x <- dots$x <- resc$X
        dots$y <- resc$Y
        dots$type <- "l"
        if(is.null(dots$main)) dots$main <- "Cniper region plot"
        if(is.null(dots$xlab)) dots$xlab <- "Dirac point"
        if(is.null(dots$ylab)) dots$ylab <- "Asymptotic Risk difference (IC1 - IC2)"

        colSet <- ltySet <- lwdSet <- FALSE
        if(!is.null(dots$col)) {colSet <- TRUE; colo <- eval(dots$col)}
        if(colSet) {
           colo <- rep(colo,length.out=2)
           dots$col <- colo[1]
        }
        if(!is.null(dots$lwd)) {lwdSet <- TRUE; lwdo <- eval(dots$lwd)}
        if(lwdSet) {
           lwdo <- rep(lwdo,length.out=2)
           dots$lwd <- lwdo[1]
        }
        if(!is.null(dots$lty)) {ltySet <- TRUE; ltyo <- eval(dots$lty)}
        if(ltySet && ((!is.numeric(ltyo) && length(ltyo)==1)||
                        is.numeric(ltyo))){
           ltyo <- list(ltyo,ltyo)
           dots$lty <- ltyo[[1]]
        }else{ if (ltySet && !is.numeric(ltyo) && length(ltyo)==2){
                   dots$lty <- ltyo[[1]]
            }
        }

        do.call(plot,dots)

        dots <- makedotsLowLevel(dots)
        dots$x <- dots$y <- NULL
        if(colSet) dots$col <- colo[2]
        if(lwdSet) dots$lwd <- lwdo[2]
        if(ltySet) dots$lty <- ltyo[[2]]

        dots$h <- if(scaleY) scaleY.fct(0) else 0
        do.call(abline, dots)

        .plotRescaledAxis(scaleX, scaleX.fct, scaleX.inv, scaleY,scaleY.fct,
                          scaleY.inv, dots$xlim, dots$ylim, x, ypts = 400)
        if(!is.null(data))
           return(.plotData(data, dots, mc, fun, L2Fam, IC1))
        invisible(NULL)
}

cniperPoint <- function(L2Fam, neighbor, risk= asMSE(),
                        lower=getdistrOption("DistrResolution"),
                        upper=1-getdistrOption("DistrResolution")){

        lower.x <- q(L2Fam)(lower)
        upper.x <- q(L2Fam)(upper)
        riskfct <- getRiskBV(risk, biastype(risk))

        psi <- optIC(model = L2Fam, risk = asCov())
        R1 <- Risks(psi)[["trAsCov"]]
        if(is.null(R1)) R1 <- getRiskIC(psi, risk = trAsCov(), L2Fam = L2Fam)
        if(length(R1) > 1) R1 <- R1$value

        robMod <- InfRobModel(center = L2Fam, neighbor = neighbor)
        eta <- optIC(model = robMod, risk = risk)
        R2 <- Risks(eta)[["trAsCov"]]
        if(is.null(R2)) R2 <- getRiskIC(eta, risk = trAsCov(), L2Fam = L2Fam)
        if(length(R2) > 1) R2 <- R2$value

        r <- neighbor@radius

        fun <- function(x){
            y1 <- evalIC(psi, x)
            y2 <- evalIC(eta, x)
            riskfct(R1,r*fct(normtype(risk))(y1))-
            riskfct(R2,r*fct(normtype(risk))(y2))
        }

        res <- uniroot(fun, lower = lower, upper = upper)$root
        names(res) <- "cniper point"
        res
}

cniperPointPlot <- function(L2Fam, data=NULL, ..., neighbor, risk= asMSE(),
                        lower=getdistrOption("DistrResolution"),
                        upper=1-getdistrOption("DistrResolution"), n = 101,
                           scaleX = FALSE, scaleX.fct, scaleX.inv,
                           scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
                           cex.pts = 1, col.pts = par("col"),
                           pch.pts = 1, jitter.fac = 1, with.lab = FALSE,
                           lab.pts = NULL, lab.font = NULL,
                           which.lbs = NULL, which.Order  = NULL,
                           return.Order = FALSE){

        mc <- match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)
        dots <- mc$"..."

        if(missing(scaleX.fct)){
           scaleX.fct <- p(L2Fam)
           scaleX.inv <- q(L2Fam)
        }

        lower.x <- q(L2Fam)(lower)
        upper.x <- q(L2Fam)(upper)
        riskfct <- getRiskBV(risk, biastype(risk))

        psi <- optIC(model = L2Fam, risk = asCov())
        R1 <- Risks(psi)[["trAsCov"]]
        if(is.null(R1)) R1 <- getRiskIC(psi, risk = trAsCov(), L2Fam = L2Fam)
        if(length(R1) > 1) R1 <- R1$value

        robMod <- InfRobModel(center = L2Fam, neighbor = neighbor)
        eta <- optIC(model = robMod, risk = risk)
        R2 <- Risks(eta)[["trAsCov"]]
        if(is.null(R2)) R2 <- getRiskIC(eta, risk = trAsCov(), L2Fam = L2Fam)
        if(length(R2) > 1) R2 <- R2$value

        r <- neighbor@radius

        fun <- function(x){
            y1 <- evalIC(psi, x)
            y2 <- evalIC(eta, x)
            riskfct(R1,r*fct(normtype(risk))(y1))-
            riskfct(R2,r*fct(normtype(risk))(y2))
        }


        x <- q(L2Fam)(seq(lower,upper,length=n))
        resc <- .rescalefct(x, function(u) sapply(u,fun), scaleX, scaleX.fct,
                     scaleX.inv, scaleY, scaleY.fct, dots$xlim, dots$ylim, dots)
        x <- dots$x <- resc$X
        dots$y <- resc$Y

        colSet <- ltySet <- lwdSet <- FALSE
        if(!is.null(dots$col)) {colSet <- TRUE; colo <- eval(dots$col)}
        if(colSet) {
           colo <- rep(colo,length.out=2)
           dots$col <- colo[1]
        }
        if(!is.null(dots$lwd)) {lwdSet <- TRUE; lwdo <- eval(dots$lwd)}
        if(lwdSet) {
           lwdo <- rep(lwdo,length.out=2)
           dots$lwd <- lwdo[1]
        }
        if(!is.null(dots$lty)) {ltySet <- TRUE; ltyo <- eval(dots$lty)}
        if(ltySet && ((!is.numeric(ltyo) && length(ltyo)==1)||
                        is.numeric(ltyo))){
           ltyo <- list(ltyo,ltyo)
           dots$lty <- ltyo[[1]]
        }else{ if (ltySet && !is.numeric(ltyo) && length(ltyo)==2){
                   dots$lty <- ltyo[[1]]
            }
        }
        if(is.null(dots$main)) dots$main <- gettext("Cniper point plot")
        if(is.null(dots$xlab)) dots$xlab <- gettext("Dirac point")
        if(is.null(dots$ylab))
            dots$ylab <- gettext("Asymptotic Risk difference (classic - robust)")
        dots$type <- "l"
        do.call(plot, dots)

        dots <- makedotsLowLevel(dots)
        dots$x <- dots$y <- NULL
        if(colSet) dots$col <- colo[2]
        if(lwdSet) dots$lwd <- lwdo[2]
        if(ltySet) dots$lty <- ltyo[[2]]

        dots$h <- if(scaleY) scaleY.fct(0) else 0
        do.call(abline, dots)
        .plotRescaledAxis(scaleX, scaleX.fct, scaleX.inv, scaleY,scaleY.fct,
                          scaleY.inv, dots$xlim, dots$ylim, x, ypts = 400)
        if(!is.null(data))
           return(.plotData(data, dots, mc, fun, L2Fam, eta))
        return(invisible(NULL))
}

}

