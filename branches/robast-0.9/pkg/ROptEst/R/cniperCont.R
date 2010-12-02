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
    function(L2Fam, neighbor, risk, lower, upper, n = 101, ...){
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
#
#cniperPointPlot(L2Fam=N0, neighbor=ContNeighborhood(radius = 0.5), risk=asMSE(),lower=-12, n =30, upper=8, lwd=c(2,4),lty=list(c(5,1),3),col=c(2,4))
 