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
        tr.invF <- sum(diag(solve(FisherInfo(L2Fam))))
        psi <- optIC(model = L2Fam, risk = asCov())
        robMod <- InfRobModel(center = L2Fam, neighbor = neighbor)
        eta <- optIC(model = robMod, risk = asMSE())
        maxMSE <- Risks(eta)$asMSE
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
    function(L2Fam, neighbor, risk, lower, upper, n = 101){
        tr.invF <- sum(diag(solve(FisherInfo(L2Fam))))
        psi <- optIC(model = L2Fam, risk = asCov())
        robMod <- InfRobModel(center = L2Fam, neighbor = neighbor)
        eta <- optIC(model = robMod, risk = asMSE())
        maxMSE <- Risks(eta)$asMSE
        fun <- function(x){
            y <- evalIC(psi, x) 
            tr.invF + as.vector(y %*% y)*neighbor@radius^2 - maxMSE
        }
        x <- seq(from = lower, to = upper, length = n)
        y <- sapply(x, fun)
        plot(x, y, type = "l", main = "Cniper point plot", 
             xlab = "Dirac point", ylab = "Asymptotic MSE difference (classic - robust)")
#        text(min(x), max(y)/2, "Robust", pos = 4)
#        text(min(x), min(y)/2, "Classic", pos = 4)
        abline(h = 0)
        invisible()
    })
