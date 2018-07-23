..makeIC.qtl <- function (IC, L2Fam){
        mc <- match.call()
        mcl <- as.list(mc)[-1]
        mcl$IC <- IC
        D1 <- L2Fam@distribution
        D1 <- as(D1,"DistributionsIntegratingByQuantiles")
        L2Fam@distribution <- D1
        L2Fam <- as(L2Fam, "L2ParamFamily")
        mcl$L2Fam <- L2Fam
        do.call("makeIC", mcl)
        }
setMethod("makeIC", signature(IC = "IC", L2Fam = "GParetoFamily"),
    ..makeIC.qtl)

setMethod("makeIC", signature(IC = "IC", L2Fam = "GEVFamilyMuUnknown"),
    ..makeIC.qtl)
setMethod("makeIC", signature(IC = "IC", L2Fam = "GEVFamily"),
    ..makeIC.qtl)
setMethod("makeIC", signature(IC = "IC", L2Fam = "ParetoFamily"),
    ..makeIC.qtl)
