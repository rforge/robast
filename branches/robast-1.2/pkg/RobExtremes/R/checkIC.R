..checkIC.qtl <- function (IC, L2Fam, out = TRUE, ...){
        mc <- match.call(expand.dots = TRUE)
        mcF <- match.call(expand.dots = FALSE)
        dots  <- mcF$"..."
        mcl <- as.list(mc)[-1]
        mcl$out <- out
        mcl$IC <- IC
        D1 <- L2Fam@distribution
        D1 <- as(D1,"DistributionsIntegratingByQuantiles")
        L2Fam@distribution <- D1
        L2Fam <- as(L2Fam, "L2ParamFamily")
        mcl$L2Fam <- L2Fam
        if(!is.null(dots)) mcl <- c(mcl,dots)
        do.call("checkIC", mcl)
        }
setMethod("checkIC", signature(IC = "IC", L2Fam = "GParetoFamily"),
    ..checkIC.qtl)

setMethod("checkIC", signature(IC = "IC", L2Fam = "GEVFamilyMuUnknown"),
    ..checkIC.qtl)
setMethod("checkIC", signature(IC = "IC", L2Fam = "GEVFamily"),
    ..checkIC.qtl)
setMethod("checkIC", signature(IC = "IC", L2Fam = "ParetoFamily"),
    ..checkIC.qtl)
