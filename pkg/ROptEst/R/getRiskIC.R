###############################################################################
## asymptotic covariance
###############################################################################
setMethod("getRiskIC", signature(IC = "HampIC", 
                                 risk = "asCov",
                                 neighbor = "missing",
                                 L2Fam = "missing"),
    function(IC, risk){
        L2Fam <- eval(IC@CallL2Fam)
        getRiskIC(IC, risk, L2Fam)
    })

setMethod("getRiskIC", signature(IC = "HampIC", 
                                 risk = "asCov",
                                 neighbor = "missing",
                                 L2Fam = "missing"),
    function(IC, risk, L2Fam){
        Cov <- IC@Risks[["asCov"]]        
        return(list(asCov = list(distribution = .getDistr(IC@L2Fam), value = Cov)))
    })


###############################################################################
## asymptotic Bias for various types
###############################################################################
setMethod("getBiasIC", signature(IC = "HampIC",
                                 neighbor = "UncondNeighborhood"),
    function(IC, neighbor, L2Fam){
        if(missing(L2Fam)) 
           {misF <- TRUE; L2Fam <- eval(IC@CallL2Fam)}

        return(list(asBias = list(distribution = .getDistr(L2Fam), 
                    neighborhood = neighbor@type, value = IC@Risks[["asBias"]])))
    })


