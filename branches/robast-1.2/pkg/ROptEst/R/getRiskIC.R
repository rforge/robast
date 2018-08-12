.checkICWithWarning <- function(IC, L2Fam, tol, ...){
          if(!missing(L2Fam)){
             prec <- checkIC(IC, L2Fam, out = FALSE, ...)
          }else{
             prec <- checkIC(IC, out = FALSE, ...)
          }
          if(prec > tol)
            warning("The maximum deviation from the exact IC properties is ", prec,
                    "\nThis is larger than the specified 'tol' ",
                    "=> the result may be wrong")
}

###############################################################################
## asymptotic covariance
###############################################################################
setMethod("getRiskIC", signature(IC = "HampIC", 
                                 risk = "asCov",
                                 neighbor = "missing",
                                 L2Fam = "missing"),
    function(IC, risk, withCheck= TRUE, ...){
        L2Fam <- force(eval(IC@CallL2Fam))
        if(missing(withCheck)) withCheck <- TRUE
        getRiskIC(IC = IC, risk = risk, L2Fam = L2Fam, withCheck= withCheck, ...)
    })

setMethod("getRiskIC", signature(IC = "HampIC", 
                                 risk = "asCov",
                                 neighbor = "missing",
                                 L2Fam = "L2ParamFamily"),
    function(IC, risk, L2Fam, withCheck= TRUE, ...){
        if(missing(withCheck)) withCheck <- TRUE
        Cov <- eval(IC@Risks[["asCov"]])
        if(is.null(Cov)){
           if(numberOfMaps(L2Fam@L2deriv)==1){ ## L2derivDim <- L2Fam@L2deriv
              L2deriv <- L2Fam@L2derivDistr[[1]]
              A <- as.vector(IC@stand)
              c0 <- IC@clip/abs(A)
              z <- IC@cent/A
              neighbor <- ContNeighborhood(1)
              Cov <- getInfV(L2deriv = L2deriv, neighbor = neighbor,
                         biastype = biastype(IC), clip = c0, cent = z, stand = A)
              if(withCheck) .checkICWithWarning(IC, L2Fam, tol=.Machine$double.eps^.25, ...)
              return(list(asCov = list(distribution = .getDistr(L2Fam),
                          value = Cov)))
            }
            return(getRiskIC(as(IC, "IC"), risk = risk, L2Fam = L2Fam, withCheck= withCheck, ...))
        }else
            return(list(asCov = list(distribution = .getDistr(L2Fam), value = Cov)))
    })

setMethod("getRiskIC", signature(IC = "TotalVarIC", 
                                 risk = "asCov",
                                 neighbor = "missing",
                                 L2Fam = "L2ParamFamily"),
    function(IC, risk, L2Fam, withCheck = TRUE, ...){
        Cov <- eval(IC@Risks[["asCov"]])
        if(missing(withCheck)) withCheck <- TRUE
        if (is.null(Cov)){
            L2deriv <- L2Fam@L2derivDistr[[1]]
            A <- as.vector(IC@stand)
            c0 <- (IC@clipUp-IC@clipLo)/abs(A)
            z <- IC@clipLo/abs(A)
            neighbor <- TotalVarNeighborhood(1)
            if(withCheck) .checkICWithWarning(IC, L2Fam, tol=.Machine$double.eps^.25, ...)
            Cov <- getInfV(L2deriv = L2deriv, neighbor = neighbor,
                       biastype = biastype(IC), clip = c0, cent = z, stand = A)
            }
        return(list(asCov = list(distribution = .getDistr(L2Fam), value = Cov)))
    })

###############################################################################
## asymptotic Bias for various types
###############################################################################
setMethod("getBiasIC", signature(IC = "HampIC",
                                 neighbor = "UncondNeighborhood"),
    function(IC, neighbor, L2Fam,...){
        if(missing(L2Fam))
            L2Fam <- force(eval(IC@CallL2Fam))

            Bias <- IC@Risks$asBias$value

        return(list(asBias = list(distribution = .getDistr(L2Fam), 
                    neighborhood = neighbor@type, value = Bias)))
    })

setMethod("getBiasIC", signature(IC = "TotalVarIC",
                                 neighbor = "UncondNeighborhood"),
    function(IC, neighbor, L2Fam,...){
        if(missing(L2Fam))
            L2Fam <- force(eval(IC@CallL2Fam))

        Bias <- IC@Risks$asBias$value
        if (is.null(Bias)){
            Bias <- if(is(neighbor,"ContNeighborhood"))
                     max(IC@clipUp,-IC@clipLo) else IC@clipUp-IC@clipLo
        }

        return(list(asBias = list(distribution = .getDistr(L2Fam), 
                    neighborhood = neighbor@type, value = Bias)))
    })

