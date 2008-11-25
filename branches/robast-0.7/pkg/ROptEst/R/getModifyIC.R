###############################################################################
## internal functions/methods to fill slot modifyIC
###############################################################################

setMethod("getModifyIC", signature(L2FamIC = "L2ParamFamily", 
                                   neighbor = "Neighborhood", risk = "asRisk"),
    function(L2FamIC, neighbor, risk){
        modIC <- function(L2Fam, IC){}
        body(modIC) <- substitute({ infMod <- InfRobModel(L2Fam, nghb)
                                    optIC(infMod, R) },
                                  list(nghb = neighbor, R = risk))
        return(modIC)
    })

setMethod("getModifyIC", signature(L2FamIC = "L2LocationFamily", 
                                   neighbor = "UncondNeighborhood", risk = "asGRisk"),
    function(L2FamIC, neighbor, risk){
        modIC <- function(L2Fam, IC){
            D <- distribution(eval(CallL2Fam(IC)))
            if(is(L2Fam, "L2LocationFamily") && is(distribution(L2Fam), class(D))){
                CallL2Fam(IC) <- fam.call(L2Fam)
                return(IC)
            }else{
                makeIC(IC, L2Fam)
            }
        }
        return(modIC)
    })

setMethod("getModifyIC", signature(L2FamIC = "L2LocationFamily", 
                                   neighbor = "UncondNeighborhood", risk = "fiUnOvShoot"),
    getMethod("getModifyIC",signature(L2FamIC = "L2LocationFamily", 
                                   neighbor = "UncondNeighborhood", risk = "asGRisk"))
    )

setMethod("getModifyIC", signature(L2FamIC = "L2ScaleFamily", 
                                   neighbor = "ContNeighborhood", risk = "asGRisk"),
    function(L2FamIC, neighbor, risk){
        modIC <- function(L2Fam, IC){
            ICL2Fam <- eval(CallL2Fam(IC))
            if(is(L2Fam, "L2ScaleFamily") && is(distribution(L2Fam), class(distribution(ICL2Fam)))){
                sdneu <- main(L2Fam)
                sdalt <- main(ICL2Fam)
                r <- neighborRadius(IC)
                w <- weight(IC)
                clip(w) <- sdneu*clip(w)/sdalt
                cent(w) <- sdalt*cent(w)/sdneu
                stand(w) <- sdneu^2*stand(w)/sdalt^2
                weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = r), 
                              biastype = biastype(IC), 
                              normW = normtype(IC))
                A <- sdneu^2*stand(IC)/sdalt^2
                b <- sdneu*clip(IC)/sdalt
                res <- list(A = A, a = sdneu*cent(IC)/sdalt, b = b, d = NULL,
                            risk = list(asMSE = A, asBias = b, asCov = A-r^2*b^2), 
                            info = Infos(IC), w = w,
                            normtype = normtype(IC), biastype = biastype(IC),
                            modifyIC = modifyIC(IC))
                IC <- generateIC(neighbor = ContNeighborhood(radius = r),
                                 L2Fam = L2Fam, res = res)
                addInfo(IC) <- c("modifyIC", "The IC has been modified")
                addInfo(IC) <- c("modifyIC", "The entries in 'Infos' may be wrong")
                return(IC)
            }else{
                makeIC(IC, L2Fam)
            }
        }
        return(modIC)
    })

setMethod("getModifyIC", signature(L2FamIC = "L2ScaleFamily", 
                                   neighbor = "TotalVarNeighborhood", risk = "asGRisk"),
    function(L2FamIC, neighbor, risk){
        modIC <- function(L2Fam, IC){
            ICL2Fam <- eval(CallL2Fam(IC))
            if(is(L2Fam, "L2ScaleFamily") && is(distribution(L2Fam), class(distribution(ICL2Fam)))){
                sdneu <- main(L2Fam)
                sdalt <- main(ICL2Fam)
                r <- neighborRadius(IC)
                w <- weight(IC)
                clip(w) <- sdneu*clip(w)/sdalt
                stand(w) <- sdneu^2*stand(w)/sdalt^2
                weight(w) <- getweight(w, neighbor = TotalVarNeighborhood(radius = r), 
                              biastype = biastype(IC), 
                              normW = normtype(IC))
                A <- sdneu^2*stand(IC)/sdalt^2
                blo <- sdneu*clipLo(IC)/sdalt
                b <- sdneu*clipUp(IC)/sdalt - blo
                res <- list(A = A, a = blo, b = b, d = NULL,
                            risk = list(asMSE = A, asBias = b, asCov = A-r^2*b^2), 
                            info = Infos(IC), w = w,
                            normtype = normtype(IC), biastype = biastype(IC),
                            modifyIC = modifyIC(IC))
                IC <- generateIC(neighbor = TotalVarNeighborhood(radius = r),
                                 L2Fam = L2Fam, res = res)
                addInfo(IC) <- c("modifyIC", "The IC has been modified")
                addInfo(IC) <- c("modifyIC", "The entries in 'Infos' may be wrong")
                return(IC)
            }else{
                makeIC(IC, L2Fam)
            }
        }
        return(modIC)
    })

setMethod("getModifyIC", signature(L2FamIC = "L2LocationScaleFamily", 
                                   neighbor = "ContNeighborhood", risk = "asGRisk"),
    function(L2FamIC, neighbor, risk){
        modIC <- function(L2Fam, IC){
            ICL2Fam <- eval(CallL2Fam(IC))
            if(is(L2Fam, "L2LocationScaleFamily") && is(distribution(L2Fam), class(distribution(ICL2Fam)))){
                sdneu <- main(L2Fam)[2]
                sdalt <- main(ICL2Fam)[2]
                r <- neighborRadius(IC)
                w <- weight(IC)
                clip(w) <- sdneu*clip(w)/sdalt
                cent(w) <- sdalt*cent(w)/sdneu
                stand(w) <- sdneu^2*stand(w)/sdalt^2
                weight(w) <- getweight(w, neighbor = ContNeighborhood(radius = r), 
                                       biastype = biastype(IC), 
                                       normW = normtype(IC))
                A <- sdneu^2*stand(IC)/sdalt^2
                b <- sdneu*clip(IC)/sdalt
                a <- sdneu*cent(IC)/sdalt
                mse <- sum(diag(A))
                Cov <- sdneu^2*Risks(IC)$asCov/sdalt^2

                res <- list(A = A, a = sdneu*cent(IC)/sdalt, b = b, d = NULL,
                            risk = list(asCov = Cov,
                                        asMSE = mse, asBias = b, 
                                        trAsCov = mse - r^2*b^2), 
                            info = Infos(IC), w = w,
                            normtype = normtype(IC), biastype = biastype(IC),
                            modifyIC = modifyIC(IC))
                IC <- generateIC(neighbor = ContNeighborhood(radius = r),
                                L2Fam = L2Fam, res = res)
                addInfo(IC) <- c("modifyIC", "The IC has been modified")
                addInfo(IC) <- c("modifyIC", "The entries in 'Infos' may be wrong")
                return(IC)
            }else{
                makeIC(IC, L2Fam)
            }
        }
        return(modIC)
    })
