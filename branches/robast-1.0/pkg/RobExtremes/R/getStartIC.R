setMethod("getStartIC",signature(model = "L2ScaleShapeUnion", risk = "interpolRisk"),
           function(model, risk, ...){

    mc <- match.call(expand.dots=TRUE)
    mc$risk <- if(type(risk)==".MBRE") asMSE() else asBias()
    mc$neighbor <- ContNeighborhood(radius=0.5)

    gridn <- type(risk)
    nam <- gsub(" ","",name(model))
    param1 <- param(model)

    scshnm <- scaleshapename(model)
    shnam <- scshnm["shape"]

    ### check whether mc[-1] is a good strategy to delete risk parameter...!!!

    nsng <- character(0)
    sng <- try(getFromNamespace(gridn, ns = "RobAStRDA"), silent=TRUE)
    if(!is(sng,"try-error")) nsng <- names(sng)
    if(length(nsng)){
       if(nam %in% nsng){
          fctN <- .versionSuff("fun")
          interpolfct <- sng[[nam]][[fctN]]
          .modifyIC0 <- function(L2Fam, IC){
                    para <- param(L2Fam)
                    if(!.is.na.Psi(para, interpolfct, shnam))
                       return(.getPsi(para, interpolfct, L2Fam, type(risk)))
                    else
                       return(do.call(getStartIC, as.list(mc[-1]),
                              envir=parent.frame(2)))
          }
          .modifyIC <- function(L2Fam,IC){
               psi.0 <- .modifyIC0(L2Fam,IC)
               psi.0@modifyIC <- .modifyIC
               return(psi.0)
          }
          if(!.is.na.Psi(param1, interpolfct, shnam)){
             IC0 <- .getPsi(param1, interpolfct, model, type(risk))
             IC0@modifyIC <- .modifyIC
             return(IC0)
          }
       }
    }
    return(do.call(getStartIC, as.list(mc[-1]), envir=parent.frame(2)))
    })

setMethod("getStartIC",signature(model = "L2LocScaleShapeUnion", risk = "interpolRisk"),
           function(model, risk, ...){

    mc <- match.call(expand.dots=TRUE)
    mc$risk <- if(type(risk)==".MBRE") asMSE() else asBias()
    mc$neighbor <- ContNeighborhood(radius=0.5)

    gridn <- type(risk)
    nam <- gsub(" ","",name(model))
    param1 <- param(model)

    scshnm <- scaleshapename(model)
    shnam <- scshnm["shape"]

    nsng <- character(0)
    sng <- try(getFromNamespace(gridn, ns = "RobAStRDA"), silent=TRUE)
    if(!is(sng,"try-error")) nsng <- names(sng)
    if(length(nsng)){
       if(nam %in% nsng){
          fctN <- .versionSuff("fun")
          interpolfct <- sng[[nam]][[fctN]]
          .modifyIC <- function(L2Fam, IC){
                    para <- param(L2Fam)
                    if(!.is.na.Psi(para, interpolfct, shnam))
                       return(.getPsi.wL(para, interpolfct, L2Fam, type(risk)))
                    else
                       return(do.call(getStartIC, as.list(mc[-1]),
                              envir=parent.frame(2)))
          }
          if(!.is.na.Psi(param1, interpolfct, shnam)){
             IC0 <- .getPsi.wL(param1, interpolfct, model, type(risk))
             IC0@modifyIC <- .modifyIC
             return(IC0)
          }
       }
    }
    return(do.call(getStartIC, as.list(mc[-1]), envir=parent.frame(2)))
    })

