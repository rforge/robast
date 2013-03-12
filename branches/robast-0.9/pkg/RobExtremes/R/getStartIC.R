setMethod("getStartIC",signature(model = "L2ScaleShapeUnion", risk = "interpolRisk"),
           function(model, risk, ...){

    mc <- match.call(expand.dots=TRUE)

    gridn <- type(risk)
    nam <- gsub(" ","",name(model))
    param1 <- param(model)
    nsng <- character(0)
    sng <- try(getFromNamespace(gridn, ns = "RobAStRDA"), silent=TRUE)
    if(!is(sng,"try-error")) nsng <- names(sng)
    if(length(nsng)){
       if(nam %in% nsng){
          fctN <- .versionSuff("fun")
          interpolfct <- sng[[nam]][[fctN]]
          .modifyIC <- function(L2Fam, IC){
                   para <- param(L2Fam)
                   .getPsi(para, interpolfct, L2Fam, type(risk))}
          IC0 <- .getPsi(param1, interpolfct, model, type(risk))
          IC0@modifyIC <- .modifyIC
          return(IC0)
       }
    }
    mc$risk <- if(type(risk)==".MBRE") asMSE() else asBias()
    mc$neighbor <- ContNeighborhood(radius=0.5)
    return(do.call(getStartIC, as.list(mc[-1]), envir=parent.frame(2)))
    })

