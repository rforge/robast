setMethod("getStartIC",signature(model = "L2ScaleShapeUnion", risk = "interpolRisk"),
           function(model, risk, ..., withMakeIC = FALSE, ..debug=FALSE,
           modifyICwarn = NULL){

    mc <- match.call(call = sys.call(sys.parent(1)))
    mc$risk <- if(type(risk)==".MBRE") asMSE() else asBias()
    mc$neighbor <- ContNeighborhood(radius=0.5)

    gridn <- gsub("\\.","",type(risk))

    nam <- paste(".",gsub("[F,f]amily","",gsub(" ","",name(model))),sep="")
    if(nam==".Weibull") nam <- ".WeibullFamily"
    if(nam==".GeneralizedPareto") nam <- ".GPareto"

    param1 <- param(model)

    scshnm <- scaleshapename(model)
    shnam <- scshnm["shape"]

    ### check whether mc[-1] is a good strategy to delete risk parameter...!!!

    nsng <- character(0)
    famg <- try(getFromNamespace(nam, ns = "RobAStRDA"), silent=TRUE)
    #sng <- try(getFromNamespace(gridn, ns = "RobAStRDA"), silent=TRUE)
    if(!is(famg,"try-error")) nsng <- names(famg)
    if(length(nsng)){
       if(gridn %in% nsng){
          interpolfct <- famg[[gridn]][[.versionSuff("fun")]]
          rm(famg, nsng, gridn)
          .modifyIC0 <- function(L2Fam, IC, withMakeIC = FALSE){
                    para <- param(L2Fam)
                    if(!.is.na.Psi(para, interpolfct, shnam))
                       return(.getPsi(para, interpolfct, L2Fam, type(risk)))
                    else{
                       IC0 <- do.call(getStartIC, as.list(mc[-1]),
                              envir=parent.frame(2))
                       return(IC0)
                    }
          }
          .modifyIC <- function(L2Fam,IC, withMakeIC = FALSE, ...){
               psi.0 <- .modifyIC0(L2Fam,IC)
               psi.0@modifyIC <- .modifyIC
               if(withMakeIC) psi.0 <- makeIC(psi.0, L2Fam, ...)
               return(psi.0)
          }

          if(!.is.na.Psi(param1, interpolfct, shnam)){
             IC0 <- .getPsi(param1, interpolfct, model, type(risk))
             IC0@modifyIC <- .modifyIC
             if(withMakeIC) IC0 <- makeIC(IC0, model, ...)
             return(IC0)
          }
          rm(mc)
       }
    }
    rm(famg, nsng,gridn)
    IC <- do.call(getStartIC, as.list(mc[-1]), envir=parent.frame(2))
    if(withMakeIC) IC <- makeIC(IC,model,...)
    return(IC)
    })

setMethod("getStartIC",signature(model = "L2LocScaleShapeUnion", risk = "interpolRisk"),
           function(model, risk, ..., withMakeIC = FALSE, ..debug=FALSE,
           modifyICwarn = NULL){

    mc <- match.call(call = sys.call(sys.parent(1)))
    mc$risk <- if(type(risk)==".MBRE") asMSE() else asBias()
    mc$neighbor <- ContNeighborhood(radius=0.5)

    gridn <- gsub("\\.","",type(risk))

    nam <- paste(".",gsub("[F,f]amily","",gsub(" ","",name(model))),sep="")
    if(nam==".GEV") nam <- ".GEVU" 

    param1 <- param(model)

    locscshnm <- locscaleshapename(model)
    shnam <- locscshnm["shape"]
    nsng <- character(0)
    famg <- try(getFromNamespace(nam, ns = "RobAStRDA"), silent=TRUE)
    if(!is(famg,"try-error")) nsng <- names(famg)
    if(length(nsng)){
       if(gridn %in% nsng){
          interpolfct <- famg[[gridn]][[.versionSuff("fun")]]
          rm(famg, nsng, gridn)
          .modifyIC0 <- function(L2Fam, IC){
                    para <- param(L2Fam)
                    if(!.is.na.Psi(para, interpolfct, shnam))
                       return(.getPsi.wL(para, interpolfct, L2Fam, type(risk)))
                    else{
                       IC0 <- do.call(getStartIC, as.list(mc[-1]),
                              envir=parent.frame(2))
                       return(IC0)
                    }
          }
          .modifyIC <- function(L2Fam,IC, withMakeIC = FALSE, ...){
               psi.0 <- .modifyIC0(L2Fam,IC)
               psi.0@modifyIC <- .modifyIC
               if(withMakeIC) psi.0 <- makeIC(psi.0, L2Fam, ...)
               return(psi.0)
          }

          if(!.is.na.Psi(param1, interpolfct, shnam)){
             IC0 <- .getPsi.wL(param1, interpolfct, model, type(risk))
             IC0@modifyIC <- .modifyIC
             if(withMakeIC) IC0 <- makeIC(IC0, model, ...)
             return(IC0)
          }
       }
    }
    rm(famg, nsng,gridn)
    IC <- do.call(getStartIC, as.list(mc[-1]), envir=parent.frame(2))
    if(withMakeIC) IC <- makeIC(IC,model,...)
    return(IC)
    })

