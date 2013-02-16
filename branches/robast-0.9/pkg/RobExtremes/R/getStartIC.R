setMethod("getStartIC",signature(model = "L2ScaleShapeUnion", risk = "interpolRisk"),
           function(model, risk, ...){

    mc <- match.call(expand.dots=TRUE)

    gridn <- type(risk)
    nam <- name(model)
    xi <- main(param(model))["shape"] #[scaleshapename(model)["shape"]]
    beta <- main(param(model))["scale"] #[scaleshapename(model)["scale"]]
    nsng <- character(0)
    sng <- try(getFromNamespace(.versionSuff(gridn), ns = "RobAStRDA"),
                                 silent=TRUE)
    if(!is(sng,"try-error")) nsng <- names(sng)
    if(length(nsng)){
       if(nam %in% nsng){
          interpolfct <- sng[[nam]]$fct
          .modifyIC <- function(L2Fam, IC){
                   para <- param(L2Fam)
                   xi0 <- main(para)["shape"]#[scaleshapename(L2Fam)["scale"]]
                   beta0 <- main(para)["scale"]#[scaleshapename(L2Fam)["scale"]]
                   .getPsi(xi0,beta0, interpolfct, L2Fam, type(risk))}
          IC0 <- .getPsi(xi, beta, interpolfct, model, type(risk))
          IC0@modifyIC <- .modifyIC
          return(IC0)
       }
    }
    mc$risk <- if(type(risk)==".MBRE") asMSE() else asBias()
    mc$neighbor <- ContNeighborhood(radius=0.5)
    return(do.call(getStartIC, as.list(mc[-1]), envir=parent.frame(2)))
    })

