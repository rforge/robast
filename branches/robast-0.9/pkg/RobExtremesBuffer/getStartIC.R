setMethod("getStartIC",signature(model = "ANY", risk = "ANY"),
           function(model, risk, ...) stop("not yet implemented"))

setMethod("getStartIC",signature(model = "L2ParamFamily", risk = "asRisk"),
           function(model, risk, ...){
    mc <- match.call(expand=TRUE)
    eps <- mc$eps
    dots <- mc$dots

    if(is.null(eps$e))){
        r.lower <- eps$sqn * eps$lower
        r.upper <- eps$sqn * eps$upper
        ICstart <- do.call(radiusMinimaxIC,
                    c(list(L2Fam = mc$L2FamStart, neighbor = mc$neighbor,
                                   risk = mc$risk,
                                   loRad = r.lower, upRad = r.upper,
                                   verbose = mc$verbose,
                                   OptOrIter = mc$OptOrIter),dots))
        if(!isTRUE(all.equal(mc$fsCor, 1, tol = 1e-3))){
            neighbor@radius <- neighborRadius(ICstart)*mc$fsCor
            infMod <- InfRobModel(center = mc$L2FamStart, neighbor = mc$neighbor)
            ICstart <- do.call(optIC, c(list( model = mc$infMod, risk = mc$risk,
                               verbose = mc$verbose, OptOrIter = mc$OptOrIter),
                               dots))
        }
    }else{
        neighbor@radius <- eps$sqn*eps$e*mc$fsCor
        infMod <- InfRobModel(center = mc$L2FamStart, neighbor = mc$neighbor)
        ICstart <- do.call(optIC, c(list(model = mc$infMod, risk = mc$risk,
                           verbose = mc$verbose, OptOrIter = mc$OptOrIter),
                           dots))
    }
  return(ICstart)
           })



setMethod("getStartIC",signature(model = "L2ScaleShapeUnion", risk = "interpolRisk"),
           function(model, risk, ...){

    mc <- match.call(expand=TRUE)

    gridn <- type(risk)
    nam <- name(model)
    xi <- main(param(model))[scaleshapename(model)["scale"]]
    nsng <- character(0)
    sng <- try(getFromNamespace(gridn, ns = "ROptEst"),silent=TRUE)
    if(!is(sng,"try-error")) nsng <- names(sng)
    if(length(nsng)){
       if(nam %in% nsng){
          interpolfct <- sng[[nam]]$fct
          return(.getPsi(xi, interpolfct, L2Fam, type(risk)))
       }
    }
    mc$risk <- if(type(risk)==".MBRE") asMSE(r=0.5) else asBias()
    return(do.call(getStartIC, mc[-1]))
    })

