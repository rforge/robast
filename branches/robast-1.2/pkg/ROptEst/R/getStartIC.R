setMethod("getStartIC",signature(model = "ANY", risk = "ANY"),
           function(model, risk, ...) stop("not yet implemented"))

setMethod("getStartIC",signature(model = "L2ParamFamily", risk = "asGRisk"),
           function(model, risk, ..., withEvalAsVar = TRUE, withMakeIC = FALSE,
           ..debug=FALSE, modifyICwarn = NULL){
    mc <- match.call(expand.dots=FALSE, call = sys.call(sys.parent(1)))
    dots <- as.list(mc$"...")
    if(missing(..debug)||!is.logical(..debug)) ..debug <- FALSE
    if("fsCor" %in% names(dots)){
        fsCor <- eval(dots[["fsCor"]])
        dots$fsCor <- NULL
    }else fsCor <- 1
    if("eps" %in% names(dots)){
       eps <- dots[["eps"]]
       dots$eps <- NULL
    }else eps <- NULL
    if("neighbor" %in% names(dots)){
       neighbor <- eval(dots[["neighbor"]])
       dots$neighbor <- NULL
    }else neighbor <- ContNeighborhood()


    if(is.null(eps[["e"]])){
        sm.rmx <- selectMethod("radiusMinimaxIC", signature(
               class(model),class(neighbor),class(risk)))
        dots.rmx <- .fix.in.defaults(dots, sm.rmx)
        dots.rmx$L2Fam <- NULL
        dots.rmx$neighbor <- NULL
        dots.rmx$risk <- NULL
        dots.rmx$modifyICwarn <- modifyICwarn
        dots.rmx[["warn"]] <- FALSE
        if(!is.null(dots[["warn"]]))if(eval(dots[["warn"]])) dots.rmx[["warn"]] <- TRUE

        infMod <- InfRobModel(center = model, neighbor = neighbor)
        dots.rmx$loRad <- eps$sqn * eps$lower
        dots.rmx$upRad <- eps$sqn * eps$upper
        arg.rmx <- c(list(L2Fam = model, neighbor = neighbor,
                                   risk = risk), dots.rmx)
        if(..debug) print(c(arg.rmx=arg.rmx))
        ICstart <- do.call(radiusMinimaxIC, args=arg.rmx, envir=parent.frame(2))
        if(..debug) print(ICstart)
        if(!isTRUE(all.equal(fsCor, 1, tol = 1e-3))){
            neighbor@radius <- neighborRadius(ICstart)*fsCor
            arg.optic <- c(list( model = infMod, risk = risk), dots.optic)
            if(..debug) print(c(arg.optic=arg.optic))
            ICstart <- do.call(optIC, args = arg.optic, envir=parent.frame(2))
        }
    }else{
        neighbor@radius <- eps$sqn * fsCor * eps$e
        infMod <- InfRobModel(center = model, neighbor = neighbor)
        sm.optic <- selectMethod("optIC", signature(
                                         class(infMod),class(risk)))
        dots.optic <- .fix.in.defaults(dots, sm.optic)
        dots.optic$model <- NULL
        dots.optic$risk <- NULL
        dots.optic$.withEvalAsVar <- withEvalAsVar
        dots.optic$withMakeIC <- withMakeIC
        dots.optic$modifyICwarn <- modifyICwarn
        dots.optic[["warn"]] <- FALSE
        if(!is.null(dots[["warn"]]))if(eval(dots[["warn"]])) dots.optic[["warn"]] <- TRUE

        arg.optic <- c(list(model = infMod, risk = risk),
                           dots.optic)
        if(..debug) print(c(arg.optic=arg.optic))
        ICstart <- do.call(optIC, args = arg.optic, envir=parent.frame(2))
    }
  return(ICstart)
  })

setMethod("getStartIC",signature(model = "L2ParamFamily", risk = "asCov"),
           function(model, risk, ..., withMakeIC = FALSE, ..debug=FALSE){
    return(optIC(model, risk, withMakeIC = withMakeIC))
  })
setMethod("getStartIC",signature(model = "L2ParamFamily", risk = "trAsCov"),
     getMethod("getStartIC", signature(model = "L2ParamFamily", risk = "asCov"))
           )

setMethod("getStartIC",signature(model = "L2ParamFamily", risk = "asBias"),
           function(model, risk, ..., withMakeIC = FALSE, ..debug=FALSE,
           modifyICwarn = NULL){
    mc <- match.call(expand.dots=FALSE, call = sys.call(sys.parent(1)))
    dots <- as.list(mc$"...")

    if("neighbor" %in% names(dots)){
       neighbor <- eval(dots[["neighbor"]])
       dots$neighbor <- NULL
    }else neighbor <- ContNeighborhood()
    if("warn" %in% names(dots)) dots$warn <- NULL

    infMod <- InfRobModel(center = model, neighbor = neighbor)
    return(do.call(optIC, c(list(infMod, risk), dots, list(warn = FALSE,
                withMakeIC = withMakeIC, modifyICwarn = modifyICwarn)),
                envir=parent.frame(2)))
           })


setMethod("getStartIC",signature(model = "L2ParamFamily", risk = "asAnscombe"),
           function(model, risk, ..., withEvalAsVar = TRUE, withMakeIC = FALSE,
           ..debug=FALSE, modifyICwarn = NULL){
    mc <- match.call(expand.dots=FALSE, call = sys.call(sys.parent(1)))
    dots <- as.list(mc$"...")
    if(missing(..debug)||!is.logical(..debug)) ..debug <- FALSE
    if("neighbor" %in% names(dots)){
       neighbor <- eval(dots[["neighbor"]])
       dots$neighbor <- NULL
    }else neighbor <- ContNeighborhood()

    infMod <- InfRobModel(center = model, neighbor = neighbor)
    sm.optic <- selectMethod("optIC", signature(
                                     class(infMod),class(risk)))
    dots.optic <- .fix.in.defaults(dots, sm.optic)
    dots.optic$model <- NULL
    dots.optic$risk <- NULL
    dots.optic$.withEvalAsVar <- withEvalAsVar
    dots.optic$withMakeIC <- withMakeIC
    dots.optic$modifyICwarn <- modifyICwarn
    dots.optic[["warn"]] <- FALSE
    if(!is.null(dots[["warn"]]))if(eval(dots[["warn"]])) dots.optic[["warn"]] <- TRUE

    arg.optic <- c(list(model = infMod, risk = risk), dots.optic)
    if(..debug) print(c(arg.optic=arg.optic))
    ICstart <- do.call(optIC, args = arg.optic, envir=parent.frame(2))
    return(ICstart)
  })
