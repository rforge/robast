## Generating function
CondTotalVarIC <- function(name, CallL2Fam = call("L2RegTypeFamily"),
                   Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x[1]*x[2]}), 
                                                    Domain = EuclideanSpace(dimension = 2))), 
                   Risks, Infos, 
                   clipUp = RealRandVariable(Map = list(function(x){ Inf }), Domain = Reals()),
                   stand = as.matrix(1), 
                   clipLo = RealRandVariable(Map = list(function(x){ -Inf }), Domain = Reals()), 
                   lowerCase = NULL, neighborRadius = 0, neighborRadiusCurve = function(x){1},
                   w = new("CondBdStWeight"),
                   normtype = NormType(), biastype = symmetricBias(),
                   modifyIC = NULL){
    if(missing(name))
        name <- "conditionally centered IC for average conditional total variation neighborhoods"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                    dimnames=list(character(0), c("method", "message")))
    
    TVIC <- new("CondTotalVarIC")
    TVIC@name <- name
    TVIC@Curve <- Curve
    TVIC@Risks <- Risks
    TVIC@Infos <- Infos
    TVIC@CallL2Fam <- CallL2Fam
    TVIC@clipLo <- clipLo
    TVIC@clipUp <- clipUp
    TVIC@stand <- stand
    TVIC@lowerCase <- lowerCase
    TVIC@neighborRadius <- neighborRadius
    TVIC@neighborRadiusCurve <- neighborRadiusCurve
    TVIC@weight <- w
    TVIC@biastype <- biastype
    TVIC@normtype <- normtype
    TVIC@modifyIC <- modifyIC

    return(TVIC)
}

## generate IC
## for internal use only!
setMethod("generateIC", signature(neighbor = "CondTotalVarNeighborhood", 
                                  L2Fam = "L2RegTypeFamily"),
    function(neighbor, L2Fam, res){
        return(CondContIC(
                name = "conditionally centered IC of contamination type",
                CallL2Fam = L2call,
                Curve = generateIC.fct(neighbor, L2Fam, res),
                clip = b,
                cent = a,
                stand = A,
                lowerCase = d,
                neighborRadius = neighbor@radius,
                neighborRadiusCurve = neighbor@radiusCurve,
                modifyIC = res$modifyIC,
                normtype = normtype,
                biastype = biastype,
                Risks = res$risk,
                Infos = matrix(res$info, ncol = 2,
                            dimnames = list(character(0), c("method", "message")))))
    })
        A <- res$A
        a <- res$a
        b <- res$b
        d <- res$d
        normtype <- res$normtype
        biastype <- res$biastype
        w <- res$w
        L2call <- L2Fam@fam.call
        L2call$trafo <- trafo(L2Fam)
        return(CondTotalVarIC(
                name = "conditionally centered IC of contamination type", 
                CallL2Fam = L2call
                Curve = generateIC.fct(neighbor, L2Fam, res),
                clipUp = b,
                clipLo = a,
                stand = A,
                lowerCase = d,
                neighborRadius = neighbor@radius,
                neighborRadiusCurve = neighbor@radiusCurve,
                modifyIC = res$modifyIC,
                normtype = normtype,
                biastype = biastype,
                Risks = res$risk,
                Infos = matrix(res$info, ncol = 2, 
                            dimnames = list(character(0), c("method", "message")))))
    })

## Access methods
setMethod("clipUp", "CondTotalVarIC", function(object) object@clipUp)
setMethod("clipLo", "CondTotalVarIC", function(object) object@clipLo)
setMethod("neighborRadiusCurve", "CondTotalVarIC", function(object) object@neighborRadiusCurve)

## replace methods
setReplaceMethod("clipUp", "CondTotalVarIC", 
    function(object, value){ 
        stopifnot(is(value, "RealRandVariable"))
        L2Fam <- eval(object@CallL2Fam)
        w <- object@weight
        clip(w) <- function(x) c(clipLo(x), value(x))
        CNB <- CondTotalVarNeighborhood(radius = object@neighborRadius,
                                    radiusCurve = object@neighborRadiusCurve)
        weight(w) <- getweight(w, neighbor = CNB,
                               biastype = object@biastype,
                               normW = object@normtype)
        res <- list(A = object@stand, a = object@clipLo, b = value, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos, w = w,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = CNB,
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("clipUp<-", "The clipping bound has been changed")
        addInfo(object) <- c("clipUp<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("clipLo", "CondTotalVarIC", 
    function(object, value){ 
        stopifnot(is(value, "RealRandVariable"))
        L2Fam <- eval(object@CallL2Fam)
        w <- object@weight
        clip(w) <- function(x) c(value(x), clipUp(x))
        CNB <- CondTotalVarNeighborhood(radius = object@neighborRadius,
                                    radiusCurve = object@neighborRadiusCurve)
        weight(w) <- getweight(w, neighbor = CNB,
                               biastype = object@biastype,
                               normW = object@normtype)
        res <- list(A = object@stand, a = value, b = object@clipUp, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos, w = w,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = CNB,
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("clipLo<-", "The centering constant has been changed")
        addInfo(object) <- c("clipLo<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("stand", "CondTotalVarIC", 
    function(object, value){ 
        stopifnot(is.matrix(value))
        L2Fam <- eval(object@CallL2Fam)
        w <- object@weight
        stand(w) <- value
        CNB <- CondTotalVarNeighborhood(radius = object@neighborRadius,
                                    radiusCurve = object@neighborRadiusCurve)
        weight(w) <- getweight(w, neighbor = CNB,
                               biastype = object@biastype,
                               normW = object@normtype)
        res <- list(A = value, a = object@clipLo, b = object@clipUp, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos, w = w,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = CNB,
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("stand<-", "The standardizing matrix has been changed")
        addInfo(object) <- c("stand<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("lowerCase", "CondTotalVarIC", 
    function(object, value){ 
        stopifnot(is.null(value)||is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = object@clipLo, b = object@clipUp, d = value,
                    risk = object@Risks, info = object@Infos, w = object@weight,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = CondTotalVarNeighborhood(radius = object@neighborRadius, 
                                                            radiusCurve = object@neighborRadiusCurve), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("lowerCase<-", "The slot 'lowerCase' has been changed")
        addInfo(object) <- c("lowerCase<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("neighborRadius", "CondTotalVarIC", 
    function(object, value){ 
        object@neighborRadius <- value
        if(any(value < 0)) # radius vector?!
            stop("'value' has to be in [0, Inf]")
        addInfo(object) <- c("neighborRadius<-", "The slot 'neighborRadius' has been changed")
        addInfo(object) <- c("neighborRadius<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("neighborRadiusCurve", "CondTotalVarIC", 
    function(object, value){ 
        object@neighborRadiusCurve <- value
        if(length(formals(value)) != 1)
            stop("'value' has to be a function of one argument")
        if(names(formals(value)) != "x")
            stop("'value' has to be a function with argument name = 'x'")
        addInfo(object) <- c("neighborRadiusCurve<-", "The slot 'neighborRadiusCurve' has been changed")
        addInfo(object) <- c("neighborRadiusCurve<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("CallL2Fam", "CondTotalVarIC",
    function(object, value){ 
        L2Fam <- eval(value)
        res <- list(A = object@stand, a = object@clipLo, b = object@clipUp, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos, w = object@weight,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = CondTotalVarNeighborhood(radius = object@neighborRadius, 
                                                            radiusCurve = object@neighborRadiusCurve), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("CallL2Fam<-", "The slot 'CallL2Fam' has been changed")
        addInfo(object) <- c("CallL2Fam<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
