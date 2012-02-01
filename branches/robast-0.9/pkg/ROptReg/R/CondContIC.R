## Generating function
CondContIC <- function(name, CallL2Fam = call("L2RegTypeFamily"),
                   Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x[1]*x[2]}), 
                                                    Domain = EuclideanSpace(dimension = 2))), 
                   Risks, Infos, 
                   clip = RealRandVariable(Map = list(function(x){ Inf }), Domain = Reals()), 
                   stand = as.matrix(1), 
                   cent = EuclRandVarList(RealRandVariable(Map = list(function(x){numeric(length(x))}),
                                                    Domain = EuclideanSpace(dimension = 2))), 
                   lowerCase = NULL, neighborRadius = 0, neighborRadiusCurve = function(x){1},
                   w = new("CondHampelWeight"),
                   normtype = NormType(), biastype = symmetricBias(),
                   modifyIC = NULL){
    if(missing(name))
        name <- "conditionally centered IC for average conditional contamination neighborhoods"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                    dimnames=list(character(0), c("method", "message")))

    contIC <- new("CondContIC")
    contIC@name <- name
    contIC@Curve <- Curve
    contIC@Risks <- Risks
    contIC@Infos <- Infos
    contIC@CallL2Fam <- CallL2Fam
    contIC@clip <- clip
    contIC@cent <- cent
    contIC@stand <- stand
    contIC@lowerCase <- lowerCase
    contIC@neighborRadius <- neighborRadius
    contIC@neighborRadiusCurve <- neighborRadiusCurve
    contIC@weight <- w
    contIC@biastype <- biastype
    contIC@normtype <- normtype
    contIC@modifyIC <- modifyIC

    return(contIC)
}

## generate IC
## for internal use only!
setMethod("generateIC", signature(neighbor = "CondContNeighborhood", 
                                  L2Fam = "L2RegTypeFamily"),
    function(neighbor, L2Fam, res){
        A <- res$A
        a <- res$a
        b <- res$b
        d <- res$d
        normtype <- res$normtype
        biastype <- res$biastype
        w <- res$w
        L2call <- L2Fam@fam.call
        L2call$trafo <- trafo(L2Fam)
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

## Access methods
setMethod("clip", "CondContIC", function(object) object@clip)
setMethod("cent", "CondContIC", function(object) object@cent)
setMethod("neighborRadiusCurve", "CondContIC", function(object) object@neighborRadiusCurve)

## replace methods
setReplaceMethod("clip", "CondContIC", 
    function(object, value){ 
        stopifnot(is(value, "RealRandVariable"))
        L2Fam <- eval(object@CallL2Fam)
        w <- object@weight
        clip(w) <- value
        CNB <- CondContNeighborhood(radius = object@neighborRadius,
                                    radiusCurve = object@neighborRadiusCurve)
        weight(w) <- getweight(w, neighbor = CNB,
                               biastype = object@biastype,
                               normW = object@normtype)
        res <- list(A = object@stand, a = object@cent, b = value, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos, w = w,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = CNB,
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("clip<-", "The clipping bound has been changed")
        addInfo(object) <- c("clip<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("cent", "CondContIC", 
    function(object, value){ 
        stopifnot(is(value, "EuclRandVarList"))
        L2Fam <- eval(object@CallL2Fam)
        w <- object@weight
        cent(w) <- as.vector(solve(object@stand) %*% value)
        CNB <- CondContNeighborhood(radius = object@neighborRadius,
                                    radiusCurve = object@neighborRadiusCurve)
        weight(w) <- getweight(w, neighbor = CNB,
                               biastype = object@biastype,
                               normW = object@normtype)
        res <- list(A = object@stand, a = value, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos, w = w,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = CNB,
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("cent<-", "The centering constant has been changed")
        addInfo(object) <- c("cent<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("CallL2Fam", "CondContIC",
    function(object, value){ 
        L2Fam <- eval(value)
        res <- list(A = object@stand, a = object@cent, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos, w = object@weight,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = CondContNeighborhood(radius = object@neighborRadius, 
                                                        radiusCurve = object@neighborRadiusCurve), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("CallL2Fam<-", "The slot 'CallL2Fam' has been changed")
        addInfo(object) <- c("CallL2Fam<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("stand", "CondContIC",
    function(object, value){
        stopifnot(is.matrix(value))
        L2Fam <- eval(object@CallL2Fam)
        w <- object@weight
        stand(w) <- value
        CNB <- CondContNeighborhood(radius = object@neighborRadius,
                                    radiusCurve = object@neighborRadiusCurve)
        weight(w) <- getweight(w, neighbor = CNB,
                               biastype = object@biastype,
                               normW = object@normtype)
        res <- list(A = value, a = object@cent, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos, w = w,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = CNB,
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("stand<-", "The standardizing matrix has been changed")
        addInfo(object) <- c("stand<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("lowerCase", "CondHampIC",
    function(object, value){
        stopifnot(is.null(value)||is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = object@cent, b = object@clip, d = value,
                    risk = object@Risks, info = object@Infos, w = object@weight,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = CondContNeighborhood(radius = object@neighborRadius,
                                                        radiusCurve = object@neighborRadiusCurve),
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("lowerCase<-", "The slot 'lowerCase' has been changed")
        addInfo(object) <- c("lowerCase<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
