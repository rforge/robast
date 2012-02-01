## Generating function
Av2CondContIC <- function(name, CallL2Fam = call("L2RegTypeFamily"),
                   Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x[1]*x[2]}), 
                                                    Domain = EuclideanSpace(dimension = 2))), 
                   Risks, Infos, clip = Inf, stand = 1, cent = 0, 
                   lowerCase = NULL, neighborRadius = 0,
                   w = new("CondHampelWeight"),
                   normtype = NormType(), biastype = symmetricBias(),
                   modifyIC = NULL, Kinv=matrix(1),D=matrix(1)){
    if(missing(name))
        name <- "conditionally centered IC for average square conditional contamination neighborhoods"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                    dimnames=list(character(0), c("method", "message")))
    
    contIC <- new("Av2CondContIC")
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
    contIC@weight <- w
    contIC@biastype <- biastype
    contIC@normtype <- normtype
    contIC@modifyIC <- modifyIC
    contIC@Kinv <- Kinv
    contIC@D <- D
    return(contIC)
}

## generate IC
## for internal use only!
setMethod("generateIC", signature(neighbor = "Av2CondContNeighborhood", 
                                  L2Fam = "L2RegTypeFamily"),
    function(neighbor, L2Fam, res){
        normtype <- res$normtype
        biastype <- res$biastype
        w <- res$w
        L2call <- L2Fam@fam.call
        res$D <- D <- L2call$trafo <- trafo(L2Fam)
        res$Kinv <- Kinv <- solve(E(L2Fam@RegDistr, fun = function(x){ x %*% t(x) }))

        return(Av2CondContIC(
                name = "conditionally centered IC of contamination type", 
                CallL2Fam = L2call,
                Curve = generateIC.fct(neighbor, L2Fam, res),
                clip = b,
                cent = z,
                stand = A,
                lowerCase = d,
                neighborRadius = neighbor@radius,
                Risks = res$risk,
                Kinv = Kinv,
                D=D,
                Infos = matrix(res$info, ncol = 2, 
                            dimnames = list(character(0), c("method", "message")))))
    })

## Access methods
setMethod("clip", "Av2CondContIC", function(object) object@clip)
setMethod("cent", "Av2CondContIC", function(object) object@cent)

## replace methods
setReplaceMethod("clip", "Av2CondContIC", 
    function(object, value){ 
        stopifnot(is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        w <- object@weight
        clip(w) <- value
        CNB <- Av2CondContNeighborhood(radius = object@neighborRadius)
        weight(w) <- getweight(w, neighbor = CNB,
                               biastype = object@biastype,
                               normW = object@normtype, Kinv =object@Kinv,
                               D = object@D)
        res <- list(A = object@stand, z = object@cent, b = value, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos, w = w,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = CNB,
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("clip<-", "The clipping bound has been changed")
        addInfo(object) <- c("clip<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("cent", "Av2CondContIC", 
    function(object, value){ 
        stopifnot(is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        w <- object@weight
        cent(w) <- as.vector(solve(object@stand) %*% value)
        CNB <- Av2CondContNeighborhood(radius = object@neighborRadius)
        weight(w) <- getweight(w, neighbor = CNB,
                               biastype = object@biastype,
                               normW = object@normtype, Kinv =object@Kinv,
                               D = object@D)
        res <- list(A = object@stand, z = value, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos, w = w,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = CNB,
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("cent<-", "The centering constant has been changed")
        addInfo(object) <- c("cent<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("stand", "Av2CondContIC", 
    function(object, value){ 
        stopifnot(is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        w <- object@weight
        stand(w) <- value
        CNB <- Av2CondContNeighborhood(radius = object@neighborRadius)
        weight(w) <- getweight(w, neighbor = CNB,
                               biastype = object@biastype,
                               normW = object@normtype, Kinv =object@Kinv,
                               D = object@D)
        res <- list(A = value, z = object@cent, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos, w = w,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = CNB,
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("stand<-", "The standardizing matrix has been changed")
        addInfo(object) <- c("stand<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("lowerCase", "Av2CondContIC", 
    function(object, value){ 
        stopifnot(is.null(value)||is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, z = object@cent, b = object@clip, d = value,
                    risk = object@Risks, info = object@Infos, w = object@weight,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = Av2CondContNeighborhood(radius = object@neighborRadius),
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("lowerCase<-", "The slot 'lowerCase' has been changed")
        addInfo(object) <- c("lowerCase<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("CallL2Fam", "Av2CondContIC",
    function(object, value){ 
        L2Fam <- eval(value)
        res <- list(A = object@stand, z = object@cent, b = object@clip, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos, w = object@weight,
                    normtype = object@normtype, biastype = object@biastype,
                    modifyIC = object@modifyIC)
        object <- generateIC(neighbor = Av2CondContNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("CallL2Fam<-", "The slot 'CallL2Fam' has been changed")
        addInfo(object) <- c("CallL2Fam<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
