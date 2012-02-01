## Generating function

Av1CondTotalVarIC <- function(name, CallL2Fam = call("L2RegTypeFamily"),
                   Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x[1]*x[2]}), 
                                                    Domain = EuclideanSpace(dimension = 2))), 
                   Risks, Infos, clipUp = Inf, stand = as.matrix(1), 
                   clipLo = RealRandVariable(Map = list(function(x){ -Inf }),
                                             Domain = EuclideanSpace(dimension = 1)), 
                   lowerCase = NULL, neighborRadius = 0,
                   w = new("CondHampelWeight"),
                   normtype = NormType(), biastype = symmetricBias(),
                   modifyIC = NULL){
    if(missing(name))
        name <- "conditionally centered IC for average conditional total variation neighborhoods"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                    dimnames=list(character(0), c("method", "message")))
    
    TVIC <- new("Av1CondTotalVarIC")
    TVIC@name <- name
    TVIC@Curve <- Curve
    TVIC@Risks <- Risks
    TVIC@Infos <- Infos
    TVIC@CallL2Fam <- CallL2Fam
    TVIC@clipUp <- clipUp
    TVIC@clipLo <- clipLo
    TVIC@stand <- stand
    TVIC@lowerCase <- lowerCase
    TVIC@neighborRadius <- neighborRadius
    TVIC@weight <- w
    TVIC@biastype <- biastype
    TVIC@normtype <- normtype
    TVIC@modifyIC <- modifyIC
    return(TVIC)
}

## generate IC
## for internal use only!
setMethod("generateIC", signature(neighbor = "Av1CondTotalVarNeighborhood", 
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
        return(Av1CondTotalVarIC(
                name = "conditionally centered IC of contamination type", 
                CallL2Fam = L2call,
                Curve = generateIC.fct(neighbor, L2Fam, res),
                clipUp = b,
                clipLo = a,
                stand = A,
                lowerCase = d,
                neighborRadius = neighbor@radius,
                Risks = res$risk,
                Infos = matrix(res$info, ncol = 2, 
                            dimnames = list(character(0), c("method", "message")))))
    })

## Access methods
setMethod("clipUp", "Av1CondTotalVarIC", function(object) object@clipUp)
setMethod("clipLo", "Av1CondTotalVarIC", function(object) object@clipLo)

## replace methods
setReplaceMethod("clipUp", "Av1CondTotalVarIC", 
    function(object, value){ 
        stopifnot(is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = object@clipLo, b = value, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = Av1CondTotalVarNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("clipUp<-", "The clipping bound has been changed")
        addInfo(object) <- c("clipUp<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("clipLo", "Av1CondTotalVarIC", 
    function(object, value){ 
        stopifnot(is(value, "RealRandVariable"))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = value, b = object@clipUp, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = Av1CondTotalVarNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("clipLo<-", "The centering constant has been changed")
        addInfo(object) <- c("clipLo<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("stand", "Av1CondTotalVarIC", 
    function(object, value){ 
        stopifnot(is.matrix(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = value, a = object@clipLo, b = object@clipUp, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = Av1CondTotalVarNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("stand<-", "The standardizing matrix has been changed")
        addInfo(object) <- c("stand<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("lowerCase", "Av1CondTotalVarIC", 
    function(object, value){ 
        stopifnot(is.null(value)||is.numeric(value))
        L2Fam <- eval(object@CallL2Fam)
        res <- list(A = object@stand, a = object@clipLo, b = object@clipUp, d = value,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = Av1CondTotalVarNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("lowerCase<-", "The slot 'lowerCase' has been changed")
        addInfo(object) <- c("lowerCase<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("neighborRadius", "Av1CondTotalVarIC", 
    function(object, value){ 
        object@neighborRadius <- value
        if(any(value < 0)) # radius vector?!
            stop("'value' has to be in [0, Inf]")
        addInfo(object) <- c("neighborRadius<-", "The slot 'neighborRadius' has been changed")
        addInfo(object) <- c("neighborRadius<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("CallL2Fam", "Av1CondTotalVarIC",
    function(object, value){ 
        L2Fam <- eval(value)
        res <- list(A = object@stand, a = object@clipLo, b = object@clipUp, d = object@lowerCase,
                    risk = object@Risks, info = object@Infos)
        object <- generateIC(neighbor = Av1CondTotalVarNeighborhood(radius = object@neighborRadius), 
                             L2Fam = L2Fam, res = res)
        addInfo(object) <- c("CallL2Fam<-", "The slot 'CallL2Fam' has been changed")
        addInfo(object) <- c("CallL2Fam<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
