## generating function
IC <- function(name, Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x}), 
               Domain = Reals())), Risks, Infos, CallL2Fam = call("L2ParamFamily"),
               modifyIC = NULL){
    if(missing(name))
        name <- "square integrable (partial) influence curve"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                     dimnames=list(character(0), c("method", "message")))

    if(!is(Domain(Curve[[1]]), "EuclideanSpace"))
        stop("The domain of 'Curve' has to be a Euclidean space")
    if(!is.character(Infos))
        stop("'Infos' contains no matrix of characters")
    for(char in names(Risks))
        if(!extends(char, "RiskType"))
            stop(paste(char, "is no valid 'RiskType'"))
    if(ncol(Infos)!=2)
        stop("'Infos' must have two columns")

    L2Fam <- eval(CallL2Fam)
    trafo <- trafo(L2Fam@param)
    if(nrow(trafo) != dimension(Curve))
        stop("wrong dimension of 'Curve'")
    if(dimension(Domain(L2Fam@L2deriv[[1]])) != dimension(Domain(Curve[[1]])))
        stop("dimension of 'Domain' of 'L2deriv' != dimension of 'Domain' of 'Curve'")

    IC1 <- new("IC")
    IC1@name <- name
    IC1@Curve <- Curve
    IC1@Risks <- Risks
    IC1@Infos <- Infos
    IC1@CallL2Fam <- CallL2Fam
    IC1@modifyIC <- modifyIC

    return(IC1)
}

# alias to generator function IC needed in functions makeIC in file CheckMakeIC.R
.IC <- IC


## access methods
setMethod("CallL2Fam", "IC", function(object) object@CallL2Fam)
setMethod("modifyIC", "IC", function(object) object@modifyIC)

## replace methods
setReplaceMethod("CallL2Fam", "IC",
    function(object, value){
        object@CallL2Fam <- value
        object
    })

## moved checkIC and makeIC methods in file CheckMakeIC.R in rev 1128

## evaluate IC
setMethod("evalIC", signature(IC = "IC", x = "numeric"), 
    function(IC, x){ 
        if(!is.null(IC@Curve[[1]]@Domain)){
            if(length(x) != IC@Curve[[1]]@Domain@dimension)
                stop("x has wrong dimension")
        }

        dimn <- dimension(IC@Curve)
        Curve <- as(diag(dimn) %*% IC@Curve, "EuclRandVariable")

        return(as.vector(evalRandVar(Curve, x)))
    })
setMethod("evalIC", signature(IC = "IC", x = "matrix"), 
    function(IC, x){ 
        if(!is.null(IC@Curve[[1]]@Domain)){
            if(ncol(x) != IC@Curve[[1]]@Domain@dimension)
                stop("x has wrong dimension")
        }

        dimn <- dimension(IC@Curve)
        Curve <- as(diag(dimn) %*% IC@Curve, "EuclRandVariable")

        if(dimn == 1)
            return(t(evalRandVar(Curve, x)[,,1]))
        else
            return(evalRandVar(Curve, x)[,,1])
    })

