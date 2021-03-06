### HampIC is only used internally; so no generating function exists;

## Access methods
setMethod("biastype", "HampIC", function(object) object@biastype)
setMethod("normtype", "HampIC", function(object) object@normtype)
setMethod("stand", "HampIC", function(object) object@stand)
setMethod("weight", "HampIC", function(object) object@weight)
setMethod("lowerCase", "HampIC", function(object) object@lowerCase)
setMethod("neighborRadius", "HampIC", function(object) object@neighborRadius)

setReplaceMethod("neighborRadius", "HampIC",
    function(object, value){
        object@neighborRadius <- value
        if(any(value < 0)) # radius vector?!
            stop("'value' has to be in [0, Inf]")
        addInfo(object) <- c("neighborRadius<-", "The slot 'neighborRadius' has been changed")
        addInfo(object) <- c("neighborRadius<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })

## evaluate IC
setMethod("evalIC.v", signature(IC = ".fastIC", x = "numeric"),
    function(IC, x){
        if(is.null(IC@.fastFct)){
           res <- setMethod("evalIC.v", signature(IC = "IC", x = "numeric"))(IC,x)
           ## cast to matrix ICdim x nobs
        }else{
           res <- IC@.fastFct(x)
        }
    })
setMethod("evalIC.v", signature(IC = ".fastIC", x = "matrix"),
    function(IC, x){
        if(is.null(IC@.fastFct)){
           res <- setMethod("evalIC.v", signature(IC = "IC", x = "matrix"))(IC,x)
           ## cast to matrix ICdim x nobs
        }else{
           res <- IC@.fastFct(x)
        }
    })
