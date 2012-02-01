### CondHampIC is only used internally; so no generating function exists;
## Access methods
setMethod("neighborRadiusCurve", "CondHampIC", function(object) object@neighborRadiusCurve)

## replace methods
setReplaceMethod("neighborRadius", "CondHampIC",
    function(object, value){ 
        object@neighborRadius <- value
        if(any(value < 0)) # radius vector?!
            stop("'value' has to be in [0, Inf]")
        addInfo(object) <- c("neighborRadius<-", "The slot 'neighborRadius' has been changed")
        addInfo(object) <- c("neighborRadius<-", "The entries in 'Risks' and 'Infos' may be wrong")
        object
    })
setReplaceMethod("neighborRadiusCurve", "CondHampIC",
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
