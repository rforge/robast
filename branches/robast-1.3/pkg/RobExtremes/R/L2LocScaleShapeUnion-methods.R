setMethod("locscaleshapename", signature(object = "L2LocScaleShapeUnion"),
           function(object) object@locscaleshapename)
setMethod("locscalename", signature(object = "L2LocScaleShapeUnion"),
           function(object) object@locscaleshapename[c("location","scale")])

setMethod("scaleshapename", signature(object = "L2LocScaleShapeUnion"),
           function(object) object@locscaleshapename[c("scale","shape")])

setMethod("scalename", signature(object = "L2LocScaleShapeUnion"),
           function(object) object@locscaleshapename["scale"])

setMethod("shapename", signature(object = "L2LocScaleShapeUnion"),
           function(object) object@scaleshapename["shape"])

setMethod("locationname", signature(object = "L2LocScaleShapeUnion"),
           function(object) object@locscaleshapename["location"])


setReplaceMethod("locscaleshapename", "L2LocScaleShapeUnion",
    function(object, value){
        if(length(value)!=3)
           stop("value of slot 'locscaleshapename' must be of length three")
        if(is.null(names(value))) names(value) <- c("location","scale","shape")
        object@locscalename <- value
        object
    })

