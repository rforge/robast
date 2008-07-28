###############################################################################
## one-step estimator
###############################################################################
setMethod("oneStepEstimator", signature(x = "numeric", 
                                        IC = "InfluenceCurve",
                                        start = "numeric"),
    function(x, IC, start){
        nrvalues <- dimension(IC@Curve)
        if(is.list(start)) start <- unlist(start)
        if(nrvalues != length(start))
            stop("dimension of 'start' != dimension of 'Curve'")

        res <- start + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)

        return(res)
    })
setMethod("oneStepEstimator", signature(x = "matrix", 
                                        IC = "InfluenceCurve",
                                        start = "numeric"),
    function(x, IC, start){
        nrvalues <- dimension(IC@Curve)
        if(is.list(start)) start <- unlist(start)
        if(nrvalues != length(start))
            stop("dimension of 'start' != dimension of 'Curve'")
        if(ncol(x) != IC@Curve[[1]]@Domain@dimension)
            stop("'x' has wrong dimension")

        res <- start + rowMeans(evalIC(IC, x), na.rm = TRUE)

        return(res)
    })
setMethod("oneStepEstimator", signature(x = "numeric", 
                                        IC = "InfluenceCurve",
                                        start = "Estimate"),
    function(x, IC, start){
        nrvalues <- dimension(IC@Curve)
        start0 <- estimate(start)
        if(is.list(start0)) start0 <- unlist(start0)
        if(nrvalues != length(start0))
            stop("dimension of slot 'estimate' of 'start' != dimension of 'Curve'")

        res <- start0 + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)

        return(res)
    })
setMethod("oneStepEstimator", signature(x = "matrix", 
                                        IC = "InfluenceCurve",
                                        start = "Estimate"),
    function(x, IC, start){
        nrvalues <- dimension(IC@Curve)
        start0 <- estimate(start)
        if(is.list(start0)) start0 <- unlist(start0)
        if(nrvalues != length(start0))
            stop("dimension of slot 'estimate' of 'start' != dimension of 'Curve'")
        if(ncol(x) != IC@Curve[[1]]@Domain@dimension)
            stop("'x' has wrong dimension")

        res <- start0 + rowMeans(evalIC(IC, x), na.rm = TRUE)

        return(res)
    })
