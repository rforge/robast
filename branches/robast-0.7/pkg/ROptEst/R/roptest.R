###############################################################################
## Optimally robust estimation
###############################################################################
roptest <- function(x, L2Fam, eps, eps.lower, eps.upper, initial.est, 
                    neighbor = ContNeighborhood(), risk = asMSE(), steps = 1, 
                    distance = CvMDist, startPar = NULL, verbose = FALSE, 
                    useLast = getRobAStBaseOption("kStepUseLast"), ...){
    es.call <- match.call()
    if(missing(x))
        stop("'x' is missing with no default")
    if(missing(L2Fam))
        stop("'L2Fam' is missing with no default")
    if(!is.numeric(x)){
        if(is.data.frame(x))
            x <- data.matrix(x)
        else
            x <- as.matrix(x)
        if(!is.matrix(x))
            stop("'x' has to be a numeric vector resp. a matrix or data.frame")
    }
    if(missing(eps) && missing(eps.lower) && missing(eps.upper)){
        eps.lower <- 0
        eps.upper <- 0.5
    }
    if(missing(eps)){
        if(!missing(eps.lower) && missing(eps.upper))
            eps.upper <- 0.5
        if(missing(eps.lower) && !missing(eps.upper))
            eps.lower <- 0
        if(length(eps.lower) != 1 || length(eps.upper) != 1)
            stop("'eps.lower' and 'eps.upper' have to be of length 1")
        if(!is.numeric(eps.lower) || !is.numeric(eps.upper) || eps.lower >= eps.upper) 
            stop("'eps.lower' < 'eps.upper' is not fulfilled")
        if((eps.lower < 0) || (eps.upper > 0.5))
            stop("'eps.lower' and 'eps.upper' have to be in [0, 0.5]")
    }else{
        if(length(eps) != 1)
            stop("'eps' has to be of length 1")
        if(eps == 0)
            stop("'eps = 0'! => use functions 'mean' and 'sd' for estimation")
        if((eps < 0) || (eps > 0.5))
            stop("'eps' has to be in (0, 0.5]")
    }
    if(!is.integer(steps))
        steps <- as.integer(steps)
    if(steps < 1){
        stop("'steps' has to be some positive integer value")
    }
    if(length(steps) != 1){
        stop("'steps' has to be of length 1")
    }

    if(missing(initial.est))
        initial.est <- MDEstimator(x = x, ParamFamily = L2Fam, distance = distance,
                                            startPar = startPar, ...)
    if(is(initial.est, "Estimate")) initial.est <- estimate(initial.est)
    newParam <- param(L2Fam)
    main(newParam)[] <- initial.est
    L2FamStart <- modifyModel(L2Fam, newParam)
    if(is.matrix(x))
        sqrtn <- sqrt(ncol(x))
    else
        sqrtn <- sqrt(length(x))
    if(missing(eps)){
        r.lower <- sqrtn*eps.lower
        r.upper <- sqrtn*eps.upper
        ICstart <- radiusMinimaxIC(L2Fam = L2FamStart, neighbor = neighbor, risk = risk, 
                                   loRad = r.lower, upRad = r.upper, verbose = verbose)
    }else{
        r <- sqrtn*eps
        neighbor@radius <- r
        infMod <- InfRobModel(center = L2FamStart, neighbor = neighbor)
        ICstart <- optIC(model = infMod, risk = risk, verbose = verbose)
    }
    res <- kStepEstimator(x, IC = ICstart, start = initial.est, steps = steps, useLast = useLast)
    res@estimate.call <- es.call
    Infos <- matrix(c("roptest", 
                      paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
                    ncol = 2)
    colnames(Infos) <- c("method", "message")
    Infos <- rbind(Infos, c("roptest", 
                            paste("computation of IC, asvar and asbias via useLast =", useLast)))
    Infos(res) <- Infos
    return(res)
}
