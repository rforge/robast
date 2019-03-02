## copied form distrEx from distrEx 2.8.0 and branch 1.2.0 on

## .qtlIntegrate is moved from RobExtremes (slightly modified) to distrEx
#   as of versions distrEx 2.8.0 and RobExtremes 1.2.0


setMethod("E", signature(object = "Pareto",
                         fun = "missing",
                         cond = "missing"),
    function(object, low = NULL, upp = NULL, ..., diagnostic = FALSE){
    if(!is.null(low)) if(low <= Min(object)) low <- NULL
    a <- shape(object); b <- Min(object)
    if(is.null(low) && is.null(upp)){
        if(a<=1) return(Inf)
        else return(b*a/(a-1))
     }
    else
        return(E(object=object,fun=function(x)x, low=low, upp=upp, ...,
                    diagnostic = diagnostic))
    })

### source http://mathworld.wolfram.com/ParetoDistribution.html


setMethod("E", signature(object = "Gumbel",
                         fun = "missing",
                         cond = "missing"),
    function(object, low = NULL, upp = NULL, ..., diagnostic = FALSE){
    a <- loc(object); b <- scale(object)
    if(is.null(low) && is.null(upp))
           return(a- EULERMASCHERONICONSTANT * b)
    else
        return(E(object=object,fun=function(x)x, low=low, upp=upp, ...,
                    diagnostic = diagnostic))
    })
## http://mathworld.wolfram.com/GumbelDistribution.html

setMethod("E", signature(object = "GPareto",
                         fun = "missing",
                         cond = "missing"),
    function(object, low = NULL, upp = NULL, ..., diagnostic = FALSE){
    if(!is.null(low)) if(low <= Min(object)) low <- NULL
    k <- shape(object); s <- scale(object); mu <- loc(object)
    if(is.null(low) && is.null(upp)){
        if(k>=1) return(Inf)
        else return(mu+s/(1-k))
     }
    else
        return(E(object=object,fun=function(x)x, low=low, upp=upp, ...,
                    diagnostic = diagnostic))
    })

### source http://en.wikipedia.org/wiki/Pareto_distribution

setMethod("E", signature(object = "DistributionsIntegratingByQuantiles",
                         fun = "function",
                         cond = "missing"),
    function(object, fun, low = NULL, upp = NULL,
             rel.tol= getdistrExOption("ErelativeTolerance"),
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"),
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"),
             IQR.fac = max(1e4,getdistrExOption("IQR.fac")), ...,
             diagnostic = FALSE){

     dots <- list(...)
     dotsI <- .filterEargs(dots)
     dotsFun <- .filterFunargs(dots,fun)
     funwD <- function(x) do.call(fun, c(list(x=x),dotsFun))

     do.call(.qtlIntegrate, c(list(object = object, fun = funwD, low = low, upp = upp,
             rel.tol= rel.tol, lowerTruncQuantile = lowerTruncQuantile,
             upperTruncQuantile = upperTruncQuantile,
             IQR.fac = IQR.fac, ...,
             .withLeftTail = FALSE, .withRightTail = TRUE,
             diagnostic = diagnostic),dotsI))
    })

setMethod("E", signature(object = "GPareto",
                         fun = "function",
                         cond = "missing"),
    function(object, fun, low = NULL, upp = NULL,
             rel.tol= getdistrExOption("ErelativeTolerance"),
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"),
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"),
             IQR.fac = max(1e4,getdistrExOption("IQR.fac")), ...,
             diagnostic = FALSE){

        dots <- list(...)
        dots.withoutUseApply <- .filterEargs(dots)
        useApply <- TRUE
        if(!is.null(dots$useApply)) useApply <- dots$useApply
        dots.withoutUseApply$useApply <- NULL
        integrand <- function(x, dfun, ...){   di <- dim(x)
                                               y <- exp(x)
                                               if(useApply){
                                                    funy <- sapply(y,fun, ...)
                                                    dim(y) <- di
                                                    dim(funy) <- di
                                               }else funy <- fun(y,...)
                                        return(funy * y * dfun(y)) }

        if(is.null(low)) low <- -Inf
        if(is.null(upp)) upp <- Inf

        Ib <- .getIntbounds(object, low, upp, lowerTruncQuantile,
              upperTruncQuantile, IQR.fac)
        low <- if(Ib["low"]<=0) -Inf else log(Ib["low"])
        upp <- log(Ib["upp"])

        return(do.call(distrExIntegrate, c(list(f = integrand,
                    lower = low,
                    upper = upp,
                    rel.tol = rel.tol,
                    distr = object, dfun = d(object)), dots.withoutUseApply,
                    diagnostic = diagnostic)))

    })


setMethod("E", signature(object = "GEV",
                         fun = "missing",
                         cond = "missing"),
    function(object, low = NULL, upp = NULL, ..., diagnostic = FALSE){
    if(!is.null(low)) if(low <= Min(object)) low <- NULL
    xi <- shape(object); sigma <- scale(object); mu <- loc(object)
    if(is.null(low) && is.null(upp)){
        if (xi==0) return(mu+sigma*EULERMASCHERONICONSTANT)
        else if(xi>=1) return(Inf)
        else return(mu+sigma*(gamma(1-xi)-1)/xi)
        }
    else
        return(E(object, low=low, upp=upp, fun = function(x)x, ...,
                 diagnostic = diagnostic))
    })

setMethod("E", signature(object = "GEV", fun = "function", cond = "missing"),
           getMethod("E",
           signature(object = "DistributionsIntegratingByQuantiles",
                     fun = "function", cond = "missing")))

## these routines are moved back to package distrEx from distrEx 2.8.0 / RobExtremes 1.2.0 on

#setMethod("E", signature(object = "Weibull", fun = "function", cond = "missing"),
#           getMethod("E",
#           signature(object = "DistributionsIntegratingByQuantiles",
#                     fun = "function", cond = "missing")))

#setMethod("E", signature(object = "Gammad", fun = "function", cond = "missing"),
#           getMethod("E",
#           signature(object = "DistributionsIntegratingByQuantiles",
#                     fun = "function", cond = "missing")))
