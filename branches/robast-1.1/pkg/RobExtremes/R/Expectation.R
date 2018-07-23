
setMethod("E", signature(object = "Pareto", 
                         fun = "missing", 
                         cond = "missing"),
    function(object, low = NULL, upp = NULL, ...){
    if(!is.null(low)) if(low <= Min(object)) low <- NULL
    a <- shape(object); b <- Min(object)
    if(is.null(low) && is.null(upp)){
        if(a<=1) return(Inf)
        else return(b*a/(a-1))
     }   
    else
        return(E(as(object,"AbscontDistribution"), low=low, upp=upp, ...))    
    })

### source http://mathworld.wolfram.com/ParetoDistribution.html


setMethod("E", signature(object = "Gumbel", 
                         fun = "missing", 
                         cond = "missing"),
    function(object, low = NULL, upp = NULL, ...){a <- loc(object); b <- scale(object)
    if(is.null(low) && is.null(upp))
           return(a- EULERMASCHERONICONSTANT * b)
    else
        return(E(as(object,"AbscontDistribution"), low=low, upp=upp, ...))    
    })
## http://mathworld.wolfram.com/GumbelDistribution.html

setMethod("E", signature(object = "GPareto", 
                         fun = "missing", 
                         cond = "missing"),
    function(object, low = NULL, upp = NULL, ...){
    if(!is.null(low)) if(low <= Min(object)) low <- NULL
    k <- shape(object); s <- scale(object); mu <- loc(object)
    if(is.null(low) && is.null(upp)){
        if(k>=1) return(Inf)
        else return(mu+s/(1-k))
     }   
    else
        return(E(as(object,"AbscontDistribution"), low=low, upp=upp, ...))    
    })

### source http://en.wikipedia.org/wiki/Pareto_distribution

setMethod("E", signature(object = "DistributionsIntegratingByQuantiles",
                         fun = "function",
                         cond = "missing"),
    function(object, fun, low = NULL, upp = NULL,
             rel.tol= getdistrExOption("ErelativeTolerance"),
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"),
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"),
             IQR.fac = max(1e4,getdistrExOption("IQR.fac")), ...
             ){

        dots <- list(...)
        dots.withoutUseApply <- dots
        useApply <- TRUE
        if(!is.null(dots$useApply)) useApply <- dots$useApply

        dots.withoutUseApply$useApply <- NULL
        dots.withoutUseApply$stop.on.error <- NULL

        integrand <- function(x, dfun, ...){   di <- dim(x)
                                               y <- q.l(object)(x)##quantile transformation
                                               if(useApply){
                                                    funy <- sapply(y,fun, ...)
                                                    dim(y) <- di
                                                    dim(funy) <- di
                                               }else funy <- fun(y,...)
                                        return(funy) }

         if(is.null(low)) low <- -Inf
         if(is.null(upp)) upp <- Inf

         Ib <- .getIntbounds(object, low, upp, lowerTruncQuantile,
               upperTruncQuantile, IQR.fac)
         low <- p(object)(Ib["low"])
         upp <- p(object)(Ib["upp"])
         if(is.nan(low)) low <- 0
         if(is.nan(upp)) upp <- 1

         if(upp < 0.98){
           int <- do.call(distrExIntegrate, c(list(f = integrand,
                    lower = low,
                    upper = upp,
                    rel.tol = rel.tol, stop.on.error = FALSE,
                    distr = object, dfun = dunif), dots.withoutUseApply))
         }else{
           int1 <- do.call(distrExIntegrate, c(list(f = integrand,
                    lower = low,
                    upper = 0.98,
                    rel.tol = rel.tol, stop.on.error = FALSE,
                    distr = object, dfun = dunif), dots.withoutUseApply))
           int2 <- do.call(distrExIntegrate, c(list(f = integrand,
                    lower = 0.98,
                    upper = upp,
                    rel.tol = rel.tol, stop.on.error = FALSE,
                    distr = object, dfun = dunif), dots.withoutUseApply))
           int <- int1+int2
         }

         return(int)

    })

setMethod("E", signature(object = "GPareto",
                         fun = "function",
                         cond = "missing"),
    function(object, fun, low = NULL, upp = NULL,
             rel.tol= getdistrExOption("ErelativeTolerance"),
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"),
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"),
             IQR.fac = max(1e4,getdistrExOption("IQR.fac")), ...
             ){

        dots <- list(...)
        dots.withoutUseApply <- dots
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
                    distr = object, dfun = d(object)), dots.withoutUseApply)))

    })


setMethod("E", signature(object = "GEV",
                         fun = "missing", 
                         cond = "missing"),
    function(object, low = NULL, upp = NULL, ...){
    if(!is.null(low)) if(low <= Min(object)) low <- NULL
    xi <- shape(object); sigma <- scale(object); mu <- loc(object)
    if(is.null(low) && is.null(upp)){
        if (xi==0) return(mu+sigma*EULERMASCHERONICONSTANT)
        else if(xi>=1) return(Inf)
        else return(mu+sigma*(gamma(1-xi)-1)/xi)
        }       
    else
        return(E(object, low=low, upp=upp, fun = function(x)x, ...))
    })

setMethod("E", signature(object = "GEV", fun = "function", cond = "missing"),
           getMethod("E",
           signature(object = "DistributionsIntegratingByQuantiles",
                     fun = "function", cond = "missing")))

setMethod("E", signature(object = "Weibull", fun = "function", cond = "missing"),
           getMethod("E",
           signature(object = "DistributionsIntegratingByQuantiles",
                     fun = "function", cond = "missing")))

setMethod("E", signature(object = "Gammad", fun = "function", cond = "missing"),
           getMethod("E",
           signature(object = "DistributionsIntegratingByQuantiles",
                     fun = "function", cond = "missing")))
