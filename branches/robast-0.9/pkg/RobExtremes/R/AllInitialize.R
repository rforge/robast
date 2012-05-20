## initialize method 
setMethod("initialize", "Gumbel",
    function(.Object, loc = 0, scale = 1) {
        .Object@img <- Reals()
        .Object@param <- new("GumbelParameter", loc = loc, scale = scale, 
                             name = gettext("parameter of a Gumbel distribution"))
        .Object@r <- function(n){}
        body(.Object@r) <- substitute({ rgumbel(n, loc = loc1, scale = scale1) },
                                     list(loc1 = loc, scale1 = scale))
        .Object@d <- function(x, log = FALSE){}
        body(.Object@d) <- substitute({  dgumbel(x, loc = loc1, scale = scale1, log = log) },
                                     list(loc1 = loc, scale1 = scale))
        .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){}
        body(.Object@p) <- substitute({p1 <- pgumbel(q, loc = loc1, scale = scale1, lower.tail = lower.tail) 
                                       return(if(log.p) log(p1) else p1)},
                                     list(loc1 = loc, scale1 = scale))
        .Object@q <- function(p, loc = loc1, scale = scale1, lower.tail = TRUE, log.p = FALSE){}
            body(.Object@q) <- substitute({
                        ## P.R.: changed to vectorized form 
                        p1 <- if(log.p) exp(p) else p
                                                                        
                        in01 <- (p1>1 | p1<0)
                        i01 <- .isEqual01(p1) 
                        i0 <- (i01 & p1<1)   
                        i1 <- (i01 & p1>0)
                        ii01 <- .isEqual01(p1) | in01
                                      
                        p0 <- p
                        p0[ii01] <- if(log.p) log(0.5) else 0.5
                                      
                        q1 <- qgumbel(p0, loc = loc1, scale = scale1, 
                                      lower.tail = lower.tail) 
                        q1[i0] <- if(lower.tail) -Inf else Inf
                        q1[i1] <- if(!lower.tail) -Inf else Inf
                        q1[in01] <- NaN
                        
                        return(q1)  
                     },  list(loc1 = loc, scale1 = scale))
        .Object@.withSim   <- FALSE
        .Object@.withArith <- FALSE
        .Object@.logExact <- FALSE
        .Object@.lowerExact <- TRUE
        .Object
    })

## Class: Pareto distribution
setMethod("initialize", "Pareto",
          function(.Object, shape = 1, Min = 1, .withArith = FALSE) {
            .Object@img <- new("Reals")
            .Object@param <- new("ParetoParameter", shape = shape, Min =  Min)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
                           { rpareto1(n, shape = shapeSub,  min = MinSub) },
                             list(shapeSub = shape,  MinSub =  Min)
                                       )
            body(.Object@d) <- substitute(
                           { dpareto1(x, shape = shapeSub,  min =  MinSub, 
                                    log = log) },
                             list(shapeSub = shape,  MinSub =  Min)
                                         )
            body(.Object@p) <- substitute(
                           { ppareto1(q, shape = shapeSub,  min =  MinSub, 
                                    lower.tail = lower.tail, log.p = log.p) },
                             list(shapeSub = shape,  MinSub =  Min)
                                         )
            body(.Object@q) <- substitute({
                        ## P.R.: changed to vectorized form 
                        p1 <- if(log.p) exp(p) else p
                                                                        
                        in01 <- (p1>1 | p1<0)
                        i01 <- .isEqual01(p1) 
                        i0 <- (i01 & p1<1)   
                        i1 <- (i01 & p1>0)
                        ii01 <- .isEqual01(p1) | in01
                                      
                        p0 <- p
                        p0[ii01] <- if(log.p) log(0.5) else 0.5
                                      
                        q1 <- qpareto1(p0, shape = shapeSub,  min =  MinSub, 
                                    lower.tail = lower.tail, log.p = log.p) 
                        q1[i0] <- if(lower.tail) -Inf else Inf
                        q1[i1] <- if(!lower.tail) -Inf else Inf
                        q1[in01] <- NaN
                        
                        return(q1)  
                     },  list(shapeSub = shape,  MinSub =  Min))
            .Object@.withArith <- .withArith
            .Object@.logExact <- TRUE
            .Object@.lowerExact <- TRUE
            .Object
          })

## Class: Generalized Pareto distribution
setMethod("initialize", "GPareto",
          function(.Object, loc = 0, scale = 1, shape = 1) {
            .Object@img <- new("Reals")
            .Object@param <- new("GParetoParameter", loc = loc, scale = scale, shape = shape)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
                           { rgpd(n, loc = locSub, scale = scaleSub,  shape = shapeSub) },
                             list(locSub = loc, scaleSub = scale, shapeSub = shape)
                                       )
            body(.Object@d) <- substitute(
                           { dgpd(x, loc = locSub, scale = scaleSub, shape = shapeSub, 
                                    log = log) },
                             list(locSub = loc, scaleSub = scale, shapeSub = shape)
                                         )
            body(.Object@p) <- substitute(
                           { if(!lower.tail && log.p){
                             q0 <- (q-locSub)/scaleSub
                             return(-log(1+shapeSub*q0)/shapeSub)
                             }else{
                             p0 <- pgpd(q, loc = locSub, scale = scaleSub, 
                                        shape = shapeSub)
                             if(!lower.tail ) p0 <- 1-p0
                             if(log.p) p0 <- log(p0)
                             return(p0)}
                           }, list(locSub = loc, scaleSub = scale, 
                                   shapeSub = shape)
                                         )
            body(.Object@q) <- substitute({
                        if(!lower.tail && log.p){
                             p1 <- p
                             p1[p<.Machine$double.eps] <- 0.5
                             q0 <- (exp(-shapeSub*p1)-1)/shapeSub*scaleSub + locSub
                             q0[p<.Machine$double.eps] <- NaN
                             return(q0)
                        }else{
                             
                        ## P.R.: changed to vectorized form 
                           p1 <- if(log.p) exp(p) else p
                                                                        
                           in01 <- (p1>1 | p1<0)
                           i01 <- .isEqual01(p1) 
                           i0 <- (i01 & p1<1)   
                           i1 <- (i01 & p1>0)
                           ii01 <- .isEqual01(p1) | in01
                                      
                           p0 <- p
                           p0[ii01] <- if(log.p) log(0.5) else 0.5
                           if(!lower.tail) p0 <- 1-p0
                                      
                           q1 <- qgpd(p0, loc = locSub, scale = scaleSub, 
                                      shape = shapeSub) 
                           q1[i0] <- if(lower.tail)  locSub else Inf
                           q1[i1] <- if(!lower.tail) locSub else Inf
                           q1[in01] <- NaN
                        
                           return(q1) 
                         }   
                     },  list(locSub = loc, scaleSub = scale, shapeSub = shape))

            .Object@.withSim   <- FALSE
            .Object@.withArith <- FALSE
            .Object@.logExact <- TRUE
            .Object@.lowerExact <- TRUE
            .Object
          })


## Class: Generalized extreme value distribution
setMethod("initialize", "GEV",
          function(.Object, loc = 0, scale = 1, shape = 1) {
            .Object@img <- new("Reals")
            .Object@param <- new("GEVParameter", loc = loc, scale = scale, shape = shape)
            .Object@r <- function(n){}
            .Object@d <- function(x, log = FALSE){}
            .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){} 
            .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){} 
            body(.Object@r) <- substitute(
                           { rgev(n, loc = locSub, scale = scaleSub,  shape = shapeSub) },
                             list(locSub = loc, scaleSub = scale, shapeSub = shape)
                                       )
            body(.Object@d) <- substitute(
                           { dgev(x, loc = locSub, scale = scaleSub, shape = shapeSub, log = log) },
                             list(locSub = loc, scaleSub = scale, shapeSub = shape)
                                         )
            body(.Object@p) <- substitute(
                           { if(lower.tail && log.p){
                             q0 <- (q-locSub)/scaleSub
                             p0 <- -(1+shapeSub*q0)^(-1/shapeSub)
                             p0[q0<(-1)] <- -Inf 
                             return(p0)
                             }else{
                             p0 <- pgev(q, loc = locSub, scale = scaleSub, shape = shapeSub,lower.tail=TRUE)
                             if(!lower.tail ) p0 <- 1-p0
                             if(log.p) p0 <- log(p0)
                             return(p0)}
                           }, list(locSub = loc, scaleSub = scale, 
                                   shapeSub = shape)
                                         )
            body(.Object@q) <- substitute({
                        if(lower.tail && log.p){
                             q0 <-((-p)^(-shapeSub)-1)/shapeSub*scaleSub+locSub  
                             #q0[p>0|p< -Inf] <- NaN
                             #q0[.isEqual01(p)& p<1] <- Inf
                             #q0[!is.finite(p)& p<0] <- locSub-scaleSub/shapeSub                             
                             p0 <- exp(p)
                             q0[p0>1|p0<0] <- NaN
                             q0[(.isEqual01(p) & p0>0)] <- Inf
                             q0[(.isEqual01(p) & p0<1)] <- locSub-scaleSub/shapeSub 
                             return(q0)
                        }else{
                           ##higher tolerance for .isEqual01
                           tol=1e-20
                           distroptions(TruncQuantile=tol)
                           p1 <- if(log.p) exp(p) else p
                           in01 <- (p1>1 | p1<0)
                           i01 <- .isEqual01(p1) 
                           i0 <- (i01 & p1<1)   
                           i1 <- (i01 & p1>0)
                           ii01 <- .isEqual01(p1) | in01
                           p0 <- p
                           p0[ii01] <- if(log.p) log(0.5) else 0.5
                           #if(!lower.tail) p0 <- 1-p0
                           q1 <- qgev(p0, loc = locSub, scale = scaleSub, shape = shapeSub, lower.tail=lower.tail) 
                           q1[i0] <- if(lower.tail)  locSub-scaleSub/shapeSub else Inf
                           q1[i1] <- if(!lower.tail) locSub-scaleSub/shapeSub else Inf
                           q1[in01] <- NaN
                           return(q1) 
                         }   
                     },  list(locSub = loc, scaleSub = scale, shapeSub = shape))

            .Object@.withSim   <- FALSE
            .Object@.withArith <- FALSE
            .Object@.logExact <- TRUE
            .Object@.lowerExact <- TRUE
            .Object
          })

