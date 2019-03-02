.locEfunLoop <- function(object, fun, useApply = TRUE, ..., diagnostic = FALSE){
        dimn <- length(fun)
        nrdim <- fun@Range@dimension
        res <- numeric(dimn)
        if(nrdim > 1)
           res <- matrix(0, nrow = dimn, ncol = nrdim)
        diagn <- NULL
        if(diagnostic)  diagn <- vector("list", dimn)
        for(i in 1:dimn){
               buf <- E(object, fun = Map(fun)[[i]], useApply = useApply, ..., diagnostic = diagnostic)
               if(diagnostic) diagn[[i]] <- attr(buf, "diagnostic")
               if(nrdim>1) res[i,] <- buf else res[i] <- buf
        }
        if(!is.null(diagn)){
           attr(res,"diagnostic") <- diagn
           class(attr(res,"diagnostic")) <- "DiagnosticClass"
        }
        return(res)
    }

.locEfunLoopCond <-     function(object, fun, cond, withCond = FALSE, useApply = TRUE, ..., diagnostic = FALSE){
        dimn <- length(fun)
        res <- numeric(dimn)
        diagn <- if(diagnostic) vector("list", dimn) else NULL
        dots <- list(...)
        dotsI <- .filterEargs(dots)
        Eargs0 <- list(object=object)
        Eargs1 <- list(cond=cond, withCond=withCond, useApply = useApply, diagnostic = diagnostic)

        for(i in 1:dimn){
            dotsFun  <- .filterFunargs(dots, fun@Map[[i]])

            funwD <- function(x)  do.call(fun@Map[[i]], c(list(x), eval.parent(dotsFun,1)))
            funwDc <- function(x,cond){ y <- c(x,cond);  do.call(fun@Map[[i]], c(list(x=y), eval.parent(dotsFun,1)))}

            Eargs <- c(Eargs0, list(fun=if(withCond)funwDc else funwD), Eargs1, dotsI)
            res[i] <- buf <- do.call(E, Eargs)
            if(diagnostic) diagn[[i]] <- attr(buf, "diagnostic")
        }
        if(!is.null(diagn)){
           attr(res,"diagnostic") <- diagn
           class(attr(res,"diagnostic")) <- "DiagnosticClass"
        }
        return(res)
   }

.locEfun <- function(object, fun, useApply = TRUE, ..., diagnostic = FALSE){
        if(!is(fun@Domain, "EuclideanSpace"))
            stop("'Domain' of the random variable is no Euclidean space")
        if(dimension(fun@Domain) != 1)
            stop("dimension of 'Domain' of the random variable has to be 1")
        if(dimension(fun@Range) != 1)
            stop("dimension of 'Range' of the random variable has to be 1")
        .locEfunLoop(object = object, fun = fun, useApply = useApply, ..., diagnostic = diagnostic)
    }

.locEmatfun <- function(object, fun, useApply = TRUE, ..., diagnostic = FALSE){
        diagn <- NULL
        res <- E(object, as(fun, "EuclRandVariable"), useApply = useApply, ..., diagnostic = diagnostic)
        if(diagnostic) diagn <- attr(res, "diagnostic")
        res <- matrix(res, nrow = nrow(fun))
        if(!is.null(diagn)){
           attr(res,"diagnostic") <- diagn
           class(attr(res,"diagnostic")) <- "DiagnosticClass"
        }
        return(res)
    }
.locElistfun <- function(object, fun, useApply = TRUE, ..., diagnostic = FALSE){
        nrvalues <- length(fun)
        res <- vector("list", nrvalues)
        diagn <- NULL
        if(diagnostic)  diagn <- vector("list", nrvalues)
        for(i in 1:nrvalues){
#               print(list(object, fun = fun[[i]], useApply = useApply, ..., diagnostic = diagnostic))
               res[[i]] <- buf <- E(object, fun = fun[[i]], useApply = useApply, ..., diagnostic = diagnostic)
               if(diagnostic) diagn[[i]] <- attr(buf, "diagnostic")
        }
        if(!is.null(diagn)){
           attr(res,"diagnostic") <- diagn
           class(attr(res,"diagnostic")) <- "DiagnosticClass"
        }
        return(res)
    }

.locEMVfun <- function(object, fun, useApply = TRUE, ..., diagnostic = FALSE){
#        print(list(object, fun, useApply, ..., diagnostic))
        if(!is(fun@Domain, "EuclideanSpace"))
            stop("'Domain' of the random variable is no Euclidean space")
        if(fun@Domain@dimension != object@img@dimension)
            stop("dimension of 'Domain' of the random variable is not equal\n",
                 "to dimension of 'img' of the distribution")
        res <- .locEfunLoop(object = object, fun = fun, useApply = useApply, ..., diagnostic = diagnostic)
        dim(res) <- c(length(fun),fun@Range@dimension)
        return(res)
    }


.locEfunCond <-
    function(object, fun, cond, withCond = FALSE, useApply = TRUE, ..., diagnostic = FALSE){
        if(!is(fun@Domain, "EuclideanSpace"))
            stop("'Domain' of the random variable has to be a Euclidean Space")
        if(withCond){
            if(fun@Domain@dimension != (1+length(cond)))
                stop("wrong dimension of 'Domain' of 'fun'")
        }else{
            if(fun@Domain@dimension != 1)
                stop("dimension of 'Domain' of 'fun' has to be 1")
        }
        if(dimension(fun@Range) != 1)
            stop("dimension of 'Range' of the random variable has to be 1")

       return(.locEfunLoopCond(object = object, fun = fun, cond = cond, withCond = withCond,
                        useApply = useApply, ..., diagnostic = diagnostic))
    }

.locEmatfunCond <- function(object, fun, cond, withCond = FALSE, useApply = TRUE, ..., diagnostic = FALSE){
        diagn <- NULL
        res <- E(object, as(fun, "EuclRandVariable"), cond = cond,
                 withCond = withCond, useApply = useApply, ..., diagnostic = diagnostic)
        if(diagnostic) diagn <- attr(res, "diagnostic")
        res <- matrix(res, nrow = nrow(fun))
        if(!is.null(diagn)){
           attr(res,"diagnostic") <- diagn
           class(attr(res,"diagnostic")) <- "DiagnosticClass"
        }
        return(res)
    }

.locElistfunCond <- function(object, fun, cond, withCond = FALSE, useApply = TRUE, ..., diagnostic = FALSE){
        nrvalues <- length(fun)
        diagn <- if(diagnostic) vector("list", nrvalues) else NULL
        res <- vector("list", nrvalues)
        for(i in 1:nrvalues){
            res[[i]] <- buf <- E(object, fun=fun[[i]], cond = cond, withCond = withCond, useApply = useApply, ..., diagnostic = diagnostic)
            if(diagnostic) diagn[[i]] <- attr(buf, "diagnostic")
        }
        if(!is.null(diagn)){
           attr(res,"diagnostic") <- diagn
           class(attr(res,"diagnostic")) <- "DiagnosticClass"
        }
        return(res)
    }



setMethod("E", signature(object = "UnivariateDistribution", 
                         fun = "EuclRandVariable", 
                         cond = "missing"),
          function(object, fun, useApply = TRUE, ..., diagnostic = FALSE){
             mc <- match.call()
             dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
             res <- do.call(.locEfun, c(list(object=object, fun= fun, useApply = useApply,
                             diagnostic = diagnostic), dotsI))
             if(diagnostic){
                diagn <- attr(res,"diagnostic")
                diagn[["call"]] <- mc
                attr(res,"diagnostic") <- diagn
                class(attr(res,"diagnostic")) <- "DiagnosticClass"
             }
             return(res)
          })
setMethod("E", signature(object = "AbscontDistribution", 
                         fun = "EuclRandVariable", 
                         cond = "missing"),
          function(object, fun, useApply = TRUE, ..., diagnostic = FALSE){
             dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
             res <- do.call(.locEfun, c(list(object=object, fun= fun, useApply = useApply,
                             diagnostic = diagnostic), dotsI))
             if(diagnostic){
                diagn <- attr(res,"diagnostic")
                diagn[["call"]] <- match.call()
                attr(res,"diagnostic") <- diagn
                class(attr(res,"diagnostic")) <- "DiagnosticClass"
             }
             return(res)
          })
setMethod("E", signature(object = "DiscreteDistribution",
                         fun = "EuclRandVariable", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
     dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
     do.call(.locEfun,c(list(object = object, fun = fun, useApply = useApply, diagnostic= FALSE), dotsI))
     })

setMethod("E", signature(object = "UnivariateDistribution", 
                         fun = "EuclRandMatrix", 
                         cond = "missing"),
          function(object, fun, useApply = TRUE, ..., diagnostic = FALSE){
             dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
             res <- do.call(.locEmatfun, c(list(object=object, fun= fun, useApply = useApply,
                             diagnostic = diagnostic), dotsI))
             if(diagnostic){
                diagn <- attr(res,"diagnostic")
                diagn[["call"]] <- match.call()
                attr(res,"diagnostic") <- diagn
                class(attr(res,"diagnostic")) <- "DiagnosticClass"
             }
             return(res)
          })

setMethod("E", signature(object = "AbscontDistribution", 
                         fun = "EuclRandMatrix", 
                         cond = "missing"),
          function(object, fun, useApply = TRUE, ..., diagnostic = FALSE){
             dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
             res <- do.call(.locEmatfun, c(list(object=object, fun= fun, useApply = useApply,
                             diagnostic = diagnostic), dotsI))
             if(diagnostic){
                diagn <- attr(res,"diagnostic")
                diagn[["call"]] <- match.call()
                attr(res,"diagnostic") <- diagn
                class(attr(res,"diagnostic")) <- "DiagnosticClass"
             }
             return(res)
          })

setMethod("E", signature(object = "DiscreteDistribution",
                         fun = "EuclRandMatrix", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
     dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
     do.call(.locEmatfun,c(list(object = object, fun = fun, useApply = useApply, diagnostic= FALSE), dotsI))
    })

setMethod("E", signature(object = "UnivariateDistribution", 
                         fun = "EuclRandVarList", 
                         cond = "missing"),
          function(object, fun, useApply = TRUE, ..., diagnostic = FALSE){
             dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
             res <- do.call(.locElistfun, c(list(object=object, fun= fun, useApply = useApply,
                             diagnostic = diagnostic), dotsI))
             if(diagnostic){
                diagn <- attr(res,"diagnostic")
                diagn[["call"]] <- match.call()
                attr(res,"diagnostic") <- diagn
                class(attr(res,"diagnostic")) <- "DiagnosticClass"
             }
             return(res)
          })

setMethod("E", signature(object = "AbscontDistribution", 
                         fun = "EuclRandVarList", 
                         cond = "missing"),
          function(object, fun, useApply = TRUE, ..., diagnostic = FALSE){
             dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
             res <- do.call(.locElistfun, c(list(object=object, fun= fun, useApply = useApply,
                             diagnostic = diagnostic), dotsI))
             if(diagnostic){
                diagn <- attr(res,"diagnostic")
                diagn[["call"]] <- match.call()
                attr(res,"diagnostic") <- diagn
                class(attr(res,"diagnostic")) <- "DiagnosticClass"
             }
             return(res)
          })

setMethod("E", signature(object = "DiscreteDistribution",
                         fun = "EuclRandVarList", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
     dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
     do.call(.locElistfun,c(list(object = object, fun = fun, useApply = useApply, diagnostic= FALSE), dotsI))
     })

setMethod("E", signature(object = "MultivariateDistribution", 
                         fun = "EuclRandVariable", 
                         cond = "missing"),
          function(object, fun, useApply = TRUE, ..., diagnostic = FALSE){
             dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
             res <- do.call(.locEMVfun, c(list(object=object, fun= fun, useApply = useApply,
                             diagnostic = diagnostic), dotsI))
             if(diagnostic){
                diagn <- attr(res,"diagnostic")
                diagn[["call"]] <- match.call()
                attr(res,"diagnostic") <- diagn
                class(attr(res,"diagnostic")) <- "DiagnosticClass"
             }
             return(res)
          })


setMethod("E", signature(object = "DiscreteMVDistribution", 
                         fun = "EuclRandVariable", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
     dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
     do.call(.locEMVfun,c(list(object = object, fun = fun, useApply = useApply, diagnostic= FALSE), dotsI))
     })

setMethod("E", signature(object = "MultivariateDistribution",
                         fun = "EuclRandMatrix", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ..., diagnostic = FALSE){
        if(!diagnostic){
           return(array(E(object, as(fun, "EuclRandVariable"), useApply = useApply, ..., diagnostic = diagnostic),
              c(nrow(fun), ncol(fun), fun@Range@dimension)))
        }else{
           res <- E(object, as(fun, "EuclRandVariable"), useApply = useApply, ..., diagnostic = diagnostic)
           diagn <- attr(res,"diagnostic")
           diagn[["call"]] <- match.call()
           res <- array(res, c(nrow(fun), ncol(fun), fun@Range@dimension))
           attr(res, "diagnostic") <- diagn
           class(attr(res,"diagnostic")) <- "DiagnosticClass"
           return(res)
        }
    })
setMethod("E", signature(object = "DiscreteMVDistribution", 
                         fun = "EuclRandMatrix", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
        array(E(object, as(fun, "EuclRandVariable"), useApply = useApply, ...),
              c(nrow(fun), ncol(fun), fun@Range@dimension))
    })
setMethod("E", signature(object = "MultivariateDistribution", 
                         fun = "EuclRandVarList", 
                         cond = "missing"),
          function(object, fun, useApply = TRUE, ..., diagnostic = FALSE){
             dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
             res <- do.call(.locElistfun, c(list(object=object, fun= fun, useApply = useApply,
                             diagnostic = diagnostic), dotsI))
             if(diagnostic){
                diagn <- attr(res,"diagnostic")
                diagn[["call"]] <- match.call()
                attr(res,"diagnostic") <- diagn
                class(attr(res,"diagnostic")) <- "DiagnosticClass"
             }
             return(res)
          })

setMethod("E", signature(object = "DiscreteMVDistribution",
                         fun = "EuclRandVarList", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
     dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
     do.call(.locElistfun,c(list(object = object, fun = fun, useApply = useApply, diagnostic= FALSE), dotsI))
     })

setMethod("E", signature(object = "UnivariateCondDistribution", 
                         fun = "EuclRandVariable", 
                         cond = "numeric"),
          function(object, fun, cond, withCond = FALSE, useApply = TRUE, ..., diagnostic = FALSE){
             dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
             res <- do.call(.locEfunCond, c(list(object=object, fun= fun, cond=cond,
                                 withCond = withCond, useApply = useApply,
                                 diagnostic = diagnostic), dotsI))
             if(diagnostic){
                diagn <- attr(res,"diagnostic")
                diagn[["call"]] <- match.call()
                attr(res,"diagnostic") <- diagn
                class(attr(res,"diagnostic")) <- "DiagnosticClass"
             }
             return(res)
          })

setMethod("E", signature(object = "AbscontCondDistribution",
                         fun = "EuclRandVariable", 
                         cond = "numeric"),
          function(object, fun, cond, withCond = FALSE, useApply = TRUE, ..., diagnostic = FALSE){
             dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
             res <- do.call(.locEfunCond, c(list(object=object, fun= fun, cond=cond,
                                 withCond = withCond, useApply = useApply,
                                 diagnostic = diagnostic), dotsI))
             if(diagnostic){
                diagn <- attr(res,"diagnostic")
                diagn[["call"]] <- match.call()
                attr(res,"diagnostic") <- diagn
                class(attr(res,"diagnostic")) <- "DiagnosticClass"
             }
             return(res)
          })

setMethod("E", signature(object = "DiscreteCondDistribution",
                         fun = "EuclRandVariable", 
                         cond = "numeric"),
    function(object, fun, cond, withCond = FALSE, useApply = TRUE, ...){
     dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
     do.call(.locEfunCond,c(list(object = object, fun = fun, cond=cond, withCond = withCond,
               useApply = useApply, diagnostic= FALSE), dotsI))
     })

setMethod("E", signature(object = "UnivariateCondDistribution",
                         fun = "EuclRandMatrix", 
                         cond = "numeric"),
          function(object, fun, cond, withCond = FALSE, useApply = TRUE, ..., diagnostic = FALSE){
             dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
             res <- do.call(.locEmatfunCond, c(list(object=object, fun= fun, cond=cond,
                                 withCond = withCond, useApply = useApply,
                                 diagnostic = diagnostic), dotsI))
             if(diagnostic){
                diagn <- attr(res,"diagnostic")
                diagn[["call"]] <- match.call()
                attr(res,"diagnostic") <- diagn
                class(attr(res,"diagnostic")) <- "DiagnosticClass"
             }
             return(res)
          })

setMethod("E", signature(object = "AbscontCondDistribution",
                         fun = "EuclRandMatrix", 
                         cond = "numeric"),
          function(object, fun, cond, withCond = FALSE, useApply = TRUE, ..., diagnostic = FALSE){
             dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
             res <- do.call(.locEmatfunCond, c(list(object=object, fun= fun, cond=cond,
                                 withCond = withCond, useApply = useApply,
                                 diagnostic = diagnostic), dotsI))
             if(diagnostic){
                diagn <- attr(res,"diagnostic")
                diagn[["call"]] <- match.call()
                attr(res,"diagnostic") <- diagn
                class(attr(res,"diagnostic")) <- "DiagnosticClass"
             }
             return(res)
          })

setMethod("E", signature(object = "DiscreteCondDistribution",
                         fun = "EuclRandMatrix", 
                         cond = "numeric"),
    function(object, fun, cond, withCond = FALSE, useApply = TRUE, ...){
     dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
     do.call(.locEmatfunCond,c(list(object = object, fun = fun, cond=cond, withCond = withCond,
               useApply = useApply, diagnostic= FALSE), dotsI))
     })

setMethod("E", signature(object = "UnivariateCondDistribution",
                         fun = "EuclRandVarList", 
                         cond = "numeric"),
          function(object, fun, cond, withCond = FALSE, useApply = TRUE, ..., diagnostic = FALSE){
             dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
             res <- do.call(.locElistfunCond, c(list(object=object, fun= fun, cond=cond,
                                 withCond = withCond, useApply = useApply,
                                 diagnostic = diagnostic), dotsI))
             if(diagnostic){
                diagn <- attr(res,"diagnostic")
                diagn[["call"]] <- match.call()
                attr(res,"diagnostic") <- diagn
                class(attr(res,"diagnostic")) <- "DiagnosticClass"
             }
             return(res)
          })

setMethod("E", signature(object = "AbscontCondDistribution", 
                         fun = "EuclRandVarList", 
                         cond = "numeric"),
          function(object, fun, cond, withCond = FALSE, useApply = TRUE, ..., diagnostic = FALSE){
             dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
             res <- do.call(.locElistfunCond, c(list(object=object, fun= fun, cond=cond,
                                 withCond = withCond, useApply = useApply,
                                 diagnostic = diagnostic), dotsI))
             if(diagnostic){
                diagn <- attr(res,"diagnostic")
                diagn[["call"]] <- match.call()
                attr(res,"diagnostic") <- diagn
                class(attr(res,"diagnostic")) <- "DiagnosticClass"
             }
             return(res)
          })

setMethod("E", signature(object = "DiscreteCondDistribution",
                         fun = "EuclRandVarList", 
                         cond = "numeric"),
    function(object, fun, cond, withCond = FALSE, useApply = TRUE, ...){
         dots <- list(...); dotsI <- .filterEargs(dots); dotsI$diagnostic <- NULL
         do.call(.locElistfunCond,c(list(object = object, fun = fun, cond=cond, withCond = withCond,
               useApply = useApply, diagnostic= FALSE), dotsI))
     })

