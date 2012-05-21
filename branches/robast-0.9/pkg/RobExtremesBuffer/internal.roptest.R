.fix.in.defaults <- function(call.list, fun){
 formals.fun <- formals(fun)
 k <- length(call.list)
 L <- length(formals.fun)
 if("..." %in% names(formals.fun)) L <- L-1
 for(i in 1:L){
     if(!is(formals.fun[[i]],"name")){
        if(!names(formals.fun)[i] %in% names(call.list)&&!is.null(formals.fun[[i]])){
           k <- k + 1
           call.list[[k]] <- formals.fun[[i]]
           names(call.list)[k] <- names(formals.fun)[i]
        }
     }
  }
 return(call.list)

}

.pretreat <- function(x, na.rm = TRUE){
    if(missing(x))
        stop("'x' is missing with no default")
    if(!is.numeric(x)){
        if(is.data.frame(x))
            x <- data.matrix(x)
        else
            x <- as.matrix(x)
        if(!is.matrix(x))
            stop("'x' has to be a numeric vector resp. a matrix or data.frame")
    }
    completecases <- complete.cases(x)
    if(na.rm) x <- na.omit(x)
}
.check.eps <- function(...){
   mc <- match.call(expand=TRUE)
   
   eps <- eps.lower <- eps.upper <- NULL
   if(is.null(mc$eps) && is.null(mc$eps.lower) && is.null(mc$eps.upper)){
        eps.lower <- 0
        eps.upper <- 0.5
    }
    if(is.null(mc$eps)){
        if(!is.null(mc$eps.lower) && is.null(mc$eps.upper))
            eps.upper <- 0.5
        if(is.null(mc$eps.lower) && !is.null(mc$eps.upper))
            eps.lower <- 0
        if(length(eps.lower) != 1 || length(eps.upper) != 1)
            stop("'eps.lower' and 'eps.upper' have to be of length 1")
        if(!is.numeric(eps.lower) || !is.numeric(eps.upper) || eps.lower >= eps.upper)
            stop("'eps.lower' < 'eps.upper' is not fulfilled")
        if((eps.lower < 0) || (eps.upper > 0.5))
            stop("'eps.lower' and 'eps.upper' have to be in [0, 0.5]")
    }else{
        eps <- mc$eps
        if(length(eps) != 1)
            stop("'eps' has to be of length 1")
        if(eps == 0)
            stop("'eps = 0'! => use functions 'mean' and 'sd' for estimation")
        if((eps < 0) || (eps > 0.5))
            stop("'eps' has to be in (0, 0.5]")
    }
    x <- mc$x
    if(is.matrix(x))
        sqrtn <- sqrt(ncol(x))
    else
        sqrtn <- sqrt(length(x))

    return(list(e=eps,lower=eps.lower, upper=eps.upper, sqn = sqrtn))
}

.isOKsteps <- function(steps){
    if(!is.integer(steps))
        steps <- as.integer(steps)
    if(steps < 1){
        stop("'steps' has to be some positive integer value")
    }
    if(length(steps) != 1){
        stop("'steps' has to be of length 1")
    }
   return(invisible(NULL))
}
.isOKfsCor <- function(fsCor){}
    if(fsCor <= 0)
        stop("'fsCor' has to be positive")
    if(length(fsCor) != 1){
        stop("'fsCor' has to be of length 1")
   return(invisible(NULL))
}


.getROptICstart <- function(...){
    mc <- match.call(expand=TRUE)
    eps <- mc$eps
    dots <- mc$dots
    
    if(is.null(eps$e))){
        r.lower <- eps$sqn * eps$lower
        r.upper <- eps$sqn * eps$upper
        ICstart <- do.call(radiusMinimaxIC,
                    c(list(L2Fam = mc$L2FamStart, neighbor = mc$neighbor,
                                   risk = mc$risk,
                                   loRad = r.lower, upRad = r.upper,
                                   verbose = mc$verbose,
                                   OptOrIter = mc$OptOrIter),dots))
        if(!isTRUE(all.equal(mc$fsCor, 1, tol = 1e-3))){
            neighbor@radius <- neighborRadius(ICstart)*mc$fsCor
            infMod <- InfRobModel(center = mc$L2FamStart, neighbor = mc$neighbor)
            ICstart <- do.call(optIC, c(list( model = mc$infMod, risk = mc$risk,
                               verbose = mc$verbose, OptOrIter = mc$OptOrIter),
                               dots))
        }
    }else{
        neighbor@radius <- eps$sqn*eps$e*mc$fsCor
        infMod <- InfRobModel(center = mc$L2FamStart, neighbor = mc$neighbor)
        ICstart <- do.call(optIC, c(list(model = mc$infMod, risk = mc$risk,
                           verbose = mc$verbose, OptOrIter = mc$OptOrIter),
                           dots))
    }
  return(ICstart)
}

genkStepCtrl <- function(useLast = getRobAStBaseOption("kStepUseLast"),
                    withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
                    IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
                    withICList = getRobAStBaseOption("withICList"),
                    withPICList = getRobAStBaseOption("withPICList"),
                    scalename = "scale", withLogScale = TRUE){
  es.call <- match.call()
  es.list <- as.list(es.call[-1])
  es.list <- .fix.in.defaults(es.list,genkStepCtrl)
 return(es.list)
}
genstartCtrl<- function(initial.est = NULL, initial.est.ArgList = NULL,
                        startPar = NULL, distance = CvMDist){
  es.call <- match.call()
  es.list <- as.list(es.call[-1])
  es.list <- .fix.in.defaults(es.list,genstartCtrl)
 return(es.list)
}
gennbCtrl <- function(neighbor = ContNeighborhood(),
                      eps, eps.lower, eps.upper){
  es.call <- match.call()
  es.list <- as.list(es.call[-1])
  es.list <- .fix.in.defaults(es.list,genstartCtrl)
 return(es.list)
}