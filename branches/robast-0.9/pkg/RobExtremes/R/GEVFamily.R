#################################
##
## Class: GEVFamily for positive shape
##
################################

## class
setClass("GEVFamily", contains="L2ParamFamily")

## methods
setMethod("validParameter",signature(object="GEVFamily"),
           function(object, param, tol =.Machine$double.eps){
             if (is(param, "ParamFamParameter")) 
                 param <- main(param)
             if (!all(is.finite(param))) 
                 return(FALSE)
             if (any(param[1] <= tol)) 
                 return(FALSE)
             if (any(param[2] <= tol))
                 return(FALSE)
             return(TRUE)
           })


## generating function 
## loc: known/fixed threshold/location parameter
## scale: scale parameter
## shape: shape parameter
## trafo: optional parameter transformation
## start0Est: startEstimator for MLE and MDE --- if NULL HybridEstimator is used;

GEVFamily <- function(loc = 0, scale = 1, shape = 0.5,trafo = NULL,start0Est = NULL){
    if(is.null(trafo)) trafo = diag(2)
    
    theta <- c(loc, scale, shape)

    ##symmetry
    distrSymm <- NoSymmetry()

    ## parameters
    names(theta) <- c("loc", "scale", "shape")

    param <- ParamFamParameter(name = "theta", main = theta[2:3], 
                               fixed = theta[1],
                               trafo = trafo)

    ## distribution
    distribution <- GEV(loc = loc, scale = scale, shape = shape)

    ## starting parameters
    startPar <- function(x,...){
        mu <- theta[1]
        
        ## Pickand estimator
        if(is.null(start0Est)){
        source("kMedMad_Qn_Estimators.R")
           e0 <- EVTEst(x,est="kMAD")
        }else{
           if(is(start0Est,"function")){
              e1 <- start0Est(x, ...)
              e0 <-  if(is(e1,"Estimate")) estimate(e1) else e1
           }
           if(!is.null(names(e0)))
               e0 <- e0[c("scale", "shape")]
        }
        if(any(x < mu-e0["scale"]/e0["shape"])) stop("some data smaller than 'loc-scale/shape' ")

        names(e0) <- NULL
        return(e0)
    }

    modifyPar <- function(theta){
        theta <- abs(theta)
        GEV(loc = loc, scale = theta[1], shape = theta[2])
    }

    ## what to do in case of leaving the parameter domain
    makeOKPar <- function(theta) {
        theta <- abs(theta)
        theta[2] <- pmin(theta[2],10)
        return(theta)
    }

    ## L2-derivative of the distribution
    L2deriv.fct <- function(param) {
        beta <- force(main(param)[1])
        xi <- force(main(param)[2])
        mu <- fixed(param)[1] 

        Lambda1 <- function(x) {
         y <- x*0
         ind <- x>(mu-beta/xi)
         y[ind] <- -1/beta - xi*(-1/xi-1)*(x[ind]-mu)/beta^2/(1+xi*(x[ind]-mu)/beta) - (x[ind]-mu)*(1+xi*(x[ind]-mu)/beta)^(-1/xi-1)/beta^2
         return(y)
        }
        Lambda2 <- function(x) {
         y <- x*0
         ind <- x>(mu-beta/xi)
         y[ind]<- log(1+xi*(x[ind]-mu)/beta)/xi^2+(-1/xi-1)*(x[ind]-mu)/beta/(1+xi*(x[ind]-mu)/beta) - (1+xi*(x[ind]-mu)/beta)^(-1/xi)*log(1+xi*(x[ind]-mu)/beta)/xi^2 + (1+xi*(x[ind]-mu)/beta)^(-1/xi-1)*(x[ind]-mu)/beta/xi
         return(y)
        }
        ## additional centering of scores to increase numerical precision!
        z1 <- E(distribution, fun=Lambda1)
        z2 <- E(distribution, fun=Lambda2)
        return(list(function(x){ Lambda1(x)-z1 },function(x){ Lambda2(x)-z2 }))
    }

    ## Fisher Information matrix as a function of parameters
    FisherInfo.fct <- function(param) {
        beta <- force(main(param)[1])
        xi <- force(main(param)[2])
        mu <- force(fixed(param)[1])
        fct <- L2deriv.fct(param)
        P <-  GEV(loc = mu, scale = beta, shape = xi)
        E11 <- E(P,function(x)fct[[1]](x)*fct[[1]](x))
        E12 <- E(P,function(x)fct[[1]](x)*fct[[2]](x))
        E22 <- E(P,function(x)fct[[2]](x)*fct[[2]](x))
        return(PosSemDefSymmMatrix(matrix(c(E11,E12,E12,E22),2,2)))
    }

    FisherInfo <- FisherInfo.fct(param)
    name <- "Generalized Extreme Value Family with positive shape parameter: Frechet Family"

    ## initializing the GPareto family with components of L2-family
    res <- L2ParamFamily(name = name, param = param, 
                         distribution = distribution, 
                         L2deriv.fct = L2deriv.fct, 
                         FisherInfo.fct = FisherInfo.fct,
                         FisherInfo = FisherInfo,
                         startPar = startPar,
                         makeOKPar = makeOKPar,
                         modifyParam = modifyPar,
                         .returnClsName = "GEVFamily")
    f.call <- substitute(GEVFamily(loc = loc0, scale = scale0, shape = shape0,trafo = trafo0), 
                         list(loc0 = loc, scale0 = scale, shape0 = shape,trafo0 = trafo))
    res@fam.call <- f.call
    return(res)
}

