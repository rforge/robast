#################################
##
## Class: GParetoFamily
##
################################


## methods
setMethod("validParameter",signature(object="GParetoFamily"),
           function(object, param, tol =.Machine$double.eps){
             if (is(param, "ParamFamParameter")) 
                 param <- main(param)
             if (!all(is.finite(param))) 
                 return(FALSE)
             #if (any(param[1] <= tol))
             #    return(FALSE)
             if(object@param@withPosRestr)
                 if (any(param[2] <= tol))
                     return(FALSE)
             return(TRUE)
           })


## generating function 
## loc: known/fixed threshold/location parameter
## scale: scale parameter
## shape: shape parameter
## of.interest: which parameters, transformations are of interest
##              posibilites are: scale, shape, quantile, expected loss, expected shortfall
##              a maximum number of two of these may be selected
## p: probability needed for quantile and expected shortfall
## N: expected frequency for expected loss
## trafo: optional parameter transformation
## start0Est: startEstimator for MLE and MDE --- if NULL HybridEstimator is used;
### now uses exp-Trafo for scale!

GParetoFamily <- function(loc = 0, scale = 1, shape = 0.5, 
                          of.interest = c("scale", "shape"), 
                          p = NULL, N = NULL, trafo = NULL,
                          start0Est = NULL, withPos = TRUE){
    theta <- c(loc, scale, shape)

    of.interest <- .pretreat.of.interest(of.interest,trafo)
             ## code .pretreat.of.interest in GEV.family.R

    ##symmetry
    distrSymm <- NoSymmetry()

    ## parameters
    names(theta) <- c("loc", "scale", "shape")
    scaleshapename <- c("scale", "shape")

    btq <- bDq <- btes <- bDes <- btel <- bDel <- NULL
    if(!is.null(p)){
       btq <- substitute({ q <- loc0 + theta[1]*((1-p0)^(-theta[2])-1)/theta[2]
                           names(q) <- "quantile"
                           }, list(loc0 = loc, p0 = p))

       bDq <- substitute({ scale <- theta[1];  shape <- theta[2]
                        D1 <- ((1-p0)^(-shape)-1)/shape
                        D2 <- -scale/shape*(D1 + log(1-p0)*(1-p0)^(-shape))
                        D <- t(c(D1, D2))
                        rownames(D) <- "quantile"; colnames(D) <- NULL
                        D }, list(p0 = p))
       btes <- substitute({ if(theta[2]>=1L) es <- NA else {
                            q <- loc0 + theta[1]*((1-p0)^(-theta[2])-1)/theta[2]
                            es <- (q + theta[1] - theta[2]*loc0)/(1-theta[2])}
                            names(es) <- "expected shortfall"
                            es }, list(loc0 = loc, p0 = p))
       bDes <- substitute({ if(theta[2]>=1L){ D1 <- D2 <- NA}else{
                            scale <- theta[1]; shape <- theta[2]
                            q <- loc0 + theta[1]*((1-p0)^(-theta[2])-1)/theta[2]
                            dq1 <- ((1-p0)^(-shape)-1)/shape
                            dq2 <- -scale/shape*(dq1 + log(1-p0)*(1-p0)^(-shape))
                            D1 <- (dq1 + 1)/(1-shape)
                            D2 <- (dq2 - loc0)/(1-shape) + (q + scale -
                                       loc0*shape)/(1-shape)^2}
                            D <- t(c(D1, D2))
                            rownames(D) <- "expected shortfall"
                            colnames(D) <- NULL
                            D }, list(loc0 = loc, p0 = p))
    }
    if(!is.null(N)){
       btel <- substitute({ if(theta[2]>=1L) el <- NA else {
                            el <- N0*(loc0 + theta[1]/(1-theta[2]))}
                            names(el) <- "expected loss"
                            el }, list(loc0 = loc,N0 = N))
       bDel <- substitute({ if(theta[2]>=1L){ D1 <- D2 <- NA}else{
                            scale <- theta[1]; shape <- theta[2]
                            D1 <- N0/(1-shape)
                            D2 <- D1*scale/(1-shape)}
                            D <- t(c(D1, D2))
                            rownames(D) <- "expected loss"
                            colnames(D) <- NULL
                            D }, list(loc0 = loc, N0 = N))
    }

    if(is.null(trafo))
       trafo <- .define.tau.Dtau(of.interest, btq, bDq, btes, bDes,
                                 btel, bDel, p, N)
    else if(is.matrix(trafo) & nrow(trafo) > 2)
           stop("number of rows of 'trafo' > 2")
           # code .define.tau.Dtau is in file GEVFamily.R

    param <- ParamFamParameter(name = "theta", main = c(theta[2],theta[3]),
                               fixed = theta[1],
                               trafo = trafo, withPosRestr = withPos,
                               .returnClsName ="ParamWithScaleAndShapeFamParameter")

    ## distribution
    distribution <- GPareto(loc = loc, scale = scale, shape = shape)

    ## starting parameters
    startPar <- function(x,...){
        tr <- theta[1]
        
        ## Pickand estimator
        if(is.null(start0Est)){
           e0 <- estimate(medkMADhybr(x, k=10, ParamFamily=GParetoFamily(loc = theta[1],
                            scale = theta[2], shape = theta[3]),
                            q.lo = 1e-3, q.up = 15))
        }else{
           if(is(start0Est,"function")){
              e1 <- start0Est(x, ...)
              e0 <-  if(is(e1,"Estimate")) estimate(e1) else e1
           }
           if(!is.null(names(e0)))
               e0 <- e0[c("scale", "shape")]
        }

        if(any(x < tr-e0["scale"]/e0["shape"]))
               stop("some data smaller than 'loc-scale/shape' ")

        names(e0) <- NULL
        return(e0)
    }


    ## what to do in case of leaving the parameter domain
    makeOKPar <- function(theta) {
        if(withPos){
           theta <- abs(theta)
        }else{
           if(!is.null(names(theta))){
              theta["scale"] <- abs(theta["scale"])
           }else{
              theta[1] <- abs(theta[1])
           }
        }
        return(theta)
    }

    modifyPar <- function(theta){
        theta <- makeOKPar(theta)

        if(!is.null(names(theta))){
            sc <- theta["scale"]
            sh <- theta["shape"]
        }else{
            theta <- abs(theta)
            sc <- theta[1]
            sh <- theta[2]
        }
        GPareto(loc = loc, scale = sc, shape = sh)
    }


    ## L2-derivative of the distribution
    L2deriv.fct <- function(param) {
        sc <- force(main(param)[1])
        k <- force(main(param)[2])
        tr <- fixed(param)[1] 

        Lambda1 <- function(x) {
            y <- x*0
            ind <- (x > tr) #
            x <- (x[ind]-tr)/sc
            x1 <- 1 + k * x
            y[ind] <- -1/sc + (1+k)/x1*x/sc
            return(y)
        }
        Lambda2 <- function(x) {
            y <- x*0
            ind <- (x > tr) #
            x <- (x[ind]-tr)/sc
            x1 <- 1 + k * x
            y[ind] <- log(x1)/k^2 - (1/k+1)*x/x1
            return(y)
        }
        ## additional centering of scores to increase numerical precision!
        z1 <- E(distribution, fun=Lambda1)
        z2 <- E(distribution, fun=Lambda2)
        return(list(function(x){ Lambda1(x)-z1 },function(x){ Lambda2(x)-z2 }))
    }

    ## Fisher Information matrix as a function of parameters
    FisherInfo.fct <- function(param) {
        sc <- force(main(param)[1])
        k <- force(main(param)[2])
#        tr <- force(fixed(param)[1])
#        fct <- L2deriv.fct(param)
#        P2 <-  GPareto(loc = tr, scale = sc, shape = k)
        E11 <- sc^-2
        E12 <- (sc*(1+k))^-1
        E22 <- 2/(1+k)
        mat <- PosSemDefSymmMatrix(matrix(c(E11,E12,E12,E22)/(1+2*k),2,2))
        dimnames(mat) <- list(scaleshapename,scaleshapename)
        return(mat)
    }

    FisherInfo <- FisherInfo.fct(param)
    name <- "Generalized Pareto Family"

    ## initializing the GPareto family with components of L2-family
    L2Fam <- new("GParetoFamily")
    L2Fam@scaleshapename <- scaleshapename
    L2Fam@name <- name
    L2Fam@param <- param
    L2Fam@distribution <- distribution
    L2Fam@L2deriv.fct <- L2deriv.fct
    L2Fam@FisherInfo.fct <- FisherInfo.fct
    L2Fam@FisherInfo <- FisherInfo
    L2Fam@startPar <- startPar
    L2Fam@makeOKPar <- makeOKPar
    L2Fam@modifyParam <- modifyPar
    L2Fam@L2derivSymm <- FunSymmList(NonSymmetric(), NonSymmetric())
    L2Fam@L2derivDistrSymm <- DistrSymmList(NoSymmetry(), NoSymmetry())

    L2deriv <- EuclRandVarList(RealRandVariable(L2deriv.fct(param),
                               Domain = Reals()))

    L2Fam@fam.call <- substitute(GParetoFamily(loc = loc0, scale = scale0,
                                 shape = shape0, of.interest = of.interest0,
                                 p = p0, N = N0, trafo = trafo0,
                                 withPos = withPos0),
                         list(loc0 = loc, scale0 = scale, shape0 = shape,
                              of.interest0 = of.interest, p0 = p, N0 = N,
                              trafo0 = trafo, withPos0 = withPos))

    L2Fam@LogDeriv <- function(x) (shape+1)/(scale+shape*(x-loc))
    L2Fam@L2deriv <- L2deriv

    suppressWarnings(
      L2Fam@L2derivDistr <- imageDistr(RandVar = L2deriv, distr = distribution)
    )

    return(L2Fam)
}

