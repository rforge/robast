#################################
##
## Class: WeibullFamily
##
################################


## methods
setMethod("validParameter",signature(object="WeibullFamily"),
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

WeibullFamily <- function(scale = 1, shape = 0.5, 
                          of.interest = c("scale", "shape"), 
                          p = NULL, N = NULL, trafo = NULL,
                          start0Est = NULL, withPos = TRUE,
                          withCentL2 = FALSE,
                          withL2derivDistr  = FALSE,
                          ..ignoreTrafo = FALSE){
    theta <- c(scale, shape)

    of.interest <- .pretreat.of.interest(of.interest,trafo)
             ## code .pretreat.of.interest in GEV.family.R

    ##symmetry
    distrSymm <- NoSymmetry()

    ## parameters
    names(theta) <- c("scale", "shape")
    scaleshapename <- c("scale"="scale", "shape"="shape")

    btq <- bDq <- btes <- bDes <- btel <- bDel <- NULL
    if(!is.null(p)){
       btq <- substitute({ q <- theta[1]*(-log(1-p0))^(1/theta[2])
                           names(q) <- "quantile"
                           q
                           }, list(p0 = p))

       bDq <- substitute({ scale <- theta[1];  shape <- theta[2]
                        lp <- -log(1-p0)
                        D1 <- lp^(1/shape)
                        D2 <- -scale/shape^2 *lp^(1/shape)*log(lp)
                        D <- t(c(D1, D2))
                        rownames(D) <- "quantile"; colnames(D) <- NULL
                        D }, list(p0 = p))
       btes <- substitute({ if(theta[2]<= (-1L)) es <- NA else {
                            s1 <- 1+1/theta[2]
                            pg <- pgamma(-log(1-p0),s1, lower.tail = FALSE)
                            g0 <- gamma(s1)
                            es <- theta[1] * g0 * pg /(1-p0)}
                            names(es) <- "expected shortfall"
                            es }, list(p0 = p))
       bDes <- substitute({ if(theta[2]<= (-1L)){ D1 <- D2 <- NA} else {
                            s1 <- 1+1/theta[2]
                            pg <- pgamma(-log(1-p0), s1, lower.tail = FALSE)
                            g0 <- gamma(s1)
                            ## dd <- ddigamma(Inf,s1)-ddigamma(-log(1-p0),s1)
                            dd <- digamma(s1)*g0 - ddigamma(-log(1-p0),s1)
                            D1 <-  g0 * pg
                            D2 <- - theta[1] * dd /theta[2]^2}
                            D <- t(c(D1, D2))/(1-p0)
                            rownames(D) <- "expected shortfall"
                            colnames(D) <- NULL
                            D }, list(p0 = p))
    }
    if(!is.null(N)){
       btel <- substitute({ el <- N0*(theta[1]*gamma(1+1/theta[2]))
                            names(el) <- "expected loss"
                            el }, list(N0 = N))
       bDel <- substitute({ scale <- theta[1]; shape <- theta[2]
                            s1 <- 1+1/shape
                            D1 <- N0*gamma(s1)
                            D2 <- -N0*theta[1]*digamma(s1)*gamma(s1)/shape^2
                            D <- t(c(D1, D2))
                            rownames(D) <- "expected loss"
                            colnames(D) <- NULL
                            D }, list(N0 = N))
    }

    fromOfInt <- FALSE
    if(is.null(trafo)||..ignoreTrafo){fromOfInt <- TRUE
       trafo <- .define.tau.Dtau(of.interest, btq, bDq, btes, bDes,
                                 btel, bDel, p, N)
    }else if(is.matrix(trafo) & nrow(trafo) > 2)
           stop("number of rows of 'trafo' > 2")
           # code .define.tau.Dtau is in file GEVFamily.R


    param <- ParamFamParameter(name = "theta", main = c(theta[1],theta[2]),
                               trafo = trafo, withPosRestr = withPos,
                               .returnClsName ="ParamWithScaleAndShapeFamParameter")

    ## distribution
    distribution <- Weibull(scale = scale, shape = shape)

    ## starting parameters
    startPar <- function(x,...){
       

        ## Pickand estimator
        if(is.null(start0Est)){
           e1 <- QuantileBCCEstimator(x)
           e0 <- estimate(e1)
        }else{
           if(is(start0Est,"function")){
              e1 <- start0Est(x, ...)
              e0 <-  if(is(e1,"Estimate")) estimate(e1) else e1
           }else stop("Argument 'start0Est' must be a function or NULL.")
           if(!is.null(names(e0)))
               e0 <- e0[c("scale", "shape")]
        }
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
        Weibull(scale = sc, shape = sh)
    }


    ## L2-derivative of the distribution
    L2deriv.fct <- function(param) {
        sc <- force(main(param)[1])
        sh <- force(main(param)[2])

        Lambda1 <- function(x) {
            y <- x*0
            x1 <- x[x>0]
            z <- x1/sc
            y[x>0] <- sh/sc*(z^sh-1)###
            return(y)
        }
        Lambda2 <- function(x) {
            y <- x*0
            x1 <- x[x>0]
            z <- x1/sc
            y[x>0] <- 1/sh-log(z)*(z^sh-1)###
            return(y)
        }
        ## additional centering of scores to increase numerical precision!
        if(withCentL2){
           dist0 <- Weibull(scale = sc, shape = sh)
           z1 <- E(dist0, fun=Lambda1)
           z2 <- E(dist0, fun=Lambda2)
        }else{z1 <- z2 <- 0}
        return(list(function(x){ Lambda1(x)-z1 },function(x){ Lambda2(x)-z2 }))
    }

    ## Fisher Information matrix as a function of parameters
    ## Fisher Information matrix as a function of parameters
    FisherInfo.fct <- function(param) {
        sc <- force(main(param)[1])
        k <- force(main(param)[2])
        g1 <- 0.42278433509847 # 1+digamma(1)
        g2 <-   1.8236806608529 # 1+trigamma(1)+digamma(1)^2+2*digamma(1)
        I11 <- k^2/sc^2
        I12 <- -g1/sc
        I22 <- g2/k^2
        mat <- PosSemDefSymmMatrix(matrix(c(I11,I12,I12,I22),2,2))
        dimnames(mat) <- list(scaleshapename,scaleshapename)
        return(mat)
    }



    FisherInfo <- FisherInfo.fct(param)
    name <- "Weibull Family"

    ## initializing the Weibull family with components of L2-family
    L2Fam <- new("WeibullFamily")
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
    L2derivDistr <- NULL
    if(withL2derivDistr){
       suppressWarnings(L2derivDistr <-
          imageDistr(RandVar = L2deriv, distr = distribution))
    }

    if(fromOfInt){
       L2Fam@fam.call <- substitute(WeibullFamily(scale = scale0,
                                 shape = shape0, of.interest = of.interest0,
                                 p = p0, N = N0,
                                 withPos = withPos0, withCentL2 = FALSE,
                                 withL2derivDistr  = FALSE, ..ignoreTrafo = TRUE),
                         list(scale0 = scale, shape0 = shape,
                              of.interest0 = of.interest, p0 = p, N0 = N,
                              withPos0 = withPos))
    }else{
       L2Fam@fam.call <- substitute(WeibullFamily(scale = scale0,
                                 shape = shape0, of.interest = NULL,
                                 p = p0, N = N0, trafo = trafo0,
                                 withPos = withPos0, withCentL2 = FALSE,
                                 withL2derivDistr  = FALSE),
                         list(scale0 = scale, shape0 = shape,
                              p0 = p, N0 = N,
                              withPos0 = withPos, trafo0 = trafo))
    }

    L2Fam@LogDeriv <- function(x){ z <- x/scale
            log(shape)-log(scale)+(shape-1)*log(z)-shape*z^(shape-1)
            }###
    L2Fam@L2deriv <- L2deriv
    L2Fam@L2derivDistr <- L2derivDistr

    L2Fam@.withMDE <- FALSE
    L2Fam@.withEvalAsVar <- FALSE
    L2Fam@.withEvalL2derivDistr <- FALSE

    return(L2Fam)
}

