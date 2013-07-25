#################################
##
## Class: GEVFamily for positive shape and mu unknown
##
################################

## methods
setMethod("validParameter",signature(object="GEVFamilyMuUnknown"),
           function(object, param, tol =.Machine$double.eps){
             if (is(param, "ParamFamParameter")) 
                 param <- main(param)
             if (!all(is.finite(param))) 
                 return(FALSE)
             if (any(param[2] <= tol))
                 return(FALSE)
             if (any(param[3] <= tol))
                 return(FALSE)
             return(TRUE)
           })


## generating function 
## loc: known/fixed threshold/location parameter
## scale: scale parameter
## shape: shape parameter
## trafo: optional parameter transformation
## start0Est: startEstimator for MLE and MDE --- if NULL HybridEstimator is used;

GEVFamilyMuUnknown <- function(loc = 0, scale = 1, shape = 0.5,
                          of.interest = c("scale", "shape"),
                          p = NULL, N = NULL, trafo = NULL,
                          start0Est = NULL, withPos = TRUE,
                          withCentL2 = FALSE,
                          withL2derivDistr  = FALSE,
                          ..ignoreTrafo = FALSE){
    theta <- c(loc, scale, shape)
    .warningGEVShapeLarge(shape)
    
    of.interest <- .pretreat.of.interest(of.interest,trafo)

    ##symmetry
    distrSymm <- NoSymmetry()

    ## parameters
    names(theta) <- c("loc", "scale", "shape")
    scaleshapename <- c("scale"="scale", "shape"="shape")


    btq <- bDq <- btes <- bDes <- btel <- bDel <- NULL
    if(!is.null(p)){
       btq <- substitute({ q <- theta[1] + theta[2]*((-log(p0))^(-theta[3])-1)/theta[3]
                           names(q) <- "quantile"
                           q
                           }, list(p0 = p))

       bDq <- substitute({ loc <- theta[1]; scale <- theta[2];  shape <- theta[3]
                        D1 <- ((-log(p0))^(-shape)-1)/shape
                        D2 <- -scale/shape*(D1 + log(-log(p0))*(-log(p0))^(-shape))
                        D <- t(c(1, D1, D2))
                        rownames(D) <- "quantile"; colnames(D) <- NULL
                        D }, list(p0 = p))
       btes <- substitute({ if(theta[3]>=1L) es <- NA else {
                            pg <- pgamma(-log(p0),1-theta[3], lower.tail = TRUE)
                            es <- theta[2] * (gamma(1-theta[3]) * pg/ (1-p0) - 1 )/
                                   theta[3]  + theta[1] }
                            names(es) <- "expected shortfall"
                            es }, list(p0 = p))
       bDes <- substitute({ if(theta[3]>=1L){ D0 <- NA, D1 <- D2 <- NA} else {
                            loc <- theta[1]; scale <- theta[2]; shape <- theta[3]
                            pg <- pgamma(-log(p0), 1-theta[3], lower.tail = TRUE)
                            dd <- ddigamma(-log(p0),1-theta[3])
                            g0 <- gamma(1-theta[3])
                            D0 <- 1
                            D1 <- (g0*pg/(1-p0)-1)/theta[3]
                            D21 <- D1/theta[2]
                            D22 <- dd/(1-p0)/theta[2]
                            D2 <- -theta[1]*(D21+D22)}
                            D <- t(c(D0,D1, D2))
                            rownames(D) <- "expected shortfall"
                            colnames(D) <- NULL
                            D }, list(p0 = p))
    }
    if(!is.null(N)){
       btel <- substitute({ if(theta[3]>=1L) el <- NA else{
                            el <- N0*(theta[1]+theta[2]*(gamma(1-theta[3])-1)/theta[3])}
                            names(el) <- "expected loss"
                            el }, list(N0 = N))
       bDel <- substitute({ if(theta[3]>=1L){ D0 <- D1 <- D2 <- NA}else{
                            loc <- theta[1]; scale <- theta[2]; shape <- theta[3]
                            ga <- gamma(1-shape)
                            D0 <- 1
                            D1 <- N0*(ga-1)/shape
                            D2 <- -N0*scale*ga*digamma(1-shape)/shape-
                                   D1*scale/shape}
                            D <- t(c(D0, D1, D2))
                            rownames(D) <- "expected loss"
                            colnames(D) <- NULL
                            D }, list(loc0 = loc, N0 = N))
    }

    fromOfInt <- FALSE
    if(is.null(trafo)||..ignoreTrafo){fromOfInt <- TRUE
       trafo <- .define.tau.Dtau(of.interest, btq, bDq, btes, bDes,
                                 btel, bDel, p, N)
    }else if(is.matrix(trafo) & nrow(trafo) > 3)
           stop("number of rows of 'trafo' > 3")
####
    param <- ParamFamParameter(name = "theta", main = c(theta[1],theta[2],theta[3]),
                               fixed = NULL,
                               trafo = trafo, withPosRestr = withPos,
                               .returnClsName ="ParamWithLocAndScaleAndShapeFamParameter")

    ## distribution
    distribution <- GEV(loc = loc, scale = scale, shape = shape)

    ## starting parameters
    startPar <- function(x,...){
        mu <- min(x)
        
        ## Pickand estimator
        if(is.null(start0Est)){
        #source("kMedMad_Qn_Estimators.R")
           PF <- GEVFamily(loc = theta[1], scale = theta[2], shape = theta[3])
           e1 <- PickandsEstimator(x,ParamFamily=PF)
           e0 <- estimate(e1)
        }else{
           if(is(start0Est,"function")){
              e1 <- start0Est(x, ...)
              e0 <-  if(is(e1,"Estimate")) estimate(e1) else e1
           }else stop("Argument 'start0Est' must be a function or NULL.")
           if(!is.null(names(e0)))
               e0 <- e0[c("scale", "shape")]
        }
#        print(e0); print(str(x)); print(head(summary(x))); print(mu)
        if(any(x < mu-e0["scale"]/e0["shape"]))
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
        .warningGEVShapeLarge(theta["shape"])
        if(!is.null(names(theta))){
            loc <- theta["loc"]
            sc <- theta["scale"]
            sh <- theta["shape"]
        }else{
            loc <- theta[1]
            theta[2:3] <- abs(theta[2:3])
            sc <- theta[2]
            sh <- theta[3]
        }
        GEV(loc = theta[1], scale = theta[2], shape = theta[3])
    }


    ## L2-derivative of the distribution
    L2deriv.fct <- function(param) {
        sc <- force(main(param)[2])
        k <- force(main(param)[3])
        tr <- force(main(param)[1])
        .warningGEVShapeLarge(k)

        k1 <- k+1
        Lambda0 <- function(x) {
         y <- x*0
         ind <- (x > tr-sc/k) # = [later] (x1>0)
         x <- (x[ind]-tr)/sc
         x1 <- 1 + k * x
         t1 <- x1^(-1/k)
         y[ind] <- (k1-t1)/x1/sc
#         xi*(-1/xi-1)*(x[ind]-mu)/beta^2/(1+xi*(x[ind]-mu)/beta) - (x[ind]-mu)*(1+xi*(x[ind]-mu)/beta)^(-1/xi-1)/beta^2
         return(y)
        }

        Lambda1 <- function(x) {
         y <- x*0
         ind <- (x > tr-sc/k) # = [later] (x1>0)
         x <- (x[ind]-tr)/sc
         x1 <- 1 + k * x
         y[ind] <- (x*(1-x1^(-1/k))-1)/x1/sc
#         xi*(-1/xi-1)*(x[ind]-mu)/beta^2/(1+xi*(x[ind]-mu)/beta) - (x[ind]-mu)*(1+xi*(x[ind]-mu)/beta)^(-1/xi-1)/beta^2
         return(y)
        }
        Lambda2 <- function(x) {
         y <- x*0
         ind <- (x > tr-sc/k) # = [later] (x1>0)
         x <- (x[ind]-tr)/sc
         x1 <- 1 + k * x
         x2 <- x / x1
         y[ind]<- (1-x1^(-1/k))/k*(log(x1)/k-x2)-x2
#         log(1+xi*(x[ind]-mu)/beta)/xi^2+(-1/xi-1)*(x[ind]-mu)/beta/(1+xi*(x[ind]-mu)/beta) - (1+xi*(x[ind]-mu)/beta)^(-1/xi)*log(1+xi*(x[ind]-mu)/beta)/xi^2 + (1+xi*(x[ind]-mu)/beta)^(-1/xi-1)*(x[ind]-mu)/beta/xi
         return(y)
        }
        ## additional centering of scores to increase numerical precision!
        if(withCentL2){
           dist0 <- GEV(scale = sc, shape = k, loc = tr)
           suppressWarnings({
             z0 <- E(dist0, fun=Lambda0)
             z1 <- E(dist0, fun=Lambda1)
             z2 <- E(dist0, fun=Lambda2)
           })
        }else{z0 <- z1 <- z2 <- 0}
        return(list(function(x){ Lambda0(x)-z0 },
                    function(x){ Lambda1(x)-z1 },function(x){ Lambda2(x)-z2 }))
    }

    ## Fisher Information matrix as a function of parameters
    FisherInfo.fct <- function(param) {
        tr <- force(main(param)[1])
        sc <- force(main(param)[2])
        k <- force(main(param)[3])
        k1 <- k+1
        .warningGEVShapeLarge(k)
        G20 <- gamma(2*k)
        G10 <- gamma(k)
        G11 <- digamma(k)*gamma(k)
        G01 <- -0.57721566490153 # digamma(1)
        G02 <- 1.9781119906559 #trigamma(1)+digamma(1)^2
        x0 <- k1^2*2*k
        I00 <- (2*k)*k1^2*G20/sc^2
        I01 <- (G10-k1*2*G20)*k1/sc^2
        I02 <- [k1*2 * gamma(2k)- k1* gamma(k) -  gamma(k)-k * G11]*k1/k
        I02 <- (2*k1*G20 -(k+2)*G10-k*G11)*k1/k/sc
        I11 <- G20*x0-2*G10*k*(k+1)+1
        I11 <- I11/sc^2/k^2
        I12 <- G20*(-x0)+ G10*(k^3+4*k^2+3*k) - k -1
        I12 <- I12 + G11*(k^3+k^2) -G01*k
        I12 <- I12/sc/k^3
        I22 <- G20*x0 +(k+1)^2 -G10*(x0+2*k*(k+1))
        I22 <- I22 - G11*2*k^2*(k+1) + G01*2*k*(1+k)+k^2 *G02
        I22 <- I22 /k^4
        mat <- PosSemDefSymmMatrix(matrix(c(I00,I01,I02,I01,I11,I12,I02,I12,I22),3,3))
        dimnames(mat) <- list(scaleshapename,scaleshapename)
        return(mat)
    }



    FisherInfo <- FisherInfo.fct(param)
    name <- "GEV Family"

    ## initializing the GPareto family with components of L2-family
    L2Fam <- new("GEVFamilyMuUnknown")
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
       L2Fam@fam.call <- substitute(GEVFamily(loc = loc0, scale = scale0,
                                 shape = shape0, of.interest = of.interest0,
                                 p = p0, N = N0,
                                 withPos = withPos0, withCentL2 = FALSE,
                                 withL2derivDistr  = FALSE, ..ignoreTrafo = TRUE),
                         list(loc0 = loc, scale0 = scale, shape0 = shape,
                              of.interest0 = of.interest, p0 = p, N0 = N,
                              withPos0 = withPos))
    }else{
       L2Fam@fam.call <- substitute(GEVFamily(loc = loc0, scale = scale0,
                                 shape = shape0, of.interest = NULL,
                                 p = p0, N = N0, trafo = trafo0,
                                 withPos = withPos0, withCentL2 = FALSE,
                                 withL2derivDistr  = FALSE),
                         list(loc0 = loc, scale0 = scale, shape0 = shape,
                              p0 = p, N0 = N,
                              withPos0 = withPos, trafo0 = trafo))
    }

    L2Fam@LogDeriv <- function(x){
                  x0 <- (x-loc)/scale
                  x1 <- 1 + x0 * shape
                  (shape+1)/scale/x1 + x1^(-1-1/shape)/scale
                  }

    L2Fam@L2deriv <- L2deriv
    L2Fam@L2derivDistr <- L2derivDistr
    L2Fam@.withMDE <- FALSE
    L2Fam@.withEvalAsVar <- FALSE
    L2Fam@.withEvalL2derivDistr <- FALSE
    return(L2Fam)
}

#ddigamma(t,s) is d/ds \int_0^t exp(-x) x^(s-1) dx

ddigamma <- function(t,s){
              int <- function(x) exp(-x)*(log(x))*x^(s-1)
              integrate(int, lower=0, upper=t)$value
              }
              