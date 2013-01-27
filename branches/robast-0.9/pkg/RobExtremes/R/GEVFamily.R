#################################
##
## Class: GEVFamily for positive shape
##
################################

### some reusable blocks of code (to avoid redundancy):

### pretreatment of of.interest
.pretreat.of.interest <- function(of.interest,trafo){
    if(is.null(trafo)){
        of.interest <- unique(of.interest)
        if(length(of.interest) > 2)
            stop("A maximum number of two parameters resp. parameter transformations may be selected.")
        if(!all(of.interest %in% c("scale", "shape", "quantile", "expected loss", "expected shortfall")))
            stop("Parameters resp. transformations of interest have to be selected from: ",
                "'scale', 'shape', 'quantile', 'expected loss', 'expected shortfall'.")

        ## reordering of of.interest
        if(("scale" %in% of.interest) && ("scale" != of.interest[1])){
            of.interest[2] <- of.interest[1]
            of.interest[1] <- "scale"
        }
        if(!("scale" %in% of.interest) && ("shape" %in% of.interest) && ("shape" != of.interest[1])){
            of.interest[2] <- of.interest[1]
            of.interest[1] <- "shape"
        }
        if(!any(c("scale", "shape") %in% of.interest) && ("quantile" %in% of.interest)
          && ("quantile" != of.interest[1])){
            of.interest[2] <- of.interest[1]
            of.interest[1] <- "quantile"
        }
        if(!any(c("scale", "shape", "quantile") %in% of.interest)
          && ("expected shortfall" %in% of.interest)
          && ("expected shortfall" != of.interest[1])){
            of.interest[2] <- of.interest[1]
            of.interest[1] <- "expected shortfall"
        }
    }
  return(of.interest)
}

.define.tau.Dtau <- function(of.interest, btq, bDq, btes,
                     bDes, btel, bDel, p, N){
        tau <- NULL
        if("scale" %in% of.interest){
            tau <- function(theta){ th <- theta[1]; names(th) <- "scale";  th}
            Dtau <- function(theta){ D <- t(c(1, 0)); rownames(D) <- "scale"; D}
        }
        if("shape" %in% of.interest){
            if(is.null(tau)){
               tau <- function(theta){th <- theta[2]; names(th) <- "shape"; th}
               Dtau <- function(theta){D <- t(c(0,1));rownames(D) <- "shape";D}
            }else{
                tau <- function(theta){th <- theta
                            names(th) <- c("scale", "shape"); th}
                Dtau <- function(theta){ D <- diag(2);
                            rownames(D) <- c("scale", "shape");D}
            }
        }
        if("quantile" %in% of.interest){
            if(is.null(p)) stop("Probability 'p' has to be specified.")
            if(is.null(tau)){
                tau <- function(theta){ }; body(tau) <- btq
                Dtau <- function(theta){ };body(Dtau) <- bDq
            }else{
                tau1 <- tau
                tau <- function(theta){ }
                body(tau) <- substitute({ btq0; c(tau0(theta), q) },
                                        list(btq0=btq, tau0 = tau1))
                Dtau1 <- Dtau
                Dtau <- function(theta){}
                body(Dtau) <- substitute({ bDq0; rbind(Dtau0(theta), D) },
                                         list(Dtau0 = Dtau1, bDq0 = bDq))
            }
        }
        if("expected shortfall" %in% of.interest){
            if(is.null(p)) stop("Probability 'p' has to be specified.")
            if(is.null(tau)){
                tau <- function(theta){ };  body(tau) <- btes
                Dtau <- function(theta){ }; body(Dtau) <- bDes
            }else{
                tau1 <- tau
                tau <- function(theta){ }
                body(tau) <- substitute({ btes0; c(tau0(theta), es) },
                                        list(tau0 = tau1, btes0=btes))
                Dtau1 <- Dtau
                Dtau <- function(theta){}
                body(Dtau) <- substitute({ bDes0; rbind(Dtau0(theta), D) },
                                         list(Dtau0 = Dtau1, bDes0=bDes))
            }
        }
        if("expected loss" %in% of.interest){
            if(is.null(N)) stop("Expected frequency 'N' has to be specified.")
            if(is.null(tau)){
                tau <- function(theta){ }; body(tau) <- btel
                Dtau <- function(theta){ }; body(Dtau) <- bDel
            }else{
                tau1 <- tau
                tau <- function(theta){ }
                body(tau) <- substitute({ btel0; c(tau0(theta), el) },
                                        list(tau0 = tau1, btel0=btel))
                Dtau1 <- Dtau
                Dtau <- function(theta){}
                body(Dtau) <- substitute({ bDel0; rbind(Dtau0(theta), D) },
                                         list(Dtau0 = Dtau1, bDel0=bDel))
            }
        }
        trafo <- function(x){ list(fval = tau(x), mat = Dtau(x)) }
        return(trafo)
}

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

GEVFamily <- function(loc = 0, scale = 1, shape = 0.5,
                          of.interest = c("scale", "shape"),
                          p = NULL, N = NULL, trafo = NULL,
                          start0Est = NULL, withPos = TRUE){
    theta <- c(loc, scale, shape)

    of.interest <- .pretreat.of.interest(of.interest,trafo)

    ##symmetry
    distrSymm <- NoSymmetry()

    ## parameters
    names(theta) <- c("loc", "scale", "shape")
    scaleshapename <- c("scale", "shape")


    btq <- bDq <- btes <- bDes <- btel <- bDel <- NULL
    if(!is.null(p)){
       btq <- substitute({ q <- loc0 + theta[1]*((-log(1-p0))^(-theta[2])-1)/theta[2]
                           names(q) <- "quantile"
                           }, list(loc0 = loc, p0 = p))

       bDq <- substitute({ scale <- theta[1];  shape <- theta[2]
                        D1 <- ((-log(1-p0))^(-shape)-1)/shape
                        D2 <- -scale/shape*(D1 + log(-log(1-p0))*(-log(1-p0))^(-shape))
                        D <- t(c(D1, D2))
                        rownames(D) <- "quantile"; colnames(D) <- NULL
                        D }, list(p0 = p))
       btes <- substitute({ if(theta[2]>=1L) es <- NA else {
                            pg <- pgamma(-log(p0),1-theta[2], lower.tail = TRUE)
                            es <- theta[1] * (gamma(1-theta[2]) * pg/ (1-p0) - 1 )/
                                   theta[2]  + loc0 }
                            names(es) <- "expected shortfall"
                            es }, list(loc0 = loc, p0 = p))
       bDes <- substitute({ if(theta[2]>=1L){ D1 <- D2 <- NA} else {
                            scale <- theta[1]; shape <- theta[2]
                            pg <- pgamma(-log(p0), 1-theta[2], lower.tail = TRUE)
                            dd <- ddigamma(-log(p0),1-theta[2])
                            g0 <- gamma(1-theta[2])
                            D1 <- (g0*pg/(1-p0)-1)/theta[2]
                            D21 <- theta[1]*D1/theta[2]
                            D22 <- theta[1]*dd/(1-p0)/theta[2]
                            D2 <- -D21+D22}
                            D <- t(c(D1, D2))
                            rownames(D) <- "expected shortfall"
                            colnames(D) <- NULL
                            D }, list(loc0 = loc, p0 = p))
    }
    if(!is.null(N)){
       btel <- substitute({ if(theta[2]>=1L) el <- NA else{
                            el <- N0*(loc0+theta[1]*gamma(1-theta[2])/theta[2])}
                            names(el) <- "expected loss"
                            el }, list(loc0 = loc,N0 = N))
       bDel <- substitute({ if(theta[2]>=1L){ D1 <- D2 <- NA}else{
                            scale <- theta[1]; shape <- theta[2]
                            D1 <- N0*gamma(1-shape)/shape
                            D2 <- -N0*theta[1]*digamma(1-theta[2])/theta[2]-
                                   D1*scale/(1-shape)}
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
####
    param <- ParamFamParameter(name = "theta", main = c(theta[2],theta[3]),
                               fixed = theta[1],
                               trafo = trafo, withPosRestr = withPos,
                               .returnClsName ="ParamWithScaleAndShapeFamParameter")

    ## distribution
    distribution <- GEV(loc = loc, scale = scale, shape = shape)

    ## starting parameters
    startPar <- function(x,...){
        mu <- theta[1]
        
        ## Pickand estimator
        if(is.null(start0Est)){
        #source("kMedMad_Qn_Estimators.R")
           e0 <- PickandsEstimator(x,ParamFamily=GParetoFamily(loc = theta[1],
                            scale = theta[2], shape = theta[3]))
        }else{
           if(is(start0Est,"function")){
              e1 <- start0Est(x, ...)
              e0 <-  if(is(e1,"Estimate")) estimate(e1) else e1
           }
           if(!is.null(names(e0)))
               e0 <- e0[c("scale", "shape")]
        }
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

        if(!is.null(names(theta))){
            sc <- theta["scale"]
            sh <- theta["shape"]
        }else{
            theta <- abs(theta)
            sc <- theta[1]
            sh <- theta[2]
        }
        GEV(loc = loc, scale = theta[1], shape = theta[2])
    }


    ## L2-derivative of the distribution
    L2deriv.fct <- function(param) {
        sc <- force(main(param)[1])
        k <- force(main(param)[2])
        tr <- fixed(param)[1]

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
        z1 <- E(distribution, fun=Lambda1)
        z2 <- E(distribution, fun=Lambda2)
        return(list(function(x){ Lambda1(x)-z1 },function(x){ Lambda2(x)-z2 }))
    }

    ## Fisher Information matrix as a function of parameters
    FisherInfo.fct <- function(param) {
        sc <- force(main(param)[1])
        k <- force(main(param)[2])
        G20 <- gamma(2*k)
        G10 <- gamma(k)
        G11 <- digamma(k)*gamma(k)
        G01 <- -0.57721566490153 # digamma(1)
        G02 <- 1.9781119906559 #trigamma(1)+digamma(1)^2
        x0 <- (k+1)^2*2*k
        I11 <- G20*x0-2*G10*k*(k+1)+1
        I11 <- I11/sc^2/k^2
        I12 <- G20*(-x0)+ G10*(k^3+4*k^2+3*k) - k -1
        I12 <- I12 + G11*(k^3+k^2) -G01*k
        I12 <- I12/sc/k^3
        I22 <- G20*x0 +(k+1)^2 -G10*(x0+2*k*(k+1))
        I22 <- I22 - G11*2*k^2*(k+1) + G01*2*k*(1+k)+k^2 *G02
        I22 <- I22 /k^4
        mat <- PosSemDefSymmMatrix(matrix(c(I11,I12,I12,I22),2,2))
        dimnames(mat) <- list(scaleshapename,scaleshapename)
        return(mat)
    }



    FisherInfo <- FisherInfo.fct(param)
    name <- "Generalized Extreme Value Family with positive shape parameter: Frechet Family"

    ## initializing the GPareto family with components of L2-family
    L2Fam <- new("GEVFamily")
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

    L2Fam@fam.call <- substitute(GEVFamily(loc = loc0, scale = scale0,
                                 shape = shape0, of.interest = of.interest0,
                                 p = p0, N = N0, trafo = trafo0,
                                 withPos = withPos0),
                         list(loc0 = loc, scale0 = scale, shape0 = shape,
                              of.interest0 = of.interest, p0 = p, N0 = N,
                              trafo0 = trafo, withPos0 = withPos))

    L2Fam@LogDeriv <- function(x){
                  x0 <- (x-loc)/scale
                  x1 <- 1 + x0 * shape
                  (shape+1)/scale/x1 + x1^(-1-1/shape)/scale
                  }

    L2Fam@L2deriv <- L2deriv

    L2Fam@L2derivDistr <- imageDistr(RandVar = L2deriv, distr = distribution)

    return(L2Fam)
}

#ddigamma(t,s) is d/ds \int_t^\infty exp(-x) x^(-s) dx

ddigamma <- function(t,s){
              int <- function(x) exp(-x)*(-log(x))*x^(-s)
              integrate(int, lower=0, upper=t)$value
              }
              