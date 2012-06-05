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
                          start0Est = NULL, withPos = TRUE){
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
    theta <- c(scale, shape)

    ##symmetry
    distrSymm <- NoSymmetry()

    ## parameters
    names(theta) <- c("scale", "shape")

    if(is.null(trafo)){
        tau <- NULL
        if("scale" %in% of.interest){
            tau <- function(theta){ 
                th <- theta[1]
                names(th) <- "scale"
                th
            }
            Dtau <- function(theta){ 
                D <- t(c(1, 0))
                rownames(D) <- "scale"
                D 
            }
        }
        if("shape" %in% of.interest){
            if(is.null(tau)){
                tau <- function(theta){ 
                    th <- theta[2] 
                    names(th) <- "shape"
                    th
                }
                Dtau <- function(theta){ 
                    D <- t(c(0, 1))
                    rownames(D) <- "shape"
                    D 
                }
            }else{
                tau <- function(theta){ 
                  th <- theta 
                  names(th) <- c("scale", "shape")
                  th
                }
                Dtau <- function(theta){ 
                    D <- diag(2) 
                    rownames(D) <- c("scale", "shape")
                    D
                }
            }
        }
        if("quantile" %in% of.interest){
            if(is.null(p)) stop("Probability 'p' has to be specified.")
            if(is.null(tau)){
                tau <- function(theta){ }
                body(tau) <- substitute({ q <- exp(shape^(-1)*log(-log(1-p))+log(scale))                                          
                                         names(q) <- "quantile"
                                          q },
                                        list(p0 = p))
                Dtau <- function(theta){ }
                body(Dtau) <- substitute({ scale <- theta[1]
                                           shape <- theta[2]
                                           D1 <- exp(log(-log(1-p))/shape^2+log(scale))/scale 
                                           D2 <- (-log(-log(1-p)))*exp(log(-log(1-p))/shape^2+log(scale))/shape^2
                                           D <- t(c(D1, D2))
                                           rownames(D) <- "quantile"
                                           colnames(D) <- NULL
                                           D },
                                         list(p0 = p))
            }else{
                tau1 <- tau
                tau <- function(theta){ }
                body(tau) <- substitute({ q <- exp(theta[2]^(-1)*log(-log(1-p))+log(theta[1]))    
                                          names(q) <- "quantile"
                                          c(tau0(theta), q) },
                                        list(tau0 = tau1, p0 = p))
                Dtau1 <- Dtau
                Dtau <- function(theta){}
                body(Dtau) <- substitute({ scale <- theta[1]
                                           shape <- theta[2]
                                           D1 <- exp(log(-log(1-p))/shape^2+log(scale))/scale
                                           D2 <- (-log(-log(1-p)))*exp(log(-log(1-p))/shape^2+log(scale))/shape^2
                                           D <- t(c(D1, D2))
                                           rownames(D) <- "quantile"
                                           colnames(D) <- NULL
                                           rbind(Dtau0(theta), D) },
                                         list(Dtau0 = Dtau1, p0 = p))
            }
        }
        if("expected shortfall" %in% of.interest){
            if(is.null(p)) stop("Probability 'p' has to be specified.")
            if(is.null(tau)){
                tau <- function(theta){ }
                body(tau) <- substitute({ q <- exp(theta[2]^(-1)*log(-log(1-p))+log(theta[1]))   
                                          es <- E(q,upp=Inf,low=p0)/(1-p0) 
                                          names(es) <- "expected shortfall"
                                          es }, 
                                        list(p0 = p))
                Dtau <- function(theta){ }
                body(Dtau) <- substitute({ scale <- theta[1]
                                           shape <- theta[2]
                                           q <- exp(shape^(-1)*log(-log(1-p))+log(scale))    
                                           dq1 <- exp(log(-log(1-p))/shape^2+log(scale))/scale
                                           dq2 <- (-log(-log(1-p)))*exp(log(-log(1-p))/shape^2+log(scale))/shape^2
                                           D1 <- E(dq1,upp=Inf,low=p0)/(1-p0)
                                           D2 <- E(dq2,upp=Inf,low=p0)/(1-p0)
                                           D <- t(c(D1, D2))
                                           rownames(D) <- "expected shortfall"
                                           colnames(D) <- NULL
                                           D },
                                         list(p0 = p))
            }else{
                tau1 <- tau
                tau <- function(theta){ }
                body(tau) <- substitute({ q <- exp(theta[2]^(-1)*log(-log(1-p))+log(theta[1]))   
                                          es <- E(q,upp=Inf,low=p0)/(1-p0) 
                                          names(es) <- "expected shortfall"
                                          c(tau0(theta), es) }, 
                                        list(tau0 = tau1, p0 = p))
                Dtau1 <- Dtau
                Dtau <- function(theta){}
                body(Dtau) <- substitute({ scale <- theta[1]
                                           shape <- theta[2]
                                           q <- exp(shape^(-1)*log(-log(1-p))+log(scale))    
                                           dq1 <- exp(log(-log(1-p))/shape^2+log(scale))/scale
                                           dq2 <- (-log(-log(1-p)))*exp(log(-log(1-p))/shape^2+log(scale))/shape^2
                                           D1 <- E(dq1,upp=Inf,low=p0)/(1-p0)
                                           D2 <- E(dq2,upp=Inf,low=p0)/(1-p0)
                                           D <- t(c(D1, D2))
                                           rownames(D) <- "expected shortfall"
                                           colnames(D) <- NULL
                                           rbind(Dtau0(theta), D) },
                                         list(Dtau0 = Dtau1, p0 = p))
            }
        }
        if("expected loss" %in% of.interest){
            if(is.null(N)) stop("Expected frequency 'N' has to be specified.")
            if(is.null(tau)){
                tau <- function(theta){ }
                body(tau) <- substitute({ el <- N0*gamma(1+1/shape)*scale
                                          names(el) <- "expected loss"
                                          el },
                                        list(N0 = N))
                Dtau <- function(theta){ }
                body(Dtau) <- substitute({ scale <- theta[1]
                                           shape <- theta[2]
                                           D1 <- N0*gamma(1/xi+1)
                                           D2 <- -N0*digamma(1/xi+1)*beta/xi^2
                                           D <- t(c(D1, D2))
                                           rownames(D) <- "expected loss"
                                           colnames(D) <- NULL
                                           D },
                                         list(N0 = N))
            }else{
                tau1 <- tau
                tau <- function(theta){ }
                body(tau) <- substitute({ el <- N0*gamma(1+1/shape)*scale
                                          names(el) <- "expected loss"
                                          c(tau0(theta), el) },
                                        list(tau0 = tau1, N0 = N))
                Dtau1 <- Dtau
                Dtau <- function(theta){}
                body(Dtau) <- substitute({ scale <- theta[1]
                                           shape <- theta[2]
                                           D1 <- N0*gamma(1/xi+1)
                                           D2 <- -N0*digamma(1/xi+1)*beta/xi^2
                                           D <- t(c(D1, D2))
                                           rownames(D) <- "expected loss"
                                           colnames(D) <- NULL
                                           rbind(Dtau0(theta), D) },
                                         list(Dtau0 = Dtau1, N0 = N))
            }
        }
        trafo <- function(x){ list(fval = tau(x), mat = Dtau(x)) }
    }else{
        if(is.matrix(trafo) & nrow(trafo) > 2) stop("number of rows of 'trafo' > 2")
    }
    param <- ParamFamParameter(name = "theta", main = c(theta[1],theta[2]),
                               trafo = trafo, withPosRestr = withPos,
                               .returnClsName ="ParamWithScaleAndShapeFamParameter")

    ## distribution
    distribution <- Weibull(scale = scale, shape = shape)

    ## starting parameters
    startPar <- function(x,...){
       

        ## Pickand estimator
        if(is.null(start0Est)){
           e0 <- estimate(medkMADhybr(x, k=10, ParamFamily=WeibullFamily(scale = theta[1], shape = theta[2]),
                            q.lo = 1e-3, q.up = 15))
        }else{
           if(is(start0Est,"function")){
              e1 <- start0Est(x, ...)
              e0 <-  if(is(e1,"Estimate")) estimate(e1) else e1
           }
           if(!is.null(names(e0)))
               e0 <- e0[c("scale", "shape")]
        }
        names(e0) <- NULL
        return(e0)
    }


    ## what to do in case of leaving the parameter domain
    makeOKPar <- function(theta) {
        if(withPos){
           if(!is.null(names(theta)))
                 theta["shape"] <- abs(theta["shape"])
           else  theta[2] <- abs(theta[2])
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
            y[x>0] <- -1/sc-(sh-1)/sc+sh*x1/sc^2###
            return(y)
        }
        Lambda2 <- function(x) {
            y <- x*0
            x1 <- x[x>0]
            y[x>0] <- 1/sh+log(x1/sc)-x1/sc###
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
        sh <- force(main(param)[2])
        Lambda <- L2deriv.fct(param)
        E11 <- E(distribution,fun=function(x)Lambda[[1]](x)^2)
        E12 <- E(distribution,fun=function(x)Lambda[[1]](x)*Lambda[[2]](x))
        E22 <- E(distribution,fun=function(x)Lambda[[2]](x)^2)
        return(PosSemDefSymmMatrix(matrix(c(E11,E12,E12,E22),2,2)))
    }

    FisherInfo <- FisherInfo.fct(param)
    name <- "Generalized Pareto Family"

    ## initializing the Weibull family with components of L2-family
    L2Fam <- new("WeibullFamily")
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

    L2Fam@fam.call <- substitute(WeibullFamily(scale = scale0,
                                 shape = shape0, of.interest = of.interest0,
                                 p = p0, N = N0, trafo = trafo0,
                                 withPos = withPos0),
                         list(scale0 = scale, shape0 = shape,
                              of.interest0 = of.interest, p0 = p, N0 = N,
                              trafo0 = trafo, withPos0 = withPos))

    L2Fam@LogDeriv <- function(x) log(shape)-log(scale)+(shape-1)*log(x/scale)-shape*x/scale###
    L2Fam@L2deriv <- L2deriv

    L2Fam@L2derivDistr <- imageDistr(RandVar = L2deriv, distr = distribution)


    return(L2Fam)
}

