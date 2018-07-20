#################################
##
## Class: GParetoFamily
##
################################


## methods
setMethod("validParameter",signature(object="ParetoFamily"),
           function(object, param, tol =.Machine$double.eps){
             if (is(param, "ParamFamParameter"))
                 param <- main(param)
             if (!all(is.finite(param)))
                 return(FALSE)
             if (any(param[1] <= tol))
                 return(FALSE)
             return(TRUE)
           })


## generating function 
## Min: known/fixed threshold/location parameter
## shape: shape parameter
## trafo: optional parameter transformation
## start0Est: startEstimator for MLE and MDE

ParetoFamily <- function(Min = 1, shape = 0.5, trafo = NULL, start0Est = NULL,
                    withCentL2 = FALSE){
    theta <- c(Min, shape)


    ##symmetry
    distrSymm <- NoSymmetry()

    ## parameters
    names(theta) <- c("Min", "shape")

    if(missing(trafo)) trafo <- matrix(1, dimnames = list("shape","shape"))
    if(is.matrix(trafo)) if(nrow(trafo) > 1)
           stop("number of rows of 'trafo' > 1")
           # code .define.tau.Dtau is in file GEVFamily.R

    param <- ParamFamParameter(name = "theta", main = theta[2],
                               fixed = theta[1],
                               trafo = trafo)

    ## distribution
    distribution <- Pareto(Min = Min, shape = shape)

    ## starting parameters
    startPar <- function(x,...){
        tr <- theta[1]
        
        ## Pickand estimator
        if(is.null(start0Est)){
           e0 <- log(2)/(log(median(x))-log(tr))
        }else{
           if(is(start0Est,"function")){
              e1 <- start0Est(x, ...)
              e0 <-  if(is(e1,"Estimate")) estimate(e1) else e1
           }else stop("Argument 'start0Est' must be a function or NULL.")
        }
        if(any(x < tr))
               stop("some data smaller than 'Min' ")
        names(e0) <- NULL
        erange <- e0*c(1/10,10)
        return(erange)
    }


    ## what to do in case of leaving the parameter domain
    makeOKPar <- function(theta) {
        theta <- abs(theta)
    }

    modifyPar <- function(theta){
        theta <- abs(theta)
        Pareto(Min = Min, shape = theta)
    }


    ## L2-derivative of the distribution
    L2deriv.fct <- function(param) {
        k <- force(main(param))
        Min0 <- fixed(param)

        Lambda <- function(x) {
            y <- x*0
            ind <- (x > Min0) #
            y[ind] <- 1/k + log(Min0/x[ind])
            return(y)
        }
        ## additional centering of scores to increase numerical precision!
        z <- if(withCentL2) E(Pareto(Min = Min0, shape = k), fun=Lambda) else 0
        fct <- function(x){Lambda(x)-z}
        return(fct)
    }

    ## Fisher Information matrix as a function of parameters
    FisherInfo.fct <- function(param) {
        k <- force(main(param))
        mat <- PosSemDefSymmMatrix(matrix(1/k^2,1,1))
        dimnames(mat) <- list("shape","shape")
        return(mat)
    }

    FisherInfo <- FisherInfo.fct(param)
    name <- "Generalized Pareto Family"

    ## initializing the GPareto family with components of L2-family
    L2Fam <- new("ParetoFamily")
    L2Fam@name <- name
    L2Fam@param <- param
    L2Fam@distribution <- distribution
    L2Fam@L2deriv.fct <- L2deriv.fct
    L2Fam@FisherInfo.fct <- FisherInfo.fct
    L2Fam@FisherInfo <- FisherInfo
    L2Fam@startPar <- startPar
    L2Fam@makeOKPar <- makeOKPar
    L2Fam@modifyParam <- modifyPar
    L2Fam@L2derivSymm <- FunSymmList(NonSymmetric())
    L2Fam@L2derivDistrSymm <- DistrSymmList(NoSymmetry())
    L2Fam@L2derivDistr <- UnivarDistrList(1/shape+log(Min)-log(distribution))

    L2deriv <- EuclRandVarList(RealRandVariable(list(L2deriv.fct(param)),
                               Domain = Reals()))

    L2Fam@fam.call <- substitute(ParetoFamily(Min = Min0, shape = shape0,
                                     trafo = trafo0),
                         list(Min0 = Min, shape0 = shape, trafo0 = trafo))

    L2Fam@L2deriv <- L2deriv
    return(L2Fam)
}

