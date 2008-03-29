###############################################################################
## get minimum bias solutions
###############################################################################
setMethod("getInfRobIC", signature(L2deriv = "UnivariateDistribution", 
                                   risk = "asBias", 
                                   neighbor = "UncondNeighborhood"),
    function(L2deriv, risk, neighbor, symm, trafo, maxiter, 
             tol, ...){
        minmaxBias(L2deriv, neighbor, biastype(risk), symm, 
                   trafo, maxiter, tol)
    })
setMethod("getInfRobIC", signature(L2deriv = "RealRandVariable", 
                                   risk = "asBias", 
                                   neighbor = "ContNeighborhood"),
    function(L2deriv, risk, neighbor, Distr, DistrSymm, L2derivSymm, 
             L2derivDistrSymm, z.start, 
             A.start, Finfo, trafo, maxiter, tol, ...){                

        normtype <- normtype(risk)
        
        FI <- solve(trafo%*%solve(Finfo)%*%t(trafo))
        if(is(normtype,"InfoNorm") || is(normtype,"SelfNorm") ) 
           {QuadForm(normtype) <- PosSemDefSymmMatrix(FI); 
            normtype(risk) <- normtype}

        comp <- .getComp(L2deriv, DistrSymm, L2derivSymm,
             L2derivDistrSymm)
             
        z.comp <- comp$"z.comp"
        A.comp <- comp$"A.comp"

        minmaxBias(L2deriv = L2deriv, neighbor = neighbor, 
                   biastype = biastype(risk), normtype = normtype(risk),
             Distr = Distr, z.start = z.start, A.start = A.start, 
             z.comp = z.comp, A.comp = A.comp, trafo = trafo,
             maxiter = maxiter, tol = tol)
    })


setMethod("minmaxBias", signature(L2deriv = "UnivariateDistribution", 
                                   neighbor = "ContNeighborhood",
                                   biastype = "BiasType"),
    function(L2deriv, neighbor, biastype, symm, 
             trafo, maxiter, tol){
        zi <- sign(as.vector(trafo))
        A <- as.matrix(zi)
        z <- q(L2deriv)(0.5)
        b <- zi*as.vector(trafo)/E(L2deriv, function(x, z){abs(x - z)}, z = z)

        if(is(L2deriv, "AbscontDistribution"))
            ws0 <- 0
        else
            ws0 <- d(L2deriv)(z)
        if(ws0 > 0)
            d <- (2*p(L2deriv)(z) - ws0 - 1)/ws0
        else 
            d <- 0
        
        info <- c("minimum asymptotic bias (lower case) solution")
        asCov <- b^2*(1-ws0) + b^2*d^2*ws0
        Risk <- list(asBias = b, asCov = asCov)

        w <- new("HampelWeight")
        cent(w) <- z
        stand(w) <- A
        clip(w) <- b
        weight(w) <- minbiasweight(w, neighbor = neighbor, biastype = biastype,
                       normW= NormType())

        return(list(A = A, a = zi*z, b = b, d = d, risk = Risk, info = info, 
                    w = w, biastype = biastype, normtype = NormType()))    
    })

setMethod("minmaxBias", signature(L2deriv = "UnivariateDistribution", 
                                   neighbor = "TotalVarNeighborhood",
                                   biastype = "BiasType"),
    function(L2deriv, neighbor, biastype,
             symm, trafo, 
             maxiter, tol){
        zi <- sign(as.vector(trafo))
        A <- as.matrix(zi)
        b <- zi*as.vector(trafo)/(-m1df(L2deriv, 0))
        p0 <- p(L2deriv)(0)
        if(is(L2deriv, "AbscontDistribution"))
            ws0 <- 0
        else
            ws0 <- d(L2deriv)(0)

        if(zi == 1)
            a <- -b*(1-p0)/(1-ws0)
        else
            a <- -b*(p0-ws0)/(1-ws0)
            
        info <- c("minimum asymptotic bias (lower case) solution")
        Risk <- list(asCov = a^2*(p0-ws0) + (zi*a+b)^2*(1-p0), asBias = b)

        w <- new("BdStWeight")
        stand(w) <- A
        clip(w) <- c(a, a+b)
        weight(w) <- minbiasweight(w, neighbor = neighbor, biastype = biastype)

        return(list(A = A, a = a, b = b, d = 0, risk = Risk, info = info, 
                    w = w, biastype = biastype, normtype = NormType()))
    })

setMethod("minmaxBias", signature(L2deriv = "RealRandVariable", 
                                   neighbor = "ContNeighborhood", 
                                   biastype = "BiasType"),
    function(L2deriv, neighbor, biastype, normtype, Distr, 
             z.start, A.start,  z.comp, A.comp, trafo, maxiter,  tol){                

        eerg <- .LowerCaseMultivariate(L2deriv, neighbor, biastype,
             normtype, Distr, trafo, z.start,
             A.start, z.comp = z.comp, A.comp = A.comp,  maxiter, tol)
        erg <- eerg$erg

        b <- 1/erg$value
        param <- erg$par
        lA.comp <- sum(A.comp)
        
        p <- nrow(trafo)
        k <- ncol(trafo)
        A <- matrix(0, ncol=k, nrow=p)
        A[A.comp] <- matrix(param[1:lA.comp], ncol=k, nrow=p)
        z <- numeric(k)
        z[z.comp] <- param[(lA.comp+1):length(param)]
        a <- as.vector(A %*% z)
        d <- numeric(p)
        # computation of 'd', in case 'L2derivDistr' not abs. cont.
        
        info <- c("minimum asymptotic bias (lower case) solution")
        Risk <- list(asBias = b)

        w <- eerg$w
        normtype <- eerg$normtype
        
        return(list(A = A, a = a, b = b, d = d, risk = Risk, info = info, 
                    w = w, biastype = biastype, normtype = normtype))
    })

setMethod("minmaxBias", signature(L2deriv = "UnivariateDistribution", 
                                   neighbor = "ContNeighborhood", 
                                   biastype = "asymmetricBias"),
    function(L2deriv, neighbor, biastype, symm, 
             trafo, maxiter, tol){                
        nu1 <- nu(biastype)[1]
        nu2 <- nu(biastype)[2]
        zi <- sign(as.vector(trafo))
        A <- as.matrix(zi)
        z <- q(L2deriv)(nu1/(nu1+nu2))
        b <- zi*as.vector(trafo)/E(L2deriv, function(x, z){(x - z)*(x>z)/nu2 +
                 (z-x)*(z>x)/nu1}, z = z)

        b1 <- b / nu1
        b2 <- b / nu2

        p <- p(L2deriv)(z)

        if(is(L2deriv, "AbscontDistribution"))
            ws0 <- 0
        else
            ws0 <- d(L2deriv)(z)
        if(ws0 > 0)
            d <- (-b2*(1-p)+b1*(p-ws0))/ws0/b
        else
            d <- 0

        info <- c("minimum asymptotic bias (lower case) solution")
        asCov <- b2^2*(1-p)+b1^2*(p-ws0) + b^2*d^2*ws0
        Risk <- list(asBias = b, asCov = asCov)

        w <- new("HampelWeight")
        cent(w) <- z
        stand(w) <- A
        clip(w) <- c(b1,b2)
        weight(w) <- minbiasweight(w, neighbor = neighbor, biastype = biastype)

        return(list(A = A, a = zi*z, b = b, d = d, risk = Risk, info = info, 
                    w = w, biastype = biastype, normtype = NormType()))
           })

setMethod("minmaxBias", signature(L2deriv = "UnivariateDistribution", 
                                   neighbor = "ContNeighborhood", 
                                   biastype = "onesidedBias"),
    function(L2deriv, neighbor, biastype, symm, 
             trafo, maxiter, tol){                
        infotxt <- c("minimum asymptotic bias (lower case) solution")
        noIC <- function(){
                warntxt <- gettext("There exists no IC attaining the infimal maxBias.")
                warning(warntxt)
                return(list(A = matrix(1), a = 1, d = 0, b = 0,
                       Risk = list(asBias = 0, asCov = 0, warning = warntxt),
                       info = infotxt, w = NULL, biastype = biastype, 
                       normtype = NormType()))}
        if(is(L2deriv, "DiscreteDistribution"))
           { if(is.finite(lowerCaseRadius(L2deriv, neighbor, risk = asMSE(), biastype)))
                {
                 sign <- sign(biastype)
                 w0 <- options("warn")
                 options(warn = -1)
        
                 l <- length(support(L2deriv))
                 if (sign>0)
                      {z0 <- support(L2deriv)[1] 
                       deltahat <- support(L2deriv)[2]-z0
                 }else{
                       z0 <- support(L2deriv)[l]
                       deltahat <- z0-support(L2deriv)[l-1]
                 }
                 p0 <- d(L2deriv)(z0)   
                 v1 <- (1-p0)/p0/z0
                 v2 <- -1/z0
                 c0 <- deltahat*p0/2
                 A0 <- abs(1/z0/c0)
                 zc <- z0+sign(biastype)*deltahat*(1-p0)/2
                 a0 <- A0*zc
                 b0 <- abs(1/z0)
                 d0  <- 0 
                 asCov <- v1^2*(1-p0)+v2^2*p0
                 Risk0 <- list(asBias = b0, asCov = asCov)
                 A0 <- matrix(A0,1,1)
                 
                 w <- new("HampelWeight")
                 cent(w) <- z0
                 stand(w) <- A0
                 clip(w) <- b0
                 weight(w) <- minbiasweight(w, neighbor = neighbor, 
                               biastype = biastype)
                               
                }else{return(noIC())}
            }else{return(noIC())}                    
        return(list(A = A0, a = a0, b = b0, d = d0, risk = Risk0, 
                    info = infotxt, w = w, biastype = biastype, 
                    normtype = NormType()))
           })

           