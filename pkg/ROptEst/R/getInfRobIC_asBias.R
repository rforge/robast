###############################################################################
## get minimum bias solutions
###############################################################################
setMethod("getInfRobIC", signature(L2deriv = "UnivariateDistribution", 
                                   risk = "asBias", 
                                   neighbor = "UncondNeighborhood"),
    function(L2deriv, risk, neighbor, symm, Finfo, trafo, upper, maxiter, 
             tol, warn){
        minmaxBias(L2deriv, neighbor, biastype(risk), symm, 
                   Finfo, trafo, upper, maxiter, tol, warn)
    })
setMethod("getInfRobIC", signature(L2deriv = "RealRandVariable", 
                                   risk = "asBias", 
                                   neighbor = "ContNeighborhood"),
    function(L2deriv, risk, neighbor, Distr, DistrSymm, L2derivSymm, 
             L2derivDistrSymm, Finfo, z.start, A.start, trafo, upper, 
             maxiter, tol, warn){                
        minmaxBias(L2deriv, neighbor, biastype(risk), normtype(risk),
             Distr, DistrSymm, L2derivSymm, L2derivDistrSymm, Finfo, 
             z.start, A.start, trafo, upper, 
             maxiter, tol, warn)
    })


setMethod("minmaxBias", signature(L2deriv = "UnivariateDistribution", 
                                   neighbor = "ContNeighborhood",
                                   biastype = "BiasType"),
    function(L2deriv, neighbor, biastype = symmetricBias(), symm, 
             Finfo, trafo, upper, maxiter, tol, warn){
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
                                   normtype = NormType())

        return(list(A = A, a = zi*z, b = b, d = d, risk = Risk, info = info, 
                    w = w, biastype = biastype, normtype = NormType()))    
    })

setMethod("minmaxBias", signature(L2deriv = "UnivariateDistribution", 
                                   neighbor = "TotalVarNeighborhood",
                                   biastype = "BiasType"),
    function(L2deriv, neighbor, biastype = symmetricBias(),
             symm, Finfo, trafo, 
             upper, maxiter, tol, warn){
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
            a <- b*(p0-ws0)/(1-ws0)
            
        info <- c("minimum asymptotic bias (lower case) solution")
        Risk <- list(asCov = a^2*(p0-ws0) + (zi*a+b)^2*(1-p0), asBias = b)

        w <- new("BdStWeight")
        stand(w) <- A
        clip(w) <- c(a, a+b)
        weight(w) <- minbiasweight(w, neighbor = neighbor, biastype = biastype, 
                               normtype = NormType())

        return(list(A = A, a = a, b = b, d = 1, risk = Risk, info = info, 
                    w = w, biastype = biastype, normtype = NormType()))
    })

setMethod("minmaxBias", signature(L2deriv = "RealRandVariable", 
                                   neighbor = "ContNeighborhood", 
                                   biastype = "BiasType"),
    function(L2deriv, neighbor, biastype = symmetricBias(), 
             normtype = NormType(), Distr, DistrSymm, L2derivSymm, 
             L2derivDistrSymm, Finfo, z.start, A.start, trafo, upper, 
             maxiter, tol, warn){                
        if(is.null(z.start)) z.start <- numeric(ncol(trafo))
        if(is.null(A.start)) A.start <- trafo


        
        abs.fct <- function(x, L2, stand, cent, norm){ 
            X <- evalRandVar(L2, as.matrix(x))[,,1] - cent
            Y <- stand %*% X
            return(fct(norm)(Y))
        }

        bmin.fct <- function(param, L2deriv, Distr, trafo, z.comp){
            p <- nrow(trafo)
            k <- ncol(trafo)
            A <- matrix(param[1:(p*k)], ncol=k, nrow=p)
            z <- numeric(k)
            z[z.comp] <- param[(p*k+1):length(param)]
            aL <- list(normtype = normtype, FI = Finfo, 
                       L2 = L2deriv, stand = A, cent = z, clip = 0, 
                       Distr = Distr, norm = fct(normtype))
            normtype <<- do.call(updateNorm, aL) 
            return(E(object = Distr, fun = abs.fct, L2 = L2deriv, stand = A, 
                     cent = z, useApply = FALSE)/sum(diag(A %*% t(trafo))))
        }
        
        nrvalues <- length(L2deriv)
        z.comp <- rep(TRUE, nrvalues)
        for(i in 1:nrvalues)
            if(is(L2derivDistrSymm[[i]], "SphericalSymmetry"))
                if(L2derivDistrSymm[[i]]@SymmCenter == 0)
                    z.comp[i] <- FALSE

        A.vec <- as.vector(A.start)
        force(normtype)
        erg <- optim(c(A.vec, z.start[z.comp]), bmin.fct, method = "Nelder-Mead", 
                    control = list(reltol = tol, maxit = 100*maxiter), 
                    L2deriv = L2deriv, Distr = Distr, trafo = trafo, z.comp = z.comp)
        b <- 1/erg$value
        param <- erg$par
        p <- nrow(trafo)
        k <- ncol(trafo)
        A <- matrix(param[1:(p*k)], ncol=k, nrow=p)
        z <- numeric(k)
        z[z.comp] <- param[(p*k+1):length(param)]
        a <- as.vector(A %*% z)
        d <- numeric(p)
        # computation of 'd', in case 'L2derivDistr' not abs. cont.
        
        info <- c("minimum asymptotic bias (lower case) solution")
        Risk <- list(asBias = b)

        w <- new("HampelWeight")
        cent(w) <- z
        stand(w) <- A
        clip(w) <- b
        weight(w) <- minbiasweight(w, neighbor = neighbor, biastype = biastype, 
                               normtype = normtype)

        return(list(A = A, a = a, b = b, d = d, risk = Risk, info = info, 
                    w = w, biastype = biastype, normtype = normtype))
    })

setMethod("minmaxBias", signature(L2deriv = "UnivariateDistribution", 
                                   neighbor = "ContNeighborhood", 
                                   biastype = "asymmetricBias"),
    function(L2deriv, neighbor, biastype, symm, 
             Finfo, trafo, upper, maxiter, tol, warn){                
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
             Finfo, trafo, upper, maxiter, tol, warn){                

        infotxt <- c("minimum asymptotic bias (lower case) solution")
        noIC <- function(){
                warntxt <- gettext("There exists no IC attaining the infimal maxBias.")
                warning(warntxt)
                return(list(A = 1, a = 1, d = 0, b = 0,
                       Risk = list(asBias = 0, asCov = 0, warning = warntxt),
                       info = infotxt, w = NULL, biastype = biastype, 
                       normtype = NormType()))}
        if(!is(L2deriv, "DiscreteDistribution"))
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

           