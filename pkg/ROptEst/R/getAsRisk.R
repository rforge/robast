###############################################################################
## asymptotic MSE
###############################################################################
setMethod("getAsRisk", signature(risk = "asMSE",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "Neighborhood", biastype = "BiasType"),
    function(risk, L2deriv, neighbor, biastype = symmetricBias(), 
             clip = NULL, cent = NULL, stand, trafo){
        if(!is.finite(neighbor@radius))
            mse <- Inf
        else
            mse <- as.vector(stand)*as.vector(trafo)
        return(list(asMSE = mse))
    })
setMethod("getAsRisk", signature(risk = "asMSE",
                                 L2deriv = "EuclRandVariable",
                                 neighbor = "Neighborhood", biastype = "BiasType"),
    function(risk, L2deriv, neighbor, biastype = symmetricBias(), 
             clip = NULL, cent = NULL, stand, trafo){
        if(!is.finite(neighbor@radius))
            mse <- Inf
        else
            mse <- sum(diag(stand %*% t(trafo)))
        return(list(asMSE = mse))
    })

###############################################################################
## minimum asymptotic Bias
###############################################################################
setMethod("getAsRisk", signature(risk = "asBias",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "ContNeighborhood", biastype = "BiasType"),
    function(risk, L2deriv, neighbor, biastype = symmetricBias(), trafo){
        z <- q(L2deriv)(0.5)                                
        bias <- abs(as.vector(trafo))/E(L2deriv, function(x, z){abs(x - z)}, 
                                        useApply = FALSE, z = z)

        return(list(asBias = bias))
    })
setMethod("getAsRisk", signature(risk = "asBias",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "TotalVarNeighborhood", biastype = "BiasType"),
    function(risk, L2deriv, neighbor, biastype = symmetricBias(), trafo){
        bias <- abs(as.vector(trafo))/(-m1df(L2deriv, 0))

        return(list(asBias = bias))
    })
setMethod("getAsRisk", signature(risk = "asBias",
                                 L2deriv = "RealRandVariable",
                                 neighbor = "ContNeighborhood", biastype = "BiasType"),
    function(risk, L2deriv, neighbor, biastype = symmetricBias(), Distr, 
             L2derivDistrSymm, trafo, 
             z.start, A.start,  maxiter, tol){                
        if(is.null(z.start)) z.start <- numeric(ncol(trafo))
        if(is.null(A.start)) A.start <- trafo

        abs.fct <- function(x, L2, stand, cent){
            X <- evalRandVar(L2, as.matrix(x))[,,1] - cent
            Y <- apply(X, 2, "%*%", t(stand)) 

            return(sqrt(colSums(Y^2)))
        }
        bmin.fct <- function(param, L2deriv, Distr, trafo, z.comp){
            p <- nrow(trafo)
            k <- ncol(trafo)
            A <- matrix(param[1:(p*k)], ncol=k, nrow=p)
            z <- numeric(k)
            z[z.comp] <- param[(p*k+1):length(param)]

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
        erg <- optim(c(A.vec, z.start[z.comp]), bmin.fct, method = "Nelder-Mead", 
                    control = list(reltol = tol, maxit = 100*maxiter), 
                    L2deriv = L2deriv, Distr = Distr, trafo = trafo, z.comp = z.comp)
        bias <- 1/erg$value
        
        return(list(asBias = bias))
    })

###############################################################################
## asymptotic covariance
###############################################################################
setMethod("getAsRisk", signature(risk = "asCov",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "ContNeighborhood", biastype = "BiasType"),
    function(risk, L2deriv, neighbor, biastype = symmetricBias(), clip, cent, stand){
        c0 <- clip/abs(as.vector(stand))
        D1 <- L2deriv - cent/as.vector(stand)
        Cov <- (clip^2*(p(D1)(-c0) + 1 - p(D1)(c0))
                + as.vector(stand)^2*(m2df(D1, c0) - m2df(D1, -c0)))

        return(list(asCov = Cov))
    })
setMethod("getAsRisk", signature(risk = "asCov",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "TotalVarNeighborhood", biastype = "BiasType"),
    function(risk, L2deriv, neighbor, biastype = symmetricBias(), clip, cent, stand){
        g0 <- cent/abs(as.vector(stand))
        c0 <- clip/abs(as.vector(stand))
        Cov <- (abs(as.vector(stand))^2*(g0^2*p(L2deriv)(g0) 
                + (g0+c0)^2*(1 - p(L2deriv)(g0+c0))
                + m2df(L2deriv, g0+c0) - m2df(L2deriv, g0)))

        return(list(asCov = Cov))
    })
setMethod("getAsRisk", signature(risk = "asCov",
                                 L2deriv = "RealRandVariable",
                                 neighbor = "ContNeighborhood", biastype = "BiasType"),
    function(risk, L2deriv, neighbor, biastype = symmetricBias(), Distr, clip, cent, 
             stand, norm = EuclideanNorm){

                 return(list(asCov = .asCovMB(L2deriv, stand, cent, clip, Distr, 
                             norm = norm)))
        })


#        Y <- as(stand %*% L2deriv - cent, "EuclRandVariable")
#        absY <- sqrt(Y %*% Y)
#        
#        nrvalues <- nrow(stand)
#        ICfct <- vector(mode = "list", length = nrvalues)
#        for(i in 1:nrvalues){
#            ICfct[[i]] <- function(x){ Yi(x)*pmin(1, b/absY(x)) }
#            body(ICfct[[i]]) <- substitute({ Yi(x)*pmin(1, b/absY(x)) },
#                                    list(Yi = Y@Map[[i]], absY = absY@Map[[1]], b = clip))
#        }
#        IC <- RealRandVariable(Map = ICfct, Domain = Y@Domain, Range = Y@Range)
#        Cov <- matrix(E(Distr, IC %*% t(IC)), ncol = nrvalues)
#
#        return(list(asCov = Cov))
#    })

### helping function

.asCovMB <- function(L2, stand, cent, clip, Distr, norm){                 
                 p <- nrow(stand)
                 idx <- matrix(1:p^2,p,p)
                 idx <- idx[col(idx)<=row(idx)]
                 Cv <- matrix(0,p,p)

                 if (clip == 0){
                     Cv[idx] <- E(object = Distr, fun = function(x){
                                X <- evalRandVar(L2, as.matrix(x))[,,1] - cent
                                Y <- stand %*% X
                                norm0 <- norm(Y)                      
                                ind <- 1-.eq(norm0)                   
                                Y0 <- Y*ind/(norm0+1-ind)
                                Y02 <- apply(Y0,2,function(x)x%*%t(x))[idx,]
                                }, useApply = FALSE)
                 }else{
                    Cv[idx] <- E(object = Distr, fun = function(x){
                               X <- evalRandVar(L2, as.matrix(x))[,,1] - cent
                               Y <- stand %*% X
                               norm0 <- norm(Y)                      
                               ind2 <- (norm0 < b/2)
                               norm1 <- ind2*clip/2 + (1-ind2)*norm0
                               ind1 <- (norm0 < b)
                               ind1 + (1-ind1)*clip/norm1
                               Y0 <- Y*ind1
                               Y02 <- apply(Y0,2,function(x)x%*%t(x))[idx,]
                       }, useApply = FALSE)
                 }
                 dCv <- diag(Cv)
                 return(PosSemDefSymmMatrix(Cv + t(Cv) - dCv))
        }


###############################################################################
## trace of asymptotic covariance
###############################################################################
setMethod("getAsRisk", signature(risk = "trAsCov",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "UncondNeighborhood", biastype = "BiasType"),
    function(risk, L2deriv, neighbor, biastype = symmetricBias(), clip, cent, stand){
        Cov <- getAsRisk(risk = asCov(), L2deriv = L2deriv, neighbor = neighbor,
                         biastype = biastype, clip = clip, cent = cent, stand = stand)$asCov

        return(list(trAsCov = as.vector(Cov)))
    })
setMethod("getAsRisk", signature(risk = "trAsCov",
                                 L2deriv = "RealRandVariable",
                                 neighbor = "ContNeighborhood", biastype = "BiasType"),
    function(risk, L2deriv, neighbor, biastype = symmetricBias(), Distr, clip, cent, stand){
        Cov <- getAsRisk(risk = asCov(), L2deriv = L2deriv, neighbor = neighbor,
                         biastype = biastype, Distr = Distr, clip = clip, 
                         cent = cent, stand = stand)$asCov

        return(list(trAsCov = sum(diag(Cov))))
    })

###############################################################################
## asymptotic under-/overshoot risk
###############################################################################
setMethod("getAsRisk", signature(risk = "asUnOvShoot",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "UncondNeighborhood", biastype = "BiasType"),
    function(risk, L2deriv, neighbor, biastype = symmetricBias(), clip, cent, stand, trafo){
        if(identical(all.equal(neighbor@radius, 0), TRUE))
            return(list(asUnOvShoot = pnorm(-risk@width/sqrt(as.vector(stand)))))
        
        g0 <- cent/abs(as.vector(stand))
        c0 <- clip/abs(as.vector(stand))
        s <- sqrt(g0^2*p(L2deriv)(g0) 
                  + (g0+c0)^2*(1 - p(L2deriv)(g0+c0))
                  + m2df(L2deriv, g0+c0) - m2df(L2deriv, g0))

        return(list(asUnOvShoot = pnorm(-risk@width*s)))
    })

###############################################################################
## asymptotic onesided bias
###############################################################################
setMethod("getAsRisk", signature(risk = "asBias",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "ContNeighborhood", biastype = "onesidedBias"),
    function(risk, L2deriv, neighbor, biastype, trafo){

        D1 <- L2deriv
        if(!is(D1, "DiscreteDistribution")) 
            return(list(asBias = 0, warn = gettext("not attained by IC")))

        sign <- sign(biastype)
        w0 <- options("warn")
        options(warn = -1)
        
        l <- length(support(L2deriv))
        if (sign>0)
           {z0 <- support(L2deriv)[1]; deltahat <- support(L2deriv)[2]-z0}
        else
           {z0 <- support(L2deriv)[l]; deltahat <- z0-support(L2deriv)[l-1]}

        bias <- abs(as.vector(trafo))/abs(z0)
        return(list(asBias = bias))
    })

###############################################################################
## asymptotic asymmetric bias
###############################################################################

setMethod("getAsRisk", signature(risk = "asBias",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "ContNeighborhood", 
                                 biastype = "asymmetricBias"),
    function(risk, L2deriv, neighbor, biastype, trafo){
        nu1 <- nu(biastype)[1]
        nu2 <- nu(biastype)[2]
        num <- nu2/(nu1+nu2)        
        z <- q(L2deriv)(num)
        Int <- E(L2deriv, function(x, m){abs(x-m)}, m = z)
        omega <- 2/(Int/nu1+Int/nu2)
        bias <- abs(as.vector(trafo))*omega
        return(list(asBias = bias))
    })

###############################################################################
## asymptotic semivariance
###############################################################################

setMethod("getAsRisk", signature(risk = "asSemivar",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "Neighborhood", 
                                 biastype = "onesidedBias"),
    function(risk, L2deriv, neighbor, biastype = positiveBias(), 
             clip, cent, stand, trafo){
        A <- as.vector(stand)*as.vector(trafo)
        r <- neighbor@radius
        b <- clip*A
        if (sign(biastype)>0)
            v <- E(L2deriv, function(x) A^2*pmin(x-cent,clip)^2)
        else
            v <- E(L2deriv, function(x) A^2*pmax(x-cent,-clip)^2)
        sv <- r*b/sqrt(v)
        if(!is.finite(r))
            semvar <- Inf
        else
            semvar <- (v+r^2*b^2)*pnorm(sv)+ r*b*sqrt(v)*dnorm(sv)
        return(list(asSemivar = semvar))
    })
