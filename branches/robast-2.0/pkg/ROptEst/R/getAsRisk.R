###############################################################################
## asymptotic MSE
###############################################################################           
setMethod("getAsRisk", signature(risk = "asMSE",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "Neighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, clip = NULL, cent = NULL, stand, trafo){
        if(!is.finite(neighbor@radius))
            mse <- Inf
        else
            mse <- as.vector(stand)*as.vector(trafo)
        return(list(asMSE = mse))
    })

setMethod("getAsRisk", signature(risk = "asMSE",
                                 L2deriv = "EuclRandVariable",
                                 neighbor = "Neighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, clip = NULL, cent = NULL, stand, trafo){
        if(!is.finite(neighbor@radius))
            mse <- Inf
        else{
            p <- nrow(trafo)
            std <- if(is(normtype(risk),"QFNorm")) 
                      QuadForm(normtype(risk)) else diag(p)
                        
            mse <- sum(diag( std %*% stand %*% t(trafo)))
        }
        return(list(asMSE = mse))
    })

###############################################################################
## minimum asymptotic Bias
###############################################################################
setMethod("getAsRisk", signature(risk = "asBias",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "ContNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, trafo){
        z <- q(L2deriv)(0.5)                                
        bias <- abs(as.vector(trafo))/E(L2deriv, function(x, z){abs(x - z)}, 
                                        useApply = FALSE, z = z)

        return(list(asBias = bias))
    })
setMethod("getAsRisk", signature(risk = "asBias",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "TotalVarNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype,  trafo){
        bias <- abs(as.vector(trafo))/(-m1df(L2deriv, 0))

        return(list(asBias = bias))
    })
setMethod("getAsRisk", signature(risk = "asBias",
                                 L2deriv = "RealRandVariable",
                                 neighbor = "ContNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, Distr, DistrSymm, L2derivSymm,
             L2derivDistrSymm, trafo, z.start, A.start,  maxiter, tol, warn){                
        
        normtype <- normtype(risk)
        biastype <- biastype(risk)

        if(is(normtype,"SelfNorm")){
                warntxt <- paste(gettext(
                "Using self-standardization, there are problems with the existence\n"
                               ),gettext(
                "of a minmax Bias IC. Instead we return a lower bound.\n"
                               ))
                if(warn) cat(warntxt)
                return(list(asBias = sqrt(nrow(trafo))))        
        }
        comp <- .getComp(L2deriv, DistrSymm, L2derivSymm, L2derivDistrSymm)
        z.comp <- comp$"z.comp"
        A.comp <- comp$"A.comp"
        
        eerg <- .LowerCaseMultivariate(L2deriv = L2deriv, neighbor = neighbor, 
             biastype = biastype, normtype = normtype, Distr = Distr, 
             trafo = trafo, z.start = z.start, A.start, z.comp = z.comp, 
             A.comp = A.comp,  maxiter = maxiter, tol = tol)
        erg <- eerg$erg
        bias <- 1/erg$value
        
        return(list(asBias = bias, normtype = eerg$normtype))
    })


###############################################################################
## asymptotic covariance
###############################################################################
setMethod("getAsRisk", signature(risk = "asCov",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "ContNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, clip, cent, stand){
#        c0 <- clip/abs(as.vector(stand))
#        D1 <- L2deriv - cent/as.vector(stand)
#        Cov <- (clip^2*(p(D1)(-c0) + 1 - p(D1)(c0))
#                + as.vector(stand)^2*(m2df(D1, c0) - m2df(D1, -c0)))
        return(list(asCov = 
        getInfV(L2deriv, neighbor, biastype(risk), clip/abs(as.vector(stand)), 
                cent/abs(as.vector(stand)), abs(as.vector(stand)))
               ))
    })
setMethod("getAsRisk", signature(risk = "asCov",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "TotalVarNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, clip, cent, stand){
#        g0 <- cent/abs(as.vector(stand))
#        c0 <- clip/abs(as.vector(stand))
#        Cov <- (abs(as.vector(stand))^2*(g0^2*p(L2deriv)(g0) 
#                + (g0+c0)^2*(1 - p(L2deriv)(g0+c0))
#                + m2df(L2deriv, g0+c0) - m2df(L2deriv, g0)))
#        return(list(asCov = Cov))
        return(list(asCov = 
        getInfV(L2deriv, neighbor, biastype(risk), clip/abs(as.vector(stand)), 
                cent/abs(as.vector(stand)), abs(as.vector(stand)))
               ))

    })
setMethod("getAsRisk", signature(risk = "asCov",
                                 L2deriv = "RealRandVariable",
                                 neighbor = "ContNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, Distr, cent, 
             stand, V.comp =  matrix(TRUE, ncol = nrow(stand), nrow = nrow(stand)), 
             w){
        return(getInfV(L2deriv = L2deriv, neighbor = neighbor, 
                       biastype = biastype(risk), Distr = Distr, 
                       V.comp = V.comp, cent = cent, 
                       stand = stand, w = w))
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



###############################################################################
## trace of asymptotic covariance
###############################################################################
setMethod("getAsRisk", signature(risk = "trAsCov",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "UncondNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, clip, cent, stand){
        Cov <- getAsRisk(risk = asCov(), L2deriv = L2deriv, neighbor = neighbor,
                         biastype = biastype(risk), clip = clip, cent = cent, stand = stand)$asCov

        return(list(trAsCov = as.vector(Cov)))
    })
setMethod("getAsRisk", signature(risk = "trAsCov",
                                 L2deriv = "RealRandVariable",
                                 neighbor = "ContNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, Distr, clip, cent, stand, normtype){
        Cov <- getAsRisk(risk = asCov(), L2deriv = L2deriv, neighbor = neighbor,
                         biastype = biastype(risk), Distr = Distr, clip = clip, 
                         cent = cent, stand = stand)$asCov

        p <- nrow(stand)
        std <- if(is(normtype,"QFNorm")) QuadForm(normtype) else diag(p)
        
        return(list(trAsCov = sum(diag(std%*%Cov))))
    })

###############################################################################
## asymptotic under-/overshoot risk
###############################################################################
setMethod("getAsRisk", signature(risk = "asUnOvShoot",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "UncondNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, clip, cent, stand, trafo){
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
                                 neighbor = "ContNeighborhood", 
                                 biastype = "onesidedBias"),
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
    function(risk, L2deriv, neighbor, biastype, 
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
