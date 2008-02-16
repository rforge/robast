###############################################################################
## centering constant for asymptotic MSE and asymptotic Hampel risk 
###############################################################################
setMethod("getInfCent", signature(L2deriv = "UnivariateDistribution",
                                  neighbor = "ContNeighborhood",
                                  biastype = "BiasType"),
    function(L2deriv, neighbor, biastype = symmetricBias(), 
             clip, cent, tol.z, symm, trafo){
        if(symm) return(0)

        z.fct <- function(z, c0, D1){
            return(c0 + (z-c0)*p(D1)(z-c0) - (z+c0)*p(D1)(z+c0) + m1df(D1, z+c0) - m1df(D1, z-c0))
        }
        lower <- q(L2deriv)(getdistrOption("TruncQuantile"))
        upper <- q(L2deriv)(1-getdistrOption("TruncQuantile"))

        return(uniroot(z.fct, lower = lower, upper = upper, tol = tol.z, 
                    c0=clip, D1=L2deriv)$root)
    })
setMethod("getInfCent", signature(L2deriv = "UnivariateDistribution",
                                  neighbor = "TotalVarNeighborhood",
                                  biastype = "BiasType"),
    function(L2deriv, neighbor, biastype = symmetricBias(), 
             clip, cent, tol.z, symm, trafo){
        if(symm) return(-clip/2)

        D1 <- sign(as.vector(trafo))*L2deriv
        g.fct <- function(g, c0, D1){
            return(g*p(D1)(g) + (g+c0)*(1-p(D1)(g+c0)) - m1df(D1, g) + m1df(D1, g+c0))
        }
        lower <- q(L2deriv)(getdistrOption("TruncQuantile"))
        upper <- q(L2deriv)(1-getdistrOption("TruncQuantile"))

        return(uniroot(g.fct, lower = lower, upper = upper, tol = tol.z, 
                    c0 = clip, D1 = D1)$root)
    })
setMethod("getInfCent", signature(L2deriv = "RealRandVariable",
                                  neighbor = "ContNeighborhood",
                                  biastype = "BiasType"),
    function(L2deriv, neighbor, biastype = symmetricBias(), 
             Distr, z.comp, stand, cent, clip){
        integrand1 <- function(x, L2, clip, cent, stand){
            X <- evalRandVar(L2, as.matrix(x))[,,1] - cent
            Y <- apply(X, 2, "%*%", t(stand)) 
            h.vct <- sqrt(colSums(Y^2))
            ind2 <- (h.vct < clip/2)
            h.vct <- ind2*clip/2 + (1-ind2)*h.vct
            ind1 <- (h.vct < clip)

            return(ind1 + (1-ind1)*clip/h.vct)
        }
        integrand2 <- function(x, L2.i, L2, clip, cent, stand){
            return(L2.i(x)*integrand1(x = x, L2 = L2, clip = clip, cent = cent, stand = stand))
        }

        res1 <- E(object = Distr, fun = integrand1, L2 = L2deriv, clip = clip, 
                  cent = cent, stand = stand, useApply = FALSE)
        nrvalues <- length(L2deriv)
        res2 <- numeric(nrvalues)
        for(i in 1:nrvalues){
            if(z.comp[i]){
                res2[i] <- E(object = Distr, fun = integrand2, L2.i = L2deriv@Map[[i]], 
                             L2 = L2deriv, clip = clip, cent = cent, stand = stand,
                             useApply = FALSE)
            }else{            
                res2[i] <- 0
            }
        }

        return(res2/res1)
    })
###############################################################################
## centering constant for asymptotic one-sided MSE and asymptotic one-sided Hampel risk
###############################################################################
setMethod("getInfCent", signature(L2deriv = "UnivariateDistribution",
                                  neighbor = "ContNeighborhood",
                                  biastype = "onesidedBias"),
    function(L2deriv, neighbor, biastype = positiveBias(), clip, cent, tol.z, symm, trafo){
        if (sign(biastype)> 0){
        z.fct <- function(z, c0, D1){
            return(c0 - (z+c0)*p(D1)(z+c0) + m1df(D1, z+c0))
        }
        lower <- q(L2deriv)(getdistrOption("TruncQuantile"))
        upper <- 0
        }else{
        z.fct <- function(z, c0, D1){
            return(- z + (z-c0)*p(D1)(z-c0) - m1df(D1, z-c0))
        }
        lower <- 0
        upper <- q(L2deriv)(1-getdistrOption("TruncQuantile"))
        }
        return(uniroot(z.fct, lower = lower, upper = upper, tol = tol.z,
                    c0=clip, D1=L2deriv)$root)
    })

setMethod("getInfCent", signature(L2deriv = "UnivariateDistribution",
                                  neighbor = "ContNeighborhood",
                                  biastype = "asymmetricBias"),
    function(L2deriv, neighbor, biastype = asymmetricBias(), clip, cent, tol.z, symm, trafo){
        nu1 <- nu(biastype)[1]
        nu2 <- nu(biastype)[2]

        z.fct <- function(z, c0, D1){
            return(c0/nu2 + (z-c0/nu1)*p(D1)(z-c0/nu1) -
                   (z+c0/nu2)*p(D1)(z+c0/nu2) + m1df(D1, z+c0/nu2) -
                   m1df(D1, z-c0/nu1))
        }
        lower <- q(L2deriv)(getdistrOption("TruncQuantile"))
        upper <- q(L2deriv)(1-getdistrOption("TruncQuantile"))

        return(uniroot(z.fct, lower = lower, upper = upper, tol = tol.z,
                    c0=clip, D1=L2deriv)$root)
    })
