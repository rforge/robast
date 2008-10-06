###############################################################################
## get classical optimal IC
###############################################################################
setMethod("getInfRobIC", signature(L2deriv = "UnivariateDistribution", 
                                   risk = "asCov", 
                                   neighbor = "ContNeighborhood"),
    function(L2deriv, risk, neighbor, Finfo, trafo, verbose = FALSE){
            info <- c("optimal IC in sense of Cramer-Rao bound")
            A <- trafo %*% solve(Finfo)
            
            b <- abs(as.vector(A))*max(abs(q(L2deriv)(1)),abs(q(L2deriv)(0)))
            
            asCov <- A %*% t(trafo)
            r <- neighbor@radius
            Risk <- list(asCov = asCov, 
                         asBias = list(value = b, biastype = symmetricBias(), 
                                       normtype = NormType(), 
                                       neighbortype = class(neighbor)), 
                         trAsCov = list(value = asCov, normtype = NormType()),
                         asMSE = list(value = asCov + r^2*b^2, 
                                      r = r,
                                      at = neighbor))

            return(list(A = A, a = 0, b = b, d = NULL, risk = Risk, info = info))
    })
setMethod("getInfRobIC", signature(L2deriv = "UnivariateDistribution", 
                                   risk = "asCov", 
                                   neighbor = "TotalVarNeighborhood"),
    function(L2deriv, risk, neighbor, Finfo, trafo, verbose = FALSE){
            info <- c("optimal IC in sense of Cramer-Rao bound")
            A <- trafo %*% solve(Finfo)
            b <- abs(as.vector(A))*(q(L2deriv)(1)-q(L2deriv)(0))
            asCov <- A %*% t(trafo)
            r <- neighbor@radius
            Risk <- list(asCov = asCov, 
                         asBias = list(value = b, biastype = symmetricBias(), 
                                       normtype = NormType(), 
                                       neighbortype = class(neighbor)), 
                         trAsCov = list(value = asCov, normtype = NormType()),
                         asMSE = list(value = asCov + r^2*b^2, 
                                      r = r,
                                      at = neighbor))

            return(list(A = A, a = -b/2, b = b, d = NULL, risk = Risk, info = info))
    })
setMethod("getInfRobIC", signature(L2deriv = "RealRandVariable", 
                                   risk = "asCov", 
                                   neighbor = "ContNeighborhood"),
    function(L2deriv, risk, neighbor, Distr, Finfo, trafo, 
             QuadForm = diag(nrow(trafo)), verbose = FALSE){
            info <- c("optimal IC in sense of Cramer-Rao bound")
            A <- trafo %*% solve(Finfo)
            IC <- A %*% L2deriv
            if(is(Distr, "UnivariateDistribution")){
                lower <- ifelse(is.finite(q(Distr)(0)), q(Distr)(1e-8), q(Distr)(0))
                upper <- ifelse(is.finite(q(Distr)(1)), q(Distr)(1-1e-8), q(Distr)(1))
                x <- seq(from = lower, to = upper, length = 1e5)
                x <- x[x!=0] # problems with NaN=log(0)!
                b <- evalRandVar(IC, as.matrix(x))^2
                b <- sqrt(max(colSums(b, na.rm = TRUE)))
            }else{
                b <- Inf # not yet implemented
            }
            
            asCov <- A %*% t(trafo)
            trAsCov <- sum(diag(QuadForm%*%asCov))
            r <- neighbor@radius
            nt <- if(identical(QuadForm,diag(nrow(trafo)))) NormType() else 
                     QFNorm(QuadForm = PosSemDefSymmMatrix(QuadForm))
            Risk <- list(asCov = asCov,  
                         asBias = list(value = b, biastype = symmetricBias(), 
                                       normtype = nt, 
                                       neighbortype = class(neighbor)), 
                         trAsCov = list(value = trAsCov, normtype = nt),
                         asMSE = list(value = trAsCov + r^2*b^2, 
                                      r = r,
                                      at = neighbor))
            return(list(A = A, a = numeric(nrow(trafo)), b = b, d = NULL, risk = Risk, info = info))
    })
