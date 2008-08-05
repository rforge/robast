###############################################################################
## IC algorithm for asymptotic Hampel risk
###############################################################################
setMethod("getInfRobIC", signature(L2deriv = "UnivariateDistribution", 
                                   risk = "asHampel", 
                                   neighbor = "UncondNeighborhood"),
    function(L2deriv, risk, neighbor, symm, Finfo, trafo, 
             upper, maxiter, tol, warn, noLow = FALSE, verbose = FALSE){
        biastype <- biastype(risk)
        normtype <- normtype(risk)

        A <- trafo / E(L2deriv, function(x){x^2})
        b <- risk@bound

        bmax <- abs(as.vector(A))*max(abs(q(L2deriv)(0)), q(L2deriv)(1))
        if(b >= bmax){
            if(warn) cat("'b >= maximum asymptotic bias' => (classical) optimal IC\n", 
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobIC(L2deriv = L2deriv, risk = asCov(), 
                               neighbor = neighbor, Finfo = Finfo, trafo = trafo,
                               verbose = verbose)
            res <- c(res, list(biastype = biastype, normtype = NormType()))
            Cov <- res$risk$asCov
            r <- neighbor@radius
            res$risk$asBias <- list(value = b, biastype = biastype, 
                                   normtype = normtype, 
                                   neighbortype = class(neighbor))
            res$risk$asMSE <- list(value = Cov + r^2*b^2, 
                                   r = r,
                                   at = neighbor)
            return(res)
        }

        if(!noLow){
            res <- getInfRobIC(L2deriv = L2deriv, risk = asBias(biastype = biastype), 
                               neighbor = neighbor, symm = symm,  
                               trafo = trafo, maxiter = maxiter, tol = tol, Finfo = Finfo,
                               warn = warn, verbose = verbose)
            bmin <- res$b
            cat("minimal bound:\t", bmin, "\n")
         }else bmin <- b/2

        if(b <= bmin){
            if(warn) cat("'b <= minimum asymptotic bias'\n",
                         "=> the minimum asymptotic bias (lower case) solution is returned\n")
            Risk <- list(asMSE = res$risk$asCov + neighbor@radius^2*bmin^2)
            res$risk <- c(Risk, res$risk)
            return(res)
        }
#        bmin <- getAsRisk(risk = asBias(biastype = biastype, normtype = normtype), 
#                          L2deriv = L2deriv, neighbor = neighbor, 
#                          biastype = biastype, trafo = trafo, Finfo = Finfo,
#                          warn = warn)$asBias
#        if(b <= bmin){
#            if(warn) cat("'b <= minimum asymptotic bias'\n",
#                         "=> the minimum asymptotic bias (lower case) solution is returned\n")
#            res <- getInfRobIC(L2deriv = L2deriv, risk = asBias(biastype = biastype), 
#                            neighbor = neighbor, symm = symm,  
#                            trafo = trafo, maxiter = maxiter, tol = tol, Finfo = Finfo,
#                            warn = warn)
#            Risk <- list(asMSE = res$risk$asCov + neighbor@radius^2*bmin^2)
#            res$risk <- c(Risk, res$risk)
#            return(res)
#        }
        c0 <- b/as.vector(A)
        if(is(symm, "SphericalSymmetry")) 
            S <- symm@SymmCenter == 0
        else
            S <- FALSE
        z <- getInfCent(L2deriv = L2deriv, neighbor = neighbor,  
                        biastype = biastype, clip = c0, cent = 0, 
                        trafo = trafo, tol.z = tol, symm = S)
        iter <- 0
        repeat{
            iter <- iter + 1
            A.old <- A
            z.old <- z
            A <- getInfStand(L2deriv = L2deriv, neighbor = neighbor, 
                         biastype = biastype,
                         clip = c0, cent = z, trafo = trafo)
            c0 <- b/as.vector(A)
            z <- getInfCent(L2deriv = L2deriv, neighbor = neighbor, 
                         biastype = biastype,
                         clip = c0, cent = z, trafo = trafo, tol.z = tol, symm = S)
            if(max(abs(as.vector(A-A.old)), abs(z-z.old)) < tol) break
            if(iter > maxiter){
                cat("maximum iterations reached!\n", "achieved precision:\t", 
                    max(abs(as.vector(A-A.old)), abs(z-z.old)), "\n")
                break
            }
        }
        info <- paste("optimally robust IC for 'asHampel' with bound =", round(b,3))
        a <- as.vector(A)*z
        Cov <- getInfV(L2deriv = L2deriv, neighbor = neighbor, 
                       biastype = biastype, clip = c0, cent = z, stand = A)

        # getAsRisk(risk = asHampel(), L2deriv = L2deriv, neighbor = neighbor, 
        #          biastype = biastype, clip = b, cent = a, stand = A)$asCov

        r <- neighbor@radius
        Risk <- list(asCov = Cov,
                     asBias = list(value = b, biastype = biastype, 
                                   normtype = normtype, 
                                   neighbortype = class(neighbor)), 
                     trAsCov = list(value = Cov, normtype = normtype),
                     asMSE = list(value = Cov + r^2*b^2, 
                                  r = r,
                                  at = neighbor))

        if(is(neighbor,"ContNeighborhood")){
            w <- new("HampelWeight")
            clip(w) <- b
            cent(w) <- z 
            stand(w) <- A
        }else{
            w <- new("BdStWeight")
            clip(w) <- c(0,b)+a
            stand(w) <- A
        } 
        weight(w) <- getweight(w, neighbor = neighbor, biastype = biastype, 
                               normW = NormType())
        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info, 
                    w = w, biastype = biastype, normtype = NormType()))
    })

setMethod("getInfRobIC", signature(L2deriv = "RealRandVariable", 
                                   risk = "asHampel", 
                                   neighbor = "ContNeighborhood"),
    function(L2deriv, risk, neighbor, Distr, DistrSymm, L2derivSymm,
             L2derivDistrSymm, Finfo, trafo, onesetLM = FALSE,
             z.start, A.start, upper, maxiter, tol, warn, verbose = FALSE){

        biastype <- biastype(risk)
        normtype <- normtype(risk)
        p <- nrow(trafo)

        FI <- solve(trafo%*%solve(Finfo)%*%t(trafo))
        if(is(normtype,"InfoNorm") || is(normtype,"SelfNorm") ) 
           {QuadForm(normtype) <- PosSemDefSymmMatrix(FI); normtype(risk) <- normtype}

        std <- if(is(normtype,"QFNorm")) QuadForm(normtype) else diag(p)

        if(is.null(z.start)) z.start <- numeric(ncol(trafo))
        if(is.null(A.start)) A.start <- trafo

        ClassIC <- trafo %*% solve(Finfo) %*% L2deriv
        lower <- q(Distr)(getdistrOption("TruncQuantile"))
        upper <- q(Distr)(1-getdistrOption("TruncQuantile"))
        x <- seq(from = lower, to = upper, by = 0.01)
        bmax <- evalRandVar(ClassIC, as.matrix(x))^2
        bmax <- sqrt(max(colSums(bmax)))
        b <- risk@bound
        cat("numerical approximation of maximal bound:\t", bmax, "\n")
        if(b >= bmax){
            if(warn) cat("'b >= maximum asymptotic bias' => (classical) optimal IC\n", 
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobIC(L2deriv = L2deriv, risk = asCov(), neighbor = neighbor, 
                                Distr = Distr, Finfo = Finfo, trafo = trafo, 
                                QuadForm = std, verbose = verbose)
            res <- c(res, list(biastype = biastype, normtype = normtype))
            trAsCov <- sum(diag(std%*%res$risk$asCov)); 
            r <- neighbor@radius
            res$risk$trAsCov <- list(value = trAsCov, normtype = normtype)
            res$risk$asBias <- list(value = b, biastype = biastype, 
                                   normtype = normtype, 
                                   neighbortype = class(neighbor))
            res$risk$asMSE <- list(value = trAsCov + r^2*b^2, 
                                   r = r,
                                   at = neighbor)
            return(res)
        }

        res <- getInfRobIC(L2deriv = L2deriv, 
                     risk = asBias(biastype = biastype, normtype = normtype), 
                     neighbor = neighbor, Distr = Distr, DistrSymm = DistrSymm, 
                     L2derivSymm = L2derivSymm, L2derivDistrSymm = L2derivDistrSymm, 
                     z.start = z.start, A.start = A.start, trafo = trafo, 
                     maxiter = maxiter, tol = tol, warn = warn, Finfo = Finfo, 
                     verbose = verbose)
        bmin <- res$b

        cat("minimal bound:\t", bmin, "\n")
        if(b <= bmin){
            if(warn) cat("'b <= minimum asymptotic bias'\n",
                         "=> the minimum asymptotic bias (lower case) solution is returned\n")

            asMSE <- sum(diag(std%*%res$risk$asCov)) + neighbor@radius^2*bmin^2
            if(!is.null(res$risk$asMSE)) res$risk$asMSE <- asMSE 
               else     res$risk <- c(list(asMSE = asMSE), res$risk)

            return(res)
        }

        comp <- .getComp(L2deriv, DistrSymm, L2derivSymm,
             L2derivDistrSymm)

        z.comp <- comp$"z.comp"
        A.comp <- comp$"A.comp"

        w <- new("HampelWeight")
        clip(w) <- b
        z <- z.start
        A <- A.start
        iter <- 0
        repeat{
            iter <- iter + 1
            z.old <- z
            A.old <- A
            cent(w) <- z 
            stand(w) <- A 

            weight(w) <- getweight(w, neighbor = neighbor, biastype = biastype, 
                                   normW = normtype)


            z <- getInfCent(L2deriv = L2deriv, neighbor = neighbor,  
                            biastype = biastype, Distr = Distr, z.comp = z.comp, 
                            w = w)
            A <- getInfStand(L2deriv = L2deriv, neighbor = neighbor, 
                         biastype = biastype, Distr = Distr, A.comp = A.comp, 
                         cent = z, trafo = trafo, w = w)

            normtype.old <- normtype
            if(is(normtype,"SelfNorm")){
               normtype(risk) <- normtype <- updateNorm(normtype = normtype,  
                   L2 = L2deriv, neighbor = neighbor, biastype = biastype,
                   Distr = Distr, V.comp = A.comp, cent = z, stand = A, w = w)
               std <- QuadForm(normtype)
            }

            prec <- max(max(abs(A-A.old)), max(abs(z-z.old)))
            if(verbose)
                cat("current precision in IC algo:\t", prec, "\n")
            if(prec < tol) break
            if(iter > maxiter){
                cat("maximum iterations reached!\n", "achieved precision:\t", prec, "\n")
                break
            }
        }
        if (onesetLM){
            cent(w) <- z 
            stand(w) <- A 

            weight(w) <- getweight(w, neighbor = neighbor, biastype = biastype, 
                                   normW = normtype)
        }
        else normtype <- normtype.old

        info <- paste("optimally robust IC for 'asHampel' with bound =", round(b,3))
        a <- as.vector(A %*% z)
        Cov <- getInfV(L2deriv = L2deriv, neighbor = neighbor, 
                       biastype = biastype, Distr = Distr, 
                       V.comp = A.comp, cent = a, 
                       stand = A, w = w)
        #getAsRisk(risk = asCov(), L2deriv = L2deriv, neighbor = neighbor, 
        #          biastype = biastype, Distr = Distr, clip = b, cent = a, 
        #          stand = A)$asCov
        trAsCov <- sum(diag(std%*%Cov)); r <- neighbor@radius
        Risk <- list(trAsCov = list(value = trAsCov, 
                                    normtype = normtype),
                     asCov = Cov,  
                     asBias = list(value = b, biastype = biastype, 
                                   normtype = normtype, 
                                   neighbortype = class(neighbor)),
                     asMSE = list(value = trAsCov + r^2*b^2, 
                                  r = r,
                                  at = neighbor))

        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info, 
                    w = w, biastype = biastype, normtype = normtype))
    })
