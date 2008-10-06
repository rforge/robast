###############################################################################
## get optimally robust IC for convex asymptotic risks
###############################################################################
setMethod("getInfRobIC", signature(L2deriv = "UnivariateDistribution", 
                                   risk = "asGRisk", 
                                   neighbor = "UncondNeighborhood"),
    function(L2deriv, risk, neighbor, symm, Finfo, trafo, upper, maxiter, tol, 
             warn, noLow = FALSE, verbose = FALSE){
        biastype <- biastype(risk)
        radius <- neighbor@radius
        if(identical(all.equal(radius, 0), TRUE)){
            if(warn) cat("'radius == 0' => (classical) optimal IC\n", 
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobIC(L2deriv = L2deriv, risk = asCov(), 
                        neighbor = neighbor, Finfo = Finfo, trafo = trafo,
                        verbose = verbose)
            res <- c(res, list(biastype = biastype, normtype = NormType()))
            Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                              biastype = biastype, clip = res$b, cent = res$a, 
                              stand = res$A, trafo = trafo)
            res$risk <- c(Risk, res$risk)
            Cov <- res$risk$asCov
            res$risk$asBias <- list(value = b, biastype = biastype, 
                                   normtype = NormType(), 
                                   neighbortype = class(neighbor))
            res$risk$asMSE <- list(value = Cov + radius^2*b^2, 
                                   r = radius,
                                   at = neighbor)
            return(res)
        }
        z <- 0
        c0 <- 0
        iter <- 0
        if(is(symm, "SphericalSymmetry")) 
            S <- symm@SymmCenter == 0
        else
            S <- FALSE

        prec <- 1
        repeat{
            iter <- iter + 1
            z.old <- z
            c0.old <- c0
            ## new
            lower0 <- getL1normL2deriv(L2deriv = L2deriv, cent = z) / 
                                      (1 + neighbor@radius^2)
            upper0 <- sqrt( as.numeric( Finfo + z^2 )/(( 1 + neighbor@radius^2)^2 - 1) )
            if (!is.null(upper)|(iter == 1)) 
                    {lower <- .Machine$double.eps^0.75
                }else{ lower <- lower0; upper <- upper0}
            ##
            c0 <- try(uniroot(getInfClip, 
                  ## new
                        lower = lower, upper = upper,
                  ##
                        tol = tol, L2deriv = L2deriv, risk = risk, 
                        neighbor = neighbor,  biastype = biastype,
                        cent = z, symm = S, 
                        trafo = trafo)$root, silent = TRUE)
            if(!is.numeric(c0)){
                if(warn) cat("The IC algorithm did not converge!\n", 
                             "'radius >= maximum radius' for the given risk?\n",
                             "=> the minimum asymptotic bias (lower case) solution is returned\n")
                res <- getInfRobIC(L2deriv = L2deriv, risk = asBias(biastype = biastype(risk),
                                                              normtype = normtype(risk)), 
                                neighbor = neighbor, Finfo = Finfo, 
                                symm = symm, trafo = trafo, upper = upper, 
                                maxiter = maxiter, tol = tol, warn = warn,
                                verbose = verbose)
                Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                                  biastype = biastype,  clip = res$b, cent = res$a, 
                                  stand = res$A, trafo = trafo)
                res$risk <- c(Risk, res$risk)
                return(res)
            }
            z <- getInfCent(L2deriv = L2deriv, neighbor = neighbor,  biastype = biastype,
                            clip = c0, cent = z, symm = S, trafo = trafo, tol.z = tol)
#            cat("c0:\t", c0, "c0.old:\t", c0.old, "z:\t", z, "z.old:\t", z.old, "\n")
            if(S) break

            prec.old <- prec
            prec <- max(abs(z - z.old), abs(c0-c0.old))
            if(verbose)
                cat("current precision in IC algo:\t", prec, "\n")
            if(prec < tol) break
            if(abs(prec.old - prec) < 1e-10){
                cat("algorithm did not converge!\n", "achieved precision:\t", prec, "\n")
                break
            }
            if(iter > maxiter){
                cat("maximum iterations reached!\n", "achieved precision:\t", prec, "\n")
                break
            }
        }
        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
        A <- getInfStand(L2deriv = L2deriv, neighbor = neighbor,
                     biastype = biastype, clip = c0, cent = z, trafo = trafo)
        a <- as.vector(A)*z
        b <- abs(as.vector(A))*c0
        Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                          biastype = biastype, clip = b, cent = a, stand = A, 
                          trafo = trafo)
        Cov <- getInfV(L2deriv = L2deriv, neighbor = neighbor, 
                       biastype = biastype, clip = c0, cent = z, stand = A)

        Risk <- c(Risk, list(asCov = Cov,  
                     asBias = list(value = b, biastype = biastype, 
                                   normtype = normtype(risk), 
                                   neighbortype = class(neighbor)), 
                     trAsCov = list(value = Cov, normtype = normtype(risk)),
                     asMSE = list(value = Cov + radius^2*b^2, 
                                  r = radius,
                                  at = neighbor)))

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

        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info, w = w,
                    biastype = biastype, normtype = normtype(risk)))
    })



################################################################################


setMethod("getInfRobIC", signature(L2deriv = "RealRandVariable", 
                                   risk = "asGRisk", 
                                   neighbor = "ContNeighborhood"),
    function(L2deriv, risk, neighbor, Distr, DistrSymm, L2derivSymm, 
             L2derivDistrSymm, Finfo, trafo, onesetLM = FALSE, 
             z.start, A.start, upper, maxiter, tol, warn, verbose = FALSE){
        biastype <- biastype(risk)
        normtype <- normtype(risk)

        FI <- solve(trafo%*%solve(Finfo)%*%t(trafo))
        if(is(normtype,"InfoNorm") || is(normtype,"SelfNorm") ) 
           {QuadForm(normtype) <- PosSemDefSymmMatrix(FI); 
            normtype(risk) <- normtype}
        QF <- if(is(normtype,"QFNorm")) QuadForm(normtype) else diag(nrow(trafo))

        if(is.null(z.start)) z.start <- numeric(ncol(trafo))
        if(is.null(A.start)) A.start <- trafo %*% solve(Finfo)

        radius <- neighbor@radius
        if(identical(all.equal(radius, 0), TRUE)){
            if(warn) cat("'radius == 0' => (classical) optimal IC\n", 
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobIC(L2deriv = L2deriv, risk = asCov(), neighbor = neighbor, 
                               Distr = Distr, Finfo = Finfo, trafo = trafo, 
                               QuadForm = QF, verbose = verbose)
            res <- c(res, list(biastype = biastype, normtype = normtype))
            Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                              biastype = biastype, cent = res$a, 
                              stand = res$A, trafo = trafo)
            res$risk <- c(Risk, res$risk)
            trAsCov <- sum(diag(QF%*%res$risk$asCov)); 
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

        comp <- .getComp(L2deriv, DistrSymm, L2derivSymm, L2derivDistrSymm)

        z.comp <- comp$"z.comp"
        A.comp <- comp$"A.comp"

        w <- new("HampelWeight")
        z <- z.start
        A <- A.start
        b <- 0
        iter <- 0
        prec <- 1
        repeat{
            iter <- iter + 1
            z.old <- z
            b.old <- b
            A.old <- A
            ##
            cent(w) <- z 
            stand(w) <- A 

            ## new
            lower0 <- getL1normL2deriv(L2deriv = L2deriv, cent = z, stand = A, 
                                       Distr = Distr, normtype = normtype)/(1+neighbor@radius^2)

            QF <- if(is(normtype,"QFNorm")) QuadForm(normtype) else diag(nrow(A))
            upper0 <- sqrt( (sum( diag(QF%*%A%*%Finfo%*%t(A))) + t(A%*%z)%*%QF%*%(A%*%z)) / 
                          ((1 + neighbor@radius^2)^2-1))

            if (!is.null(upper)|(iter == 1)) 
                    {lower <- .Machine$double.eps^0.75; 
                     if(is.null(upper)) upper <- 10*upper0
                }else{ lower <- lower0; upper <- upper0}

            ##
            b <- try(uniroot(getInfClip, 
                  ## new
                         lower = lower, upper = upper,
                  ##
                         tol = tol, L2deriv = L2deriv, risk = risk, 
                         biastype = biastype, Distr = Distr, neighbor = neighbor, 
                         stand = A, cent = z, trafo = trafo)$root, silent = TRUE)
            if(!is.numeric(b)){
                if(warn) cat("Could not determine optimal clipping bound!\n", 
                             "'radius >= maximum radius' for the given risk?\n",
                             "=> the minimum asymptotic bias (lower case) solution is returned\n",
                             "If 'no' => Try again with modified starting values ",
                             "'z.start' and 'A.start'\n")
                             res <- getInfRobIC(L2deriv = L2deriv, 
                                        risk =  asBias(biastype = biastype(risk),
                                                       normtype = normtype(risk)), 
                                    neighbor = neighbor, Distr = Distr, DistrSymm = DistrSymm,
                                    L2derivSymm = L2derivSymm, L2derivDistrSymm = L2derivDistrSymm,
                                    z.start = z.start, A.start = A.start, trafo = trafo, 
                                    maxiter = maxiter, tol = tol, warn = warn, Finfo = Finfo,
                                    verbose = verbose)
                             normtype(risk) <- res$normtype
                             Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                                  biastype = biastype, clip = NULL,   
                                  cent = res$a, stand = res$A, trafo = trafo)
                res$risk <- c(Risk, res$risk)
                return(res)
            }
            clip(w) <- b


            weight(w) <- getweight(w, neighbor = neighbor, biastype = biastype, 
                                   normW = normtype)

            z <- getInfCent(L2deriv = L2deriv, neighbor = neighbor,  
                            biastype = biastype, Distr = Distr, z.comp = z.comp, 
                            w = w)
            A <- getInfStand(L2deriv = L2deriv, neighbor = neighbor, 
                         biastype = biastype, Distr = Distr, A.comp = A.comp, 
                         cent = z, trafo = trafo, w = w)

            normtype.old <- normtype

            if (is(normtype,"SelfNorm"))
                {normtype(risk) <- normtype <- updateNorm(normtype = normtype, 
                   L2 = L2deriv, neighbor = neighbor, biastype = biastype,
                   Distr = Distr, V.comp = A.comp, cent = as.vector(A %*% z), 
                   stand = A, w = w)}

            prec.old <- prec
            prec <- max(abs(b-b.old), max(abs(A-A.old)), max(abs(z-z.old)))
            if(verbose)
                cat("current precision in IC algo:\t", prec, "\n")
            if(prec < tol) break
            if(abs(prec.old - prec) < 1e-10){
                cat("algorithm did not converge!\n", "achieved precision:\t", prec, "\n")
                break
            }
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
        }else normtype <- normtype.old

        a <- as.vector(A %*% z)
        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
        Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                          biastype = biastype, clip = b, cent = a, stand = A, 
                          trafo = trafo)
        Cov <- getInfV(L2deriv = L2deriv, neighbor = neighbor, 
                       biastype = biastype, Distr = Distr, 
                       V.comp = A.comp, cent = a, 
                       stand = A, w = w)

        QF <- if(is(normtype,"QFNorm")) QuadForm(normtype) else diag(nrow(A))

        trAsCov <- sum(diag(QF%*%Cov))
        Risk <- c(Risk, list(asCov = Cov,  
                     asBias = list(value = b, biastype = biastype, 
                                   normtype = normtype, 
                                   neighbortype = class(neighbor)), 
                     trAsCov = list(value = trAsCov, 
                                   normtype = normtype),
                     asMSE = list(value = trAsCov + radius^2*b^2, 
                                  r = radius,
                                  at = neighbor)))

        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info, w = w, 
                    biastype = biastype, normtype = normtype))
    })
