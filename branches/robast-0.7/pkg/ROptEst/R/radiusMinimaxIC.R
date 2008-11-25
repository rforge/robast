###############################################################################
## radius minimax optimally robust IC 
## for L2ParamFamily and asymptotic risks
###############################################################################
setMethod("radiusMinimaxIC", signature(L2Fam = "L2ParamFamily", 
                                       neighbor = "UncondNeighborhood",
                                       risk = "asGRisk"),
    function(L2Fam, neighbor, risk, loRad, upRad, z.start = NULL, 
             A.start = NULL, upper = 1e5, maxiter = 50, 
             tol = .Machine$double.eps^0.4, warn = FALSE, verbose = FALSE){
        ow <- options("warn")
        on.exit(options(ow))
        if(length(loRad) != 1)
            stop("'loRad' is not of length == 1")
        if(length(upRad) != 1)
            stop("'upRad' is not of length == 1")
        if(loRad >= upRad)
            stop("'upRad < loRad' is not fulfilled")
        biastype <- biastype(risk)
        L2derivDim <- numberOfMaps(L2Fam@L2deriv)

        if(is(normtype(risk),"SelfNorm")||is(normtype(risk),"InfoNorm"))
           upRad <- min(upRad,10) 

        if(L2derivDim == 1){
            options(warn = -1)
            upper.b <- upper
            lower <- ifelse(identical(all.equal(loRad, 0), TRUE), 1e-4, loRad)
            upper <- ifelse(upRad == Inf, max(loRad+1, 2), upRad)

            if(identical(all.equal(loRad, 0), TRUE)){
                loRad <- 0
                loRisk <- 1/as.vector(L2Fam@FisherInfo)
            }else{
                neighbor@radius <- loRad
                resLo <- getInfRobIC(L2deriv = L2Fam@L2derivDistr[[1]], neighbor = neighbor, 
                            risk = risk, symm = L2Fam@L2derivDistrSymm[[1]],
                            Finfo = L2Fam@FisherInfo, upper = upper.b,
                            trafo = L2Fam@param@trafo, maxiter = maxiter, tol = tol, 
                            warn = warn, verbose = verbose)
                loRisk <- getAsRisk(risk = risk, L2deriv = L2Fam@L2derivDistr[[1]], 
                                    neighbor = neighbor, biastype = biastype,
                                    clip = resLo$b, cent = resLo$a, 
                                    stand = resLo$A, trafo = L2Fam@param@trafo)[[1]]
            }

            if(upRad == Inf){
                bmin <- getAsRisk(risk = asBias(biastype = biastype), 
                                  L2deriv = L2Fam@L2derivDistr[[1]], 
                                  neighbor = neighbor, biastype = biastype, 
                                  trafo = L2Fam@param@trafo)$asBias
                upRisk <- bmin^2
            }else{
                neighbor@radius <- upRad
                resUp <- getInfRobIC(L2deriv = L2Fam@L2derivDistr[[1]], neighbor = neighbor, 
                            risk = risk, symm = L2Fam@L2derivDistrSymm[[1]],
                            Finfo = L2Fam@FisherInfo, upper = upper.b,
                            trafo = L2Fam@param@trafo, maxiter = maxiter, tol = tol, 
                            warn = warn, verbose = verbose)
                upRisk <- getAsRisk(risk = risk, L2deriv = L2Fam@L2derivDistr[[1]], 
                                    neighbor = neighbor, biastype = biastype, 
                                    clip = resUp$b, cent = resUp$a, 
                                    stand = resUp$A, trafo = L2Fam@param@trafo)[[1]]
            }

            loNorm<- upNorm <- NormType()
            leastFavR <- uniroot(getIneffDiff, lower = lower, upper = upper, 
                            tol = .Machine$double.eps^0.25, L2Fam = L2Fam, neighbor = neighbor, 
                            upper.b = upper.b, risk = risk, loRad = loRad, upRad = upRad, 
                            loRisk = loRisk, upRisk = upRisk, eps = tol, 
                            MaxIter = maxiter, warn = warn, 
                            loNorm = loNorm, upNorm = upNorm)$root
            neighbor@radius <- leastFavR
            res <- getInfRobIC(L2deriv = L2Fam@L2derivDistr[[1]], neighbor = neighbor, 
                        risk = risk, symm = L2Fam@L2derivSymm[[1]],
                        Finfo = L2Fam@FisherInfo, upper = upper.b,
                        trafo = L2Fam@param@trafo, maxiter = maxiter, tol = tol, 
                        warn = warn, verbose = verbose)
            options(ow)
            res$info <- c("radiusMinimaxIC", paste("radius minimax IC for radius interval [", 
                            round(loRad, 3), ", ", round(upRad, 3), "]", sep=""))
            res$info <- rbind(res$info, c("radiusMinimaxIC", 
                            paste("least favorable radius: ", round(leastFavR, 3), sep="")))
            res$info <- rbind(res$info, c("radiusMinimaxIC", 
                            paste("maximum ", sQuote(class(risk)[1]), "-inefficiency: ",
                            round(ineff, 3), sep="")))
            res <- c(res, modifyIC = getModifyIC(L2FamIC = L2Fam, 
                                                 neighbor = neighbor, 
                                                 risk = risk))
            return(generateIC(neighbor, L2Fam, res))
        }else{
            if(is(L2Fam@distribution, "UnivariateDistribution")){
                if((length(L2Fam@L2deriv) == 1) & is(L2Fam@L2deriv[[1]], "RealRandVariable")){
                    L2deriv <- L2Fam@L2deriv[[1]]
                    L2derivSymm <- L2Fam@L2derivSymm
                    L2derivDistrSymm <- L2Fam@L2derivDistrSymm
                }else{
                    L2deriv <- diag(dimension(L2Fam@L2deriv)) %*% L2Fam@L2deriv
                    L2deriv <- RealRandVariable(Map = L2deriv@Map, Domain = L2deriv@Domain)
                    nrvalues <- numberOfMaps(L2deriv)
                    if(numberOfMaps(L2Fam@L2deriv) != nrvalues){
                        L1 <- vector("list", nrvalues)
                        L2 <- vector("list", nrvalues)
                        for(i in 1:nrvalues){
                            L1[[i]] <- NonSymmetric()
                            L2[[i]] <- NoSymmetry()
                        }
                        L2derivSymm <- new("FunSymmList", L1)
                        L2derivDistrSymm <- new("DistrSymmList", L2)
                    }
                }
                normtype <- normtype(risk)

                Finfo <- L2Fam@FisherInfo
                trafo <- L2Fam@param@trafo

                p <- nrow(trafo)
                FI0 <- trafo%*%solve(Finfo)%*%t(trafo)
                FI <- solve(FI0)

                if(is(normtype,"InfoNorm") || is(normtype,"SelfNorm") ) 
                     {QuadForm(normtype) <- PosSemDefSymmMatrix(FI); 
                      normtype(risk) <- normtype}
                std <- if(is(normtype,"QFNorm")) QuadForm(normtype) else diag(p)

                options(warn = -1)
                upper.b <- upper
                lower <- ifelse(identical(all.equal(loRad, 0), TRUE), 1e-4, loRad)
                upper <- ifelse(upRad == Inf, max(loRad+1, 2), upRad)

                if(identical(all.equal(loRad, 0), TRUE)){
                    loRad <- 0
                    loRisk <- sum(diag(std%*%FI0))
                    loNorm <- normtype
                }else{
                    neighbor@radius <- loRad
                    resLo <- getInfRobIC(L2deriv = L2deriv, neighbor = neighbor, risk = risk, 
                                Distr = L2Fam@distribution, DistrSymm = L2Fam@distrSymm, 
                                L2derivSymm = L2derivSymm, L2derivDistrSymm = L2derivDistrSymm, 
                                Finfo = L2Fam@FisherInfo, trafo = trafo, z.start = z.start, 
                                A.start = A.start, upper = upper.b, maxiter = maxiter, 
                                tol = tol, warn = warn, verbose = verbose)
                    riskLo <- risk
                    normtype(riskLo) <- resLo$normtype
                    loRisk <- getAsRisk(risk = riskLo, L2deriv = L2deriv, 
                                        neighbor = neighbor, biastype = biastype, 
                                        clip = resLo$b, cent = resLo$a, 
                                        stand = resLo$A, trafo = trafo)[[1]]
                    loNorm <- resLo$normtype
                }

                if(upRad == Inf){
                    biasR <- getAsRisk(risk = asBias(biastype = biastype(risk), 
                                      normtype = normtype), L2deriv = L2deriv, 
                                      neighbor = neighbor, biastype = biastype, 
                                      Distr = L2Fam@distribution, 
                                      DistrSymm = L2Fam@distrSymm, 
                                      L2derivSymm = L2derivSymm, 
                                      L2derivDistrSymm= L2derivDistrSymm,
                                trafo = trafo, z.start = z.start, 
                                A.start = A.start, 
                                maxiter = maxiter, tol = tol,
                                warn = warn)
                    bmin <- biasR$asBias
                    upNorm <- biasR$normtype
                    upRisk <- bmin^2
                }else{
                    neighbor@radius <- upRad
                    resUp <- getInfRobIC(L2deriv = L2deriv, neighbor = neighbor, risk = risk, 
                                Distr = L2Fam@distribution, DistrSymm = L2Fam@distrSymm, 
                                L2derivSymm = L2derivSymm, L2derivDistrSymm = L2derivDistrSymm, 
                                Finfo = L2Fam@FisherInfo, trafo = trafo, z.start = z.start, 
                                A.start = A.start, upper = upper.b, maxiter = maxiter, 
                                tol = tol, warn = warn, verbose = verbose)
                    riskUp <- risk
                    normtype(riskUp) <- resUp$normtype
                    upRisk <- getAsRisk(risk = riskUp, L2deriv = L2deriv, neighbor = neighbor, 
                                biastype = biastype, clip = resUp$b, cent = resUp$a, stand = resUp$A, trafo = trafo)[[1]]
                    upNorm <- resUp$normtype
                }
                leastFavR <- uniroot(getIneffDiff, lower = lower, upper = upper, 
                                tol = .Machine$double.eps^0.25, L2Fam = L2Fam, neighbor = neighbor, 
                                z.start = z.start, A.start = A.start, upper.b = upper.b, risk = risk, 
                                loRad = loRad, upRad = upRad, loRisk = loRisk, upRisk = upRisk, 
                                eps = tol, MaxIter = maxiter, warn = warn, 
                                loNorm = loNorm, upNorm = upNorm, verbose = verbose)$root
                neighbor@radius <- leastFavR
                res <- getInfRobIC(L2deriv = L2deriv, neighbor = neighbor, risk = risk, 
                            Distr = L2Fam@distribution, DistrSymm = L2Fam@distrSymm, 
                            L2derivSymm = L2derivSymm, L2derivDistrSymm = L2derivDistrSymm, 
                            Finfo = L2Fam@FisherInfo, trafo = trafo, z.start = z.start, 
                            A.start = A.start, upper = upper.b, maxiter = maxiter, 
                            tol = tol, warn = warn, verbose = verbose)
                res$info <- c("radiusMinimaxIC", paste("radius minimax IC for radius interval [", 
                                round(loRad, 3), ", ", round(upRad, 3), "]", sep=""))
                res$info <- rbind(res$info, c("radiusMinimaxIC", 
                                paste("least favorable radius: ", round(leastFavR, 3), sep="")))
                res$info <- rbind(res$info, c("radiusMinimaxIC", 
                                paste("maximum ", sQuote(class(risk)[1]), "-inefficiency: ",
                            round(ineff, 3), sep="")))
                res <- c(res, modifyIC = getModifyIC(L2FamIC = L2Fam, 
                                                 neighbor = neighbor, 
                                                 risk = risk))
                return(generateIC(neighbor, L2Fam, res))
            }else{
                stop("not yet implemented")
            }
        }
    })
