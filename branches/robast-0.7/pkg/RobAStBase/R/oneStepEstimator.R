###############################################################################
## one-step estimator
###############################################################################
setMethod("oneStepEstimator", signature(x = "numeric", 
                                        IC = "InfluenceCurve",
                                        start = "numeric"),
    function(x, IC, start, useLast = getRobAStBaseOption("kStepUseLast")){
        es.call <- match.call()
        es.call[[1]] <- as.name("oneStepEstimator")
        nrvalues <- dimension(IC@Curve)
        if(is.list(start)) start <- unlist(start)
        if(nrvalues != length(start))
            stop("dimension of 'start' != dimension of 'Curve'")

        res <- start + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)

        if(is(IC, "IC")){
            L2Fam <- eval(CallL2Fam(IC))
            Infos <- matrix(c("oneStepEstimator", 
                            paste("1-step estimate for", name(L2Fam))),
                            ncol = 2)
            colnames(Infos) <- c("method", "message")
            if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
            if(useLast && !is(modifyIC(IC), "NULL") ){
                newParam <- param(L2Fam)
                main(newParam)[] <- res
                newL2Fam <- modifyModel(L2Fam, newParam)
                IC <- modifyIC(IC)(newL2Fam, IC)
                Infos <- rbind(Infos, c("oneStepEstimator", 
                                        "computation of IC, asvar and asbias via useLast = TRUE"))
            }else{
                if(useLast && is(modifyIC(IC), "NULL")){
                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
                             is filled with some function!")
                    Infos <- rbind(Infos, c("oneStepEstimator", 
                                            "slot 'modifyIC' of 'IC' was not filled!"))
                }
                Infos <- rbind(Infos, c("oneStepEstimator", 
                                        "computation of IC, asvar and asbias via useLast = FALSE"))
            }
            if("asCov" %in% names(Risks(IC)))
                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
                    asVar <- Risks(IC)$asCov
                else
                    asVar <- Risks(IC)$asCov$value
            else
                asVar <- getRiskIC(IC, risk = asCov())$asCov$value

            if("asBias" %in% names(Risks(IC))){
                if(length(Risks(IC)$asBias) == 1)
                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
                else
                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
            }else{
                if(is(IC, "HampIC")){
                    r <- neighborRadius(IC)
                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
                }else{
                    asBias <- NULL
                }
            }
        }else{
            Infos <- matrix(c("oneStepEstimator", "1-step estimate"), ncol = 2)
            colnames(Infos) <- c("method", "message")
            asVar <- NULL
            asBias <- NULL
        }

        new("kStepEstimate", name = "1-step estimate", estimate = res, 
            estimate.call = es.call, samplesize = length(x), asvar = asVar, 
            asbias = asBias, pIC = IC, steps = 1L, Infos = Infos)
    })
setMethod("oneStepEstimator", signature(x = "matrix", 
                                        IC = "InfluenceCurve",
                                        start = "numeric"),
    function(x, IC, start, useLast = getRobAStBaseOption("kStepUseLast")){
        es.call <- match.call()
        es.call[[1]] <- as.name("oneStepEstimator")
        nrvalues <- dimension(IC@Curve)
        if(is.list(start)) start <- unlist(start)
        if(nrvalues != length(start))
            stop("dimension of 'start' != dimension of 'Curve'")
        if(ncol(x) != IC@Curve[[1]]@Domain@dimension)
            stop("'x' has wrong dimension")

        res <- start + rowMeans(evalIC(IC, x), na.rm = TRUE)

        if(is(IC, "IC")){
            L2Fam <- eval(CallL2Fam(IC))
            Infos <- matrix(c("oneStepEstimator", 
                            paste("1-step estimate for", name(L2Fam))),
                            ncol = 2)
            colnames(Infos) <- c("method", "message")
            if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
            if(useLast && !is(modifyIC(IC), "NULL") ){
                newParam <- param(L2Fam)
                main(newParam)[] <- res
                newL2Fam <- modifyModel(L2Fam, newParam)
                IC <- modifyIC(IC)(newL2Fam, IC)
                Infos <- rbind(Infos, c("oneStepEstimator", 
                                        "computation of IC, asvar and asbias via useLast = TRUE"))
            }else{
                if(useLast && is(modifyIC(IC), "NULL")){
                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
                             is filled with some function!")
                    Infos <- rbind(Infos, c("oneStepEstimator", 
                                            "slot 'modifyIC' of 'IC' was not filled!"))
                }
                Infos <- rbind(Infos, c("oneStepEstimator", 
                                        "computation of IC, asvar and asbias via useLast = FALSE"))
            }
            if("asCov" %in% names(Risks(IC)))
                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
                    asVar <- Risks(IC)$asCov
                else
                    asVar <- Risks(IC)$asCov$value
            else
                asVar <- getRiskIC(IC, risk = asCov())$asCov$value

            if("asBias" %in% names(Risks(IC))){
                if(length(Risks(IC)$asBias) == 1)
                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
                else
                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
            }else{
                if(is(IC, "HampIC")){
                    r <- neighborRadius(IC)
                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
                }else{
                    asBias <- NULL
                }
            }
        }else{
            Infos <- matrix(c("oneStepEstimator", "1-step estimate"), ncol = 2)
            colnames(Infos) <- c("method", "message")
            asVar <- NULL
            asBias <- NULL
        }

        new("kStepEstimate", name = "1-step estimate", estimate = res, 
            estimate.call = es.call, samplesize = ncol(x), asvar = asVar, 
            asbias = asBias, pIC = IC, steps = 1L, Infos = Infos)
    })
setMethod("oneStepEstimator", signature(x = "numeric", 
                                        IC = "InfluenceCurve",
                                        start = "Estimate"),
    function(x, IC, start, useLast = getRobAStBaseOption("kStepUseLast")){
        es.call <- match.call()
        es.call[[1]] <- as.name("oneStepEstimator")
        nrvalues <- dimension(IC@Curve)
        start0 <- estimate(start)
        if(is.list(start0)) start0 <- unlist(start0)
        if(nrvalues != length(start0))
            stop("dimension of slot 'estimate' of 'start' != dimension of 'Curve'")

        res <- start0 + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)

        if(is(IC, "IC")){
            L2Fam <- eval(CallL2Fam(IC))
            Infos <- matrix(c("oneStepEstimator", 
                            paste("1-step estimate for", name(L2Fam))),
                            ncol = 2)
            colnames(Infos) <- c("method", "message")
            if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
            if(useLast && !is(modifyIC(IC), "NULL") ){
                newParam <- param(L2Fam)
                main(newParam)[] <- res
                newL2Fam <- modifyModel(L2Fam, newParam)
                IC <- modifyIC(IC)(newL2Fam, IC)
                Infos <- rbind(Infos, c("oneStepEstimator", 
                                        "computation of IC, asvar and asbias via useLast = TRUE"))
            }else{
                if(useLast && is(modifyIC(IC), "NULL")){
                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC' 
                             is filled with some function!")
                    Infos <- rbind(Infos, c("oneStepEstimator", 
                                            "slot 'modifyIC' of 'IC' was not filled!"))
                }
                Infos <- rbind(Infos, c("oneStepEstimator", 
                                        "computation of IC, asvar and asbias via useLast = FALSE"))
            }
            if("asCov" %in% names(Risks(IC)))
                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
                    asVar <- Risks(IC)$asCov
                else
                    asVar <- Risks(IC)$asCov$value
            else
                asVar <- getRiskIC(IC, risk = asCov())$asCov$value

            if("asBias" %in% names(Risks(IC))){
                if(length(Risks(IC)$asBias) == 1)
                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
                else
                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
            }else{
                if(is(IC, "HampIC")){
                    r <- neighborRadius(IC)
                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
                }else{
                    asBias <- NULL
                }
            }
        }else{
            Infos <- matrix(c("oneStepEstimator", "1-step estimate"), ncol = 2)
            colnames(Infos) <- c("method", "message")
            asVar <- NULL
            asBias <- NULL
        }

        new("kStepEstimate", name = "1-step estimate", estimate = res, 
            estimate.call = es.call, samplesize = length(x), asvar = asVar, 
            asbias = asBias, pIC = IC, steps = 1L, Infos = Infos)
    })
setMethod("oneStepEstimator", signature(x = "matrix", 
                                        IC = "InfluenceCurve",
                                        start = "Estimate"),
    function(x, IC, start, useLast = getRobAStBaseOption("kStepUseLast")){
        es.call <- match.call()
        es.call[[1]] <- as.name("oneStepEstimator")
        nrvalues <- dimension(IC@Curve)
        start0 <- estimate(start)
        if(is.list(start0)) start0 <- unlist(start0)
        if(nrvalues != length(start0))
            stop("dimension of slot 'estimate' of 'start' != dimension of 'Curve'")
        if(ncol(x) != IC@Curve[[1]]@Domain@dimension)
            stop("'x' has wrong dimension")

        res <- start0 + rowMeans(evalIC(IC, x), na.rm = TRUE)

        if(is(IC, "IC")){
            L2Fam <- eval(CallL2Fam(IC))
            Infos <- matrix(c("oneStepEstimator", 
                            paste("1-step estimate for", name(L2Fam))),
                            ncol = 2)
            colnames(Infos) <- c("method", "message")
            if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
            if(useLast && !is(modifyIC(IC), "NULL") ){
                newParam <- param(L2Fam)
                main(newParam)[] <- res
                newL2Fam <- modifyModel(L2Fam, newParam)
                IC <- modifyIC(IC)(newL2Fam, IC)
                Infos <- rbind(Infos, c("oneStepEstimator", 
                                        "computation of IC, asvar and asbias via useLast = TRUE"))
            }else{
                if(useLast && is(modifyIC(IC), "NULL")){
                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
                             is filled with some function!")
                    Infos <- rbind(Infos, c("oneStepEstimator", 
                                            "slot 'modifyIC' of 'IC' was not filled!"))
                }
                Infos <- rbind(Infos, c("oneStepEstimator", 
                                        "computation of IC, asvar and asbias via useLast = FALSE"))
            }
            if("asCov" %in% names(Risks(IC)))
                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
                    asVar <- Risks(IC)$asCov
                else
                    asVar <- Risks(IC)$asCov$value
            else
                asVar <- getRiskIC(IC, risk = asCov())$asCov$value

            if("asBias" %in% names(Risks(IC))){
                if(length(Risks(IC)$asBias) == 1)
                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
                else
                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
            }else{
                if(is(IC, "HampIC")){
                    r <- neighborRadius(IC)
                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
                }else{
                    asBias <- NULL
                }
            }
        }else{
            Infos <- matrix(c("oneStepEstimator", "1-step estimate"), ncol = 2)
            colnames(Infos) <- c("method", "message")
            asVar <- NULL
            asBias <- NULL
        }

        new("kStepEstimate", name = "1-step estimate", estimate = res, 
            estimate.call = es.call, samplesize = ncol(x), asvar = asVar, 
            asbias = asBias, pIC = IC, steps = 1L, Infos = Infos)
    })
