###############################################################################
## one-step estimator
###############################################################################
setMethod("kStepEstimator", signature(x = "numeric", 
                                      IC = "IC",
                                      start = "numeric"),
    function(x, IC, start, steps = 1L, useLast = getRobAStBaseOption("kStepUseLast")){
        es.call <- match.call()
        es.call[[1]] <- as.name("kStepEstimator")
        if(!is.integer(steps))
          steps <- as.integer(steps)
        if(steps < 1)
            stop("steps needs to be a positive integer")

        nrvalues <- dimension(IC@Curve)
        if(is.list(start)) start <- unlist(start)
        if(nrvalues != length(start))
            stop("dimension of 'start' != dimension of 'Curve'")

        res <- start + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)

        L2Fam <- eval(CallL2Fam(IC))
        Infos <- matrix(c("kStepEstimator", 
                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
                        ncol = 2)
        colnames(Infos) <- c("method", "message")
        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE

        if(steps == 1){
            if(useLast && !is(modifyIC(IC), "NULL") ){
                newParam <- param(L2Fam)
                main(newParam) <- res
                newL2Fam <- modifyModel(L2Fam, newParam)
                IC <- modifyIC(IC)(newL2Fam, IC)
                Infos <- rbind(Infos, c("kStepEstimator", 
                                        "computation of IC, asvar and asbias via useLast = TRUE"))
            }else{
                if(useLast && is(modifyIC(IC), "NULL")){
                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
                             is filled with some function!")
                    Infos <- rbind(Infos, c("kStepEstimator", 
                                            "slot 'modifyIC' of 'IC' was not filled!"))
                }
                Infos <- rbind(Infos, c("kStepEstimator", 
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
            return(new("kStepEstimate", estimate.call = es.call,
                       name = paste(steps, "-step estimate", sep = ""),
                       estimate = res, samplesize = length(x), asvar = asVar,
                       asbias = asBias, pIC = IC, steps = 1L, Infos = Infos))
        }else{
            if(is(modifyIC(IC), "NULL"))
                stop("slot 'modifyIC' of 'IC' is 'NULL'!")
            for(i in 2:steps){
                start <- res
                newL2Fam <- eval(CallL2Fam(IC))
                newParam <- param(newL2Fam)
                main(newParam) <- start
                newL2Fam <- modifyModel(newL2Fam, newParam)
                IC <- modifyIC(IC)(newL2Fam, IC)
                res <- start + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)
            }
            if(useLast){
                newL2Fam <- eval(CallL2Fam(IC))
                newParam <- param(newL2Fam)
                main(newParam) <- res
                newL2Fam <- modifyModel(newL2Fam, newParam)
                IC <- modifyIC(IC)(newL2Fam, IC)
                Infos <- rbind(Infos, c("kStepEstimator", 
                                        "computation of IC, asvar and asbias via useLast = TRUE"))
            }else{
                Infos <- rbind(Infos, c("kStepEstimator", 
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
            return(new("kStepEstimate", estimate.call = es.call,
                       name = paste(steps, "-step estimate", sep = ""),
                       estimate = res, samplesize = length(x), asvar = asVar,
                       asbias = asBias, pIC = IC, steps = steps, Infos = Infos))
        }
    })

setMethod("kStepEstimator", signature(x = "matrix", 
                                      IC = "IC",
                                      start = "numeric"),
    function(x, IC, start, steps = 1, useLast = getRobAStBaseOption("kStepUseLast")){
        es.call <- match.call()
        es.call[[1]] <- as.name("kStepEstimator")
        if(!is.integer(steps))
          steps <- as.integer(steps)
        if(steps < 1)
            stop("steps needs to be a positive integer")

        nrvalues <- dimension(IC@Curve)
        if(is.list(start)) start <- unlist(start)
        if(nrvalues != length(start))
            stop("dimension of 'start' != dimension of 'Curve'")
        if(ncol(x) != IC@Curve[[1]]@Domain@dimension)
            stop("'x' has wrong dimension")

        res <- start + rowMeans(evalIC(IC, x), na.rm = TRUE)

        L2Fam <- eval(CallL2Fam(IC))
        Infos <- matrix(c("kStepEstimator", 
                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
                        ncol = 2)
        colnames(Infos) <- c("method", "message")
        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE

        if(steps == 1){
            if(useLast && !is(modifyIC(IC), "NULL") ){
                newParam <- param(L2Fam)
                main(newParam) <- res
                newL2Fam <- modifyModel(L2Fam, newParam)
                IC <- modifyIC(IC)(newL2Fam, IC)
                Infos <- rbind(Infos, c("kStepEstimator", 
                                        "computation of IC, asvar and asbias via useLast = TRUE"))
            }else{
                if(useLast && is(modifyIC(IC), "NULL")){
                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
                             is filled with some function!")
                    Infos <- rbind(Infos, c("kStepEstimator", 
                                            "slot 'modifyIC' of 'IC' was not filled!"))
                }
                Infos <- rbind(Infos, c("kStepEstimator", 
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
            return(new("kStepEstimate", estimate.call = es.call,
                       name = paste(steps, "-step estimate", sep = ""),
                       estimate = res, samplesize = ncol(x), asvar = asVar,
                       asbias = asBias, pIC = IC, steps = 1L, Infos = Infos))
        }else{
            if(is(modifyIC(IC), "NULL"))
                stop("slot 'modifyIC' of 'IC' is 'NULL'!")
            for(i in 2:steps){
                start <- res
                newL2Fam <- eval(CallL2Fam(IC))
                newParam <- param(newL2Fam)
                main(newParam) <- start
                newL2Fam <- modifyModel(newL2Fam, newParam)
                IC <- modifyIC(IC)(newL2Fam, IC)
                res <- start + rowMeans(evalIC(IC, x), na.rm = TRUE)
            }
            if(useLast){
                newL2Fam <- eval(CallL2Fam(IC))
                newParam <- param(newL2Fam)
                main(newParam) <- res
                newL2Fam <- modifyModel(newL2Fam, newParam)
                IC <- modifyIC(IC)(newL2Fam, IC)
                Infos <- rbind(Infos, c("kStepEstimator", 
                                        "computation of IC, asvar and asbias via useLast = TRUE"))
            }else{
                Infos <- rbind(Infos, c("kStepEstimator", 
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
            return(new("kStepEstimate", estimate.call = es.call,
                       name = paste(steps, "-step estimate", sep = ""),
                       estimate = res, samplesize = ncol(x), asvar = asVar,
                       asbias = asBias, pIC = IC, steps = steps, Infos = Infos))
        }
    })
setMethod("kStepEstimator", signature(x = "numeric", 
                                      IC = "IC",
                                      start = "Estimate"),
    function(x, IC, start, steps = 1, useLast = getRobAStBaseOption("kStepUseLast")){
        es.call <- match.call()
        es.call[[1]] <- as.name("kStepEstimator")
        if(!is.integer(steps))
          steps <- as.integer(steps)
        if(steps < 1)
            stop("steps needs to be a positive integer")

        nrvalues <- dimension(IC@Curve)
        start0 <- estimate(start)
        if(is.list(start0)) start0 <- unlist(start0)
        if(nrvalues != length(start0))
            stop("dimension of slot 'estimate' of 'start' != dimension of 'Curve'")

        res <- start0 + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)

        L2Fam <- eval(CallL2Fam(IC))
        Infos <- matrix(c("kStepEstimator", 
                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
                        ncol = 2)
        colnames(Infos) <- c("method", "message")
        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE

        if(steps == 1){
            if(useLast && !is(modifyIC(IC), "NULL") ){
                newParam <- param(L2Fam)
                main(newParam) <- res
                newL2Fam <- modifyModel(L2Fam, newParam)
                IC <- modifyIC(IC)(newL2Fam, IC)
                Infos <- rbind(Infos, c("kStepEstimator", 
                                        "computation of IC, asvar and asbias via useLast = TRUE"))
            }else{
                if(useLast && is(modifyIC(IC), "NULL")){
                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
                             is filled with some function!")
                    Infos <- rbind(Infos, c("kStepEstimator", 
                                            "slot 'modifyIC' of 'IC' was not filled!"))
                }
                Infos <- rbind(Infos, c("kStepEstimator", 
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
            return(new("kStepEstimate", estimate.call = es.call,
                       name = paste(steps, "-step estimate", sep = ""),
                       estimate = res, samplesize = length(x), asvar = asVar,
                       asbias = asBias, pIC = IC, steps = 1L, Infos = Infos))
        }else{
            if(is(modifyIC(IC), "NULL"))
                stop("slot 'modifyIC' of 'IC' is 'NULL'!")
            for(i in 2:steps){
                start0 <- res
                newL2Fam <- eval(CallL2Fam(IC))
                newParam <- param(newL2Fam)
                main(newParam) <- start0
                newL2Fam <- modifyModel(newL2Fam, newParam)
                IC <- modifyIC(IC)(newL2Fam, IC)
                res <- start0 + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)
            }
            if(useLast){
                newL2Fam <- eval(CallL2Fam(IC))
                newParam <- param(newL2Fam)
                main(newParam) <- res
                newL2Fam <- modifyModel(newL2Fam, newParam)
                IC <- modifyIC(IC)(newL2Fam, IC)
                Infos <- rbind(Infos, c("kStepEstimator", 
                                        "computation of IC, asvar and asbias via useLast = TRUE"))
            }else{
                Infos <- rbind(Infos, c("kStepEstimator", 
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
            return(new("kStepEstimate", estimate.call = es.call,
                       name = paste(steps, "-step estimate", sep = ""),
                       estimate = res, samplesize = length(x), asvar = asVar,
                       asbias = asBias, pIC = IC, steps = steps, Infos = Infos))
        }
    })
setMethod("kStepEstimator", signature(x = "matrix", 
                                      IC = "IC",
                                      start = "Estimate"),
    function(x, IC, start, steps = 1, useLast = getRobAStBaseOption("kStepUseLast")){
        es.call <- match.call()
        es.call[[1]] <- as.name("kStepEstimator")
        if(!is.integer(steps))
          steps <- as.integer(steps)
        if(steps < 1)
            stop("steps needs to be a positive integer")

        nrvalues <- dimension(IC@Curve)
        start0 <- estimate(start)
        if(is.list(start0)) start0 <- unlist(start0)
        if(nrvalues != length(start0))
            stop("dimension of slot 'estimate' of 'start' != dimension of 'Curve'")
        if(ncol(x) != IC@Curve[[1]]@Domain@dimension)
            stop("'x' has wrong dimension")

        res <- start0 + rowMeans(evalIC(IC, x), na.rm = TRUE)

        L2Fam <- eval(CallL2Fam(IC))
        Infos <- matrix(c("kStepEstimator", 
                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
                        ncol = 2)
        colnames(Infos) <- c("method", "message")
        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE

        if(steps == 1){
            if(useLast && !is(modifyIC(IC), "NULL") ){
                newParam <- param(L2Fam)
                main(newParam) <- res
                newL2Fam <- modifyModel(L2Fam, newParam)
                IC <- modifyIC(IC)(newL2Fam, IC)
                Infos <- rbind(Infos, c("kStepEstimator", 
                                        "computation of IC, asvar and asbias via useLast = TRUE"))
            }else{
                if(useLast && is(modifyIC(IC), "NULL")){
                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
                             is filled with some function!")
                    Infos <- rbind(Infos, c("kStepEstimator", 
                                            "slot 'modifyIC' of 'IC' was not filled!"))
                }
                Infos <- rbind(Infos, c("kStepEstimator", 
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
            return(new("kStepEstimate", estimate.call = es.call,
                       name = paste(steps, "-step estimate", sep = ""),
                       estimate = res, samplesize = ncol(x), asvar = asVar,
                       asbias = asBias, pIC = IC, steps = 1L, Infos = Infos))
        }else{
            if(is(modifyIC(IC), "NULL"))
                stop("slot 'modifyIC' of 'IC' is 'NULL'!")
            for(i in 2:steps){
                start0 <- res
                newL2Fam <- eval(CallL2Fam(IC))
                newParam <- param(newL2Fam)
                main(newParam) <- start0
                newL2Fam <- modifyModel(newL2Fam, newParam)
                IC <- modifyIC(IC)(newL2Fam, IC)
                res <- start0 + rowMeans(evalIC(IC, x), na.rm = TRUE)
            }
            if(useLast){
                newL2Fam <- eval(CallL2Fam(IC))
                newParam <- param(newL2Fam)
                main(newParam) <- res
                newL2Fam <- modifyModel(newL2Fam, newParam)
                IC <- modifyIC(IC)(newL2Fam, IC)
                Infos <- rbind(Infos, c("kStepEstimator", 
                                        "computation of IC, asvar and asbias via useLast = TRUE"))
            }else{
                Infos <- rbind(Infos, c("kStepEstimator", 
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
            return(new("kStepEstimate", estimate.call = es.call,
                       name = paste(steps, "-step estimate", sep = ""),
                       estimate = res, samplesize = ncol(x), asvar = asVar,
                       asbias = asBias, pIC = IC, steps = steps, Infos = Infos))
        }
    })
