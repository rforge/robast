###############################################################################
## k-step estimator
###############################################################################

### helper function from distrMod:
.isUnitMatrix <- distrMod:::.isUnitMatrix


### no dispatch on top layer -> keep product structure of dependence
kStepEstimator <- function(x, IC, start, steps = 1L,
                           useLast = getRobAStBaseOption("kStepUseLast"),
                           withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
                           IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
                           na.rm = TRUE, ...){
## save call

        es.call <- match.call()
        es.call[[1]] <- as.name("kStepEstimator")

## get some dimensions
        L2Fam <- eval(CallL2Fam(IC))
        Param <- param(L2Fam)
        tf <- trafo(L2Fam,Param)
        Dtau <- tf$mat
        trafoF <- tf$fct
        p <- nrow(Dtau)
        k <- ncol(Dtau)

## check input
        if(!is.integer(steps))
          steps <- as.integer(steps)
        if(steps < 1)
            stop("steps needs to be a positive integer")
        if(! is(IC, "IC"))
           stop("Argument 'IC' must be of class 'IC'")

### transform if necessary
        x0 <- x
        x0 <- if(is.numeric(x) && ! is.matrix(x)) {
                x0 <- as.matrix(x)
                }
        completecases <- complete.cases(x0)
        if(na.rm) x0 <- na.omit(x0)

### use dispatch here  (dispatch only on start)

        start.val <- kStepEstimator.start(start, x=x0, nrvalues = k, na.rm = na.rm, ...)


### a starting value in k-space
        u.theta <- start.val
        theta <- if(is(start.val,"Estimate")) estimate(start.val)
                 else trafoF(u.theta)$fval

        ### update - function
        updateStep <- function(u.theta, theta, IC, L2Fam, Param, withModif = TRUE){

                IC.c <- as(diag(p) %*% IC@Curve, "EuclRandVariable")


                theta <- theta + rowMeans(evalRandVar(IC.c, x0),
                                          na.rm = na.rm )

                tf <- trafo(L2Fam, Param)
                Dtau <- tf$mat

                if(!.isUnitMatrix(Dtau)){
                     Dminus <- solve(Dtau, generalized = TRUE)
                     projtau <- diag(k) - Dminus %*% Dtau

                     IC.tot1 <- Dminus %*% IC.c
                     IC.tot2 <- 0 * IC.tot1

                     if(sum(diag(projtau))>0.5 && ### is EM-D^-D != 0 (i.e. rk D<p)
                        withUpdateInKer){
                            if(!is.null(IC.UpdateInKer)&&!is(IC.UpdateInKer,"IC"))
                               warning("'IC.UpdateInKer' is not of class 'IC'; we use default instead.")
                            IC.tot2 <- if(is.null(IC.UpdateInKer))
                                 getBoundedIC(L2Fam, D = projtau) else
                                 as(projtau %*% IC.UpdateInKer@Curve,
                                    "EuclRandVariable")
                     }
                     IC.tot <- IC.tot1 + IC.tot2

                     u.theta <- u.theta + rowMeans(evalRandVar(IC.tot, x0),
                                                   na.rm = na.rm)
                }else{
                     IC.tot <- IC.c
                     u.theta <- theta
                }

                cnms <-  if(is.null(names(u.theta))) colnames(Dtau) else names(u.theta)
                u.var <- matrix(E(L2Fam, IC.tot %*% t(IC.tot)),
                                  k,k, dimnames = list(cnms,cnms))

                if(withModif){
                   main(Param)[] <- as.numeric(u.theta)
                   L2Fam <- modifyModel(L2Fam, Param)
                   IC <- modifyIC(IC)(L2Fam, IC)
                }

                return(list(IC = IC, Param = Param, L2Fam = L2Fam,
                            theta = theta, u.theta = u.theta, u.var = u.var))
        }

        Infos <- matrix(c("kStepEstimator",
                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
                        ncol = 2)
        colnames(Infos) <- c("method", "message")
        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE

        ### iteration


        if(!is(modifyIC(IC), "NULL") ){
           for(i in 1:steps){
               if(i>1){
                  IC <- upd$IC
                  L2Fam <- upd$L2Fam
                  Param <- upd$Param
                  tf <- trafo(L2Fam, Param)
               }
               upd <- updateStep(u.theta,theta,IC, L2Fam, Param,
                                 withModif = (steps>1) | useLast)
               u.theta <- upd$u.theta
               theta <- upd$theta
               u.var <- upd$u.var
           }
           if(useLast){
              IC <- upd$IC
              L2Fam <- upd$L2Fam
              Param <- upd$Param
              tf <- trafo(L2Fam, Param)
              Infos <- rbind(Infos, c("kStepEstimator",
                  "computation of IC, trafo, asvar and asbias via useLast = TRUE"))
           }else{
              Infos <- rbind(Infos, c("kStepEstimator",
                       "computation of IC, trafo, asvar and asbias via useLast = FALSE"))
           }
        }else{
           if(steps > 1)
              stop("slot 'modifyIC' of 'IC' is 'NULL'!")
           upd <- updateStep(u.theta,theta,IC, L2Fam, Param, withModif = FALSE)
           u.theta <- upd$u.theta
           theta <- upd$theta
           u.var <- upd$u.var
           if(useLast){
              warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
                     is filled with some function!")
              Infos <- rbind(Infos, c("kStepEstimator",
                          "slot 'modifyIC' of 'IC' was not filled!"))
           }else
            Infos <- rbind(Infos, c("kStepEstimator",
            "computation of IC, asvar and asbias via useLast = FALSE"))
        }

        if(! distrMod:::.isUnitMatrix(trafo(L2Fam)))
             Infos <- rbind(Infos, c("kStepEstimator",
                            paste("computation of IC",
                                   ifelse(withUpdateInKer,"with","without") ,
                                   "modification in ker(trafo)")))

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

        if(is.null(names(theta))) names(theta) <- rownames(Dtau)
        if(is.null(names(u.theta))) names(u.theta) <- colnames(u.theta)
        dim(theta) <- NULL
        dim(u.theta) <- NULL
        
        nuis.idx <- if(is(start,"Estimate")) start@nuis.idx else NULL
        fixed <- if(is(start,"Estimate")) start@fixed else NULL

        return(new("kStepEstimate", estimate.call = es.call,
                       name = paste(steps, "-step estimate", sep = ""),
                       estimate = theta, samplesize = nrow(x0), asvar = asVar,
                       trafo = tf, fixed = fixed,
                       nuis.idx = nuis.idx, untransformed.estimate = u.theta,
                       completecases= completecases,
                       untransformed.asvar = u.var,
                       asbias = asBias, pIC = IC, steps = steps, Infos = Infos))
}

#### old method:

# setMethod("kStepEstimator", signature(x = "numeric",
#                                      IC = "IC",
#                                      start = "numeric"),
#    function(x, IC, start, steps = 1L, useLast = getRobAStBaseOption("kStepUseLast")){
#        es.call <- match.call()
#        es.call[[1]] <- as.name("kStepEstimator")
#        if(!is.integer(steps))
#          steps <- as.integer(steps)
#        if(steps < 1)
#            stop("steps needs to be a positive integer")
#
#        nrvalues <- dimension(IC@Curve)
#        if(is.list(start)) start <- unlist(start)
#        if(nrvalues != length(start))
#            stop("dimension of 'start' != dimension of 'Curve'")
#
#        res <- start + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)
#
#        L2Fam <- eval(CallL2Fam(IC))
#        Infos <- matrix(c("kStepEstimator",
#                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
#                        ncol = 2)
#        colnames(Infos) <- c("method", "message")
#        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
#
#        if(steps == 1){
#            if(useLast && !is(modifyIC(IC), "NULL") ){
#                newParam <- param(L2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(L2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                if(useLast && is(modifyIC(IC), "NULL")){
#                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
#                             is filled with some function!")
#                    Infos <- rbind(Infos, c("kStepEstimator",
#                                            "slot 'modifyIC' of 'IC' was not filled!"))
#                }
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = length(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = 1L, Infos = Infos))
#        }else{
#            if(is(modifyIC(IC), "NULL"))
#                stop("slot 'modifyIC' of 'IC' is 'NULL'!")
#            for(i in 2:steps){
#                start <- res
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- start
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                res <- start + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)
#            }
#            if(useLast){
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = length(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = steps, Infos = Infos))
#        }
#    })
#
#setMethod("kStepEstimator", signature(x = "matrix",
#                                      IC = "IC",
#                                      start = "numeric"),
#    function(x, IC, start, steps = 1, useLast = getRobAStBaseOption("kStepUseLast")){
#        es.call <- match.call()
#        es.call[[1]] <- as.name("kStepEstimator")
#        if(!is.integer(steps))
#          steps <- as.integer(steps)
#        if(steps < 1)
#            stop("steps needs to be a positive integer")
#
#        nrvalues <- dimension(IC@Curve)
#        if(is.list(start)) start <- unlist(start)
#        if(nrvalues != length(start))
#            stop("dimension of 'start' != dimension of 'Curve'")
#        if(ncol(x) != IC@Curve[[1]]@Domain@dimension)
#            stop("'x' has wrong dimension")
#
#        res <- start + rowMeans(evalIC(IC, x), na.rm = TRUE)
#
#        L2Fam <- eval(CallL2Fam(IC))
#        Infos <- matrix(c("kStepEstimator",
#                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
#                        ncol = 2)
#        colnames(Infos) <- c("method", "message")
#        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
#
#        if(steps == 1){
#            if(useLast && !is(modifyIC(IC), "NULL") ){
#                newParam <- param(L2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(L2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                if(useLast && is(modifyIC(IC), "NULL")){
#                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
#                             is filled with some function!")
#                    Infos <- rbind(Infos, c("kStepEstimator",
#                                            "slot 'modifyIC' of 'IC' was not filled!"))
#                }
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#             }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = ncol(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = 1L, Infos = Infos))
#        }else{
#            if(is(modifyIC(IC), "NULL"))
#                stop("slot 'modifyIC' of 'IC' is 'NULL'!")
#            for(i in 2:steps){
#                start <- res
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- start
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                res <- start + rowMeans(evalIC(IC, x), na.rm = TRUE)
#            }
#            if(useLast){
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = ncol(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = steps, Infos = Infos))
#        }
#    })
#setMethod("kStepEstimator", signature(x = "numeric",
#                                      IC = "IC",
#                                      start = "Estimate"),
#    function(x, IC, start, steps = 1, useLast = getRobAStBaseOption("kStepUseLast")){
#        es.call <- match.call()
#        es.call[[1]] <- as.name("kStepEstimator")
#        if(!is.integer(steps))
#          steps <- as.integer(steps)
#        if(steps < 1)
#            stop("steps needs to be a positive integer")
#
#        nrvalues <- dimension(IC@Curve)
#        start0 <- estimate(start)
#        if(is.list(start0)) start0 <- unlist(start0)
#        if(nrvalues != length(start0))
#            stop("dimension of slot 'estimate' of 'start' != dimension of 'Curve'")
#
#        res <- start0 + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)
#
#        L2Fam <- eval(CallL2Fam(IC))
#        Infos <- matrix(c("kStepEstimator",
#                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
#                        ncol = 2)
#        colnames(Infos) <- c("method", "message")
#        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
#
#        if(steps == 1){
#            if(useLast && !is(modifyIC(IC), "NULL") ){
#                newParam <- param(L2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(L2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                if(useLast && is(modifyIC(IC), "NULL")){
#                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
#                             is filled with some function!")
#                    Infos <- rbind(Infos, c("kStepEstimator",
#                                            "slot 'modifyIC' of 'IC' was not filled!"))
#                }
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = length(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = 1L, Infos = Infos))
#        }else{
#            if(is(modifyIC(IC), "NULL"))
#                stop("slot 'modifyIC' of 'IC' is 'NULL'!")
#            for(i in 2:steps){
#                start0 <- res
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- start0
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                res <- start0 + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)
#            }
#            if(useLast){
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                 Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = length(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = steps, Infos = Infos))
#        }
#    })
#setMethod("kStepEstimator", signature(x = "matrix",
#                                      IC = "IC",
#                                      start = "Estimate"),
#    function(x, IC, start, steps = 1, useLast = getRobAStBaseOption("kStepUseLast")){
#        es.call <- match.call()
#        es.call[[1]] <- as.name("kStepEstimator")
#        if(!is.integer(steps))
#          steps <- as.integer(steps)
#        if(steps < 1)
#            stop("steps needs to be a positive integer")
#
#        nrvalues <- dimension(IC@Curve)
#        start0 <- estimate(start)
#        if(is.list(start0)) start0 <- unlist(start0)
#        if(nrvalues != length(start0))
#            stop("dimension of slot 'estimate' of 'start' != dimension of 'Curve'")
#        if(ncol(x) != IC@Curve[[1]]@Domain@dimension)
#            stop("'x' has wrong dimension")
#
#        res <- start0 + rowMeans(evalIC(IC, x), na.rm = TRUE)
#
#        L2Fam <- eval(CallL2Fam(IC))
#        Infos <- matrix(c("kStepEstimator",
#                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
#                        ncol = 2)
#        colnames(Infos) <- c("method", "message")
#        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
#
#        if(steps == 1){
#            if(useLast && !is(modifyIC(IC), "NULL") ){
#                newParam <- param(L2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(L2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                if(useLast && is(modifyIC(IC), "NULL")){
#                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
#                             is filled with some function!")
#                    Infos <- rbind(Infos, c("kStepEstimator",
#                                            "slot 'modifyIC' of 'IC' was not filled!"))
#                }
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = ncol(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = 1L, Infos = Infos))
#        }else{
#            if(is(modifyIC(IC), "NULL"))
#                stop("slot 'modifyIC' of 'IC' is 'NULL'!")
#            for(i in 2:steps){
#                start0 <- res
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- start0
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                res <- start0 + rowMeans(evalIC(IC, x), na.rm = TRUE)
#            }
#            if(useLast){
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = ncol(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = steps, Infos = Infos))
#        }
#    })

