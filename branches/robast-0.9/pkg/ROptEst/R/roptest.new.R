###############################################################################
## Optimally robust estimation
###############################################################################

roptest.n <- function(x, L2Fam, eps, eps.lower, eps.upper, fsCor = 1, initial.est,
                    neighbor = ContNeighborhood(), risk = asMSE(), steps = 1L,
                    distance = CvMDist, startPar = NULL, verbose = NULL,
                    OptOrIter = "iterate",
                    useLast = getRobAStBaseOption("kStepUseLast"),
                    withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
                    IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
                    withICList = getRobAStBaseOption("withICList"),
                    withPICList = getRobAStBaseOption("withPICList"),
                    na.rm = TRUE, initial.est.ArgList, ...,
                    withLogScale = TRUE,..withCheck=FALSE){
    es.call <- match.call()
    es.call0 <- match.call(expand.dots=FALSE)
    dots <- es.call0[["..."]]
    es.call0$"..." <- NULL
    es.call1 <- .constructArg.list(roptest.n,es.call0, onlyFormal=FALSE,
                            debug = ..withCheck)$mc

    res <- .constructArg.list(gennbCtrl,es.call1, onlyFormal=TRUE,
                    debug = ..withCheck)
     nbCtrl0 <-  as.list(res$mc[-1])
     es.call1 <- res$esc
    res <- .constructArg.list(genstartCtrl,es.call1, onlyFormal=TRUE,
                   debug = ..withCheck)
     startCtrl0 <-  as.list(res$mc[-1])
     es.call1 <- res$esc
    res <- .constructArg.list(genkStepCtrl,es.call1, onlyFormal=TRUE,
                   debug = ..withCheck)
     kStepCtrl0 <-  as.list(res$mc[-1])
     es.call1 <- res$esc

    list1 <- as.list(es.call1[-1])
    res <- do.call(gennbCtrl,nbCtrl0)
      if(length(res)) list1$nbCtrl <- res
    res <- do.call(genstartCtrl,startCtrl0)
      if(length(res)) list1$startCtrl <- res
    res <- do.call(genkStepCtrl,kStepCtrl0)
      if(length(res)) list1$kStepCtrl <- res
    list1 <- c(list1,dots)
    list1$..withCheck <- NULL
    if(..withCheck) print(list1)
    if(..withCheck) return(substitute(do.call(roptEst, lis), list(lis=list1)))
    else{
    res <- do.call(roptEst, list1)
    res@estimate.call <- es.call
    return(res)}
}
#roptest.n(x=1:10,L2Fam=GammaFamily(),also=3,..withCheck=TRUE)

#roptest.n(x=1:10,L2Fam=GammaFamily(),also=3)

roptEst <- function(x, L2Fam,  fsCor = 1,
                     risk = asMSE(), steps = 1L,
                      verbose = NULL,
                    OptOrIter = "iterate",
                    nbCtrl = gennbCtrl(),
                    startCtrl = genstartCtrl(),
                    kStepCtrl = genkStepCtrl(),
                    na.rm = TRUE, ..., debug = FALSE){

    es.call <- match.call()
    es.call0 <- match.call(expand.dots=FALSE)
    dots <- es.call0[["..."]]
    es.call0$"..." <- NULL
#    if(!is.null(dots[["escall"]])) es.call <- dots[["escall"]]

    if(missing(nbCtrl)||is.null(nbCtrl))          nbCtrl <- gennbCtrl()
    if(missing(startCtrl)||is.null(startCtrl)) startCtrl <- genstartCtrl()
    if(missing(kStepCtrl)||is.null(kStepCtrl)) kStepCtrl <- genkStepCtrl()
    nbCtrl    <- .fix.in.defaults(nbCtrl, gennbCtrl)
    startCtrl <- .fix.in.defaults(startCtrl, genstartCtrl)
    kStepCtrl <- .fix.in.defaults(kStepCtrl, genkStepCtrl)
        
    es.list <- as.list(es.call0[-1])
    es.list <- c(es.list,nbCtrl)
    es.list$nbCtrl <- NULL
    es.list0 <-  c(es.list,dots)

    if(missing(verbose)|| is.null(verbose))
           es.list$verbose <- getRobAStBaseOption("all.verbose")

    if(missing(L2Fam))
        stop("'L2Fam' is missing with no default")
    if(!is(L2Fam, "L2ParamFamily"))
        stop("'L2Fam' must be of class 'L2ParamFamily'.")


    res.x <- .pretreat(x,na.rm)
    x <- res.x$x
    completecases <- res.x$completecases

    es.list[["x"]] <- x
    if(debug){   cat("\nes.list:\n");print(es.list);cat("\n\n")}

    .isOKfsCor(fsCor)

    .isOKsteps(steps)

    if(debug){
      if(is.null(startCtrl$initial.est)){
       print(substitute(MDEstimator(x = x0, ParamFamily = L2Fam0,
                        distance = distance0,
                        startPar = startPar0,dots0),
                        list(x0=x,L2Fam0=L2Fam, distance0=startCtrl$distance,
                             startPar0=startCtrl$startPar0, dots0=dots)))
        startCtrl$initial.est <- "BLUB"
      }
    }else{
      if(is.null(startCtrl$initial.est))
         startCtrl$initial.est <- MDEstimator(x = x, ParamFamily = L2Fam,
                                  distance = startCtrl$distance,
                                  startPar = startCtrl$startPar, ...)
    }
    nrvalues <-  length(L2Fam@param)

    if(debug){
      if(!is.null(startCtrl$initial.est.ArgList)){
       initial.est <- print(substitute({kStepEstimator.start(initial.est0, x = x0,
                                        nrvalues = nrvalues0, na.rm = na.rm0,
                                        L2Fam = L2Fam0,
                                        startList = startCtrl0)},
                                        list(x0=x,L2Fam0=L2Fam,
                                        initial.est0=startCtrl$initial.est,
                                        nrvalues0=nrvalues,na.rm0=na.rm,
                                     startCtrl0 = startCtrl$initial.est.ArgList)
                                 ))
      }else{
       print(substitute(kStepEstimator.start(initial.est0, x = x0,
                                        nrvalues = nrvalues0, na.rm = na.rm0,
                                        L2Fam = L2Fam0),list(x0=x,L2Fam0=L2Fam,
                                        initial.est0=startCtrl$initial.est,
                                        nrvalues0=nrvalues,na.rm0=na.rm)
                                 ))
      }
    }else{
    sy <- system.time({
      initial.est <-  kStepEstimator.start(eval(startCtrl$initial.est), x = x,
                                        nrvalues = nrvalues, na.rm = na.rm,
                                        L2Fam = L2Fam)
     })
     print(sy)
    }

    newParam <- param(L2Fam)

    if(!debug){
       main(newParam)[] <- as.numeric(initial.est)
       L2FamStart <- modifyModel(L2Fam, newParam)
    }
    if(debug) print(risk)

    if(!is(risk,"interpolRisk"))
       es.list0[["eps"]] <- do.call(.check.eps, args=c(nbCtrl,list("x"=x)))
    es.list0$risk <- NULL
    es.list0$L2Fam <- NULL
    neighbor <- eval(es.list0$neighbor)
    es.list0$neighbor <- NULL
    es.list0$kStepCtrl <- NULL
    es.list0$startCtrl <- NULL
    es.list0$distance <- NULL
    es.list0$x <- NULL
    es.list0$steps <- NULL
    es.list0$eps.lower <- NULL
    es.list0$eps.upper <- NULL
    es.list0$na.rm <- NULL
    if(debug) {cat("\n\n\n::::\n\n")
    argList <- c(list(model=L2Fam,risk=risk,neighbor=neighbor),
                                             es.list0)
    print(argList)
    cat("\n\n\n")
    }
    if(!debug){
      sy <-  system.time({
       ICstart <- do.call(getStartIC, args=c(list(model=L2FamStart,risk=risk,
                              neighbor=neighbor),es.list0))
     })
     print(sy)
     }
      if(debug){
         argList <- list(x, IC = ICstart, start = initial.est, steps = steps,
                            useLast = kStepCtrl$useLast,
                            withUpdateInKer = kStepCtrl$withUpdateInKer,
                            IC.UpdateInKer = kStepCtrl$IC.UpdateInKer,
                            withICList = kStepCtrl$withICList,
                            withPICList = kStepCtrl$withPICList,
                            na.rm = na.rm,
                            scalename = kStepCtrl$scalename,
                            withLogScale = kStepCtrl$withLogScale)
         print(argList) }
      sy <- system.time({
         res <- kStepEstimator(x, IC = ICstart, start = initial.est, steps = steps,
                            useLast = kStepCtrl$useLast,
                            withUpdateInKer = kStepCtrl$withUpdateInKer,
                            IC.UpdateInKer = kStepCtrl$IC.UpdateInKer,
                            withICList = kStepCtrl$withICList,
                            withPICList = kStepCtrl$withPICList,
                            na.rm = na.rm,
                            scalename = kStepCtrl$scalename,
                            withLogScale = kStepCtrl$withLogScale)
                            })
       print(sy)

    if(!debug) res@estimate.call <- es.call
    Infos <- matrix(c("roptest", 
                      paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
                     ncol = 2)
    colnames(Infos) <- c("method", "message")


    if(! distrMod:::.isUnitMatrix(trafo(L2Fam)))
       Infos <- rbind(Infos, c("method"="roptest",
                    "message"=paste("computation of IC",
                               ifelse(kStepCtrl$withUpdateInKer,"with","without") ,
                                   "modification in ker(trafo)")))

    Infos <- rbind(Infos, c("method"="roptest",
     "message"=paste("computation of IC, asvar and asbias via useLast =",
                  kStepCtrl$useLast)))

    if(debug) print(Infos)
    if(debug) return(NULL)

    Infos(res) <- Infos
    res@completecases <- completecases
    res@start <- initial.est
    return(res)
}
