###############################################################################
## Optimally robust estimation
###############################################################################
roptest <- function(x, L2Fam,  fsCor = 1,
                     risk = asMSE(), steps = 1L,
                      verbose = NULL,
                    OptOrIter = "iterate",
                    nbCtrl = gennbCtrl(neighbor = ContNeighborhood(),
                                       eps, eps.lower, eps.upper)
                    startCtrl = genstartCtrl(initial.est = NULL,
                                             initial.est.ArgList = NULL,
                                       startPar = NULL, distance = CvMDist),
                    kstepCtrl = genkstepCtrl(
                         useLast = getRobAStBaseOption("kStepUseLast"),
                         withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
                         IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
                         withICList = getRobAStBaseOption("withICList"),
                         withPICList = getRobAStBaseOption("withPICList"),
                         withLogScale = TRUE),
                    na.rm = TRUE, ...){

    es.call <- match.call()
    dots <- match.call(expand=FALSE)$dots
    es.list <- as.list(es.call[-1])
    es.list <- .fix.in.defaults(es.list,roptest)
    es.list <- c(es.list,nbCtrl)
    es.list$dots <- dots

    if(missing(verbose)|| is.null(verbose))
           es.list$verbose <- getRobAStBaseOption("all.verbose")

    if(missing(L2Fam))
        stop("'L2Fam' is missing with no default")

    x <- .pretreat(x,na.rm)
    
    es.list$eps <- do.call(.check.eps, args=es.list)
    
    .isOKfsCor(fsCor)

    .isOKsteps(steps)

    if(is.null(startCtrl$initial.est))
        startCtrl$initial.est <- MDEstimator(x = x, ParamFamily = L2Fam,
                                            distance = startCtrl$distance,
                                            startPar = startCtrl$startPar, ...)
    nrvalues <-  length(L2Fam@param)
    initial.est <- kStepEstimator.start(initial.est, x = x,
                                        nrvalues = nrvalues, na.rm = na.rm,
                                        L2Fam = L2Fam,
                                        startList = startCtrl$initial.est.ArgList)


    newParam <- param(L2Fam)
    main(newParam)[] <- as.numeric(initial.est)
    L2FamStart <- modifyModel(L2Fam, newParam)

    es.list0 <- es.list
    es.list$risk <- NULL
    es.list$L2Fam <- NULL
    ICstart <- do.call(getstartIC, args=c(list(model=L2Fam,risk=risk),es.list))

    res <- kStepEstimator(x, IC = ICstart, start = initial.est, steps = steps,
                          useLast = kStepCtrl$useLast,
                          withUpdateInKer = kStepCtrl$withUpdateInKer,
                          IC.UpdateInKer = kStepCtrl$IC.UpdateInKer,
                          withICList = kStepCtrl$withICList,
                          withPICList = kStepCtrl$withPICList,
                          na.rm = na.rm,
                          scalename = kstepCtrl$scalename,
                          withLogScale = kstepCtrl$withLogScale)


    res@estimate.call <- es.call
    Infos <- matrix(c("roptest", 
                      paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
                    ncol = 2)
    colnames(Infos) <- c("method", "message")

    if(! distrMod:::.isUnitMatrix(trafo(L2Fam)))
       Infos <- rbind(Infos, c("roptest",
                            paste("computation of IC",
                                   ifelse(withUpdateInKer,"with","without") ,
                                   "modification in ker(trafo)")))

    Infos <- rbind(Infos, c("roptest",
                            paste("computation of IC, asvar and asbias via useLast =", useLast)))
    Infos(res) <- Infos
    res@completecases <- completecases
    res@start <- initial.est
    return(res)
}
