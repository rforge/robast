RMXEstimator <- function(x, L2Fam, fsCor = 1, initial.est,
                    neighbor = ContNeighborhood(), steps = 1L,
                    distance = CvMDist, startPar = NULL, verbose = NULL,
                    OptOrIter = "iterate",
                    useLast = getRobAStBaseOption("kStepUseLast"),
                    withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
                    IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
                    withICList = getRobAStBaseOption("withICList"),
                    withPICList = getRobAStBaseOption("withPICList"),
                    na.rm = TRUE, initial.est.ArgList, ..., withLogScale = TRUE,
                    ..withCheck=FALSE, withTimings = FALSE, withMDE = NULL,
                    withEvalAsVar = NULL, withMakeIC = FALSE,
                    modifyICwarn = NULL, E.argList = NULL,
                    diagnostic = FALSE){

   mc <- match.call(expand.dots=FALSE)
   dots <- mc$"..."

   gsANY <- selectMethod("getStartIC", c(model="ANY",risk="ANY"))@defined
   clsL2Fam <- c(class(L2Fam))
   gsCUR <- selectMethod("getStartIC", c(model=clsL2Fam, risk="interpolRisk"))@defined
   risk0 <- asMSE()
   if(!all(all.equal(gsANY,gsCUR)==TRUE)) risk0 <- RMXRRisk()

   roptestArgList <- list(x = x, L2Fam = L2Fam, fsCor = fsCor,
                       neighbor = neighbor, risk = risk0, steps = steps,
                       distance = distance, startPar = startPar, verbose = verbose,
                       OptOrIter = OptOrIter, useLast = useLast,
                       withUpdateInKer = withUpdateInKer, IC.UpdateInKer = IC.UpdateInKer,
                       withICList = withICList, withPICList = withPICList, na.rm = na.rm,
                       withLogScale = withLogScale, ..withCheck = ..withCheck,
                       withTimings = withTimings, withMDE = withMDE,
                       withEvalAsVar = withEvalAsVar, withMakeIC = withMakeIC,
                       modifyICwarn = modifyICwarn, E.argList = E.argList,
                       diagnostic = diagnostic)

   if(!is.null(dots)) roptestArgList <- c(roptestArgList, dots)
   if(!missing(initial.est)) roptestArgList$initial.est <- initial.est
   if(!missing(initial.est.ArgList)) roptestArgList$initial.est.ArgList <- initial.est

   res <- do.call(roptest, roptestArgList)
   res@roptestCall <- res@estimate.call
   res@estimate.call <- mc
   return(res)
}

OMSEstimator <- function(x, L2Fam, eps =0.5, fsCor = 1, initial.est,
                    neighbor = ContNeighborhood(), steps = 1L,
                    distance = CvMDist, startPar = NULL, verbose = NULL,
                    OptOrIter = "iterate",
                    useLast = getRobAStBaseOption("kStepUseLast"),
                    withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
                    IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
                    withICList = getRobAStBaseOption("withICList"),
                    withPICList = getRobAStBaseOption("withPICList"),
                    na.rm = TRUE, initial.est.ArgList, ..., withLogScale = TRUE,
                    ..withCheck=FALSE, withTimings = FALSE, withMDE = NULL,
                    withEvalAsVar = NULL, withMakeIC = FALSE,
                    modifyICwarn = NULL, E.argList = NULL,
                    diagnostic = FALSE){

   if(!is.numeric(eps)||length(eps)>1||any(eps<0))
      stop("Radius 'eps' must be given, of length 1 and non-negative.")
   mc <- match.call(expand.dots=FALSE)
   dots <- mc$"..."

   gsANY <- selectMethod("getStartIC", c(model="ANY",risk="ANY"))@defined
   clsL2Fam <- c(class(L2Fam))
   gsCUR <- selectMethod("getStartIC", c(model=clsL2Fam, risk="interpolRisk"))@defined
   risk0 <- asMSE()
   if(!all(all.equal(gsANY,gsCUR)==TRUE)&& abs(eps-0.5)<1e-3) risk0 <- OMSRRisk()

   roptestArgList <- list(x = x, L2Fam = L2Fam, eps = 0.5, fsCor = fsCor,
                       neighbor = neighbor, risk = risk0, steps = steps,
                       distance = distance, startPar = startPar, verbose = verbose,
                       OptOrIter = OptOrIter, useLast = useLast,
                       withUpdateInKer = withUpdateInKer, IC.UpdateInKer = IC.UpdateInKer,
                       withICList = withICList, withPICList = withPICList, na.rm = na.rm,
                       withLogScale = withLogScale, ..withCheck = ..withCheck,
                       withTimings = withTimings, withMDE = withMDE,
                       withEvalAsVar = withEvalAsVar, withMakeIC = withMakeIC,
                       modifyICwarn = modifyICwarn, E.argList = E.argList,
                       diagnostic = diagnostic)

   if(!is.null(dots)) roptestArgList <- c(roptestArgList, dots)
   if(!missing(initial.est)) roptestArgList$initial.est <- initial.est
   if(!missing(initial.est.ArgList)) roptestArgList$initial.est.ArgList <- initial.est

   res <- do.call(roptest, roptestArgList)
   res@roptestCall <- res@estimate.call
   res@estimate.call <- mc
   return(res)
}

OBREstimator <- function(x, L2Fam, eff=0.95, fsCor = 1, initial.est,
                    neighbor = ContNeighborhood(), steps = 1L,
                    distance = CvMDist, startPar = NULL, verbose = NULL,
                    OptOrIter = "iterate",
                    useLast = getRobAStBaseOption("kStepUseLast"),
                    withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
                    IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
                    withICList = getRobAStBaseOption("withICList"),
                    withPICList = getRobAStBaseOption("withPICList"),
                    na.rm = TRUE, initial.est.ArgList, ..., withLogScale = TRUE,
                    ..withCheck=FALSE, withTimings = FALSE, withMDE = NULL,
                    withEvalAsVar = NULL, withMakeIC = FALSE,
                    modifyICwarn = NULL, E.argList = NULL,
                    diagnostic = FALSE){

   if(!is.numeric(eff)||length(eff)>1||any(eff<0|eff>1))
      stop("Efficiency loss (in the ideal model) 'eff' must be given, of length 1 and in [0,1].")
   mc <- match.call(expand.dots=FALSE)
   dots <- mc$"..."

   risk0 <- asAnscombe(eff)

   roptestArgList <- list(x = x, L2Fam = L2Fam, fsCor = fsCor,
                       neighbor = neighbor, risk = risk0, steps = steps,
                       distance = distance, startPar = startPar, verbose = verbose,
                       OptOrIter = OptOrIter, useLast = useLast,
                       withUpdateInKer = withUpdateInKer, IC.UpdateInKer = IC.UpdateInKer,
                       withICList = withICList, withPICList = withPICList, na.rm = na.rm,
                       withLogScale = withLogScale, ..withCheck = ..withCheck,
                       withTimings = withTimings, withMDE = withMDE,
                       withEvalAsVar = withEvalAsVar, withMakeIC = withMakeIC,
                       modifyICwarn = modifyICwarn, E.argList = E.argList,
                       diagnostic = diagnostic)

   if(!is.null(dots)) roptestArgList <- c(roptestArgList, dots)
   if(!missing(initial.est)) roptestArgList$initial.est <- initial.est
   if(!missing(initial.est.ArgList)) roptestArgList$initial.est.ArgList <- initial.est

   res <- do.call(roptest, roptestArgList)
   res@roptestCall <- res@estimate.call
   res@estimate.call <- mc
   return(res)
}

MBREstimator <- function(x, L2Fam, fsCor = 1, initial.est,
                    neighbor = ContNeighborhood(), steps = 1L,
                    distance = CvMDist, startPar = NULL, verbose = NULL,
                    OptOrIter = "iterate",
                    useLast = getRobAStBaseOption("kStepUseLast"),
                    withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
                    IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
                    withICList = getRobAStBaseOption("withICList"),
                    withPICList = getRobAStBaseOption("withPICList"),
                    na.rm = TRUE, initial.est.ArgList, ..., withLogScale = TRUE,
                    ..withCheck=FALSE, withTimings = FALSE, withMDE = NULL,
                    withEvalAsVar = NULL, withMakeIC = FALSE,
                    modifyICwarn = NULL, E.argList = NULL,
                    diagnostic = FALSE){

   mc <- match.call(expand.dots=FALSE)
   dots <- mc$"..."

   gsANY <- selectMethod("getStartIC", c(model="ANY",risk="ANY"))@defined
   clsL2Fam <- c(class(L2Fam))
   gsCUR <- selectMethod("getStartIC", c(model=clsL2Fam, risk="interpolRisk"))@defined
   risk0 <- asBias()
   if(!all(all.equal(gsANY,gsCUR)==TRUE)) risk0 <- MBRRisk()

   roptestArgList <- list(x = x, L2Fam = L2Fam, fsCor = fsCor,
                       neighbor = neighbor, risk = risk0, steps = steps,
                       distance = distance, startPar = startPar, verbose = verbose,
                       OptOrIter = OptOrIter, useLast = useLast,
                       withUpdateInKer = withUpdateInKer, IC.UpdateInKer = IC.UpdateInKer,
                       withICList = withICList, withPICList = withPICList, na.rm = na.rm,
                       withLogScale = withLogScale, ..withCheck = ..withCheck,
                       withTimings = withTimings, withMDE = withMDE,
                       withEvalAsVar = withEvalAsVar, withMakeIC = withMakeIC,
                       modifyICwarn = modifyICwarn, E.argList = E.argList,
                       diagnostic = diagnostic)

   if(!is.null(dots)) roptestArgList <- c(roptestArgList, dots)
   if(!missing(initial.est)) roptestArgList$initial.est <- initial.est
   if(!missing(initial.est.ArgList)) roptestArgList$initial.est.ArgList <- initial.est

   res <- do.call(roptest, roptestArgList)
   res@roptestCall <- res@estimate.call
   res@estimate.call <- mc
   return(res)

}

