.getpsi <- function(nameInSysdata =".OMSE",
                    PFam = GParetoFamily(),
                    with.correct=FALSE, r=.5){

   sng <- try(getFromNamespace(nameInSysdata, ns = "ROptEst"),silent=TRUE)
   if(is(sng,"try-error"))
      stop("Grid of L.M.s for scale-shape family not yet availabe.")
   sngr <- sng[[name(PFam)]][["grid"]]
   snf <- sng[[name(PFam)]][["fct"]]
   ncgr <- ncol(sngr)
   b <- function(xi0) sng(xi0,1)
   a <- function(xi0) c(sng(xi0,2),sng(xi0,3))
   aw <- function(xi0) c(sng(xi0,4),sng(xi0,5))
   A <- function(xi0) {am <- mean(c(sng(xi0,7),sng(xi0,8)))
                       matrix(c(sng(xi0,6),am,am,sng(xi0,9)),2,2)}
   Aw <- function(xi0) {am <- mean(c(sng(xi0,11),sng(xi0,12)))
                       matrix(c(sng(xi0,10),am,am,sng(xi0,13)),2,2)}

   normt <- NormType()
   biast <- symmetricBias()
   ICT <- paste("optimally robust IC for", switch(nameInSysdata,
                      c(".OMSE"="maxMSE",".RMXE"="RMX", ".MBRE"="maxBias")))
   riskT <- if(nameInSysdata!=".MBRE") "asGRisk" else "asBias"

   res.xi <- function(xi0){
         neighbor <- ContNeighborhood(radius = r.0(xi0))
         w <- new("HampelWeight")
         stand(w) <- Aw(xi0)
         cent(w) <- aw(xi0)
         clip(w) <- b(xi0)
         if(nameInSysdata!=".MBRE")
            weight(w) <- getweight(w, neighbor = neighbor, biastype = biast,
                                normW = normt)
         else weight(w) <- minbiasweight(w, neighbor = neighbor, biastype = biast,
                                normW = normt)
         res = c(a = a(xi0), A = A(xi0), b = b(xi0), d = 0,
                    normtype = normt, biastype = biast, w = w,
                    info = c("optIC", ICT), risk = riskT,
                    modifyIC = NULL
                )
          return(res)
     }
   IC.fct.xi.0 <- function(xi0){
         param0 <- ParamFamParameter(name = "theta", main = c(1,xi0),
                               fixed = 0,
                               trafo =  param(PFam)@trafo)
         PFam@param <- param0
         Lx <- PFam@L2deriv.fct(param0)
         IC.fctx <- function(x){
                   Lx.v <- t(cbind(Lx[[1]](x),Lx[[2]](x)))
                   Y <- A(xi0)%*%(Lx.v-a(xi0))*weight(w)(x)
                   return(Y)
         }
         return(IC.fctx)
   }
   IC.xi.0 <- function(xi0){
         param0 <- ParamFamParameter(name = "theta", main = c(1,xi0),
                               fixed = 0,
                               trafo =  param(PFam)@trafo)
         PFam@param <- param0
         res0 <- res(xi0)
         res0$modifyIC <- function(L2Fam, IC){
             param1 <- param(L2Fam)
             xi <- main(param1)["shape"]
             beta <- main(param1)["scale"]
             IC <- IC.xi(xi)
             IC1M <- IC[[1]]@Map
             IC2M <- list(function(x)IC[[1]]@Map[[1]](x/beta),
                          function(x)IC[[1]]@Map[[2]](x/beta))
             IC[[1]]@Map <- IC2M
             return(IC)
         }

         IC <- generateIC(neighbor= neighbor,
                                     L2Fam = PFam, res = res0)
         return(IC)
   }
   IC.fct.xi <- IC.fct.xi.0
   IC.xi <- IC.xi.0
   if(with.correct){
      IC.fct.xi <- function(xi0){
            res0 <- IC.fct.xi.0(xi0)
            distr <- distribution(PFam)

            param1 <- param(distribution)
            param1["shape"] <- xi0
            param1["scale"] <- 1
            distr@param <- param1
            
            I.w <- E(distr, fun=weight(w))
            I.1 <- E(distr, fun=function(x) weight(w)(x)*Lx[[1](x))
            I.2 <- E(distr, fun=function(x) weight(w)(x)*Lx[[2](x))
            z.c <- c(I.1,I.2)/I.w

            L1.0 <- function(x) Lx[[1]](x)-z.c[1]
            L2.0 <- function(x) Lx[[2]](x)-z.c[2]
            I.11 <- E(distr, fun=function(x) weight(w)(x)*L1.0(x)^2)
            I.21 <- E(distr, fun=function(x) weight(w)(x)*L1.0(x)*L2.0(x))
            I.22 <- E(distr, fun=function(x) weight(w)(x)*L2.0(x)^2)

            A.c <- solve(matrix(c(I.11,I.21,I.21,I.22),2,2))
            IC.fct <- function(x) A.c%*%(res0$IC.fct(x)-z.c)
            return(IC.fct)
         }
      IC.xi <- function(xi0){
            res0 <- IC.xi.0(xi0)
            param0 <- ParamFamParameter(name = "theta", main = c(1,xi0),
                               fixed = 0,
                               trafo =  param(PFam)@trafo)
            PFam@param <- param0
            IC <- makeIC(res0$IC,PFam)
            return(IC)
            }
   }
 return(list(IC.fct.xi=IC.fct.xi,IC.xi=IC.xi))
}

roptestQuick <- function(x, PFam=GParetoFamily(), type=c("RMXE","OMSE","MBRE")){
   theta0 <- estimate(medkMADhybr(x,PFam=PFam))
   type <- match.arg(type)
   nametype <- paste(".",type,sep="")
   ICfct <- .getpsi(nameInSysdata = nametype,
                    PFam = PFam)$IC.fct.xi
   IC <- ICfct(xi=theta0["shape"])
   return(theta0+ rowMeans(IC(x)))
}

OMSE.5 <- function(x, PFam = GParetoFamily()) roptestQuick(x,PFam,"OMSE")
RMXE <- function(x, PFam = GParetoFamily())   roptestQuick(x,PFam,"RMXE")
MBRE <- function(x, PFam = GParetoFamily())   roptestQuick(x,PFam,"MBRE")

roptestScSh <- function(x, PFam=GParetoFamily(), type=c("RMXE","OMSE","MBRE"),
                        initial.est, steps = 1L, verbose = NULL,
                        useLast = getRobAStBaseOption("kStepUseLast"),
                    withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
                    IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
                    withICList = getRobAStBaseOption("withICList"),
                    withPICList = getRobAStBaseOption("withPICList"),
                    na.rm = TRUE, initial.est.ArgList, ...,
                    scalename = "scale", withLogScale = TRUE){
    if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")
    es.call <- match.call()

    type <- match.arg(type)
    nametype <- paste(".",type,sep="")

    if(missing(x))
        stop("'x' is missing with no default")
    if(missing(PFam))
        stop("'PFam' is missing with no default")
    if(!is.numeric(x)){
        if(is.data.frame(x))
            x <- data.matrix(x)
        else
            x <- as.matrix(x)
        if(!is.matrix(x))
            stop("'x' has to be a numeric vector resp. a matrix or data.frame")
    }
    completecases <- complete.cases(x)
    if(na.rm) x <- na.omit(x)

        if(!is.integer(steps))
        steps <- as.integer(steps)
    if(steps < 1){
        stop("'steps' has to be some positive integer value")
    }
    if(length(steps) != 1){
        stop("'steps' has to be of length 1")
    }

    if(missing(initial.est))
       initial.est <- medkMADhybr(x,PFam=PFam)
        
    nrvalues <-  length(PFam@param)
    initial.est <- kStepEstimator.start(initial.est, x = x,
                                        nrvalues = nrvalues, na.rm = na.rm,
                                        L2Fam = PFam,
                                        startList = initial.est.ArgList)

    newParam <- param(PFam)
    main(newParam)[] <- as.numeric(initial.est)
    L2FamStart <- modifyModel(PFam, newParam)
    ICStart <- .getpsi(nameInSysdata = nametype,
                    PFam = PFam, with.correct=TRUE)$IC.xi(main(newParam)[2])

    res <- kStepEstimator(x, IC = ICstart, start = initial.est, steps = steps,
                          useLast = useLast, withUpdateInKer = withUpdateInKer,
                          IC.UpdateInKer = IC.UpdateInKer,
                          withICList = withICList, withPICList = withPICList,
                          na.rm = na.rm, scalename = "scale", withLogScale = TRUE)

    res@estimate.call <- es.call
    Infos <- matrix(c("roptestScSh",
                      paste(steps, "-step estimate for ", name(PFam), sep = "")),
                    ncol = 2)
    colnames(Infos) <- c("method", "message")

    if(! distrMod:::.isUnitMatrix(trafo(PFam)))
       Infos <- rbind(Infos, c("roptestScSh",
                            paste("computation of IC",
                                   ifelse(withUpdateInKer,"with","without") ,
                                   "modification in ker(trafo)")))

    Infos <- rbind(Infos, c("roptestScSh",
                            paste("computation of IC, asvar and asbias via useLast =", useLast)))
    Infos(res) <- Infos
    res@completecases <- completecases
    res@start <- initial.est
    return(res)
}
