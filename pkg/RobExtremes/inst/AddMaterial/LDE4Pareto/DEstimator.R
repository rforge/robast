## optional file if one wants analoga to LDEstimators in ParetoCase.

.DMatch <- function(x.0,disp.est.0, disp.fctal.0, ParamFamily.0,
                    disp.est.ctrl.0 = NULL, disp.fctal.ctrl.0=NULL,
                    q.lo.0 =0, q.up.0=Inf, log.q.0 =TRUE, ..., vdbg=FALSE
                        ){
    dots <- list(...)

    disp.emp <- do.call(disp.est.0, args = .prepend(x.0,disp.est.ctrl.0, dots))
    q.f <- function(xi){
       distr.new <- ParamFamily.0@modifyParam(xi)
       sc.th <- do.call(disp.fctal.0, args = .prepend(distr.new,disp.fctal.ctrl.0, dots))
       val <- if(log.q.0) log(sc.th) - log(disp.emp) else
                        sc.th-disp.emp
       if(vdbg) print(val)
       return(val)
    }
    xi.01 <- try(uniroot(q.f,lower=q.lo.0,upper=q.up.0), silent=TRUE)
    if(is(xi.01, "try-error")) stop("Error in calculating LD-estimator: 'uniroot' did not converge.")
    xi.0 <- xi.01$root
    names(xi.0) <- c("shape")
    distr.new.0 <- ParamFamily.0@modifyParam(xi)
    val <-   c(xi.0, disp.emp)
    names(val) <- c("shape", "disp")
    return(val)
}

DEstimator <- function(x, disp.est, disp.fctal, ParamFamily,
                        disp.est.ctrl = NULL, disp.fctal.ctrl=NULL,
                        q.lo =1e-3, q.up=15, log.q =TRUE,
                        name, Infos, asvar = NULL, nuis.idx = NULL,
                        trafo = NULL, fixed = NULL, asvar.fct  = NULL, na.rm = TRUE,
                        ..., .withEvalAsVar = FALSE, vdbg = FALSE){
    param0 <- main(param(ParamFamily))
    name.est <- "DEstimator"
    es.call <- match.call()
    if(missing(name))
        name <- "Some estimator"
    Dname <- paste("Dispersion:",
                     paste(deparse(substitute(disp.fctal))))

    DMval <- NULL
    estimator <- function(x,...){
         DMval <<- .LDMatch(x.0= x,
                         disp.est.0 =  disp.est,
                         disp.fctal.0 =  disp.fctal,
                         ParamFamily.0 = ParamFamily,
                         disp.est.ctrl.0 = disp.est.ctrl,
                         disp.fctal.ctrl.0 = disp.fctal.ctrl,
                         q.lo.0 = q.lo,
                         q.up.0 = q.up,
                         log.q.0 = log.q, vdbg = vdbg)
         return(LDMval[1])
    }

    asvar.fct0 <- asvar.fct
    asvar.0 <- asvar
    nuis.idx.0 <- nuis.idx
    trafo.0 <- trafo
    if(is.null(fixed)) fixed <- fixed(ParamFamily)
    fixed.0 <- fixed
    na.rm.0 <- na.rm

    estimate <- Estimator(x, estimator, name, Infos,
                      asvar = asvar.0, nuis.idx = nuis.idx.0,
                      trafo = trafo.0, fixed = fixed.0,
                      asvar.fct = asvar.fct0,
                      na.rm = na.rm.0, ...,
                      .withEvalAsVar = .withEvalAsVar,
                      ParamFamily = ParamFamily)

    estimate@estimate.call <- es.call

    if(missing(Infos))
        Infos <- matrix(c("DEstimator", Dname),
                           ncol=2, dimnames=list(character(0), c("method", "message")))
    else{
        Infos <- matrix(c(rep("DEstimator", length(Infos)+1), c(Dname,Infos)),
                          ncol = 2)
        colnames(Infos) <- c("method", "message")
    }
    estimate@Infos <- Infos

    estim <- new("DEstimate")

    sln <- names(getSlots(class(estimate)))
    for( i in 1:length(sln))
        slot(estim, sln[i]) <- slot(estimate, sln[i])
    rm(estimate)

    estim@dispersion <- DMval["disp"]

    return(.checkEstClassForParamFamily(ParamFamily,estim))
}

kMADShapeEstimator <- function(x, ParamFamily, k=1, q.lo =1e-3, q.up=15, nuis.idx = NULL,
                        trafo = NULL, fixed = NULL, asvar.fct = NULL, na.rm = TRUE,
                        ..., .withEvalAsVar = FALSE, vdbg = FALSE){
      es.call <- match.call()
      if(missing(k)) k <- 1

      if (is.null(asvar.fct)){asvar.fct <- asvarkMAD
                              asvar <- asvarkMAD(ParamFamily, k=k)}


      es <- LDEstimator(x, disp.est = kMAD, disp.fctal = kMAD,
                     ParamFamily = ParamFamily,
                     disp.est.ctrl = list(k=k, na.rm = na.rm),
                     disp.fctal.ctrl=list(k=k),
                     q.lo =q.lo, q.up=q.up, log.q=TRUE,
                     name = "kMADShEst", Infos="kMADShEst",
                     asvar = asvar, nuis.idx = nuis.idx, trafo = trafo, fixed = fixed,
                     asvar.fct = asvar.fct, na.rm = na.rm, ...,
                      .withEvalAsVar = .withEvalAsVar, vdbg = vdbg)
      es@estimate.call <- es.call

      return(.checkEstClassForParamFamily(ParamFamily,es))
                     }

QnShapeEstimator <- function(x,  ParamFamily, q.lo =1e-3, q.up=15, nuis.idx = NULL,
                        trafo = NULL, fixed = NULL, asvar.fct = NULL, na.rm = TRUE,
                        ..., .withEvalAsVar = FALSE){
    es.call <- match.call()
    es <- LDEstimator(x, disp.est = Qn, disp.fctal = Qn,
                     ParamFamily = ParamFamily,
                     disp.est.ctrl = list(constant=1,na.rm = na.rm),
                     disp.fctal.ctrl = NULL,
                     q.lo =q.lo, q.up=q.up, log.q=TRUE,
                     name = "QnShEst", Infos="QnShEst",
                     asvar = NULL, nuis.idx = nuis.idx, trafo = trafo, fixed = fixed,
                     asvar.fct = asvar.fct, na.rm = na.rm, ...,
                      .withEvalAsVar = .withEvalAsVar)
      es@estimate.call <- es.call
      return(.checkEstClassForParamFamily(ParamFamily,es))
                     }

SnShapeEstimator <- function(x, ParamFamily, q.lo =1e-3, q.up=10, nuis.idx  = NULL,
                        trafo = NULL, fixed = NULL, asvar.fct = NULL, na.rm = TRUE,
                        accuracy = 100, ..., .withEvalAsVar = FALSE){
      es.call <- match.call()
      es <- LDEstimator(x, disp.est = Sn, disp.fctal = Sn,
                     ParamFamily = ParamFamily,
                     disp.est.ctrl = list(constant=1,na.rm = na.rm),
                     disp.fctal.ctrl = list(accuracy=accuracy),
                     q.lo =q.lo, q.up=q.up, log.q=TRUE,
                     name = "SnShEst", Infos="SnShEst",
                     asvar = NULL, nuis.idx = nuis.idx, trafo = trafo, fixed = fixed,
                     asvar.fct = asvar.fct, na.rm = na.rm, ...,
                      .withEvalAsVar = .withEvalAsVar)
      es@estimate.call <- es.call
      return(.checkEstClassForParamFamily(ParamFamily,es))
      }

kMADhybrShapeEstimator <- function(x, ParamFamily, k=1, q.lo =1e-3, q.up=15,
                        KK=20, nuis.idx = NULL,
                        trafo = NULL, fixed = NULL, asvar.fct = NULL, na.rm = TRUE,
                        ..., .withEvalAsVar = FALSE){
 i <- 1
 es <- try(kMADShapeEstimator(x, ParamFamily = ParamFamily, k = k,
                            q.lo = q.lo, q.up = q.up,
                            nuis.idx = nuis.idx, trafo = trafo,
                            fixed = fixed, asvar.fct = asvar.fct, na.rm = na.rm,
                             ..., .withEvalAsVar = FALSE),
                             silent=TRUE)
 if(! any(is.na(estimate(es))) && !is(es,"try-error"))
   {return(.checkEstClassForParamFamily(ParamFamily,es))}

 k1 <- 3.23
 while(i<KK){
      i <- i + 1
      es <- try(kMADShapeEstimator(x, k = k1, ParamFamily = ParamFamily,
                            q.lo = q.lo, q.up = q.up,
                            nuis.idx = nuis.idx, trafo = trafo,
                            fixed = fixed, asvar.fct = asvar.fct, na.rm = na.rm,
                             ..., .withEvalAsVar = FALSE), silent=TRUE)
      k1 <- k1 * 3
      if(! any(is.na(es)) && !is(es,"try-error"))
         {if(!missing(asvar.fct)) if(!is.null(asvar.fct)) if(.withEvalAsVar){
             if(is.call(es@asvar)) es@asvar <- eval(es@asvar)
             if(is.call(es@untransformed.asvar))
                es@untransformed.asvar <- eval(es@untransformed.asvar)
             }
          return(.checkEstClassForParamFamily(ParamFamily,es))}
      }
 return(c("scale"=NA,"shape"=NA))
}

setMethod("dispersion", "DEstimate", function(object) object@dispersion)
setMethod("location", "LEstimate", function(object) object@location)

medShapeEstimator <- function(x, ParamFamily, q.lo =1e-3, q.up=10, nuis.idx  = NULL,
                        trafo = NULL, fixed = NULL, asvar.fct = NULL, na.rm = TRUE,
                        accuracy = 100, ..., .withEvalAsVar = FALSE){
      es.call <- match.call()
      if (is.null(asvar.fct)){asvar.fct <- asvarMedShape
                              asvar <- asvarMedShape(ParamFamily)}
    Lname <- paste("Location: Median")
    LMval <- NULL
    estimator <- function(x,...){
         LMval <<- .DMatch(x.0= x,
                         disp.est.0 =  median,
                         disp.fctal.0 =  median,
                         ParamFamily.0 = ParamFamily,
                         disp.est.ctrl.0 = list(na.rm = na.rm),,
                         disp.fctal.ctrl.0 = NULL,
                         q.lo.0 = q.lo,
                         q.up.0 = q.up,
                         log.q.0 = log.q, vdbg = vdbg)
         return(LMval[1])
    }

    asvar.fct0 <- asvar.fct
    asvar.0 <- asvar
    nuis.idx.0 <- nuis.idx
    trafo.0 <- trafo
    if(is.null(fixed)) fixed <- fixed(ParamFamily)
    fixed.0 <- fixed
    na.rm.0 <- na.rm

    estimate <- Estimator(x, estimator, name, Infos,
                      asvar = asvar.0, nuis.idx = nuis.idx.0,
                      trafo = trafo.0, fixed = fixed.0,
                      asvar.fct = asvar.fct0,
                      na.rm = na.rm.0, ...,
                      .withEvalAsVar = .withEvalAsVar,
                      ParamFamily = ParamFamily)

    estimate@estimate.call <- es.call

    if(missing(Infos))
        Infos <- matrix(c("DEstimator", Lname),
                           ncol=2, dimnames=list(character(0), c("method", "message")))
    else{
        Infos <- matrix(c(rep("DEstimator", length(Infos)+1), c(Lname,Infos)),
                          ncol = 2)
        colnames(Infos) <- c("method", "message")
    }
    estimate@Infos <- Infos

    estim <- new("LEstimate")

    sln <- names(getSlots(class(estimate)))
    for( i in 1:length(sln))
        slot(estim, sln[i]) <- slot(estimate, sln[i])
    rm(estimate)

    estim@location <- LMval[2]
    estim@estimate.call <- es.call

    return(.checkEstClassForParamFamily(ParamFamily,estim))
    }

setClass("DEstimate",
         representation(dispersion = "numeric"
                        ),
         prototype(name = "Disp estimate",
                   estimate = numeric(0),
                   samplesize = numeric(0),
                   completecases = logical(0),
                   asvar = NULL,
                   estimate.call = call("{}"),
                   location = 0,
                   dispersion = 1,
                   Infos = matrix(c(character(0),character(0)), ncol=2,
                                  dimnames=list(character(0), c("method", "message"))),
                   nuis.idx = NULL,
                   trafo = list(fct = function(x){
                                      list(fval = x, mat = matrix(1))},
                                mat = matrix(1))
                   ),
         contains = "Estimate")

setClass("LEstimate",
         representation(location = "numeric",
                        ),
         prototype(name = "Loc estimate",
                   estimate = numeric(0),
                   samplesize = numeric(0),
                   completecases = logical(0),
                   asvar = NULL,
                   estimate.call = call("{}"),
                   location = 0,
                   dispersion = 1,
                   Infos = matrix(c(character(0),character(0)), ncol=2,
                                  dimnames=list(character(0), c("method", "message"))),
                   nuis.idx = NULL,
                   trafo = list(fct = function(x){
                                      list(fval = x, mat = matrix(1))},
                                mat = matrix(1))
                   ),
         contains = "Estimate")

setClass("ParetoEstimate", contains="Estimate")
setClass("ParetoDEstimate", contains=c("DEstimate", "ParetoEstimate"))
setClass("ParetoLEstimate", contains=c("LEstimate", "ParetoEstimate"))
setClass("ParetokStepEstimate", contains=c("kStepEstimate", "ParetoEstimate"))
setClass("ParetoMCEstimate", contains=c("MCEstimate", "ParetoEstimate"))


asvarkMADShape <- function(model, k=1){

  xi <- main(model@param)
  beta <- fixed(model@param)

  M <- kMAD(model@distribution, k=k)
  m <- q.l(model)(.5)

  x1.0 <- m - M
  x2.0 <- m + k * M

  ## joint Variance of median and kMAD, see Serfling Mazumder
  dmm <- d(model)(x1.0)
  dmp <- d(model)(x2.0)
  dm <-  d(model)(m)
  alpha <- p(model)(m-M)+p(model)(m+k*M)
  betA <- dmm-dmp
  ceta <- dmm+k*dmp
  eta <- betA^2 + 4*(1-alpha)*betA*dm

  V <- 1/4/ceta^2*(1+eta/dm^2)

   psi_kMad <- function(x){
       cp <- k*dmp+dmm
       cm = dmp-dmm
       return((0.5-((x<=m+k*M)&(x>=m-M)))/cp + cm/cp*((x<=m)-0.5)/dm)
   }

   L2d <- model@L2deriv[[1]]
   L_xi.f = function(x) evalRandVar(L2d,x)[2,]

   E12 <- E(distribution(model),fun=function(x) psi_kMad(x) * L_xi.f(x))

  ASV_Med <- PosSemDefSymmMatrix(V/E12^2)
  dimnames(ASV_Med) <- list(shapename(model),shapename(model))
  return(ASV_Med)
}

asvarMedShape <- function(model){

  m <- q.l(model)(.5)

  dm <-  d(model)(m)

  V <- 1/4/dm^2

   psi_med <- function(x) (0.5-(x<=m))/dm

   L2d <- model@L2deriv[[1]]
   L_xi.f = function(x) evalRandVar(L2d,x)[2,]

   E12 <- E(distribution(model),fun=function(x) psi_med(x) * L_xi.f(x))
   D <- 1/E12

  ASV_Med <- PosSemDefSymmMatrix(V/E12^2)
  dimnames(ASV_Med) <- list(shapename(model),shapename(model))
  return(ASV_Med)
}
