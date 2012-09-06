.prepend <- function(prep, list0, dots = NULL){
   if(length(list0)+length(dots)==0) return(list(prep))
   n <- length(list0) + 1
   list1 <- vector("list",n)
   list1[[1]] <- prep
   names(list1)[1] <- "x"
   if(n>1) for(i in 2:n) {
           list1[[i]] <- list0[[i-1]]
           names(list1)[i] <- names(list0)[i-1]}
   ldots <- length(dots)
   l1 <- length(list1)
   if(ldots) {
       for( i in 1:ldots){
            list1[[l1+i]] <- dots[[i]]
            names(list1)[l1+i] <- names(dots)[i]
       }
   }
   return(list1)
}

.LDMatch <- function(x.0, loc.est.0,disp.est.0,
                          loc.fctal.0, disp.fctal.0, ParamFamily.0,
                        loc.est.ctrl.0 = NULL, loc.fctal.ctrl.0=NULL,
                        disp.est.ctrl.0 = NULL, disp.fctal.ctrl.0=NULL,
                        q.lo.0 =0, q.up.0=Inf, log.q.0 =TRUE, ...
                        ){
    dots <- list(...)
    loc.emp <- do.call(loc.est.0, args = .prepend(x.0,loc.est.ctrl.0, dots))
    disp.emp <- do.call(disp.est.0, args = .prepend(x.0,disp.est.ctrl.0, dots))
    q.emp <- if(log.q.0) log(loc.emp)-log(disp.emp) else loc.emp/disp.emp
    q.f <- function(xi){
       distr.new <- ParamFamily.0@modifyParam(theta=c("scale"=1,"shape"=xi))
       loc.th <- do.call(loc.fctal.0, args = .prepend(distr.new,loc.fctal.ctrl.0, dots))
       sc.th <- do.call(disp.fctal.0, args = .prepend(distr.new,disp.fctal.ctrl.0, dots))
       val <- if(log.q.0) log(loc.th)-log(sc.th) - q.emp else
                        loc.th/sc.th-q.emp
       return(val)
    }
    xi.0 <- uniroot(q.f,lower=q.lo.0,upper=q.up.0)$root
    distr.new.0 <- ParamFamily.0@modifyParam(theta=c("scale"=1,"shape"=xi.0))
    m1xi <- do.call(loc.fctal.0, args = .prepend(distr.new.0,loc.fctal.ctrl.0, dots))
    val <-   c("shape"=xi.0,"scale"=loc.emp/m1xi, "loc"=loc.emp,"disp"=disp.emp)
    return(val)
}

LDEstimator <- function(x, loc.est, disp.est,
                        loc.fctal, disp.fctal, ParamFamily,
                        loc.est.ctrl = NULL, loc.fctal.ctrl=NULL,
                        disp.est.ctrl = NULL, disp.fctal.ctrl=NULL,
                        q.lo =1e-3, q.up=15, log.q =TRUE,
                        name, Infos, asvar = NULL, nuis.idx = NULL,
                        trafo = NULL, fixed = NULL, asvar.fct  = NULL, na.rm = TRUE,
                        ...){
    param0 <- main(param(ParamFamily))
    if(!all(c("shape","scale") %in% names(param0)))
        stop("LDEstimators expect shape-scale models.")
    name.est <- "LDEstimator"
    es.call <- match.call()
    if(missing(name))
        name <- "Some estimator"
    LDnames <- paste("Location:",
                     paste(deparse(substitute(loc.fctal))),
                     " ","Dispersion:",
                     paste(deparse(substitute(disp.fctal))))

    LDMval <- NULL
    estimator <- function(x,...){
         LDMval <<- .LDMatch(x.0= x,
                         loc.est.0 = loc.est, 
                         disp.est.0 =  disp.est,
                         loc.fctal.0 = loc.fctal, 
                         disp.fctal.0 =  disp.fctal,
                         ParamFamily.0 = ParamFamily,
                         loc.est.ctrl.0 = loc.est.ctrl,
                         loc.fctal.ctrl.0 = loc.fctal.ctrl,
                         disp.est.ctrl.0 = disp.est.ctrl,
                         disp.fctal.ctrl.0 = disp.fctal.ctrl,
                         q.lo.0 = q.lo, 
                         q.up.0 = q.up, 
                         log.q.0 = log.q)
         return(LDMval[1:2])
    }

    asvar.0 <- asvar
    nuis.idx.0 <- nuis.idx
    trafo.0 <- trafo
    fixed.0 <- fixed
    na.rm.0 <- na.rm

    print(nuis.idx.0)
    print(trafo.0)
    print(fixed.0)
    

    estimate <- Estimator(x, estimator, name, Infos,
                      asvar = asvar.0, nuis.idx = nuis.idx.0,
                      trafo = trafo.0, fixed = fixed.0,
                      na.rm = na.rm.0, ...)

    print(estimate)
    #print(estimate@untransformed.estimate)
    print(estimate@untransformed.asvar)
    cat("\n asvar",estimate@asvar,"\n")


    if(missing(asvar)) asvar <- NULL

    if((is.null(asvar))&&(!missing(asvar.fct))&&(!is.null(asvar.fct)))
          asvar <- asvar.fct(ParamFamily, estimate, ...)

    estimate@untransformed.asvar <- asvar
    #print(estimate)

    l.e <- length(estimate@untransformed.estimate)
    idx <- NULL
    idm <- 1:l.e
    if(!is.null(nuis.idx))
        {idx <- nuis.idx
         idm <- idm[-idx]
         mat <- diag(length(idm))}

    if(!.isUnitMatrix(estimate@trafo$mat)){
       estimate@estimate <- estimate@trafo$fct(estimate)
       
       if(!is.null(asvar))
           estimate@asvar <- estimate@trafo$mat%*%asvar[idm,idm]%*%t(estimate@trafo$mat)
    }

     print(estimate@asvar)

    estimate@estimate.call <- es.call

    if(missing(Infos))
        Infos <- matrix(c("LDEstimator", LDnames),
                           ncol=2, dimnames=list(character(0), c("method", "message")))
    else{
        Infos <- matrix(c(rep("LDEstimator", length(Infos)+1), c(LDnames,Infos)),
                          ncol = 2)
        colnames(Infos) <- c("method", "message")
    }
    estimate@Infos <- Infos

    estim <- new("LDEstimate")

    sln <- names(getSlots(class(estimate)))
    for( i in 1:length(sln))
        slot(estim, sln[i]) <- slot(estimate, sln[i])
    rm(estimate)
    
    estim@dispersion <- LDMval["disp"]
    estim@location <- LDMval["loc"]

    return(estim)
}


medkMAD <- function(x, k=1, ParamFamily, q.lo =1e-3, q.up=15, nuis.idx = NULL,
                        trafo = NULL, fixed = NULL, asvar.fct = NULL, na.rm = TRUE,
                        ...){
      es.call <- match.call()
      if(missing(k)) k <- 1

      if (is.null(asvar.fct)){asvar.fct <- asvarMedkMAD
                              asvar <- asvarMedkMAD(ParamFamily, k=k)}


      es <- LDEstimator(x, loc.est = median, disp.est = kMAD,
                     loc.fctal = median, disp.fctal = kMAD,
                     ParamFamily = ParamFamily,
                     loc.est.ctrl = list(na.rm = na.rm), loc.fctal.ctrl = NULL,
                     disp.est.ctrl = list(k=k, na.rm = na.rm),
                     disp.fctal.ctrl=list(k=k),
                     q.lo =q.lo, q.up=q.up, log.q=TRUE,
                     name = "medkMAD", Infos="medkMAD",
                     asvar = asvar, nuis.idx = nuis.idx, trafo = trafo, fixed = fixed,
                     asvar.fct = asvar.fct, na.rm = na.rm, ...)
      es@estimate.call <- es.call
     
      return(es)
                     }
                        
medQn <- function(x,  ParamFamily, q.lo =1e-3, q.up=15, nuis.idx = NULL,
                        trafo = NULL, fixed = NULL, asvar.fct = NULL, na.rm = TRUE,
                        ...){
    es.call <- match.call()
    es <- LDEstimator(x, loc.est = median, disp.est = Qn,
                     loc.fctal = median, disp.fctal = Qn,
                     ParamFamily = ParamFamily,
                     loc.est.ctrl = list(na.rm = na.rm), loc.fctal.ctrl = NULL,
                     disp.est.ctrl = list(constant=1,na.rm = na.rm),
                     disp.fctal.ctrl = NULL,
                     q.lo =q.lo, q.up=q.up, log.q=TRUE,
                     name = "medQn", Infos="medQn",
                     asvar = NULL, nuis.idx = nuis.idx, trafo = trafo, fixed = fixed,
                     asvar.fct = asvar.fct, na.rm = na.rm, ...)
      es@estimate.call <- es.call
      return(es)
                     }

medSn <- function(x, ParamFamily, q.lo =1e-3, q.up=10, nuis.idx  = NULL,
                        trafo = NULL, fixed = NULL, asvar.fct = NULL, na.rm = TRUE,
                        accuracy = 100, ...){
      es.call <- match.call()
      es <- LDEstimator(x, loc.est = median, disp.est = Sn,
                     loc.fctal = median, disp.fctal = Sn,
                     ParamFamily = ParamFamily,
                     loc.est.ctrl = list(na.rm = na.rm), loc.fctal.ctrl = NULL,
                     disp.est.ctrl = list(constant=1,na.rm = na.rm),
                     disp.fctal.ctrl = list(accuracy=accuracy),
                     q.lo =q.lo, q.up=q.up, log.q=TRUE,
                     name = "medSn", Infos="medSn",
                     asvar = NULL, nuis.idx = nuis.idx, trafo = trafo, fixed = fixed,
                     asvar.fct = asvar.fct, na.rm = na.rm, ...)
      es@estimate.call <- es.call
      return(es)
      }

medkMADhybr <- function(x, k=1, ParamFamily, q.lo =1e-3, q.up=15,
                        KK=20, nuis.idx = NULL,
                        trafo = NULL, fixed = NULL, asvar.fct = NULL, na.rm = TRUE,
                        ...){
 i <- 1
 es <- try(medkMAD(x, k = k, ParamFamily = ParamFamily,
                            q.lo = q.lo, q.up = q.up,
                            nuis.idx = nuis.idx, trafo = trafo,
                            fixed = fixed, asvar.fct = asvar.fct, na.rm = na.rm,
                             ...), silent=TRUE)
 if(! any(is.na(es)) && !is(es,"try-error"))
   {return(es)}

 k1 <- 3.23
 while(i<KK){
      i <- i + 1
      es <- try(medkMAD(x, k = k1, ParamFamily = ParamFamily,
                            q.lo = q.lo, q.up = q.up,
                            nuis.idx = nuis.idx, trafo = trafo,
                            fixed = fixed, asvar.fct = asvar.fct, na.rm = na.rm,
                             ...), silent=TRUE)
      k1 <- k1 * 3
      if(! any(is.na(es)) && !is(es,"try-error"))
         {return(es)}
      }
 return(c("scale"=NA,"shape"=NA))
}

setMethod("location", "LDEstimate", function(object) object@location)
setMethod("dispersion", "LDEstimate", function(object) object@dispersion)