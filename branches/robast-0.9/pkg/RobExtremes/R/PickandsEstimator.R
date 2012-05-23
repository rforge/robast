.PickandsEstimator <- function(x, alpha = 2){
 a1 <- 1-1/alpha
 a2 <- 1-1/alpha^2

 ms <- quantile(x,c(a1,a2))
 names(ms) <- NULL
 I <- ms[2]
 m <- ms[1]
 xi <- abs( log((I-m)/m)/log(1-a1))
 beta <- xi*m^2/abs(I-2*m)
 theta <- c(beta,xi)
 names(theta) <- c("scale","shape")
 return(theta)
}

PickandsEstimator <- function(x, alpha = 2, ParamFamily=GParetoFamily(),
                        name, Infos, asvar = NULL, nuis.idx = NULL,
                        trafo = NULL, fixed = NULL, asvar.fct  = NULL, na.rm = TRUE,
                        ...){
    if(!is(ParamFamily,"GParetoFamily"))
         stop("Pickands estimator only available for GPD.")
    name.est <- "PickandsEstimator"
    es.call <- match.call()
    error <- FALSE
    if(missing(alpha)) alpha <- 2
    if(length(alpha)>1 || any(!is.finite(alpha)) || any(alpha<=1))
       stop("'alpha' has to be a numeric > 1 of length 1.")

    if(missing(name))
        name <- "Some estimator"


    asvar.fct.0 <- function(L2Fam=ParamFamily, param){
                       asvarPickands(model=L2Fam, alpha = alpha)}
    asvar <- asvarPickands(model=ParamFamily, alpha = alpha)
    nuis.idx.0 <- nuis.idx
    trafo.0 <- trafo
    fixed.0 <- fixed
    na.rm.0 <- na.rm

    estimate <- Estimator(x, .PickandsEstimator, name, Infos,
                      asvar = asvar, nuis.idx = nuis.idx.0,
                      trafo = trafo.0, fixed = fixed.0,
                      na.rm = na.rm.0, alpha = alpha, ...)

    estimate@untransformed.asvar <- asvar(estimate)
    estimate@asvar <- asvar


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

    estimate@estimate.call <- es.call

    if(missing(Infos))
        Infos <- matrix(c("PickandsEstimator", ""),
                           ncol=2, dimnames=list(character(0), c("method", "message")))
    else{
        Infos <- matrix(c(rep("PickandsEstimator", length(Infos)+1), c("",Infos)),
                          ncol = 2)
        colnames(Infos) <- c("method", "message")
    }
    estimate@Infos <- Infos

    return(estimate)
}
