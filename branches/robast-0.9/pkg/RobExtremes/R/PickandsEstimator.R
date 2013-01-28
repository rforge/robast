.PickandsEstimator <- function(x, alpha = 2, GPD.l=TRUE){
 a1 <- 1-1/alpha
 a2 <- 1-1/alpha^2
 if(!GPD.l){
   a1 <- exp(-1/alpha)
   a2 <- exp(-1/alpha^2)
 }
 ms <- quantile(x,c(a1,a2))
 names(ms) <- NULL
 I <- ms[2]
 m <- ms[1]
 xi <- log((I-m)/m)/log(alpha)

 ##############
 ###
 ### the scale estimate we use, i.e. beta <- xi*m/(alpha^xi-1)
 ### differs from the one given in the original reference
 ### ---which was: beta <- xi * m^2 /(I-2*m) leading to xi-dependent bdp---
 ### the one chosen here avoids taking differences I - 2m hence does not
 ### require I > 2m; this leads to (functional) breakdown point
 ###                 min(a1,1-a2,a2-a1)
 ### note that this value is independent of xi !!
 ### for GPD the optimal choice of alpha is 2 leading to bdp 1/4
 ### for GEVD the optimal choice of alpha is 2.248 leading to bdp 0.180
 ###          (and the standard choice alpha = 2 leads to bdp 0.172)
 ### for comparison: with the original scale estimate, at xi=0.7, this
 ###      gives optimal bdp's 0.060 (GEVD) and 0.070 (GPD),
 ###      resp. bdp's 0.048 (GEVD) and 0.064 (GPD) for alpha = 2
 ###
 #############

 beta <- xi*m/(alpha^xi-1)
 ###
 theta <- c(beta,xi)
 names(theta) <- c("scale","shape")
 return(theta)
}

PickandsEstimator <- function(x, ParamFamily=GParetoFamily(), alpha = 2,
                        name, Infos, nuis.idx = NULL,
                        trafo = NULL, fixed = NULL,  na.rm = TRUE,
                        ...){
    force(ParamFamily)
    isGP <- is(ParamFamily,"GParetoFamily")
    if(!(isGP|is(ParamFamily,"GEVFamily")))
         stop("Pickands estimator only available for GPD and GEVD.")
    es.call <- match.call()
    if(missing(alpha)) alpha <- if(isGP) 2 else 2.248
    if(length(alpha)>1 || any(!is.finite(alpha)) || any(alpha<=1))
       stop("'alpha' has to be a numeric > 1 of length 1.")

    if(missing(name))
        name <- "PickandsEstimator"

    asvar.fct.0 <- function(L2Fam, param){
                       asvarPickands(model=L2Fam, alpha = alpha)}
    nuis.idx.0 <- nuis.idx
    trafo.0 <- trafo
    if(is.null(fixed)) fixed <- fixed(ParamFamily)
    fixed.0 <- fixed
    na.rm.0 <- na.rm

    .mPick <- function(x) .PickandsEstimator(x,alpha=alpha, GPD.l=isGP)
    estimate <- Estimator(x, .mPick, name, Infos,
                          asvar.fct = asvar.fct.0, asvar = NULL,
                          nuis.idx = nuis.idx.0, trafo = trafo.0,
                          fixed = fixed.0, na.rm = na.rm.0, ...,
                          ParamFamily = ParamFamily)
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
