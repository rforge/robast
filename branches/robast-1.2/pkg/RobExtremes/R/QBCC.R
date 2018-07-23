#####################################################################
# Robust starting estimator for Weibull model
# acc. Boudt Caliskan Croux (2011)
# for p1 < p2
# has bdp min(p1,1-p2,p2-p1), maximal for p1=p2=1/3
#
#####################################################################
.QBCC <- function(x, p1 = 1/3, p2 = 2/3){
 if(p1>=p2) {p<-p1; p1 <- p2; p2 <- p}
 l1 <- -log(1-p1); l2 <- -log(1-p2)
 ms <- quantile(x,c(p1,p2))
 names(ms) <- NULL
 Q2 <- ms[2]
 Q1 <- ms[1]
 xi <- (log(l2)-log(l1))/(log(Q2)-log(Q1))
 beta <- Q2*l2^(-1/xi)
 ###
 theta <- c(beta,xi)
 names(theta) <- c("scale","shape")
 return(theta)
}

QuantileBCCEstimator <- function(x, p1=1/3, p2=2/3,
                        name, Infos, nuis.idx = NULL,
                        trafo = NULL, fixed = NULL,  na.rm = TRUE,
                        ...){
    es.call <- match.call()
    force(p1); force(p2)
    if(length(p1)>1 || any(!is.finite(p1)) || p1<=0 || p1>=1)
       stop("'p1' has to be in [0,1] and of length 1.")
    if(length(p2)>1 || any(!is.finite(p2)) || p2<=0 || p2>=1 || abs(p1-p2)< 1e-8)
       stop("'p2' has to be in [0,1] and of length 1 and distinct of 'p1'.")

    if(missing(name))
        name <- "QuantileBCCEstimator"

    ParamFamily <- WeibullFamily()
    asvar.fct.0 <- function(L2Fam=ParamFamily, param){
                       asvarQBCC(model=L2Fam, p1 = p1, p2 = p2)}
    nuis.idx.0 <- nuis.idx
    trafo.0 <- trafo
    if(is.null(fixed)&!is.null( fixed(ParamFamily))) fixed <- fixed(ParamFamily)
    fixed.0 <- fixed
    na.rm.0 <- na.rm

    .mQBCC <- function(x) .QBCC(x,p1=p1,p2=p2)
    estimate <- Estimator(x, .mQBCC, name, Infos,
                          asvar.fct = asvar.fct.0, asvar = NULL,
                          nuis.idx = nuis.idx.0, trafo = trafo.0,
                          fixed = fixed.0, na.rm = na.rm.0, ...,
                          ParamFamily = ParamFamily)
    estimate@estimate.call <- es.call

    if(missing(Infos))
        Infos <- matrix(c(name, ""),
                           ncol=2, dimnames=list(character(0), c("method", "message")))
    else{
        Infos <- matrix(c(rep(name, length(Infos)+1), c("",Infos)),
                          ncol = 2)
        colnames(Infos) <- c("method", "message")
    }
    estimate@Infos <- Infos

    return(estimate)
}
