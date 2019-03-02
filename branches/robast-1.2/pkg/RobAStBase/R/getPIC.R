setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="ANY",estimator="MCEstimate"),
              function(PFam, estimator)estimator)


#setMethod(".checkEstClassForParamFamily",
#              signature=signature(PFam="ANY",estimator="MCEstimate"),
#              function(PFam, estimator) .extendbyPIC(PFam, estimator, "MCALEstimate"))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="ANY",estimator="MLEstimate"),
              function(PFam, estimator) .extendbyPIC(PFam, estimator, "ML.ALEstimate"))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="ANY",estimator="CvMMDEstimate"),
              function(PFam, estimator) .extendbyPIC(PFam, estimator, "CvMMD.ALEstimate"))

.extendbyPIC <- function(PFam, estimator, toClass){
                 fromSlotNames <- slotNames(class(estimator))
                 to <- new(toClass)
                 for(item in fromSlotNames) slot(to, item) <- slot(estimator,item)
                 to@pIC <- substitute(getPIC(estimator0), list(estimator0=estimator))
                 to
              }

.getPIC <- function(object){
       if(is.null(object@pIC)) return(NULL)
       pIC0 <- object@pIC
       if(is(pIC0, "InfluenceCurve")) return(pIC0)
       if(is.call(pIC0)) pIC0 <- eval(pIC0)
       return(pIC0)
}

.getL2Fam <- function(estimator){
       ecl <- as.list(estimator@estimate.call)[-1]
       L2Fam0 <- eval(ecl[["ParamFamily"]])
       param.0 <- param(L2Fam0)
       theta <- untransformed.estimate(estimator)
       idx <- idx.m <- seq(length(theta))
       if(!is.null(nuisance(param.0))){
          lnx <- length(nuisance(param.0))
          idx.n <- rev(rev(idx)[1:lnx])
          idx.m <- idx[-idx.n]
          th.nuis <- theta[idx.n]
          names(th.nuis) <- names(nuisance(param.0))
          param.0@nuisance <- th.nuis
       }
       th.main <- theta[idx.m]
       names(th.main)<-  names(main(param.0))
       param.0@main <- th.main
       param.0@trafo <- trafo(estimator)$mat
       L2Fam <- modifyModel(L2Fam0, param.0)
       return(L2Fam)
}


setMethod("getPIC","ANY", function(estimator)NULL)

setMethod("getPIC","MLEstimate", function(estimator){
       L2Fam <- .getL2Fam(estimator)
       pIC <- optIC(L2Fam, risk=asCov())
       return(pIC)
    })

setMethod("getPIC","CvMMDEstimate", function(estimator){
       L2Fam <- .getL2Fam(estimator)
       param.0 <- param(L2Fam)
       ecl <- as.list(estimator@estimate.call)[-1]
#       print(system.time({
       if(grepl("mu = model distr",name(estimator))){
          res <- .CvMMDCovariance(L2Fam=L2Fam, param=param.0,withpreIC=TRUE, N = 2000)
       }else{
          if(grepl("mu = emp\\. cdf",name(estimator))){
             x <- eval(ecl$x)
             res <- .CvMMDCovarianceWithMux(L2Fam = L2Fam, param=param.0,x=x,withpreIC=TRUE, N = 2000)
          }else{
             mu <- eval(ecl$mu)
             res <- .CvMMDCovariance(L2Fam=L2Fam, param=param.0,x=x,withpreIC=TRUE, mu=mu, N = 2000)
          }
       }
#       }))
       ICCurve <- res$preIC
       ICname <- "IC of CvM MDE"
       ICCallL2Fam <- L2Fam@fam.call
       ICRisks <- list(asCov = estimator@asvar)
       ICInfos = matrix(c("pIC-CvM-MDE","computed by .CvMMDCovariance[WithMux]"), ncol=2,
                                dimnames=list(character(0), c("method", "message")))
       pIC <- IC(name = ICname, Curve = ICCurve, Risks=ICRisks,
                 Infos = ICInfos, CallL2Fam = ICCallL2Fam, modifyIC = NULL)
       return(pIC)
    })
