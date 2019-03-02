setMethod("gev.diag",  "gev.fit", function(z) ismev::gev.diag(z))
setMethod("gpd.diag",  "gpd.fit", function(z) ismev::gpd.diag(z))
setMethod("gev.diag",  "GEVEstimate", function(z){
           z0 <- ..gevgpd.diag.fct(z) 
           oK <- ..checkIsNa.gevgpd.diag(z0)
           if(!oK) stop(paste("Not all covariances / se's were specified.\n",
                              "Try using 'returnlevelplot' / 'qqplot' instead.",
                              sep=""))
           ismev::gev.diag(z0)})
setMethod("gpd.diag",  "GPDEstimate", function(z, npy=365){ 
           z0 <- ..gevgpd.diag.fct(z, GPD=TRUE, npy=npy) 
           oK <- ..checkIsNa.gevgpd.diag(z0)
           if(!oK) stop(paste("Not all covariances / se's were specified.\n",
                              "Try using 'returnlevelplot' / 'qqplot' instead.",
                              sep=""))
           ismev::gpd.diag(z0)
           })
setMethod("gev.prof",  "gev.fit", function(z, m, xlow, xup, conf = 0.95, nint = 100)
           ismev::gev.prof(z, m, xlow, xup, conf = conf, nint = nint))
setMethod("gpd.prof",  "gpd.fit", function(z, m, xlow, xup, npy = 365, 
           conf = 0.95, nint = 100) ismev::gpd.prof(z, m, xlow, xup, npy = npy, 
                           conf = conf, nint = nint))
setMethod("gev.prof",  "GEVEstimate", 
           function(z, m, xlow, xup, conf = 0.95, nint = 100){
           z0 <- ..gevgpd.diag.fct(z) 
           oK <- ..checkIsNa.gevgpd.diag(z0)
           if(!oK) stop(paste("Not all covariances / se's were specified.\n",
                              "Profiling is not available for this estimator.",
                              sep=""))           
           ismev::gev.prof(z0, m, xlow, xup, 
                           conf = conf, nint = nint)
           })
setMethod("gpd.prof",  "GPDEstimate", function(z, m, xlow, xup, npy = 365, 
           conf = 0.95, nint = 100){
           z0 <- ..gevgpd.diag.fct(z, GPD=TRUE, npy = npy) 
           oK <- ..checkIsNa.gevgpd.diag(z0)
           if(!oK) stop(paste("Not all covariances / se's were specified.\n",
                              "Profiling is not available for this estimator.",
                              sep=""))           
           ismev::gpd.prof(z0, m, xlow, xup, npy = npy, conf = conf, 
                           nint = nint)
           })
setMethod("gev.profxi",  "gev.fit", function(z, xlow, xup, conf = 0.95, 
           nint = 100) ismev::gev.profxi(z, xlow, xup, conf = conf, 
                          nint = nint))
setMethod("gpd.profxi",  "gpd.fit", function(z,  xlow, xup, conf = 0.95, 
           nint = 100) ismev::gpd.profxi(z,  xlow, xup, conf = conf, 
                          nint = nint))
setMethod("gev.profxi",  "GEVEstimate", function(z, xlow, xup, conf = 0.95, 
           nint = 100){ 
           z0 <- ..gevgpd.diag.fct(z) 
           oK <- ..checkIsNa.gevgpd.diag(z0)
           if(!oK) stop(paste("Not all covariances / se's were specified.\n",
                              "Profiling is not available for this estimator.",
                              sep=""))           
           ismev::gev.profxi(z0, xlow, xup, conf = 0.95, nint = 100)
           })
setMethod("gpd.profxi",  "GPDEstimate", function(z,  xlow, xup, npy=365, 
           conf = 0.95, nint = 100){ 
           z0 <- ..gevgpd.diag.fct(z, GPD=TRUE, npy = npy) 
           oK <- ..checkIsNa.gevgpd.diag(z0)
           if(!oK) stop(paste("Not all covariances / se's were specified.\n",
                              "Profiling is not available for this estimator.",
                              sep=""))           
           ismev::gpd.profxi(z0,  xlow, xup, conf = 0.95, nint = 100)
           })

..gevgpd.diag.fct <- function(z, GPD=FALSE, npy=365){
            utfe <- untransformed.estimate(z)
            if(length(utfe)==2)
               param <- ParamFamParameter(main=utfe, nuisance=nuisance(z),
                               fixed=fixed(z))
            else  param <- ParamFamParameter(main=utfe, nuisance=nuisance(z))
            es.call <- z@estimate.call
            nm.call <- names(es.call)
            if("pIC" %in% names(getSlots(class(z)))){
               PFam0 <- eval(pIC(z)@CallL2Fam)
            }else{
                  PFam <- NULL
                  if("ParamFamily" %in% nm.call)
                     PFam <- eval(as.list(es.call)[["ParamFamily"]])
                  if(is.null(PFam))
                     stop("There is no object of class 'ProbFamily' in the call of 'z'")
                  PFam0 <- modifyModel(PFam, param)
            }
            thresh <- if(GPD)  fixed(param) else NULL
            x <- eval(es.call$x)
            n0 <- length(x)
            z0 <- list()
            if(GPD){ 
               z0$threshold <- thresh
               z0$npy <- npy
               z0$xdata <- x
               x <- x[x > thresh]
            }

            z0$data <- x
            z0$trans <- FALSE
            z0$model <- if(GPD) list(NULL,NULL) else list(NULL,NULL,NULL)
            z0$link <- rep("identity",if(GPD) 2 else 3)
            z0$conv <- 0
            z0$nllh <- -sum(d(PFam0)(x,log=TRUE))
            
            ml0 <- c(untransformed.estimate(z))
            
            if(!GPD){ 
               if(!is.null(fixed(z))) ml0 <- c(fixed(z),ml0)
               nam <- c("loc", "scale", "shape")
             }else  nam <- c("scale", "shape") 

            names(ml0) <- nam 
            z0$mle <- ml0
            
            if (GPD){  
                n <- length(x)
                z0$nexc <- n 
                z0$n <- n0
                z0$rate <- n/n0
            }else n <- n0

            cova <- asvar(z)
            if(is.null(cova)){ 
               cova <- matrix(NA,nrow=length(nam),ncol=length(nam))
               se <- rep(NA,length(nam))
            }else{
               se <- diag(cova)^.5
            }
            dimnames(cova) <- list(nam,nam)          
            z0$cov <- cova/n
            
            names(se) <- nam
            z0$se <- se

            if(GPD){
               vals <-  cbind(t(matrix(z0$mle,2,n)),rep(thresh,n))
               dimnames(vals) <- list(NULL, c("","","u"))
               z0$vals <- vals
            }else
               z0$vals <- t(matrix(z0$mle,3,n))
            if(GPD) z0 <- z0[c(5:7,1,11,4,8,9,16,10,13:15,12,2,3)]
            else z0 <- z0[c(2:6, 1,7:10)]
            class(z0) <- if(GPD) "gpd.fit" else "gev.fit"
            return(z0)
            #ismev::gev.fit(z0)
}

..checkIsNa.gevgpd.diag <- function(x){
  oK <- TRUE
  if(any(is.na(x$se))) oK <- FALSE
  if(any(is.na(x$cov))) oK <- FALSE
  return(oK)
}
