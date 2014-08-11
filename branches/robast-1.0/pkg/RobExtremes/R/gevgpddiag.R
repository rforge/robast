setMethod("gev.diag",  "gev.fit", function(z) ismev::gev.diag(z))
setMethod("gpd.diag",  "gpd.fit", function(z) ismev::gpd.diag(z))
setMethod("gev.diag",  "GEVEstimate", function(z)ismev::gev.diag(..gevgpd.diag.fct(z)))
setMethod("gpd.diag",  "GPDEstimate", function(z)ismev::gev.diag(..gevgpd.diag.fct(z)))
setMethod("gev.prof",  "gev.fit", function(z, m, xlow, xup, conf = 0.95, nint = 100)
           ismev::gev.prof(z, m, xlow, xup, conf = 0.95, nint = 100))
setMethod("gpd.prof",  "gpd.fit", function(z, m, xlow, xup, npy = 365, conf = 0.95, nint = 100)
           ismev::gpd.prof(z, m, xlow, xup, npy, conf, nint))
setMethod("gev.prof",  "GEVEstimate", function(z, m, xlow, xup, conf = 0.95, nint = 100)
           ismev::gev.prof(..gevgpd.diag.fct(z), m, xlow, xup, conf = 0.95, nint = 100))
setMethod("gpd.prof",  "GPDEstimate", function(z, m, xlow, xup, npy = 365, conf = 0.95, nint = 100)
           ismev::gev.prof(..gevgpd.diag.fct(z), m, xlow, xup, npy = 365, conf = 0.95, nint = 100))
setMethod("gev.profxi",  "gev.fit", function(z, xlow, xup, conf = 0.95, nint = 100)
           ismev::gev.profxi(z, xlow, xup, conf = 0.95, nint = 100))
setMethod("gpd.profxi",  "gpd.fit", function(z,  xlow, xup, conf = 0.95, nint = 100)
           ismev::gpd.profxi(z,  xlow, xup, conf = 0.95, nint = 100))
setMethod("gev.profxi",  "GEVEstimate", function(z, xlow, xup, conf = 0.95, nint = 100)
           ismev::gev.profxi(..gevgpd.diag.fct(z), xlow, xup, conf = 0.95, nint = 100))
setMethod("gpd.profxi",  "GPDEstimate", function(z,  xlow, xup, conf = 0.95, nint = 100)
           ismev::gev.profxi(..gevgpd.diag.fct(z),  xlow, xup, conf = 0.95, nint = 100))

..gevgpd.diag.fct <- function(z){
            utfe <- untransformed.estimate(z)
            if(length(utfe)==2)
               param <- ParamFamParameter(main=utfe, nuisance=nuisance(z),
                               fixed=fixed(z))
            else  param <- ParamFamParameter(main=utfe, nuisance=nuisance(z))
            es.call <- z@estimate.call
            nm.call <- names(es.call)
            PFam <- NULL
            if("ParamFamily" %in% nm.call)
                PFam <- eval(as.list(es.call)[["ParamFamily"]])
            if(is.null(PFam))
               stop("There is no object of class 'ProbFamily' in the call of 'z'")
            PFam0 <- modifyModel(PFam, param)
            x <- eval(es.call$x)
            z0 <- list()
            z0$trans <- FALSE
            z0$model <- list(NULL,NULL,NULL)
            z0$link <- rep("identity",3)
            z0$conv <- 0
            z0$nllh <- sum(d(PFam0)(x,log=TRUE))
            z0$data <- x
            ml0 <- c(untransformed.estimate(z))
            if(!is.null(fixed(z))) ml0 <- c(fixed(z),ml0)
            z0$mle <- ml0
            z0$cov <- asvar(z)/length(x)
            z0$se <- diag(z0$cov)^.5
            z0$vals <- t(matrix(z0$mle,3,length(x)))
            return(z0)
            #ismev::gev.fit(z0)
}
