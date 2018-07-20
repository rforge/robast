setMethod("Qn", signature(x = "ANY"),
    function(x, ...){
        dots <- list(...)
        constant <- ifelse(hasArg(constant), dots$"constant", 2.21914)
        finite.corr <- ifelse(hasArg(finite.corr), dots$"finite.corr",
                              !hasArg(constant))
        if(!is.null(dots$"na.rm")) if(dots$"na.rm") x <- x[!is.na(x)]
        robustbase::Qn(x, constant=constant, finite.corr=finite.corr)
    })

setMethod("Sn", signature(x = "ANY"),
    function(x, ...){
        dots <- list(...)
        constant <- ifelse(hasArg(constant), dots$"constant", 1.1926)
        finite.corr <- ifelse(hasArg(finite.corr), dots$"finite.corr",
                              !hasArg(constant))
        if(!is.null(dots$"na.rm")) if(dots$"na.rm") x <- x[!is.na(x)]
        robustbase::Sn(x, constant=constant, finite.corr=finite.corr)
    })

setMethod("Qn", signature(x = "UnivariateDistribution"),
    function(x, q00 = NULL,  ...){
         if(is.null(q00)) q00 <- 10*q.l(x)(3/4)

         intv <- function(xx,q0=q00) p(x)(q0+q.l(x)(xx))-5/8
         intq <- function(q){
                    sapply(q, function(q1){
                                 integrate(intv,lower=0,upper=1,q0=q1)$value})
                 }
         qc = try(uniroot(intq,lower=-q00,upper=q00)$root,silent=TRUE)
         if(is(qc,"try-error")) {
               print("error")
            return(NA)
         }else  return(qc)
    })

setMethod("Qn", signature(x = "DiscreteDistribution"),
    function(x,  ...){
         x2 <- x-x
         q.l(x2)(5/8)
    })


setMethod("Sn", signature(x = "UnivariateDistribution"),
    function(x, low=0,upp=1.01, accuracy = 1000, ...){
          m0 <- median(x)
          M0 <- mad(x)
          g <- function(xx){
               fct <- function(m) p(x)(m+xx)-p(x)(-m+xx)-0.5
               up0 <- upp*(M0+abs(m0-xx))
               m <- try(uniroot(fct, lower = low,
                        upper = up0,
                        f.lower=if(low<1e-12) -0.5 else fct(low),
                        f.upper=max(fct(up0),1e-8))$root,
                        silent = TRUE)
               if(is(m,"try-error")) {
#                        print("error")
                        return(NA)
               }else{   return(m)    }
          }

          x0 <- q.l(x)(seq(.5/accuracy,1-.5/accuracy,length=accuracy))
          y  <- sapply(x0,g)
          c0 <- median(y,na.rm=TRUE)
          return(c0)
    })

setMethod("Sn", signature(x = "DiscreteDistribution"),
    function(x, ...){

          g <- function(xx){
               median(abs(x-xx))
          }

          pr <- prob(x)
          sx <- support(x)
          y  <- sapply(sx,g)
          o <- order(y)
          yo <- y[o]
          pro <- pr[o]
          cpro <- cumsum(pro)
          ws <- min(sum(which(cpro<0.5))+1,length(o))
          return(yo[ws])
    })


setMethod("Sn", signature(x = "Norm"),
    function(x, ...){
           return(sd(x)*  0.838504603)
    })

setMethod("Qn", signature(x = "Norm"),
    function(x, ...){
           return(sd(x)*0.45062411)
    })

setMethod("Sn", signature(x = "AffLinDistribution"),
    function(x, ...){
           return(abs(x@a) * Sn(x@X0,...))
    })

setMethod("Qn", signature(x = "AffLinDistribution"),
    function(x, ...){
           return(abs(x@a) * Qn(x@X0,...))
    })

.Sn.intp <- function(x, nam){
    if(abs(scale(x)-1)< 1e-12){
       famg <- try(getFromNamespace(nam, ns = "RobAStRDA"), silent =TRUE)
#       sng <- try(getFromNamespace(".Sn", ns = "RobAStRDA"), silent =TRUE)
       if(is(famg,"try-error")) return(Sn(as(x,"AbscontDistribution")))
#       if(!nam %in% names(sng)) return(Sn(as(x,"AbscontDistribution")))
#       if(!.versionSuff("fun")=="fun.N") return(Sn(as(x,"AbscontDistribution")))
       snf <- famg[["Sn"]][[.versionSuff("fun")]]
       ret <- snf(shape(x))
       if(is.na(ret)) return(Sn(as(x,"AbscontDistribution")))
    }else ret <- scale(x)*Sn(x=x/scale(x))
    return(ret)
}

setMethod("Sn", signature(x = "GPareto"),
    function(x, ...).Sn.intp(x,".GPareto") )

setMethod("Sn", signature(x = "GEV"),
    function(x, ...).Sn.intp(x,".GEV") )

setMethod("Sn", signature(x = "Gammad"),
    function(x, ...).Sn.intp(x,".Gamma") )

setMethod("Sn", signature(x = "Weibull"),
    function(x, ...).Sn.intp(x,".Weibull") )

setMethod("Sn", signature(x = "Pareto"),
    function(x, ...).Sn.intp(x,".Pareto") )
