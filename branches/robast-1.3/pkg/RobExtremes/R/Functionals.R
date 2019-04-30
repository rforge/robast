
setMethod("var", signature(x = "Pareto"),
    function(x, propagate.names=getdistrExOption("propagate.names.functionals"), ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
        return(var(as(x,"AbscontDistribution"),...))
    else{a <- shape(x); b <- Min(x)
        if(a<=2) return(NA)
        ret.v <- b^2 * a/(a-1)^2/(a-2)
        if(!propagate.names){names(ret.v) <- NULL}
        return(ret.v)
    }})
### source http://mathworld.wolfram.com/ParetoDistribution.html

setMethod("var", signature(x = "Gumbel"),
    function(x, propagate.names=getdistrExOption("propagate.names.functionals"), ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
        return(var(as(x,"AbscontDistribution"),...))
    else{  b <- scale(x)
           ret.v <- (b^2 * pi^2/6)
           if(!propagate.names){names(ret.v) <- NULL}
           return(ret.v)
    }})
## http://mathworld.wolfram.com/GumbelDistribution.html

setMethod("var", signature(x = "GPareto"),
    function(x, propagate.names=getdistrExOption("propagate.names.functionals"), ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
        return(var(as(x,"AbscontDistribution"),...))
    else{ k <- shape(x); s <- scale(x)
        if(k>=1/2) return(NA)
        ret.v <- s^2/(1-k)^2/(1-2*k)
        if(!propagate.names){names(ret.v) <- NULL}
        return(ret.v)
    }})
### source http://en.wikipedia.org/wiki/Pareto_distribution


setMethod("var", signature(x = "GEV"),
    function(x, propagate.names=getdistrExOption("propagate.names.functionals"), ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
        return(var(as(x,"AbscontDistribution"),...))
    else{ xi <- shape(x); sigma <- scale(x)
        if(xi>=1/2) return(NA)
        if(xi==0) ret.v <- sigma^2*pi^2/6
        if((xi!=0)&&(xi<1/2)) ret.v <- sigma^2*(gamma(1-2*xi)-gamma(1-xi)^2)/xi^2
        if(!propagate.names){names(ret.v) <- NULL}
        return(ret.v)
    }})
### http://en.wikipedia.org/wiki/Generalized_extreme_value_distribution

#################################################################
# some exact medians
#################################################################


setMethod("median", signature(x = "Pareto"),
    function(x, propagate.names=getdistrExOption("propagate.names.functionals")) {
    a <- shape(x); b<- Min(x)
    ret.v <- b*2^(1/a)
    if(!propagate.names){names(ret.v) <- NULL}
    return(ret.v)
    })
setMethod("median", signature(x = "Gumbel"),
    function(x, propagate.names=getdistrExOption("propagate.names.functionals")) {
    a <- loc(x); b <- scale(x)
    ret.v <- a - b *log(log(2))
    if(!propagate.names){names(ret.v) <- NULL}
    return(ret.v)
    })
setMethod("median", signature(x = "GPareto"),
    function(x, propagate.names=getdistrExOption("propagate.names.functionals")) {
    k <- shape(x); mu <- loc(x); s <- scale(x)
    ret.v <- mu + s*(2^k-1)/k
    if(!propagate.names){names(ret.v) <- NULL}
    return(ret.v)
    })
setMethod("median", signature(x = "GEV"),
    function(x, propagate.names=getdistrExOption("propagate.names.functionals")) {
    xi <- shape(x); mu <- loc(x); sigma <- scale(x)
    if (xi != 0) ret.v <- (mu + sigma*(log(2)^(-xi)-1)/xi)
       else ret.v <- (mu-sigma*log(log(2)))
    if(!propagate.names){names(ret.v) <- NULL}
    return(ret.v)
    })

#################################################################
# some exact IQRs
#################################################################


setMethod("IQR", signature(x = "Pareto"),
    function(x, propagate.names=getdistrExOption("propagate.names.functionals")) {
    a <- shape(x); b<- Min(x)
    ret.v <- (b*(4^(1/a)-(4/3)^(1/a)))
    if(!propagate.names){names(ret.v) <- NULL}
    return(ret.v)
    })
setMethod("IQR", signature(x = "Gumbel"),
    function(x, propagate.names=getdistrExOption("propagate.names.functionals")) {
    b <- scale(x)
    ret.v <- (b * (log(log(4))-log(log(4/3))))
    if(!propagate.names){names(ret.v) <- NULL}
    return(ret.v)
    })
setMethod("IQR", signature(x = "GPareto"),
    function(x, propagate.names=getdistrExOption("propagate.names.functionals")) {
    k <- shape(x); s<- scale(x)
    ret.v <- (s/k*4^k*(1-3^(-k)))
    if(!propagate.names){names(ret.v) <- NULL}
    return(ret.v)
    })
setMethod("IQR", signature(x = "GEV"),
    function(x, propagate.names=getdistrExOption("propagate.names.functionals")) {
    xi <- shape(x); sigma<- scale(x)
    if (xi != 0) ret.v <- (sigma*((log(4/3))^(-xi)-(log(4))^(-xi))/xi)
            else ret.v <- (sigma*(log(log(4))-log(log(4/3))))
    if(!propagate.names){names(ret.v) <- NULL}
    return(ret.v)
    })
