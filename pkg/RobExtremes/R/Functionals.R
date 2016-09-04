
setMethod("var", signature(x = "Pareto"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
        return(var(as(x,"AbscontDistribution"),...))
    else{ a <- shape(x); b <- Min(x)
        if(a<=2) return(NA)
        return(b^2 * a/(a-1)^2/(a-2))
    }})
### source http://mathworld.wolfram.com/ParetoDistribution.html

setMethod("var", signature(x = "Gumbel"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
        return(var(as(x,"AbscontDistribution"),...))
    else{  b <- scale(x)
            return(b^2 * pi^2/6)
    }})
## http://mathworld.wolfram.com/GumbelDistribution.html

setMethod("var", signature(x = "GPareto"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
        return(var(as(x,"AbscontDistribution"),...))
    else{ k <- shape(x); s <- scale(x)
        if(k>=1/2) return(NA)
        return(s^2/(1-k)^2/(1-2*k))
    }})
### source http://en.wikipedia.org/wiki/Pareto_distribution


setMethod("var", signature(x = "GEV"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
        return(var(as(x,"AbscontDistribution"),...))
    else{ xi <- shape(x); sigma <- scale(x)
        if(xi>=1/2) return(NA)
        if(xi==0) return(sigma^2*pi^2/6)
        if((xi!=0)&&(xi<1/2))return(sigma^2*(gamma(1-2*xi)-gamma(1-xi)^2)/xi^2)
    }})
### http://en.wikipedia.org/wiki/Generalized_extreme_value_distribution

#################################################################
# some exact medians
#################################################################


setMethod("median", signature(x = "Pareto"),
    function(x) {a <- shape(x); b<- Min(x)
              return(b*2^(1/a))
    })
setMethod("median", signature(x = "Gumbel"),
    function(x) {a <- loc(x); b <- scale(x)
              return(a - b *log(log(2)))
    })
setMethod("median", signature(x = "GPareto"),
    function(x) {k <- shape(x); mu <- loc(x); s <- scale(x)
              return(mu + s*(2^k-1)/k)
    })
setMethod("median", signature(x = "GEV"),
    function(x) {xi <- shape(x); mu <- loc(x); sigma <- scale(x)
              if (xi != 0) return(mu + sigma*(log(2)^(-xi)-1)/xi)
              else return(mu-sigma*log(log(2)))
    })

#################################################################
# some exact IQRs
#################################################################


setMethod("IQR", signature(x = "Pareto"),
    function(x) {a <- shape(x); b<- Min(x)
              return(b*(4^(1/a)-(4/3)^(1/a)))
    })
setMethod("IQR", signature(x = "Gumbel"),
    function(x) { b <- scale(x)
              return(b * (log(log(4))-log(log(4/3))))
    })
setMethod("IQR", signature(x = "GPareto"),
    function(x) {k <- shape(x); s<- scale(x)
              return(s/k*4^k*(1-3^(-k)))
    })
setMethod("IQR", signature(x = "GEV"),
    function(x) {xi <- shape(x); sigma<- scale(x)
             if (xi != 0) return(sigma*((log(4/3))^(-xi)-(log(4))^(-xi))/xi)
             else return(sigma*(log(log(4))-log(log(4/3))))
    })
