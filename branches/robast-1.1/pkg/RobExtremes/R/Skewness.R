
#    
setMethod("skewness", signature(x = "Pareto"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp))  
        return(skewness(as(x,"AbscontDistribution"),...))
    else{
         a <- shape(x)
         if(a<=3) return(NA)
         else
         return( 2*(a+1)/(a-3)*sqrt(1-2/a) ) 
    }
})
### source http://mathworld.wolfram.com/ParetoDistribution.html

setMethod("skewness", signature(x = "Gumbel"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp))  
        return(skewness(as(x,"AbscontDistribution"),...))
    else{
         return( -12 * sqrt(6) * APERYCONSTANT / pi^3 ) 
# http://mathworld.wolfram.com/GumbelDistribution.html         
    }
})

setMethod("skewness", signature(x = "GPareto"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp))  
        return(skewness(as(x,"AbscontDistribution"),...))
    else{
         k <- shape(x)
         if(k>=1/3) return(NA)
         else
         return( 2*(1+k)*sqrt(1-2*k)/(1-3*k) ) 
    }
})
### source Maple...

setMethod("skewness", signature(x = "GEV"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp))  
        return(skewness(as(x,"AbscontDistribution"),...))
    else{
         xi <- shape(x)
         if(xi>=1/3) return(NA)
         if(xi==0) return(APERYCONSTANT/pi^3*12*6^.5)
         else
         return((gamma(1-3*xi)-3*gamma(1-xi)*gamma(1-2*xi) + 2*gamma(1-xi)^3)/(gamma(1-2*xi)-gamma(1-xi)^2)^(3/2)) 
    }
})

### source http://en.wikipedia.org/wiki/Generalized_extreme_value_distribution
###        http://en.wikipedia.org/wiki/Gumbel_distribution
###        http://en.wikipedia.org/wiki/Riemann_zeta_function 
