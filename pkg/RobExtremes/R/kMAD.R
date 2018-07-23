############################################################
##                                                       ###
## Median of absolute deviations for skewed distribution ###
##                                                       ###  
############################################################


setMethod("kMAD", signature(x = "numeric", k = "numeric"),
    function(x, k=1, na.rm = TRUE, eps = .Machine$double.eps, ...){
       if(na.rm) x <- x[!is.na(x)]
       if(! length(k)==1) stop ("k has to be a numeric of length 1")
       if(k<=0) stop ("k has to be strictly positive")
       eps1 <- min(diff(unique(sort(x))))          
       erg  <- .C(C_kMad, as.double(x),
                as.integer(length(x)),
                as.integer(k),
                d = double(1),
                eps = as.double(min(eps1,eps)))
      return(erg$d)
    })

setMethod("kMAD", signature(x = "UnivariateDistribution", k = "numeric"),
    function(x, k=1, up = NULL, ...){
       if(! length(k)==1) stop ("k has to be a numeric of length 1")
       if(k<=0) stop ("k has to be strictly positive")
       m <- median(x)
       if(is.null(up)) up <- min(3*IQR(x),q.l(x)(1)-m,m-q.l(x)(0))
       fun <- function(t)
           {p(x)(m+k*t)-p(x)(m-t)-.5}
       return(uniroot(fun,lower=0,upper=up)$root)
 })                                

