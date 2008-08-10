###############################################################################
## Evaluate roblox on columns of a matrix
###############################################################################
colRoblox <- function(x, mean, sd, eps, eps.lower, eps.upper, initial.est, k = 1){
    call.est <- match.call()
    if(missing(x))
        stop("'x' is missing with no default")
    if(is.data.frame(x))
        x <- data.matrix(x)
    else
        x <- as.matrix(x)
    if(!is.matrix(x))
        stop("'x' has to be a matrix resp. convertable to a matrix by 'as.matrix'
              or 'data.matrix'")

    res <- rowRoblox(x = t(x), mean = mean, sd = sd, eps = eps, eps.lower = eps.lower,
                     eps.upper = eps.upper, initial.est = initial.est, k = k)
    res@estimate.call <- call.est
    return(res)
}
