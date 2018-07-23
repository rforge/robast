###############################################################################
## finite-sample under-/overshoot risk
###############################################################################

## taken from ROptEstOld getFiRisk.R svn rev 874 2016-04-25

# cdf of truncated normal distribution
ptnorm <- function(x, mu, A, B){
    ((A <= x)*(x <= B)*(pnorm(x-mu)-pnorm(A-mu))/(pnorm(B-mu)-pnorm(A-mu))
    + (x > B))
}

# n-fold convolution for truncated normal distributions
conv.tnorm <- function(z, A, B, mu, n, m){
    if(n == 1) return(ptnorm(z, mu = mu, A = A, B = B))
    if(z <= n*A) return(0)
    if(z >= n*B) return(1)

    M <- 2^m
    h <- (B-A)/M
    x <- seq(from = A, to = B, by = h)
    p1 <- ptnorm(x, mu = mu, A = A, B = B)
    p1 <- p1[2:(M + 1)] - p1[1:M]

    ## FFT
    pn <- c(p1, numeric((n-1)*M))

    ## convolution theorem for DFTs
    pn <- Re(fft(fft(pn)^n, inverse = TRUE)) / (n*M)
    pn <- (abs(pn) >= .Machine$double.eps)*pn
    i.max <- n*M-(n-2)
    pn <- c(0,pn[1:i.max])
    pn <- cumsum(pn)

    ## cdf with continuity correction h/2
    x <- c(n*A,seq(from = n*A+n/2*h, to = n*B-n/2*h, by=h),n*B)
    pnfun1 <- approxfun(x = x+0.5*h, y = pn, yleft = 0, yright = pn[i.max+1])
    pnfun2 <- function(x) pnfun1(x) / pn[i.max+1]

    return(pnfun2(z))
}
