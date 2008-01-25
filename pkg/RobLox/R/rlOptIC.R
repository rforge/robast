###############################################################################
# optimally robust IC for normal location
###############################################################################
rlOptIC <- function(r, mean = 0, sd = 1, bUp = 1000, computeIC = TRUE){
    c0 <- uniroot(function(c0, r){return(r^2*c0 - 2*(dnorm(c0) - c0*pnorm(-c0)))}, 
             lower = 1e-4, upper = bUp, tol = .Machine$double.eps^0.5, r = r)$root
    
    A1 <- 1/(2*pnorm(c0) - 1)
    b <- sd*A1*c0
    A <- sd^2*A1
    
    if(computeIC){
        return(generateIC(neighbor = ContNeighborhood(radius = r), 
                    L2Fam = NormLocationFamily(mean = mean, sd = sd), 
                    res = list(A = as.matrix(A), a = 0, b = b, d = NULL, 
                               risk = list(asMSE = A, asBias = b, asCov = A - r^2*b^2), 
                               info = c("rlOptIC", "optimally robust IC for AL estimators and 'asMSE'"))))
    }else{
        return(list(A = A, a = 0, b = b))
    }
}
