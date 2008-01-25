###############################################################################
# weight function
###############################################################################
.ALrlsGetw <- function(x, b, a1, a2, a3){
    hvkt <- sqrt(a3^2*x^4 + (a1^2 - 2*a2*a3^2)*x^2 + a2^2*a3^2)
    ind1 <- (hvkt < b)
    
    return(ind1 + (1-ind1)*b/hvkt)
}

###############################################################################
# computation of r
###############################################################################
.ALrlsGetr <- function(b, r, a1, a2, a3){
    integrandr <- function(x, b, a1, a2, a3){
        hvkt <- sqrt(a3^2*x^4 + (a1^2 - 2*a2*a3^2)*x^2 + a2^2*a3^2)/b - 1
        return((hvkt > 0)*hvkt*dnorm(x))
    }
    Int <- integrate(integrandr, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, a1 = a1, a2 = a2, 
                a3 = a3, b = b)$value

    return(r-sqrt(2*Int))
}


###############################################################################
# computation of a1, a2 and a3
###############################################################################
.ALrlsGeta1a2a3 <- function(b, a1, a2, a3, inteps=1e-12){
    integrand1 <- function(x, b, a1, a2, a3){ 
        x^2*.ALrlsGetw(x, b, a1, a2, a3)*dnorm(x)
    }
    Int1 <- 2*integrate(integrand1, lower = 0, upper = Inf, 
                    rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                    a2 = a2, a3 = a3)$value
    a1 <- 1/Int1

    integrand2 <- function(x, b, a1, a2, a3){
        .ALrlsGetw(x, b, a1, a2, a3)*dnorm(x)
    }
    Int2 <- 2*integrate(integrand2, lower = 0, upper = Inf, 
                    rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                    a2 = a2, a3 = a3)$value
    a2 <- Int1/Int2
    
    integrand3 <- function(x, b, a1, a2, a3){
        (x^2 - a2)^2*.ALrlsGetw(x, b, a1, a2, a3)*dnorm(x)
    }
    Int3 <- 2*integrate(integrand3, lower = 0, upper = Inf, 
                    rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                    a2 = a2, a3 = a3)$value
    a3 <- 1/Int3

    return(list(a1=a1, a2=a2, a3=a3))
}


###############################################################################
# optimal IC
###############################################################################
rlsOptIC.AL <- function(r, mean = 0, sd = 1, A.loc.start = 1, a.sc.start = 0, 
                        A.sc.start = 0.5, bUp = 1000, delta = 1e-6, itmax = 100, 
                        check = FALSE, computeIC = TRUE){
    a1 <- A.loc.start; a2 <- 1+a.sc.start; a3 <- A.sc.start
    b <- uniroot(.ALrlsGetr, lower = 1e-4, upper = bUp, 
            tol = .Machine$double.eps^0.5, r = r, a1 = a1, a2 = a2, 
            a3 = a3)$root

    iter <- 0
    repeat{    
        iter <- iter + 1
        if(iter > itmax){
            stop("Algorithm did not converge!\n", 
                 "=> increase itmax or try different starting values",
                 "'A.loc.start', 'a.sc.start' and 'A.sc.start'\n")
        }
        a1.old <- a1; a2.old <- a2; a3.old <- a3; b.old <- b
        
        a1a2a3 <- .ALrlsGeta1a2a3(b = b, a1 = a1, a2 = a2, a3 = a3)
        a1 <- a1a2a3$a1; a2 <- a1a2a3$a2; a3 <- a1a2a3$a3

        b <- uniroot(.ALrlsGetr, lower = 1e-4, upper = bUp, 
                tol = .Machine$double.eps^0.5, r = r, a1 = a1, a2 = a2, 
                a3 = a3)$root
        if(max(abs(a1.old-a1), abs(a2.old-a2), abs(a3.old-a3), abs(b.old-b))<delta)
            break
    }

    if(check){
        integrand1 <- function(x, b, a1, a2, a3){
            x^2*.ALrlsGetw(x, b, a1, a2, a3)*dnorm(x)
        }
        Int1 <- 2*integrate(integrand1, lower = 0, upper = Inf, 
                        rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                        a2 = a2, a3 = a3)$value
        ch1 <- a1*Int1

        integrand2 <- function(x, b, a1, a2, a3){
            (x^2 - a2)^2*.ALrlsGetw(x, b, a1, a2, a3)*dnorm(x)
        }
        Int2 <- 2*integrate(integrand2, lower = 0, upper = Inf, 
                        rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                        a2 = a2, a3 = a3)$value
        ch2 <- a3*Int2
    
        integrand3 <- function(x, b, a1, a2, a3){
            (x^2 - a2)*.ALrlsGetw(x, b, a1, a2, a3)*dnorm(x)
        }
        Int3 <- 2*integrate(integrand3, lower=0, upper=Inf, 
                        rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                        a2 = a2, a3 = a3)$value
        ch3 <- a3*Int3

        ch4 <- .ALrlsGetr(b = b, r = r, a1 = a1, a2 = a2, a3 = a3)
        
        cat("Fisher consistency of eta.loc:\t", ch1-1, "\n")
        cat("centering of eta.sc:\t", ch3, "\n")
        cat("Fisher consistency of eta.sc:\t", ch2-1, "\n")
        cat("MSE equation:\t", ch4, "\n")
    }
    
    A <- sd^2*diag(c(a1, a3))
    a <- sd*c(0, a3*(a2-1))
    b <- sd*b
    mse <- sd^2*(a1 + a3)
    

    if(computeIC){
        return(generateIC(neighbor = ContNeighborhood(radius = r), 
                    L2Fam = NormLocationScaleFamily(mean = mean, sd = sd), 
                    res = list(A = A, a = a, b = b, d = NULL, 
                               risk = list(asMSE = mse, asBias = b, asCov = mse - r^2*b^2), 
                               info = c("rlsOptIC.AL", "optimally robust IC for AL estimators and 'asMSE'"))))
    }else{
        return(list(A = A, a = a, b = b))
    }
}
