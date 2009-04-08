## small util if imageDistr fails
.getImageDistr <- function(f, distr){ 
    if (is(distr, "DiscreteDistribution"))
        return(DiscreteDistribution(prob=d(distr)(support(distr)), supp=f(support(distr))))
    if(is(try(return(f(distr)), silent = TRUE), "try-error")){
        rl <- function(n){ 
            xr <- r(distr)(n) 
            f(xr) 
        }
        if(length(unique(rl(10000)))!=10000)
           return(AbscontDistribution(r = rl, .withArith = TRUE, .withSim = TRUE))
        else
           return(UnivarLebDecDistribution(r = rl))
    }
}