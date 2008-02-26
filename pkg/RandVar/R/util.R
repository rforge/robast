## small util if imageDistr fails
.getImageDistr <- function(f, distr)
{ if (is(try(return(f(distr)), silent = TRUE),
         "try-error")){
  rl <- function(n) { xr <- r(distr)(n); f(xr) }
  return(new("AbscontDistribution", r = rl, .withArith = TRUE, .withSim = TRUE))}
}         