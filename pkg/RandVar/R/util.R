## small util if imageDistr fails
.getImageDistr <- function(f, distr)
{ if (is(try(return(f(distr)), silent = TRUE),
         "try-error")){
  rl <- function(n) { xr <- r(distr)(n); f(xr) }
  dr <- RtoDPQ(r = rl)
  return(new("AbscontDistribution", d = dr$dfun, 
       r = rl, p = dr$pfun, q = dr$qfun, .withSim = TRUE))}
}       