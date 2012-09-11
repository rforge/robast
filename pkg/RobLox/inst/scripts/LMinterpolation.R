###############################################################################
## Interpolated functions to speed up computation of Lagrange Multipliers
## 
## regarding a change to .C calls and approxfun in R 2.16.0, we need to make
## distinction between version before 2.16.0 and afterwards
###############################################################################

library(RobLox)
radius <- c(1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, seq(1e-4, 0.01, by = 0.001),
            seq(0.02, 5, by = 0.01), seq(5.05, 10, by = 0.05))
location <- sapply(radius, rlOptIC, computeIC = FALSE)
scale <- sapply(radius, rsOptIC, computeIC = FALSE)

fun <- function(radius){
  print(radius)
  rlsOptIC.AL(radius, computeIC = FALSE)
}
locationScale <- sapply(radius, fun)
#locationScale <- sapply(radius, rlsOptIC.AL, computeIC = FALSE)

## location
A.loc <- unlist(location[1,])
b.loc <- unlist(location[3,])
if(getRversion() < "2.16.0"){
  .getA.loc.old <- approxfun(radius, A.loc, yleft = 1)
  .getb.loc.old <- approxfun(radius, b.loc, yleft = Inf)
  .getA.loc.old <- approxfun(radius, A.loc, yleft = 1)
  .getb.loc.old <- approxfun(radius, b.loc, yleft = Inf)
}else{
  .getA.loc.new <- approxfun(radius, A.loc, yleft = 1)
  .getb.loc.new <- approxfun(radius, b.loc, yleft = Inf)
  .getA.loc.new <- approxfun(radius, A.loc, yleft = 1)
  .getb.loc.new <- approxfun(radius, b.loc, yleft = Inf)
}

## scale
A.sc <- unlist(scale[1,])
a.sc <- unlist(scale[2,])
b.sc <- unlist(scale[3,])
if(getRversion() < "2.16.0"){
  .getA.sc.old <- approxfun(radius, A.sc, yleft = 0.5)
  .geta.sc.old <- approxfun(radius, a.sc, yleft = 0)
  .getb.sc.old <- approxfun(radius, b.sc, yleft = Inf)
}else{
  .getA.sc.new <- approxfun(radius, A.sc, yleft = 0.5)
  .geta.sc.new <- approxfun(radius, a.sc, yleft = 0)
  .getb.sc.new <- approxfun(radius, b.sc, yleft = Inf)
}

## location and scale
n <- length(radius)
A1.locsc <- unlist(locationScale[1,])[seq(1, 4*n-3, by = 4)]
A2.locsc <- unlist(locationScale[1,])[seq(4, 4*n, by = 4)]
a.locsc <- unlist(locationScale[2,])[seq(2, 2*n, by = 2)]
b.locsc <- unlist(locationScale[3,])
if(getRversion() < "2.16.0"){
  .getA1.locsc.old <- approxfun(radius, A1.locsc, yleft = 1)
  .getA2.locsc.old <- approxfun(radius, A2.locsc, yleft = 0.5)
  .geta.locsc.old <- approxfun(radius, a.locsc, yleft = 0)
  .getb.locsc.old <- approxfun(radius, b.locsc, yleft = Inf)
}else{
  .getA1.locsc.new <- approxfun(radius, A1.locsc, yleft = 1)
  .getA2.locsc.new <- approxfun(radius, A2.locsc, yleft = 0.5)
  .geta.locsc.new <- approxfun(radius, a.locsc, yleft = 0)
  .getb.locsc.new <- approxfun(radius, b.locsc, yleft = Inf)
}

if(getRversion() < "2.16.0"){
  save(.getA.loc.old, .getb.loc.old, .getA.sc.old, .geta.sc.old, .getb.sc.old, 
       .getA1.locsc.old, .getA2.locsc.old, .geta.locsc.old, .getb.locsc.old, 
       file = "savedataOld.rda")
}else{
  save(.getA.loc.new, .getb.loc.new, .getA.sc.new, .geta.sc.new, .getb.sc.new, 
       .getA1.locsc.new, .getA2.locsc.new, .geta.locsc.new, .getb.locsc.new, 
       file = "savedataNew.rda")
}

## Saving the results in sysdata.rda
#load("sysdata.rda")
#load("savedataOld.rda")
#load("savedataNew.rda")
#save(.finiteSampleRadius.loc, .finiteSampleRadius.locsc, .finiteSampleRadius.sc, 
#     .getA1.locsc.new, .getA1.locsc.old, .getA2.locsc.new, .getA2.locsc.old, 
#     .getA.loc.new, .getA.loc.old, .geta.locsc.new, .geta.locsc.old, 
#     .geta.sc.new, .getA.sc.new, .geta.sc.old, .getA.sc.old, 
#     .getb.loc.new, .getb.loc.old, .getb.locsc.new, .getb.locsc.old, 
#     .getb.sc.new, .getb.sc.old, file = "sysdata.rda")

