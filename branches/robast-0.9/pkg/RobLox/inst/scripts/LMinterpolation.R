###############################################################################
## Interpolated functions to speed up computation of Lagrange Multipliers
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

A.loc <- unlist(location[1,])
b.loc <- unlist(location[3,])
.getA.loc <- approxfun(radius, A.loc, yleft = 1)
.getb.loc <- approxfun(radius, b.loc, yleft = Inf)

A.sc <- unlist(scale[1,])
a.sc <- unlist(scale[2,])
b.sc <- unlist(scale[3,])
.getA.sc <- approxfun(radius, A.sc, yleft = 0.5)
.geta.sc <- approxfun(radius, a.sc, yleft = 0)
.getb.sc <- approxfun(radius, b.sc, yleft = Inf)

n <- length(radius)
A1.locsc <- unlist(locationScale[1,])[seq(1, 4*n-3, by = 4)]
A2.locsc <- unlist(locationScale[1,])[seq(4, 4*n, by = 4)]
a.locsc <- unlist(locationScale[2,])[seq(2, 2*n, by = 2)]
b.locsc <- unlist(locationScale[3,])
.getA1.locsc <- approxfun(radius, A1.locsc, yleft = 1)
.getA2.locsc <- approxfun(radius, A2.locsc, yleft = 0.5)
.geta.locsc <- approxfun(radius, a.locsc, yleft = 0)
.getb.locsc <- approxfun(radius, b.locsc, yleft = Inf)

save(.getA.loc, .getb.loc, .getA.sc, .geta.sc, .getb.sc, .getA1.locsc, .getA2.locsc,
     .geta.locsc, .getb.locsc, file = "savedata.rda")

