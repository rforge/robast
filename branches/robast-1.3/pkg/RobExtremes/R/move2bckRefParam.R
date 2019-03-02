setMethod("moveL2Fam2RefParam", signature(L2Fam = "L2ScaleShapeUnion"),
          function(L2Fam, ...){ param <- param(L2Fam)
                                scale <- 1
                                shape <- main(param)[2]
                                main(param) <- c(scale, shape)
                                fixed(param) <- 0
                                modifyModel(L2Fam, param)})

################################################################################

### remains to be done: Risk trafo !!!


setMethod("moveICBackFromRefParam", signature(IC = "IC",
           L2Fam = "L2ScaleShapeUnion"), function(IC, L2Fam, ...){
              param <- param(L2Fam)
              L2call <- L2Fam@fam.call
              param <- param(L2Fam)
              loc <- .loc(L2Fam)
              scale <- main(param)[1]
              IC.cf1 <- IC@Curve[[1]]@Map[[1]]
              IC@Curve[[1]]@Map[[1]] <- function(x) scale*IC.cf1((x-loc)/scale)
              IC.cf2 <- IC@Curve[[1]]@Map[[2]]
              IC@Curve[[1]]@Map[[2]] <- function(x) IC.cf2((x-loc)/scale)
              CallL2Fam(IC) <- L2call
              return(IC)})

setMethod("moveICBackFromRefParam", signature(IC = "IC",
           L2Fam = "L2LocScaleShapeUnion"), function(IC, L2Fam, ...){
              param <- param(L2Fam)
              L2call <- L2Fam@fam.call
              param <- param(L2Fam)
              loc <- main(param)[1]
              scale <- main(param)[2]
              IC.cf0 <- IC@Curve[[1]]@Map[[1]]
              IC@Curve[[1]]@Map[[1]] <- function(x) scale*IC.cf0((x-loc)/scale)
              IC.cf1 <- IC@Curve[[1]]@Map[[2]]
              IC@Curve[[1]]@Map[[2]] <- function(x) scale*IC.cf1((x-loc)/scale)
              IC.cf2 <- IC@Curve[[1]]@Map[[3]]
              IC@Curve[[1]]@Map[[3]] <- function(x) IC.cf2((x-loc)/scale)
              CallL2Fam(IC) <- L2call
              return(IC)})

