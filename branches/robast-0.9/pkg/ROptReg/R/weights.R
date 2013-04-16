setMethod("clip", "CondBoundedWeight", function(x1) x1@clip)
setReplaceMethod("clip", "CondBoundedWeight",
    function(object, value){
        object@clip <- value
        object
    })


setMethod("stand", "CondBdStWeight", function(object) object@stand)
setReplaceMethod("stand", "CondBdStWeight",
    function(object, value){
        object@stand <- value
        object
    })

setMethod("cent", "CondHampelWeight", function(object) object@cent)

setReplaceMethod("cent", "CondHampelWeight",
    function(object, value){
        object@cent <- value
        object
    })


.fac <- function(norm,D,Kinv){
   if(is(norm,"SelfNorm")||is(norm,"InfoNorm")){
      1/sqrt(nrow(D))
   }else{ if(norm,"QFNorm"){
             B <- QuadForm(norm)
             sum(diag(t(D)%*%B%*%D%*%Kinv))^-.5
          }else{
             sum(diag(t(D)%*%D%*%Kinv))^-.5
          }
   }
}

setMethod("getweight",
          signature(Weight = "CondHampelWeight", neighbor = "ContNeighborhood",
                    biastype = "BiasType"),# normtype = "NormType"),
          getMethod("getweight", signature=signature(Weight = "HampelWeight",
                    neighbor = "ContNeighborhood", biastype = "BiasType")))

setMethod("getweight",
          signature(Weight = "CondHampelWeight", neighbor = "TotalVarNeighborhood",
                    biastype = "BiasType"),
          getMethod("getweight", signature=signature(Weight = "HampelWeight",
                    neighbor = "TotalVarNeighborhood", biastype = "BiasType")))

setMethod("getweight",
          signature(Weight = "CondHampelWeight", neighbor = "Av1CondContNeighborhood",
                    biastype = "BiasType"),#  norm = "missing"),
          function(Weight, neighbor, biastype, normW, ...)
               {A <- stand(Weight)
                b <- clip(Weight)
                z <- cent(Weight)
                function(x,X){
                   y <- as.numeric(as.matrix(A)%*%(x-z(X)))
                   norm0 <- fct(normW)(y)
                   ind2 <- (norm0 < b/2)
                   norm1 <- ind2*b/2 + (1-ind2)*norm0
                   ind1 <- (norm0 < b)
                   ind1 + (1-ind1)*b/norm1
                   }
                }
          )

setMethod("getweight",
          signature(Weight = "CondHampelWeight", neighbor = "Av1CondContNeighborhood",
                    biastype = "onesidedBias"),#  norm = "missing"),
          function(Weight, neighbor, biastype, ...)
               {A <- stand(Weight)
                b <- clip(Weight)
                z <- cent(Weight)
                function(x,X){
                   y <- as.numeric(as.matrix(A)%*%(x-z(X)))*sign(biastype)
                   norm1 <- pmax(y,b/2)
                   pmin(1,b/norm1)
                   }
                }
          )

setMethod("getweight",
          signature(Weight = "CondHampelWeight", neighbor = "Av1CondContNeighborhood",
                    biastype = "asymmetricBias"),# norm = "missing"),
          function(Weight, neighbor, biastype, ...)
               {A <- stand(Weight)
                b <- clip(Weight)
                b1 <- b/nu(biastype)[1]
                b2 <- b/nu(biastype)[2]
                z <- cent(Weight)
                function(x,X){
                   y <- as.numeric(as.matrix(A)%*%(x-z(X)))
                   norm1 <- pmax(-y,b1/2)
                   norm2 <- pmax(y,b2/2)
                   pmin(1,b1/norm1,b2/norm2)
                   }
                }
          )


setMethod("getweight",
          signature(Weight = "CondBdStWeight", neighbor = "Av1CondTotalVarNeighborhood",
                    biastype = "BiasType"),#  norm = "missing"),
          function(Weight, neighbor, biastype, ...)
               {A <- stand(Weight)
                b <- clip(Weight)
                b1 <- -b[1]
                b2 <- b[2]
                function(x,X){
                   y <- as.numeric(as.matrix(A)%*%x)
                   norm1 <- pmax(-y,b1/2)
                   norm2 <- pmax(y,b2/2)
                   pmin(1,b1/norm1,b2/norm2)
                   }
                }
          )

setMethod("getweight",
          signature(Weight = "HampelWeight", neighbor = "Av2CondContNeighborhood",
                    biastype = "BiasType"),# normtype = "NormType"),
          function(Weight, neighbor, biastype, normW, Kinv, D, ...)
               {fac <- .fac(normW,D,Kinv)
                A <- stand(Weight)
                b <- clip(Weight)/fac
                z <- cent(Weight)
                function(x){
                   y <- A%*%(x-z)
                   norm0 <- fct(normW)(y)
                   ind2 <- (norm0 < b/2)
                   norm1 <- ind2*b/2 + (1-ind2)*norm0
                   ind1 <- (norm0 < b)
                   ind1 + (1-ind1)*b/norm1
                   }
                }
          )


setMethod("getweight",
          signature(Weight = "HampelWeight", neighbor = "Av2CondContNeighborhood",
                    biastype = "onesidedBias"),#  norm = "missing"),
          function(Weight, neighbor, biastype, normW, Kinv, D, ...)
               {A <- stand(Weight)
                fac <- .fac(normW,D,Kinv)
                b <- clip(Weight)/fac
                z <- cent(Weight)
                function(x){
                   y <- as.numeric(as.matrix(A)%*%(x-z))*sign(biastype)
                   norm1 <- pmax(y,b/2)
                   pmin(1,b/norm1)
                   }
                }
          )

setMethod("getweight",
          signature(Weight = "HampelWeight", neighbor = "Av2CondContNeighborhood",
                    biastype = "asymmetricBias"),# norm = "missing"),
          function(Weight, neighbor, biastype, normW, Kinv, D, ...)
               {A <- stand(Weight)
                fac <- .fac(normW,D,Kinv)
                b <- clip(Weight)/fac
                b1 <- b/nu(biastype)[1]
                b2 <- b/nu(biastype)[2]
                z <- cent(Weight)
                function(x){
                   y <- as.numeric(as.matrix(A)%*%(x-z))
                   norm1 <- pmax(-y,b1/2)
                   norm2 <- pmax(y,b2/2)
                   pmin(1,b1/norm1,b2/norm2)
                   }
                }
          )


setMethod("getweight",
          signature(Weight = "BdStWeight", neighbor = "Av2CondTotalVarNeighborhood",
                    biastype = "BiasType"),#  norm = "missing"),
          function(Weight, neighbor, biastype, normW, Kinv, D, ...)
               {A <- stand(Weight)
                fac <- .fac(normW,D,Kinv)
                b <- clip(Weight)/fac
                b1 <- -b[1]
                b2 <- b[2]
                function(x){
                   y <- as.numeric(as.matrix(A)%*%x)
                   norm1 <- pmax(-y,b1/2)
                   norm2 <- pmax(y,b2/2)
                   pmin(1,b1/norm1,b2/norm2)
                   }
                }
          )

setMethod("minbiasweight",
          signature(Weight = "CondHampelWeight", neighbor = "ContNeighborhood",
                    biastype = "BiasType"),
getMethod("minbiasweight", signature=signature(Weight = "CondHampelWeight",
                    neighbor = "ContNeighborhood",
                    biastype = "BiasType")))

setMethod("minbiasweight",
          signature(Weight = "CondHampelWeight", neighbor = "TotalVarNeighborhood",
                    biastype = "BiasType"),
getMethod("minbiasweight", signature=signature(Weight = "CondHampelWeight",
                    neighbor = "TotalVarNeighborhood",
                    biastype = "BiasType")))

setMethod("minbiasweight",
          signature(Weight = "CondHampelWeight", neighbor = "Av1CondContNeighborhood",
                    biastype = "BiasType"),#  norm = "missing"),
          function(Weight, neighbor, biastype, normW)
               {A <- stand(Weight)
                b <- clip(Weight)
                z <- cent(Weight)
                function(x,X){
                   y <- A%*%(x-z(X))
                   norm0 <- fct(normW)(y)
                   ind <- 1-.eq(norm0)
                   ind*b(X)/(norm0+1-ind)
                   }
                }
          )


setMethod("minbiasweight",
          signature(Weight = "CondHampelWeight", neighbor = "Av1CondContNeighborhood",
                    biastype = "asymmetricBias"),#  norm = "missing"),
          function(Weight, neighbor, biastype, ...)
               {A <- stand(Weight)
                b <- clip(Weight)
                b1 <- function(X) -b(X)[1]
                b2 <- function(X) b(X)[2]
                z <- cent(Weight)
                function(x,X){
                   y <- as.numeric(as.matrix(A)%*%(x-z(X)))
                   indp <- (y>0)
                   ind0 <- .eq(y)
                   indm <- (y<0)
                   indm*b1(X)/(y+ind0) + indp*b2(X)/(y+ind0)
                   }
                }
          )

setMethod("minbiasweight",
          signature(Weight = "CondHampelWeight", neighbor = "Av1CondContNeighborhood",
                    biastype = "onesidedBias"),#  norm = "missing"),
          function(Weight, neighbor, biastype, ...)
               {A <- stand(Weight)
                b <- clip(Weight)
                z <- cent(Weight)
                function(x,X){
                   y <- as.numeric(as.matrix(A)%*%(x-z(X)))
                   ind <- (y*sign(biastype) >0)
                   ind0 <- .eq(y)
                   ind*b(X)/(y+ind0)+(1-ind)
                   }
                }
          )


setMethod("minbiasweight",
          signature(Weight = "CondBdStWeight", neighbor = "Av1CondTotalVarNeighborhood",
                    biastype = "BiasType"),
          function(Weight, neighbor, biastype, ...)
               {A <- stand(Weight)
                b <- clip(Weight)
                b1 <- function(X) b(X)[1]
                b2 <- function(X) b(X)[2]
                function(x,X){
                   y <- as.numeric(as.matrix(A)%*%(x))
                   indp <- (y>0)
                   ind0 <- .eq(y)
                   indm <- (y<0)
                   indm*b1(X)/(y+ind0) + indp*b2(X)/(y+ind0)
                   }
                }
          )

setMethod("minbiasweight",
          signature(Weight = "HampelWeight", neighbor = "Av2CondContNeighborhood",
                    biastype = "BiasType"),#  norm = "NormType"),
          function(Weight, neighbor, biastype, normW, Kinv, D, ...)
               {A <- stand(Weight)
                fac <- .fac(normW,D,Kinv)
                b <- clip(Weight)*fac
                z <- cent(Weight)
                function(x){
                   y <- A%*%(x-z)
                   norm0 <- fct(normW)(y)
                   ind <- 1-.eq(norm0)
                   ind*b/(norm0+1-ind)
                   }
                }
          )

setMethod("minbiasweight",
          signature(Weight = "HampelWeight", neighbor = "Av2CondContNeighborhood",
                    biastype = "asymmetricBias"),#  norm = "missing"),
          function(Weight, neighbor, biastype, normW, Kinv, D, ...)
               {A <- stand(Weight)
                fac <- .fac(normW,D,Kinv)
                b <- clip(Weight)*fac
                b1 <- -b[1]
                b2 <- b[2]
                z <- cent(Weight)
                function(x){
                   y <- as.numeric(as.matrix(A)%*%(x-z))
                   indp <- (y>0)
                   ind0 <- .eq(y)
                   indm <- (y<0)
                   indm*b1/(y+ind0) + indp*b2/(y+ind0)
                   }
                }
          )

setMethod("minbiasweight",
          signature(Weight = "HampelWeight", neighbor = "Av2CondContNeighborhood",
                    biastype = "onesidedBias"),#  norm = "missing"),
          function(Weight, neighbor, biastype, normW, Kinv, D, ...)
               {A <- stand(Weight)
                fac <- .fac(normW,D,Kinv)
                b <- clip(Weight)*fac
                z <- cent(Weight)
                function(x){
                   y <- as.numeric(as.matrix(A)%*%(x-z))
                   ind <- (y*sign(biastype) >0)
                   ind0 <- .eq(y)
                   ind*b/(y+ind0)+(1-ind)
                   }
                }
          )


setMethod("minbiasweight",
          signature(Weight = "BdStWeight", neighbor = "Av2CondTotalVarNeighborhood",
                    biastype = "BiasType"),
          function(Weight, neighbor, biastype, normW, Kinv, D, ...)
               {A <- stand(Weight)
                fac <- .fac(normW,D,Kinv)
                b <- clip(Weight)*fac
                b1 <- b[1]
                b2 <- b[2]
                function(x){
                   y <- as.numeric(as.matrix(A)%*%(x))
                   indp <- (y>0)
                   ind0 <- .eq(y)
                   indm <- (y<0)
                   indm*b1/(y+ind0) + indp*b2/(y+ind0)
                   }
                }
          )


