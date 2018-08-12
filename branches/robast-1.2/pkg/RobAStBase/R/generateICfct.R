## generate IC
## for internal use only!
setMethod("generateIC.fct", signature(neighbor = "UncondNeighborhood", L2Fam = "L2ParamFamily"),
    function(neighbor, L2Fam, res){
        A <- as.matrix(res$A)
        a <- if(is(neighbor,"TotalVarNeighborhood")) 0 else res$a
        b <- res$b
        d <- if(!is.null(res$d)) res$d else 0
        w <- weight(res$w)
        nrvalues <- nrow(A)
        dims <- ncol(A)
        ICfct <- vector(mode = "list", length = nrvalues)
        L <- as(diag(dims)%*%L2Fam@L2deriv, "EuclRandVariable")
        distr <- distribution(L2Fam)

        L.fct <- function(x) evalRandVar(L,as.matrix(x))[,,1]
        if(nrvalues == 1){
            if(!is.null(res$d)){
                ICfct[[1]] <- function(x){}
                if(dims==1L){
                    body(ICfct[[1]]) <- substitute(
                                            { Lx <- L(x); wx <- w(Lx)
                                              Yx <- A %*% Lx - a
                                              ifelse(1-.eq(Yx),as.numeric(Yx*w(Lx)),zi*d*b) },
                                            list(L = L.fct, w = w, b = b, d = d, A = A, a = a,
                                                zi = sign(trafo(L2Fam@param)), .eq = .eq))
                }else{
                    body(ICfct[[1]]) <- substitute(
                                            { Lx <- L(x); wx <- w(Lx)
                                              Yx <- A %*% Lx - a
                                              ifelse(1-.eq(Yx), as.numeric(Yx*w(Lx)), NA) },
                                            list(L = L.fct, w = w, b = b, d = d, A = A, a = a,
                                                 .eq = .eq))
                }
            }else{
                ICfct[[1]] <- function(x){}
                body(ICfct[[1]]) <- substitute({ Lx <- L(x); wx <- w(Lx); #Lx <- as.matrix(Lx)
                                                 Yx <- A %*% Lx - a
                                                 as.numeric(Yx*wx) },
                                                 list(L = L.fct, A = A, a = a, w = w))
            }
        }else{
            if(!is.null(res$d))
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){}
                    body(ICfct[[i]]) <- substitute({Lx <- L(x)
                                                    Yix <- Ai %*% Lx - ai ; # print(dim(Yix)); print(head(Yix[,1:10]));
                                                    as.numeric(Yix*w(Lx) + .eq(Yix)*di)
                                                    },
                                                 list(L = L.fct, Ai = A[i,,drop=FALSE], ai = a[i], w = w,
                                                      di = d[i]))#,  .eq = .eq))
                }
            else
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){}
                    body(ICfct[[i]]) <- substitute({Lx <- L(x)
                                                    Yix <- Ai %*% Lx - ai
                                                    as.numeric(Yix*w(Lx))  },
                                                 list(L = L.fct, Ai = A[i,,drop=FALSE], ai = a[i], w = w))
                }
        }
        return(EuclRandVarList(EuclRandVariable(Map = ICfct, Domain = L@Domain,
                                         Range = Reals()))) # EuclideanSpace(dimension = nrvalues))))
    })

## comment 20180809: reverted changes in rev 1110 as to generate.fast.fc:
## generate fast IC fct
## for internal use only!
if(FALSE){
generateIC.fast.fct <- function(neighbor, L2Fam, res){
        A <- as.matrix(res$A)
        a <- if(is(neighbor,"TotalVarNeighborhood")) 0 else res$a
        b <- res$b
        d <- res$d
        w <- weight(res$w)
        nrvalues <- nrow(A)
        dims <- ncol(A)
        L <- as(diag(dims)%*%L2Fam@L2deriv, "EuclRandVariable")
        distr <- distribution(L2Fam)
        L.fct <- function(x) evalRandVar(L,as.matrix(x))[,,1]
        fastFct <- function(x){}
        if(nrvalues==1L){
           d0 <- if(dims==1L) d else NA
        }else{
           d0 <- if(!is.null(d)) d else 0
        }
        zi0 <- if(nrvalues==1L && dims==1L) sign(trafo(L2Fam@param)) else 1
        b0 <- if(nrvalues==1L) b else 1
        body(fastFct) <- substitute({ indS <- liesInSupport(Di,x,checkFin=TRUE)
                                      Lx <- L(x)
                                      Yx <- A %*% Lx - a
                                      ind <- 1-.eq(Yx)
                                      ifelse(ind,Yx*w(Lx), zi*d*b)*indS
                                      },
                                      list(L = L.fct, w = w, b = b0,
                                             d = d0 , A = A, a = a, zi = zi0,
                                             .eq = .eq, Di = distr))
        return(fastFct)
    }
}

.fixInLiesInSupport<- function(IC, distr){
   MapL <- IC@Curve[[1]]@Map
   for(i in 1:length(MapL))
      body(IC@Curve[[1]]@Map[[i]]) <- substitute({
         liesInSupport(distr,x,checkFin=TRUE)*fct(x)
      }, list(fct = MapL[[i]], distr=distr))
   return(IC)
}
