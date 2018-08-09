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
        dim <- ncol(A)
        ICfct <- vector(mode = "list", length = nrvalues)
        L <- as(diag(dim)%*%L2Fam@L2deriv, "EuclRandVariable")
        distr <- distribution(L2Fam)
        L.fct <- function(x) evalRandVar(L,x)
        if(nrvalues == 1){
            if(!is.null(res$d)){
                ICfct[[1]] <- function(x){}
                if(all(dim(trafo(L2Fam@param)) == c(1, 1))){
                    body(ICfct[[1]]) <- substitute(
                                            { indS <- liesInSupport(Di,x,checkFin=TRUE)
                                              Lx <- L(x)
                                              Yx <- A %*% Lx - a
                                              ind <- 1-.eq(Yx)
                                              (Yx*w(Lx) + zi*(1-ind)*d*b)*indS },
                                            list(L = L.fct, w = w, b = b, d = d, A = A, a = a,
                                                zi = sign(trafo(L2Fam@param)), .eq = .eq, Di = distr))
                }else{
                    body(ICfct[[1]]) <- substitute(
                                            { indS <- liesInSupport(Di,x,checkFin=TRUE)
                                              Lx <- L(x)
                                              Yx <- A %*% Lx - a
                                              ind <- 1-.eq(Yx)
                                              ifelse(ind, Yx*w(Lx), NA)*indS },
                                            list(L = L.fct, w = w, b = b, d = d, A = A, a = a,
                                                 .eq = .eq, Di = distr))
                }
            }else{
                ICfct[[1]] <- function(x){}
                body(ICfct[[1]]) <- substitute({ indS <- liesInSupport(Di,x,checkFin=TRUE)
                                                 Lx <- L(x)
                                                 Yx <- A %*% Lx - a
                                                 Yx*w(Lx)*indS },
                                                 list(L = L.fct, A = A, a = a, w = w, Di = distr))
            }
        }else{
            if(!is.null(res$d))
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){}
                    body(ICfct[[i]]) <- substitute({indS <- liesInSupport(Di,x,checkFin=TRUE)
                                                    Lx <- L(x)
                                                    Yix <- Ai %*% Lx - ai
                                                    ind <- 1-.eq(Yix)
                                                    (ind*Yix*w(Lx) + (1-ind)*di)*indS
                                                    },
                                                 list(L = L.fct, Ai = A[i,,drop=FALSE], ai = a[i], w = w,
                                                      di = d[i], Di = distr))#,  .eq = .eq))
                }
            else
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){}
                    body(ICfct[[i]]) <- substitute({indS <- liesInSupport(Di,x,checkFin=TRUE)
                                                    Lx <- L(x)
                                                    Yix <- Ai %*% Lx - ai
                                                    Yix*w(Lx)*indS  },
                                                 list(L = L.fct, Ai = A[i,,drop=FALSE], ai = a[i], w = w, Di = distr))
                }
        }
        return(EuclRandVarList(EuclRandVariable(Map = ICfct, Domain = L@Domain,
                                         Range = Reals()))) # EuclideanSpace(dimension = nrvalues))))
    })

## generate fast IC fct
## for internal use only!
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
        L.fct <- function(x) evalRandVar(L,x)
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

