## generate IC
## for internal use only!
setMethod("generateIC.fct", signature(neighbor = "CondContNeighborhood", L2Fam = "L2ParamFamily"),
    function(neighbor, L2Fam, res){
        A <- as.matrix(res$A)
        a <- res$a
        b <- res$b
        d <- res$d
        w <- weight(res$w)

        nrvalues <- nrow(A)
        dim <- ncol(A)

        ICfct <- vector(mode = "list", length = nrvalues)

        a1 <- as(diag(nrvalues) %*% a, "EuclRandVariable")
        Y <- as(A %*% L2Fam@L2deriv, "EuclRandVariable") - a1

        L <- as(diag(dim)%*%L2Fam@L2deriv, "EuclRandVariable")
        L.fct <- function(x) evalRandVar(L,x)

        k <- dimension(img(L2Fam@RegDistr))
        if(nrvalues == 1){
            if(!is.null(d)){
                ICfct[[1]] <- function(x){}
                if(all(dim(trafo(L2Fam@param)) == c(1, 1))){
                    body(ICfct[[1]]) <- substitute(
                                            { ind <- 1-.eq(Y(x))
                                              Y(x)*w(L(x),x[1:k]) + zi*(1-ind)*d*b(x[1:k])},
                                            list(Y = Y@Map[[1]], L = L.fct, w = w, b = b, d = d,
                                                zi = sign(trafo(L2Fam@param)), .eq = .eq, k = k))
                }else{
                    body(ICfct[[1]]) <- substitute(
                                            { ind <- 1-.eq(Y(x))
                                              ifelse(ind, Y(x)*w(L(x),x[1:k]), NA) },
                                            list(Y = Y@Map[[1]], L = L.fct, w = w, b = b, d = d,
                                                 .eq = .eq, k = k))
                }
            }else{
                ICfct[[1]] <- function(x){}
                body(ICfct[[1]]) <- substitute({ Y(x)*w(L(x),x[1:k]) },
                                                 list(Y = Y@Map[[1]], L = L.fct, w = w, k = k))
            }
        }else{
            absY <- sqrt(Y %*% Y)
            if(!is.null(d))
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){}
                    body(ICfct[[i]]) <- substitute({ind <- 1-.eq(Yi(x))
                                                    ind*Yi(x)*w(L(x),x[1:k]) + (1-ind)*d
                                                    },
                                                 list(Yi = Y@Map[[i]], L = L.fct, w = w,
                                                      b = b, d = d[i], k = k))
                }
            else
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){}
                    body(ICfct[[i]]) <- substitute({  Yi(x)*w(L(x),x[1:k])  },
                                                 list(Yi = Y@Map[[i]], L = L.fct,
                                                      w = w, k=k))
                }
        }

       return(EuclRandVarList(EuclRandVariable(Map = ICfct, Domain = Y@Domain,
                                         Range = Y@Range)))
    })

setMethod("generateIC.fct", signature(neighbor = "CondTotalVarNeighborhood", L2Fam = "L2ParamFamily"),
    function(neighbor, L2Fam, res){
        A <- as.matrix(res$A)
        a <- res$a
        b <- res$b
        d <- res$d
        w <- weight(res$w)

        nrvalues <- nrow(A)
        dim <- ncol(A)

        ICfct <- vector(mode = "list", length = nrvalues)

        a1 <- as(diag(nrvalues) %*% a, "EuclRandVariable")
        Y <- as(A %*% L2Fam@L2deriv, "EuclRandVariable") - a1

        L <- as(diag(dim)%*%L2Fam@L2deriv, "EuclRandVariable")
        L.fct <- function(x) evalRandVar(L,x)

        k <- dimension(img(L2Fam@RegDistr))
        if(nrvalues == 1){
            if(!is.null(d)){
                ICfct[[1]] <- function(x){}
                if(all(dim(trafo(L2Fam@param)) == c(1, 1))){
                    body(ICfct[[1]]) <- substitute(
                                            { ind <- 1-.eq(Y(x))
                                              Y(x)*w(L(x),x[1:k]) + zi*(1-ind)*d*b(x[1:k])},
                                            list(Y = Y@Map[[1]], L = L.fct, w = w, b = b, d = d,
                                                zi = sign(trafo(L2Fam@param)), .eq = .eq, k = k))
                }else{
                    body(ICfct[[1]]) <- substitute(
                                            { ind <- 1-.eq(Y(x))
                                              ifelse(ind, Y(x)*w(L(x),x[1:k]), NA) },
                                            list(Y = Y@Map[[1]], L = L.fct, w = w, b = b, d = d,
                                                 .eq = .eq, k = k))
                }
            }else{
                ICfct[[1]] <- function(x){}
                body(ICfct[[1]]) <- substitute({ Y(x)*w(L(x),x[1:k]) },
                                                 list(Y = Y@Map[[1]], L = L.fct, w = w, k = k))
            }
        }else{
            absY <- sqrt(Y %*% Y)
            if(!is.null(d))
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){}
                    body(ICfct[[i]]) <- substitute({ind <- 1-.eq(Yi(x))
                                                    ind*Yi(x)*w(L(x),x[1:k]) + (1-ind)*d
                                                    },
                                                 list(Yi = Y@Map[[i]], L = L.fct, w = w,
                                                      b = b, d = d[i], k = k))
                }
            else
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){}
                    body(ICfct[[i]]) <- substitute({  Yi(x)*w(L(x),x[1:k])  },
                                                 list(Yi = Y@Map[[i]], L = L.fct,
                                                      w = w, k=k))
                }
        }

       return(EuclRandVarList(EuclRandVariable(Map = ICfct, Domain = Y@Domain,
                                         Range = Y@Range)))
    })






setMethod("generateIC.fct", signature(neighbor = "Av1CondTotalVarNeighborhood", L2Fam = "L2ParamFamily"),
    function(neighbor, L2Fam, res){

        A <- as.matrix(res$A)
        normtype <- res$normtype
        biastype <- res$biastype
        w <- weight(res$w)

        nrvalues <- nrow(A)
        dim <- ncol(A)
        ICfct <- vector(mode = "list", length = 1)
        L2 <- L2Fam@ErrorL2deriv[[1]]
        k <- dimension(img(L2Fam@RegDistr))
        if(!is.null(d)){
            ICfct[[1]] <- function(x){}
            body(ICfct[[1]]) <- substitute({ ind1 <- (L2(x[k+1]) > 0); ind2 <- (L2(x[k+1]) < 0)
                                             A <- matrix(A.vec, ncol = k)
                                             Y <- as.vector(A %*% x[1:k])
                                             v <- sqrt(sum(Y^2))
                                             ax <- a(x[1:k])
                                             Y/v*((ax+b)*ind1 + ax*ind2) },
                                           list(A.vec = as.vector(A), L2 = L2@Map[[1]], a = a@Map[[1]],
                                                b = b, k = k))
        }else{
            if(b == Inf){
                ICfct[[1]]<- function(x){}
                body(ICfct[[1]]) <- substitute({ A <- matrix(A.vec, ncol = k)
                                                 v <- as.vector(sqrt(sum((A %*% x[1:k])^2)))
                                                 ax <- a(x[1:k])
                                                 if(ax == -Inf)
                                                     as.vector(A %*% x[1:k])*L2(x[k+1])
                                                 else
                                                     as.vector(A %*% x[1:k])*max(a(x[1:k])/v, L2(x[k+1])) },
                                               list(A.vec = as.vector(A), L2 = L2@Map[[1]], a = a@Map[[1]],
                                                    b = b, k = k))
            }else{
                ICfct[[1]] <- function(x){}
                body(ICfct[[1]]) <- substitute({ A <- matrix(A.vec, ncol = k)
                                                 v <- as.vector(sqrt(sum((A %*% x[1:k])^2)))
                                                 ax <- a(x[1:k])
                                                 as.vector(A %*% x[1:k])*min(max(ax/v, L2(x[k+1])), (ax+b)/v) },
                                               list(A.vec = as.vector(A), L2 = L2@Map[[1]], a = a@Map[[1]],
                                                    b = b, k = k))
            }
        }
       return(EuclRandVarList(EuclRandVariable(Map = ICfct, Domain = Y@Domain,
                                         Range = Y@Range)))
    })

setMethod("generateIC.fct", signature(neighbor = "Av2CondContNeighborhood", L2Fam = "L2ParamFamily"),
    function(neighbor, L2Fam, res){
        A <- res$A
        z <- res$z
        b <- res$b
        d <- res$d
        ICfct <- vector(mode = "list", length = 1)
        L2 <- L2Fam@ErrorL2deriv[[1]]
        k <- dimension(img(L2Fam@RegDistr))
        K.inv <- solve(E(L2Fam@RegDistr, fun = function(x){ x %*% t(x) }))
        trafo <- L2Fam@param@trafo


            if(!is.null(d)){
                #b0 <- b/sqrt(sum(diag(trafo %*% K.inv %*% t(trafo))))  ### to be changed for other norms
                ICfct[[1]] <- function(x){}
                if(all(dim(trafo(L2Fam@param)) == c(1, 1))){
                    body(ICfct[[1]]) <- substitute(
                                            { L2x <- L2(x[k+1])-z
                                              ind <- .eq(L2x)
                                              D <- matrix(D.vec, ncol = k)
                                              K.inv <- matrix(K.vec, ncol = k)
                                              D %*% K.inv %*% x[1:k]*(L2x*w(L2x) + zi*ind*d) },
                                            list(L2 = L2@Map[[1]], D.vec = as.vector(trafo),
                                                 z = z, w = w, d = d, K.vec = as.vector(K.inv),
                                                 zi = sign(trafo(L2Fam@param)), .eq = .eq, k=k))
                }else{
                    body(ICfct[[1]]) <- substitute(
                                            { L2x <- L2(x[k+1])-z
                                              ind <- .eq(L2x)
                                              D <- matrix(D.vec, ncol = k)
                                              K.inv <- matrix(K.vec, ncol = k)
                                              D %*% K.inv %*% x[1:k]*ifelse(ind, NA, L2x*w(L2x)) },
                                            list(L2 = L2@Map[[1]], D.vec = as.vector(trafo),
                                                 z = z, w = w, K.vec = as.vector(K.inv),
                                                 .eq = .eq, k = k))
                }
            }else{
                ICfct[[1]] <- function(x){}
####            c0 <- b/(A*sqrt(sum(diag(K.inv)))) ## allgemeiner:
                body(ICfct[[1]]) <- substitute({ L2x <- L2(x[k+1])-z
                                                 D <- matrix(D.vec, ncol = k)
                                                 K.inv <- matrix(K.vec, ncol = k)
                                                 A*D %*% K.inv %*% x[1:k]*L2x*w(L2x) },
                                            list(L2 = L2@Map[[1]], D.vec = as.vector(trafo),
                                                 z = z, w = w, K.vec = as.vector(K.inv),
                                                 .eq = .eq, k = k))
            }

       return(EuclRandVarList(EuclRandVariable(Map = ICfct, Domain = Y@Domain,
                                         Range = Y@Range)))
    })
##                Curve = EuclRandVarList(EuclRandVariable(Map = ICfct, Domain = L2Fam@L2deriv[[1]]@Domain,
##                                         dimension = trunc(nrow(trafo)))),
