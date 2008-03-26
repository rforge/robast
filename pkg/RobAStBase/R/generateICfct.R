## generate IC
## for internal use only!
setMethod("generateIC.fct", signature(neighbor = "UncondNeighborhood", L2Fam = "L2ParamFamily"),
    function(neighbor, L2Fam, res){
        A <- res$A
        a <- res$a
        b <- res$b
        d <- res$d
        w <- weight(res$w)
        nrvalues <- nrow(A)
        ICfct <- vector(mode = "list", length = nrvalues)
        Y <- as(A %*% L2Fam@L2deriv - a, "EuclRandVariable")
        L <- as(L2Fam@L2deriv, "EuclRandVariable")
        if(nrvalues == 1){
            if(!is.null(d)){
                ICfct[[1]] <- function(x){}
                body(ICfct[[1]]) <- substitute(
                                        { ind <- 1-.eq(Y(x))
                                          b*(Y(x)*w(L(x)) + zi*(1-ind)*d) },
                                        list(Y = Y@Map[[1]], L = L@Map[[1]], w = w, b = b, d = d,
                                             zi = sign(L2Fam@param@trafo), .eq = .eq))
            }else{
                ICfct[[1]] <- function(x){}
                body(ICfct[[1]]) <- substitute({ Y(x)*w(L(x)) },
                                                 list(Y = Y@Map[[1]], L = L@Map[[1]], w = w))
            }
        }else{
            if(!is.null(d))
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){}
                    body(ICfct[[i]]) <- substitute({ind <- 1-.eq(Y(x))
                                                    ind*b*Yi(x)*w(L(x)) + (1-ind)*d
                                                    },
                                                 list(Yi = Y@Map[[i]], L = L@Map[[1]], w = w,
                                                      b = b, d = d[i],  .eq = .eq))
                }
            else
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){}
                    body(ICfct[[i]]) <- substitute({  Yi(x)*w(L(x))  },
                                                 list(Yi = Y@Map[[i]], L = L@Map[[1]], w = w))
                }
        }
        return(EuclRandVarList(EuclRandVariable(Map = ICfct, Domain = Y@Domain,
                                         Range = Y@Range)))
    })

