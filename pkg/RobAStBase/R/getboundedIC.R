getBoundedIC <- function(L2Fam, D=trafo(L2Fam@param),..., diagnostic = FALSE){

        dotsI <- .filterEargsWEargList(list(...))
        if(is.null(dotsI$useApply)) dotsI$useApply <- FALSE

        FI <- FisherInfo(L2Fam)
        bm <- sum(diag(distr::solve(FI)))
        w <- new("BoundedWeight", clip = bm, weight = function(x){
                   norm0 <- EuclideanNorm(as.matrix(x))
                   ind2 <- (norm0 < bm/2)
                   norm1 <- ind2*bm/2 + (1-ind2)*norm0
                   ind1 <- (norm0 < bm)
                   return(ind1 + (1-ind1)*bm/norm1)})

        dims <- length(L2Fam@param)

        L2deriv <- as(diag(dims) %*% L2Fam@L2deriv, "EuclRandVariable")

        ICfct <- vector(mode = "list", length = dims)
        L.fct <- function(x) evalRandVar(L2deriv,as.matrix(x))[,,1]

        for(i in 1:dims){
                ICfct[[i]] <- function(x){}
                body(ICfct[[i]]) <- substitute({ Yi(x)*w(L(x)) },
                                                 list(Yi = L2deriv@Map[[i]],
                                                      L = L.fct,
                                                      w = weight(w)))
            }
        L2w <- EuclRandVariable(Map = ICfct, Domain = L2deriv@Domain,
                                         Range = L2deriv@Range)
        D1 <- L2Fam@distribution

        cent <- numeric(dims)
        stand.0 <- matrix(0,dims,dims)

        diagn <- if(diagnostic) vector("list", dims*(dims+3)/2) else NULL
        k <- 0
        if(diagnostic) dotsI$diagnostic <- TRUE

        for(i in 1:dims){
            fun <- function(x) {Lx <- L.fct(x); wx <- weight(w)(Lx); return(Lx[i,]*wx)}
            Eargs <- c(list(object=D1, fun=fun), dotsI)
            cent[i] <- buf <- do.call(E,Eargs)
            if(diagnostic){
               k <- k + 1
               diagn[[k]] <- attr(buf,"diagnostic")
            }
        }
        for(i in 1:dims)
           for(j in i:dims){
            fun <- function(x) {Lx <- L.fct(x); wx <- weight(w)(Lx)
                                return((Lx[i,]-cent[i])*(Lx[j,]-cent[j])*wx)}
            Eargs <- c(list(object=D1, fun=fun), dotsI)
            stand.0[i,j] <- buf <- do.call(E,Eargs)
            if(diagnostic){
               k <- k + 1
               diagn[[k]] <- attr(buf,"diagnostic")
            }
           }
        stand.0[row(stand.0)>col(stand.0)] <- t(stand.0)[row(stand.0)>col(stand.0)]

        stand <- as.matrix(D %*% distr::solve(stand.0, generalized = TRUE))
        L2w0 <- L2w - cent
        res <- as(stand %*% L2w0, "EuclRandVariable")
        if(diagnostic){
           attr(res,"diagnostic") <- diagn
           if(!is.null(attr(res,"diagnostic")))
              class(attr(res,"diagnostic")) <- "DiagnosticClass"
        }
        return(res)
}
