.LowerCaseMultivariate <- function(L2deriv, neighbor, biastype,
             normtype, Distr, trafo, z.start,
             A.start, z.comp, A.comp, maxiter, tol){

        w <- new("HampelWeight")

        if(is.null(z.start)) z.start <- numeric(ncol(trafo))
        if(is.null(A.start)) A.start <- trafo
        if(is.null(A.comp)) 
           A.comp <- matrix(TRUE, nrow = nrow(trafo), ncol = ncol(trafo))
        if(is.null(z.comp)) 
           z.comp <- rep(TRUE, nrow(trafo))

        lA.comp <- sum(A.comp)
        
        abs.fct <- function(x, L2, stand, cent, normtype){
            X <- evalRandVar(L2, as.matrix(x))[,,1] - cent
            Y <- stand %*% X
            return(fct(normtype)(Y))
        }

        bmin.fct <- function(param, L2deriv, Distr, trafo){
            p <- nrow(trafo)
            k <- ncol(trafo)
            A <- matrix(0, ncol = k, nrow = p)
            A[A.comp] <- param[1:lA.comp]
            z <- numeric(k)
            z[z.comp] <- param[(lA.comp+1):length(param)]
            
#            if(is(normtype,"SelfNorm")) 
#               A <- A/max(A)
            
            w0 <- w
            cent(w0) <- z
            stand(w0) <- A
            weight(w0) <- minbiasweight(w0, neighbor = neighbor,
                                           biastype = biastype,
                                           normW = normtype)
            w <<- w0
            if (is(normtype,"SelfNorm")){
               normtype  <<- updateNorm(normtype = normtype, L2 = L2deriv,
                                        neighbor = neighbor, biastype = biastype,
                                        Distr = Distr, V.comp = A.comp,
                                        cent = z, stand = A,  w = w0)
               }

            E1 <- E(object = Distr, fun = abs.fct, L2 = L2deriv, stand = A,
                     cent = z, normtype = normtype, useApply = FALSE)
            stA <- if (is(normtype,"QFnorm"))
                       QuadForm(normtype)%*%A else A
            erg <- E1/sum(diag(stA %*% t(trafo)))
            clip(w0) <- 1/erg
            w <<- w0
            return(erg)
        }

        A.vec <- as.vector(A.start[A.comp])
        p.vec <- c(A.vec, z.start[z.comp])
        force(normtype)
        erg <- optim(p.vec, bmin.fct, method = "Nelder-Mead",
                    control = list(reltol = tol, maxit = 100*maxiter),
                    L2deriv = L2deriv, Distr = Distr, trafo = trafo)

        return(list(erg=erg, w=w, normtype = normtype, z.comp = z.comp))
    }


