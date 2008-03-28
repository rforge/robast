.LowerCaseMultivariate <- function(L2deriv, neighbor, biastype,
             normtype, Distr, L2derivDistrSymm, trafo, z.start,
             A.start, maxiter, tol){

        w <- new("HampelWeight")


        if(is.null(z.start)) z.start <- numeric(ncol(trafo))
        if(is.null(A.start)) A.start <- trafo

        abs.fct <- function(x, L2, stand, cent, normtype){
            X <- evalRandVar(L2, as.matrix(x))[,,1] - cent
            Y <- stand %*% X
            return(fct(normtype)(Y))
        }

        bmin.fct <- function(param, L2deriv, Distr, trafo, z.comp){
            p <- nrow(trafo)
            k <- ncol(trafo)
            A <- matrix(param[1:(p*k)], ncol=k, nrow=p)
            z <- numeric(k)
            z[z.comp] <- param[(p*k+1):length(param)]

            if (is(normtype,"SelfNorm")){
               w0 <- w
               cent(w0) <- z
               stand(w0) <- A
               weight(w0) <- minbiasweight(w0, neighbor = neighbor,
                                           biastype = biastype,
                                           normtype = normtype)
               w <<- w0
               normtype  <<- updateNorm(normtype = normtype, L2 = L2deriv,
                                        neighbor = neighbor, biastype = biastype,
                                        Distr = Distr, V.comp = matrix(TRUE, p,p),
                                        cent = z, stand = A,  w = w,  ...)

               }

            E1 <- E(object = Distr, fun = abs.fct, L2 = L2deriv, stand = A,
                     cent = z, normtype = normtype, useApply = FALSE)
            stA <- if (is(normtype,"QFnorm"))
                       QuadForm(normtype)%*%A else A

            return(E1/sum(diag(stA %*% t(trafo))))
        }

        nrvalues <- length(L2deriv)
        z.comp <- rep(TRUE, nrvalues)
        for(i in 1:nrvalues)
            if(is(L2derivDistrSymm[[i]], "SphericalSymmetry"))
                if(L2derivDistrSymm[[i]]@SymmCenter == 0)
                    z.comp[i] <- FALSE

        A.vec <- as.vector(A.start)
        force(normtype)

        erg <- optim(c(A.vec, z.start[z.comp]), bmin.fct, method = "Nelder-Mead",
                    control = list(reltol = tol, maxit = 100*maxiter),
                    L2deriv = L2deriv, Distr = Distr, trafo = trafo, z.comp = z.comp)

        return(list(erg=erg, w=w, normtype = normtype))
    }



