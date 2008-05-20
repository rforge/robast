setMethod("comparePlot", signature("IC","IC"),
    function(obj1,obj2, obj3 = NULL, obj4 = NULL, ...){
        L2Fam <- eval(obj1@CallL2Fam)
        L2Fam1c <- obj1@CallL2Fam
        L2Fam2c <- obj2@CallL2Fam
        if(!identical(L2Fam1c,L2Fam2c))
            stop("ICs need to be defined for the same model")

        e1 <- L2Fam@distribution
        if(!is(e1, "UnivariateDistribution")) stop("not yet implemented")

        if(is(e1, "AbscontDistribution")){
            lower <- ifelse(is.finite(q(e1)(0)), q(e1)(0), q(e1)(getdistrOption("TruncQuantile")))
            upper <- ifelse(is.finite(q(e1)(1)), q(e1)(1), q(e1)(1 - getdistrOption("TruncQuantile")))
            h <- upper - lower
            x.vec <- seq(from = lower - 0.1*h, to = upper + 0.1*h, length = 1000)
            plty <- "l"
            lty <- "solid"
        }else{
            if(is(e1, "DiscreteDistribution")){
                x.vec <- support(e1)
                plty <- "p"
                lty <- "dotted"
            }else{
                x.vec <- r(e1)(1000)
                x.vec <- sort(unique(x.vec))
                plty <- "p"
                lty <- "dotted"
            }
        }

        dims <- nrow(L2Fam@param@trafo)
        IC1 <- as(diag(dims) %*% obj1@Curve, "EuclRandVariable")
        IC2 <- as(diag(dims) %*% obj2@Curve, "EuclRandVariable")

        obj <- obj3
        if(is(obj, "IC"))
           {
           if(!identical(L2Fam1c,obj@CallL2Fam))
               stop("ICs need to be defined for the same model")
           IC3 <- as(diag(dims) %*% obj3@Curve, "EuclRandVariable")
           }

        obj <- obj4
        if(is(obj, "IC"))
           {
           if(!identical(L2Fam1c,obj@CallL2Fam))
               stop("ICs need to be defined for the same model")
           IC4 <- as(diag(dims) %*% obj4@Curve, "EuclRandVariable")
           }

        w0 <- options("warn")
        options(warn = -1)
        opar <- par()
        nrows <- trunc(sqrt(dims))
        ncols <- ceiling(dims/nrows)
        par(mfrow = c(nrows, ncols))

        if(is(e1, "DiscreteDistribution"))
                x.vec1 <- seq(from = min(x.vec), to = max(x.vec), length = 1000)

        for(i in 1:dims){
            matp  <- cbind(sapply(x.vec, IC1@Map[[i]]),sapply(x.vec, IC2@Map[[i]]))
            if(is(obj3, "IC"))
                matp  <- cbind(matp,sapply(x.vec, IC3@Map[[i]]))
            if(is(obj4, "IC"))
                matp  <- cbind(matp,sapply(x.vec, IC4@Map[[i]]))

            matplot(x.vec, matp,
                 type = plty, lty = lty,
                 xlab = "x", ylab = "(partial) IC")
            if(is(e1, "DiscreteDistribution")){
                 matp1 <- cbind(sapply(x.vec1, IC1@Map[[i]]),sapply(x.vec1, IC2@Map[[i]]))
                 if(is(obj3, "IC"))
                    matp1  <- cbind(matp1,sapply(x.vec1, IC3@Map[[i]]))
                 if(is(obj4, "IC"))
                    matp1  <- cbind(matp1,sapply(x.vec1, IC4@Map[[i]]))
                 matlines(x.vec1, matp1, lty = "dotted")
                 }

            if(is.null(L2Fam@param@nuisance))
                title(paste("Component", i, "of (partial) ICs\nfor", name(L2Fam)[1],
                            "\nwith main parameter (", paste(round(L2Fam@param@main, 3), collapse = ", "), ")"), cex.main = 0.8)
            else
                title(paste("Component", i, "of (partial) ICs\nfor", name(L2Fam)[1],
                            "\nwith main parameter (", paste(round(L2Fam@param@main, 3), collapse = ", "),
                            ")\nand nuisance parameter (", paste(round(L2Fam@param@nuisance, 3), collapse = ", "), ")"), cex.main = 0.8)
        }
        par(opar)
        options(w0)
        invisible()
    })
