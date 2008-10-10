setMethod("infoPlot", "IC",
    function(object, ...){
        dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."

        L2Fam <- eval(object@CallL2Fam)
       
        if(!is.null(dots[["lty"]]))  dots["lty"] <- NULL
        if(!is.null(dots[["type"]])) dots["type"] <- NULL
        if(!is.null(dots[["main"]])) dots["main"] <- NULL
        if(!is.null(dots[["sub"]]))  dots["sub"] <- NULL
        if(!is.null(dots[["xlab"]])) dots["xlab"] <- NULL
        if(!is.null(dots[["ylab"]])) dots["ylab"] <- NULL
        if(!is.null(dots[["xlim"]])) dots["xlim"] <- NULL
        if(!is.null(dots[["ylim"]])) dots["ylim"] <- NULL
        if(is.null(dots[["colA"]]))     dots$"colA" <- grey(0.5)
        if(is.null(dots[["colB"]]))     dots$"colB" <- par("col")
        if(is.null(dots[["lwdA"]]))     dots$"lwdA" <- par("lwd")
        if(is.null(dots[["lwdB"]]))     dots$"lwdB" <- 2
        if(is.null(dots[["cex.main"]])) dots$"cex.main" <- 0.8
        


        e1 <- L2Fam@distribution
        if(!is(e1, "UnivariateDistribution") | is(e1, "CondDistribution"))
            stop("not yet implemented")

        if(is(e1, "UnivariateDistribution")){
            if(is(e1, "AbscontDistribution")){
                ifelse(is.finite(q(e1)(0)), lower <- q(e1)(0), lower <- q(e1)(getdistrOption("TruncQuantile")))
                ifelse(is.finite(q(e1)(1)), upper <- q(e1)(1), upper <- q(e1)(1 - getdistrOption("TruncQuantile")))
                h <- upper - lower
                x.vec <- seq(from = lower - 0.1*h, to = upper + 0.1*h, length = 1000)
                plty <- "l"
                lty <- "solid"
            }
            if(is(e1, "DiscreteDistribution")){
                x.vec <- support(e1)
                plty <- "o"
                lty <- "dotted"
            }

            trafo <- L2Fam@param@trafo
            dims <- nrow(trafo)

            QFc <- diag(dims)
            if(is(object,"ContIC") & dims>1 )
               {if (is(normtype(object),"QFNorm")) QFc <- QuadForm(normtype(object))
                QFc0 <- solve( trafo %*% solve(L2Fam@FisherInfo) %*% t(trafo ))
                if (is(normtype(object),"SelfNorm")|is(normtype(object),"InfoNorm")) 
                    QFc <- QFc0
               }
            QFc.5 <- sqrt(PosSemDefSymmMatrix(QFc))

            classIC <- as(trafo %*% solve(L2Fam@FisherInfo) %*% L2Fam@L2deriv, "EuclRandVariable")
            absInfoClass <- t(classIC) %*% QFc %*% classIC
            absInfoClass <- sapply(x.vec, absInfoClass@Map[[1]])

            QF <- diag(dims)
            if(is(object,"ContIC") & dims>1 )
               {if (is(normtype(object),"QFNorm")) QF <- QuadForm(normtype(object))}
            QF.5 <- sqrt(PosSemDefSymmMatrix(QF))

            IC1 <- as(diag(dims) %*% object@Curve, "EuclRandVariable")
            absInfo <- t(IC1) %*% QF %*% IC1
            absInfo <- sapply(x.vec, absInfo@Map[[1]])

            dots["col"] <- dots[["colA"]]
            dots["lwd"] <- dots[["lwdA"]]
            do.call(plot, args=c(list(x.vec, absInfoClass, type = plty, 
                 lty = "dashed", 
                 ylim = c(0, 2*max(absInfo, na.rm = TRUE)), xlab = "x", 
                 ylab = "absolute information"), dots))
            dots["col"] <- dots[["colB"]]
            dots["lwd"] <- dots[["lwdB"]]
            do.call(lines, args=c(list(x.vec, absInfo, type = plty, lty = lty), 
                   dots))
            legend(max(x.vec), 0, xjust = 1, yjust = 0,
                   legend = c("class. opt. IC"), lty = "dashed", 
                               col = c(dots[["colA"]]), cex=0.75)

            if(is.null(L2Fam@param@nuisance))
                do.call(title, args=c(list(paste("Absolute information of (partial) IC for", 
                        name(L2Fam)[1], "\nwith main parameter (", 
                        paste(round(L2Fam@param@main, 3), collapse = ", "), ")")), 
                              dots))
            else
                do.call(title, args=c(list(paste("Absolute information of (partial) IC for", 
                        name(L2Fam)[1], "\nwith main parameter (", 
                        paste(round(L2Fam@param@main, 3), collapse = ", "),
                            ")\nand nuisance parameter (", 
                        paste(round(L2Fam@param@nuisance, 3), collapse = ", "), 
                        ")")), dots)) 

            if(dims > 1){
                nrows <- trunc(sqrt(dims))
                ncols <- ceiling(dims/nrows)
                w0 <- options("warn")
                options(warn = -1)
                opar <- par()
                devNew()
                par(mfrow = c(nrows, ncols))
                IC1.i.5 <- QF.5%*%IC1
                classIC.i.5 <- QFc.5%*%classIC
                for(i in 1:dims){
                    y.vec <- sapply(x.vec, IC1.i.5@Map[[i]])^2/absInfo
                    dots["col"] <- dots[["colB"]]
                    dots["lwd"] <- dots[["lwdA"]]
                    do.call(plot, args=c(list(x.vec, y.vec, type = plty, 
                                  lty = lty, xlab = "x", 
                                  ylab = "relative information", 
                                  ylim = c(0, 1.1)), dots))

                    yc.vec <- sapply(x.vec, classIC.i.5@Map[[i]])^2/absInfoClass
                    dots["col"] <- dots[["colA"]]
                    dots["lwd"] <- dots[["lwdB"]]
                    do.call(lines, args=c(list(x.vec, yc.vec, type = plty, 
                          lty = "dashed"),dots))
                    legend(max(x.vec), 1.1, xjust = 1, cex = 0.6, 
                           legend = c("class. opt. IC"), lty = "dashed", 
                                    col = c(dots[["colA"]]))
                    if(is.null(L2Fam@param@nuisance))
                        do.call(title, args=c(list(paste("Relative information of\ncomponent", 
                                i, "of (partial) IC\nfor", name(L2Fam)[1], 
                                "\nwith main parameter (", 
                                paste(round(L2Fam@param@main, 3), 
                                collapse = ", "), ")")),
                                dots))
                    else
                        do.call(title, args=c(list(paste("Relative information of\ncomponent", 
                                i, "of (partial) IC\nfor", name(L2Fam)[1], 
                                    "\nwith main parameter (", 
                                    paste(round(L2Fam@param@main, 3), 
                                    collapse = ", "),
                                    ")\nand nuisance parameter (", 
                                    paste(round(L2Fam@param@nuisance, 3), 
                                    collapse = ", "), ")")),
                                dots)) 
                }
            }
            par(opar)
            options(w0)
        }
    })
