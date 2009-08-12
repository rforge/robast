###############################################################################
## get optimally robust IC for convex asymptotic risks
###############################################################################
setMethod("getInfRobIC", signature(L2deriv = "UnivariateDistribution", 
                                   risk = "asGRisk", 
                                   neighbor = "UncondNeighborhood"),
    function(L2deriv, risk, neighbor, symm, Finfo, trafo, upper = NULL,
             lower = NULL, maxiter, tol,
             warn, noLow = FALSE, verbose = FALSE){
        biastype <- biastype(risk)
        radius <- neighbor@radius
        if(identical(all.equal(radius, 0), TRUE)){
            if(warn) cat("'radius == 0' => (classical) optimal IC\n", 
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobIC(L2deriv = L2deriv, risk = asCov(), 
                        neighbor = neighbor, Finfo = Finfo, trafo = trafo,
                        verbose = verbose)
            res <- c(res, list(biastype = biastype, normtype = NormType()))
            if(!is(risk, "asMSE")){
                Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                                  biastype = biastype, clip = res$b, cent = res$a, 
                                  stand = res$A, trafo = trafo)
                res$risk <- c(Risk, res$risk)
            }
            Cov <- res$risk$asCov
            res$risk$asBias <- list(value = b, biastype = biastype, 
                                   normtype = NormType(), 
                                   neighbortype = class(neighbor))
            res$risk$asMSE <- list(value = Cov + radius^2*b^2, 
                                   r = radius,
                                   at = neighbor)
            return(res)
        }
        z <- 0
        c0 <- 0
        iter <- 0
        if(is(symm, "SphericalSymmetry")) 
            S <- symm@SymmCenter == 0
        else
            S <- FALSE

        prec <- 1
        repeat{
            iter <- iter + 1
            z.old <- z
            c0.old <- c0
            ## new
            L1n <- getL1normL2deriv(L2deriv = L2deriv, cent = z)
            lower0 <-  L1n/(1 + radius^2)
            if(is(neighbor,"TotalVarNeighborhood")) {
                   lower0 <- (L1n-z)/(1 + radius^2)/2}
            upper0 <- max(L1n/radius,
                 sqrt( as.numeric( Finfo + z^2 )/(( 1 + radius^2)^2 - 1) ))
            if (is.null(lower)|(iter == 1))
                  lower <- .Machine$double.eps^0.6
            else {if(iter>1) lower <- max(lower0,lower)}
            if (is.null(upper)|(iter == 1))
                  upper <- 5* max(abs(trafo))*max(Finfo)
            else {if(iter>1) upper <- min(upper,upper0)}
            ##
            c0 <- try(uniroot(getInfClip, 
                  ## new
                        lower = lower, upper = upper,
                  ##
                        tol = tol, L2deriv = L2deriv, risk = risk, 
                        neighbor = neighbor,  biastype = biastype,
                        cent = z, symm = S, 
                        trafo = trafo)$root, silent = TRUE)
            if(!is.numeric(c0)){
                if(warn) cat("The IC algorithm did not converge!\n", 
                             "'radius >= maximum radius' for the given risk?\n",
                             "=> the minimum asymptotic bias (lower case) solution is returned\n")
                res <- getInfRobIC(L2deriv = L2deriv, risk = asBias(biastype = biastype(risk),
                                                              normtype = normtype(risk)), 
                                neighbor = neighbor, Finfo = Finfo, 
                                symm = symm, trafo = trafo, upper = upper, 
                                maxiter = maxiter, tol = tol, warn = warn,
                                verbose = verbose)
                Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                                  biastype = biastype,  clip = res$b, cent = res$a, 
                                  stand = res$A, trafo = trafo)
                res$risk <- c(Risk, res$risk)
                return(res)
            }
            z <- getInfCent(L2deriv = L2deriv, neighbor = neighbor,  biastype = biastype,
                            clip = c0, cent = z, symm = S, trafo = trafo, tol.z = tol)
#            cat("c0:\t", c0, "c0.old:\t", c0.old, "z:\t", z, "z.old:\t", z.old, "\n")
            if(S) break

            prec.old <- prec
            prec <- max(abs(z - z.old), abs(c0-c0.old))
            if(verbose)
                cat("current precision in IC algo:\t", prec, "\n")
            if(prec < tol) break
            if(abs(prec.old - prec) < 1e-10){
                cat("algorithm did not converge!\n", "achieved precision:\t", prec, "\n")
                break
            }
            if(iter > maxiter){
                cat("maximum iterations reached!\n", "achieved precision:\t", prec, "\n")
                break
            }
        }
        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
        A <- getInfStand(L2deriv = L2deriv, neighbor = neighbor,
                     biastype = biastype, clip = c0, cent = z, trafo = trafo)
        a <- as.vector(A)*z
        b <- abs(as.vector(A))*c0
        if(!is(risk, "asMSE")){
            Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor, 
                              biastype = biastype, clip = b, cent = a, stand = A, 
                              trafo = trafo)
        }else{
            Risk <- NULL
        }
        Cov <- getInfV(L2deriv = L2deriv, neighbor = neighbor, 
                       biastype = biastype, clip = c0, cent = z, stand = A)

        Risk <- c(Risk, list(asCov = Cov,  
                     asBias = list(value = b, biastype = biastype, 
                                   normtype = normtype(risk), 
                                   neighbortype = class(neighbor)), 
                     trAsCov = list(value = Cov, normtype = normtype(risk)),
                     asMSE = list(value = Cov + radius^2*b^2, 
                                  r = radius,
                                  at = neighbor)))

        if(is(neighbor,"ContNeighborhood")){
            w <- new("HampelWeight")
            clip(w) <- b
            cent(w) <- as.numeric(z)
            stand(w) <- A
        }else if (is(neighbor,"TotalVarNeighborhood")){
            w <- new("BdStWeight")
            clip(w) <- c(0,b)+a
            stand(w) <- A
        }

        weight(w) <- getweight(w, neighbor = neighbor, biastype = biastype, 
                               normW = NormType())

        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info, w = w,
                    biastype = biastype, normtype = normtype(risk)))
    })



###################################################################################
# multivariate solution G-Risk   --- new 10-08-09
###################################################################################

setMethod("getInfRobIC", signature(L2deriv = "RealRandVariable",
                                   risk = "asGRisk",
                                   neighbor = "UncondNeighborhood"),
    function(L2deriv, risk, neighbor, Distr, DistrSymm, L2derivSymm,
             L2derivDistrSymm, Finfo, trafo, onesetLM = FALSE,
             z.start, A.start, upper = NULL, lower = NULL,
             OptOrIter = "iterate",
             maxiter, tol, warn, verbose = FALSE, ...){

        mc <- match.call()

        ## some abbreviations / checks
        radius <- neighbor@radius
        biastype <- biastype(risk)
        normtype <- normtype(risk)

        p <- nrow(trafo)
        k <- ncol(trafo)

        if( is(neighbor,"TotalVarNeighborhood") && p>1)
           stop("Not yet implemented")

        ## non-standard norms
        FI <- solve(trafo%*%solve(Finfo)%*%t(trafo))
        if(is(normtype,"InfoNorm") || is(normtype,"SelfNorm") ){
           QuadForm(normtype) <- PosSemDefSymmMatrix(FI)
           normtype(risk) <- normtype
        }
        std <- if(is(normtype,"QFNorm"))
                  QuadForm(normtype) else diag(p)

        ## starting values
        if(is.null(z.start)) z.start <- numeric(k)
        if(is.null(A.start)) A.start <- trafo %*% solve(Finfo)


        ## sort out upper solution if radius = 0
        if(identical(all.equal(radius, 0), TRUE))
           return(.getUpperSol(L2deriv = L2deriv, b = b, radius = radius,
                               risk = risk, neighbor = neighbor,
                               biastype = biastype, normtype = normtype,
                               Distr = Distr, Finfo = Finfo, trafo = trafo,
                               QF = std, verbose = verbose, warn = warn))

        ## determine which entries must be computed
        comp <- .getComp(L2deriv, DistrSymm, L2derivSymm, L2derivDistrSymm)
        z.comp <- comp$"z.comp"
        A.comp <- comp$"A.comp"

        ## selection of the algorithm
        pM <- pmatch(tolower(OptOrIter),c("optimize","iterate", "doubleiterate"))
        OptOrIter <- pM
        if (is.na(pM)) OptOrIter <- 1

        ## initialize
        if(is(neighbor,"ContNeighborhood")){
            w <- new("HampelWeight")
        }else{ if (is(neighbor,"TotalVarNeighborhood"))
            w <- new("BdStWeight")
        }
        z <- z.start
        A <- A.start
        b <- 0
        a <- 0
        iter <- 0
        prec <- 1
        iter.In <- 0

        if(is(neighbor,"ContNeighborhood")){
               w <- new("HampelWeight")
        }else if(is(neighbor,"TotalVarNeighborhood")){
               w <- new("BdStWeight")
        }

        ## determining A,a,b with either optimization of iteration:
        if(OptOrIter == 1){
            if(is.null(lower)){
               lowBerg <- .getLowerSol(L2deriv = L2deriv, risk = risk,
                                   neighbor = neighbor, Distr = Distr,
                                   DistrSymm = DistrSymm,
                                   L2derivSymm = L2derivSymm,
                                   L2derivDistrSymm = L2derivDistrSymm,
                                   z.start = z.start, A.start = A.start,
                                   trafo = trafo, maxiter = maxiter, tol = tol,
                                   warn = FALSE, Finfo = Finfo, verbose = verbose)
               lower <- lowBerg$b}
            #if(is.null(upper))
               upper <- 5*max(solve(Finfo))

            iter.In <- 0
            prec.In <- 0
            OptIterCall <- numeric(1)
            Cov <- 0
            Risk <- 1e10
            normtype.old <- normtype
            z.opt <- z
            A.opt <- A
            a.opt <- as.numeric(A.opt %*% z.opt)
            w.opt <- w
            std.opt <- std
            normtype.opt <- normtype

            asGRiskb <- function(b0){
               iter <<- iter + 1
               erg <- getLagrangeMultByOptim(b = b0, L2deriv = L2deriv, risk = risk,
                         FI = Finfo, trafo = trafo, neighbor = neighbor,
                         biastype = biastype, normtype = normtype.opt, Distr = Distr,
                         z.start = z.opt, A.start = A.opt, w.start = w.opt, std = std.opt,
                         z.comp = z.comp, A.comp = A.comp,
                         maxiter = round(maxiter/50*iter^5), tol = tol^(iter^5/40),
                         onesetLM = onesetLM, verbose = verbose, ...)

               w0 <- erg$w
               A0 <- erg$A
               a0 <- erg$a
               z0 <- erg$z
               std0 <- erg$std
               biastype0 <- erg$biastype
               normtype.old0 <- normtype
               normtype0 <- erg$normtype
               risk0 <- erg$risk
               iter.In <<- iter.In + erg$iter
#               cat("Iteration ",iter  ,"--innen:", erg$iter,"\n")
#               print(round(c(maxiter = maxiter/30*iter, tol = tol^(iter/15)),4))
               Cov0 <- getInfV(L2deriv = L2deriv, neighbor = neighbor,
                              biastype = biastype, Distr = Distr,
                              V.comp = A.comp, cent = a0,
                              stand = A0, w = w0)

               if(!is(risk, "asMSE")){
                   Risk0 <- getAsRisk(risk = risk0, L2deriv = L2deriv, neighbor = neighbor,
                                     biastype = biastype0, clip = b0, cent = a0, stand = A0,
                                     trafo = trafo)
               }else{
                   Risk0 <- sum(diag(std0%*%Cov0)) + radius^2 * b0^2
               }
               Risk <<- Risk0
#               print("A")
#               print(list(Risk=Risk0,A=A0,a=a0,z=z0,b=b0))
#               print("...")
               w.opt <<- w <<- w0
               A.opt <<- A <<- A0
               a.opt <<- a <<- a0
               z.opt <<- z <<- z0
               b <<- b0
               std.opt <<- std <<- std0
               biastype <<- biastype0
               normtype.old <<- normtype.old0
               normtype.opt <<- normtype <<- normtype0
               risk <<- risk0
               prec.In <<- erg$prec
               OptIterCall <<- erg$call
               Cov <<- Cov0

               return(Risk0-sum(diag(std0%*%A0%*%t(trafo))))
            }
            tol0 <- tol^.5
            f.l <- asGRiskb(lower)
            f.u <- asGRiskb(upper)

            du <- uniroot(asGRiskb, interval = c(lower,upper), tol = tol0,
                           f.lower=f.l, f.upper=f.u)
#            du <- optimize(asGRiskb, interval = c(lower,upper), tol = tol0,
#                           A.o = A.opt, z.o = z.opt, a.o = a.opt, w.o = w.opt,
#                           std.o = std.opt, normtype.o = normtype.opt)
        }else{
            repeat{
                iter <- iter + 1
                z.old <- z
                b.old <- b
                A.old <- A
                ##
                ## get interval for b-search:
                LUB <- .getLowUpB(L2deriv = L2deriv, Finfo = Finfo, Distr = Distr,
                                  normtype = normtype, z = z, A = A, radius = radius,
                                  iter = iter)
                lower <- LUB$lower
                upper <- if(is.null(upper)) LUB$upper else min(upper,LUB$upper)

#                print(c(lower,upper))
                ## solve for b
                b <- try(uniroot(getInfClip,
                             lower = lower, upper = upper,
                             tol = tol, L2deriv = L2deriv, risk = risk,
                             biastype = biastype, Distr = Distr, neighbor = neighbor,
                             stand = A, cent = z, trafo = trafo)$root, silent = TRUE)

                ## if no solution return lower case:
                if(!is.numeric(b))
                   return(.getLowerSol(L2deriv = L2deriv, risk = risk,
                                       neighbor = neighbor, Distr = Distr,
                                       DistrSymm = DistrSymm,
                                       L2derivSymm = L2derivSymm,
                                       L2derivDistrSymm = L2derivDistrSymm,
                                       z.start = z.start, A.start = A.start,
                                       trafo = trafo, maxiter = maxiter, tol = tol,
                                       warn = warn, Finfo = Finfo, verbose = verbose)
                           )

                maxit2 <- if(OptOrIter==2) 0 else maxiter
                erg <- getLagrangeMultByIter(b = b, L2deriv = L2deriv, risk = risk,
                              trafo = trafo, neighbor = neighbor, biastype = biastype,
                              normtype = normtype, Distr = Distr,
                              z.start = z, A.start = A, w.start = w,
                              std = std, z.comp = z.comp,
                              A.comp = A.comp, maxiter = maxit2, tol = tol,
                              onesetLM = onesetLM, verbose = verbose, warnit = (OptOrIter!=2))

                 ## read out solution
                 w <- erg$w
                 A <- erg$A
                 a <- erg$a
                 z <- erg$z
                 biastype <- erg$biastype
                 normtype.old <- normtype
                 normtype <- erg$normtype
                 risk <- erg$risk
                 iter.In <- iter.In + erg$iter
                 prec.In <- erg$prec
                 OptIterCall <- erg$call
                 std <- erg$std
#                 print(list(z=z,A=A,b=b))

                 if (onesetLM&&maxiter>1){
                     if(is(neighbor,"ContNeighborhood"))
                           cent(w) <- as.numeric(z)
                     if(is(neighbor,"TotalVarNeighborhood"))
                           clip(w) <- c(0,b)+a
                     stand(w) <- A
                     weight(w) <- getweight(w, neighbor = neighbor, biastype = biastype,
                                            normW = normtype)
                     }
                 else normtype <- normtype.old

                 ## check precision and number of iterations in outer b-loop
                 prec.old <- prec
                 prec <- max(abs(b-b.old), max(abs(A-A.old)), max(abs(z-z.old)))
                 if(verbose)
                     cat("current precision in IC algo:\t", prec, "\n")
                 if(prec < tol) break
                 if(abs(prec.old - prec) < 1e-10){
                     cat("algorithm did not converge!\n", "achieved precision:\t", prec, "\n")
                     break
                 }
                 if(iter > maxiter){
                     cat("maximum iterations reached!\n", "achieved precision:\t", prec, "\n")
                     break
                 }
            }
        ### issue some diagnostics if wanted
          if(verbose){
             cat("Iterations needed: outer (b-loop):",
                  iter," inner (A,a-loop):", iter.In,"\n")
             cat("Precision achieved: all in all (b+A,a-loop):",
                 prec," inner (A,a-loop):", prec.In,"\n")
          }
        ### determine Covariance of pIC
          Cov <- getInfV(L2deriv = L2deriv, neighbor = neighbor,
                       biastype = biastype, Distr = Distr,
                       V.comp = A.comp, cent = a,
                       stand = A, w = w)
          if(!is(risk, "asMSE")){
              Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor,
                                biastype = biastype, clip = b, cent = a, stand = A,
                                trafo = trafo)
          }else{
              Risk <- NULL
          }
        }


        ### add some further informations for the pIC-slots info and risk
        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))

        trAsCov <- sum(diag(std%*%Cov))
        Risk <- c(Risk, list(asCov = Cov,
                     asBias = list(value = b, biastype = biastype,
                                   normtype = normtype,
                                   neighbortype = class(neighbor)),
                     trAsCov = list(value = trAsCov,
                                   normtype = normtype),
                     asMSE = list(value = trAsCov + radius^2*b^2,
                                  r = radius,
                                  at = neighbor)))

        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info, w = w,
                    biastype = biastype, normtype = normtype,
                    call = mc, iter = iter, prec = prec, OIcall = OptIterCall,
                    iter.In = iter.In, prec.In = prec.In))
    })


### helper function to return the upper case solution if r=0
.getUpperSol <- function(L2deriv, b, radius, risk, neighbor, biastype,
                       normtype, Distr, Finfo, trafo,
                       QuadForm, verbose, warn){

            if(warn) cat("'radius == 0' => (classical) optimal IC\n",
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobIC(L2deriv = L2deriv, risk = asCov(), neighbor = neighbor,
                               Distr = Distr, Finfo = Finfo, trafo = trafo,
                               QuadForm = QuadForm, verbose = verbose)
            res <- c(res, list(biastype = biastype, normtype = normtype))
            if(!is(risk, "asMSE")){
                Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor,
                                  biastype = biastype, cent = res$a,
                                  stand = res$A, trafo = trafo)
                res$risk <- c(Risk, res$risk)
            }
            trAsCov <- sum(diag(QuadForm%*%res$risk$asCov));
            res$risk$trAsCov <- list(value = trAsCov, normtype = normtype)
            res$risk$asBias <- list(value = b, biastype = biastype,
                                   normtype = normtype,
                                   neighbortype = class(neighbor))
            res$risk$asMSE <- list(value = trAsCov + radius^2*b^2,
                                   r = radius,
                                   at = neighbor)
            return(res)
}

### helper function to return the lower case solution if b-search was not successful
.getLowerSol  <- function(L2deriv, risk, neighbor, Distr, DistrSymm,
                         L2derivSymm, L2derivDistrSymm,
                         z.start, A.start, trafo,
                         maxiter, tol, warn, Finfo, verbose){
                if(warn) cat("Could not determine optimal clipping bound!\n",
                             "'radius >= maximum radius' for the given risk?\n",
                             "=> the minimum asymptotic bias (lower case) solution is returned\n",
                             "If 'no' => Try again with modified starting values ",
                             "'z.start' and 'A.start'\n")
                res <- getInfRobIC(L2deriv = L2deriv,
                                   risk =  asBias(biastype = biastype(risk),
                                                  normtype = normtype(risk)),
                                   neighbor = neighbor, Distr = Distr, DistrSymm = DistrSymm,
                                   L2derivSymm = L2derivSymm, L2derivDistrSymm = L2derivDistrSymm,
                                   z.start = z.start, A.start = A.start, trafo = trafo,
                                   maxiter = maxiter, tol = tol, warn = warn, Finfo = Finfo,
                                   verbose = verbose)
                normtype(risk) <- res$normtype
                if(!is(risk, "asMSE")){
                    Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor,
                                      biastype = biastype, clip = NULL,
                                      cent = res$a, stand = res$A, trafo = trafo)
                    res$risk <- c(Risk, res$risk)
                }
                return(res)

}


### helper function to return upper & lower bounds for b for b-search
.getLowUpB <- function(L2deriv, Finfo, Distr, normtype, z, A, radius, iter){
            L1n <- getL1normL2deriv(L2deriv = L2deriv, cent = z, stand = A,
                                       Distr = Distr, normtype = normtype)
            lower0 <- L1n/(1+radius^2)
            if(is(neighbor,"TotalVarNeighborhood")) {
                   lower0 <- (L1n-A%*%z)/(1 + radius^2)/2}

            QF <- if(is(normtype,"QFNorm")) QuadForm(normtype) else diag(nrow(A))
            upper0 <- max(L1n/radius,
                    sqrt( (sum( diag(QF%*%A%*%Finfo%*%t(A))) + t(A%*%z)%*%QF%*%(A%*%z)) /
                          ((1 + radius^2)^2-1)))

            if (iter == 1){
                lower <- .Machine$double.eps^.6
                upper <- 10*upper0
            }else{
                lower <- lower0
                upper <- 2*upper0
            }

            ##
            return(list(lower=lower, upper=upper))
}

################################################################################
## old routine
################################################################################
#
#setMethod("getInfRobIC", signature(L2deriv = "RealRandVariable",
#                                   risk = "asGRisk",
#                                   neighbor = "UncondNeighborhood"),
#    function(L2deriv, risk, neighbor, Distr, DistrSymm, L2derivSymm,
#             L2derivDistrSymm, Finfo, trafo, onesetLM = FALSE,
#             z.start, A.start, upper = NULL, maxiter, tol, warn, verbose = FALSE){
#        biastype <- biastype(risk)
#        normtype <- normtype(risk)
#        p <- nrow(trafo)
#        if( is(neighbor,"TotalVarNeighborhood") && p>1)
#           stop("Not yet implemented")
#
#        FI <- solve(trafo%*%solve(Finfo)%*%t(trafo))
#        if(is(normtype,"InfoNorm") || is(normtype,"SelfNorm") )
#           {QuadForm(normtype) <- PosSemDefSymmMatrix(FI);
#            normtype(risk) <- normtype}
#        QF <- if(is(normtype,"QFNorm")) QuadForm(normtype) else diag(nrow(trafo))
#
#        if(is.null(z.start)) z.start <- numeric(ncol(trafo))
#        if(is.null(A.start)) A.start <- trafo %*% solve(Finfo)
#
#        radius <- neighbor@radius
#        if(identical(all.equal(radius, 0), TRUE)){
#            if(warn) cat("'radius == 0' => (classical) optimal IC\n",
#                         "in sense of Cramer-Rao bound is returned\n")
#            res <- getInfRobIC(L2deriv = L2deriv, risk = asCov(), neighbor = neighbor,
#                               Distr = Distr, Finfo = Finfo, trafo = trafo,
#                               QuadForm = QF, verbose = verbose)
#            res <- c(res, list(biastype = biastype, normtype = normtype))
#            if(!is(risk, "asMSE")){
#                Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor,
#                                  biastype = biastype, cent = res$a,
#                                  stand = res$A, trafo = trafo)
#                res$risk <- c(Risk, res$risk)
#            }
#            trAsCov <- sum(diag(QF%*%res$risk$asCov));
#            res$risk$trAsCov <- list(value = trAsCov, normtype = normtype)
#            res$risk$asBias <- list(value = b, biastype = biastype,
#                                   normtype = normtype,
#                                   neighbortype = class(neighbor))
#            res$risk$asMSE <- list(value = trAsCov + radius^2*b^2,
#                                   r = radius,
#                                   at = neighbor)
#            return(res)
#        }
#
#        comp <- .getComp(L2deriv, DistrSymm, L2derivSymm, L2derivDistrSymm)
#
#        z.comp <- comp$"z.comp"
#        A.comp <- comp$"A.comp"
#
#        if(is(neighbor,"ContNeighborhood")){
#            w <- new("HampelWeight")
#        }else{ if (is(neighbor,"TotalVarNeighborhood"))
#            w <- new("BdStWeight")
#        }
#        z <- z.start
#        A <- A.start
#        b <- 0
#        a <- 0
#        iter <- 0
#        prec <- 1
#        repeat{
#            iter <- iter + 1
#            z.old <- z
#            b.old <- b
#            A.old <- A
#            ##
#            if(is(neighbor,"ContNeighborhood"))
#               cent(w) <- z
#            stand(w) <- A
#
#            ## new
#            L1n <- getL1normL2deriv(L2deriv = L2deriv, cent = z, stand = A,
#                                       Distr = Distr, normtype = normtype)
#            lower0 <- L1n/(1+radius^2)
#            if(is(neighbor,"TotalVarNeighborhood")) {
#                   lower0 <- (L1n-A%*%z)/(1 + radius^2)/2}
#
#            QF <- if(is(normtype,"QFNorm")) QuadForm(normtype) else diag(nrow(A))
#            upper0 <- max(L1n/radius,
#                    sqrt( (sum( diag(QF%*%A%*%Finfo%*%t(A))) + t(A%*%z)%*%QF%*%(A%*%z)) /
#                          ((1 + radius^2)^2-1)))
#
#            if (!is.null(upper)|(iter == 1))
#                    {lower <- .Machine$double.eps^0.6;
#                     if(is.null(upper)) upper <- 10*upper0
#                }else{ lower <- lower0; upper <- upper0}
#
#            ##
#
#            b <- try(uniroot(getInfClip,
#                  ## new
#                         lower = lower, upper = upper,
#                  ##
#                         tol = tol, L2deriv = L2deriv, risk = risk,
#                         biastype = biastype, Distr = Distr, neighbor = neighbor,
#                         stand = A, cent = z, trafo = trafo)$root, silent = TRUE)
#            if(!is.numeric(b)){
#                if(warn) cat("Could not determine optimal clipping bound!\n",
#                             "'radius >= maximum radius' for the given risk?\n",
#                             "=> the minimum asymptotic bias (lower case) solution is returned\n",
#                             "If 'no' => Try again with modified starting values ",
#                             "'z.start' and 'A.start'\n")
#                res <- getInfRobIC(L2deriv = L2deriv,
#                                   risk =  asBias(biastype = biastype(risk),
#                                                  normtype = normtype(risk)),
#                                   neighbor = neighbor, Distr = Distr, DistrSymm = DistrSymm,
#                                   L2derivSymm = L2derivSymm, L2derivDistrSymm = L2derivDistrSymm,
#                                   z.start = z.start, A.start = A.start, trafo = trafo,
#                                   maxiter = maxiter, tol = tol, warn = warn, Finfo = Finfo,
#                                   verbose = verbose)
#                normtype(risk) <- res$normtype
#                if(!is(risk, "asMSE")){
#                    Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor,
#                                      biastype = biastype, clip = NULL,
#                                      cent = res$a, stand = res$A, trafo = trafo)
#                    res$risk <- c(Risk, res$risk)
#                }
#                return(res)
#            }
#            if(is(neighbor,"ContNeighborhood")){
#                clip(w) <- b
#            }else if(is(neighbor,"TotalVarNeighborhood")){
#                clip(w) <- if(iter==1) c(-1,1)*b/2 else c(0,b)+a
#
#            }
#
#
#            weight(w) <- getweight(w, neighbor = neighbor, biastype = biastype,
#                                   normW = normtype)
#
#            z <- getInfCent(L2deriv = L2deriv, neighbor = neighbor,
#                            biastype = biastype, Distr = Distr, z.comp = z.comp,
#                            w = w)
#
#            if(is(neighbor,"TotalVarNeighborhood")){
#                  a <- z
#                  z <- solve(A,a)
#                  zc <- numeric(ncol(trafo))
#            }else if(is(neighbor,"ContNeighborhood")) {
#                  zc <- z
#            }
#
#            A <- getInfStand(L2deriv = L2deriv, neighbor = neighbor,
#                         biastype = biastype, Distr = Distr, A.comp = A.comp,
#                         cent = zc, trafo = trafo, w = w)
#
#            normtype.old <- normtype
#
#            if (is(normtype,"SelfNorm"))
#                {normtype(risk) <- normtype <- updateNorm(normtype = normtype,
#                   L2 = L2deriv, neighbor = neighbor, biastype = biastype,
#                   Distr = Distr, V.comp = A.comp, cent = as.vector(A %*% z),
#                   stand = A, w = w)}
#
#            prec.old <- prec
#            prec <- max(abs(b-b.old), max(abs(A-A.old)), max(abs(z-z.old)))
#            if(verbose)
#                cat("current precision in IC algo:\t", prec, "\n")
#            if(prec < tol) break
#            if(abs(prec.old - prec) < 1e-10){
#                cat("algorithm did not converge!\n", "achieved precision:\t", prec, "\n")
#                break
#            }
#            if(iter > maxiter){
#                cat("maximum iterations reached!\n", "achieved precision:\t", prec, "\n")
#                break
#            }
#        }
#        if (onesetLM){
#            if(is(neighbor,"ContNeighborhood"))
#               cent(w) <- z
#            if(is(neighbor,"TotalVarNeighborhood"))
#               clip(w) <- c(0,b)+a
#            stand(w) <- A
#            weight(w) <- getweight(w, neighbor = neighbor, biastype = biastype,
#                                   normW = normtype)
#        }else normtype <- normtype.old
#
#        if(is(neighbor,"ContNeighborhood"))
#           a <- as.vector(A %*% z)
#
#        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
#        if(!is(risk, "asMSE")){
#            Risk <- getAsRisk(risk = risk, L2deriv = L2deriv, neighbor = neighbor,
#                              biastype = biastype, clip = b, cent = a, stand = A,
#                              trafo = trafo)
#        }else{
#            Risk <- NULL
#        }
#        Cov <- getInfV(L2deriv = L2deriv, neighbor = neighbor,
#                       biastype = biastype, Distr = Distr,
#                       V.comp = A.comp, cent = a,
#                       stand = A, w = w)
#
#        QF <- if(is(normtype,"QFNorm")) QuadForm(normtype) else diag(nrow(A))
#
#        trAsCov <- sum(diag(QF%*%Cov))
#        Risk <- c(Risk, list(asCov = Cov,
#                     asBias = list(value = b, biastype = biastype,
#                                   normtype = normtype,
#                                   neighbortype = class(neighbor)),
#                     trAsCov = list(value = trAsCov,
#                                   normtype = normtype),
#                     asMSE = list(value = trAsCov + radius^2*b^2,
#                                  r = radius,
#                                  at = neighbor)))
#
#        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info, w = w,
#                    biastype = biastype, normtype = normtype))
#    })























