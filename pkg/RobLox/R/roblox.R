###############################################################################
## computation of radius-minimax IC
## using predefined functions included in "sysdata.rda"
###############################################################################
.getlsInterval <- function(r, rlo, rup, mean, sd, delta, A.loc.start, 
                           a.sc.start, A.sc.start, bUp, itmax){
    if(r > 80){
        Ab <- rlsOptIC.AL(r = r, mean = mean, sd = sd, A.loc.start = A.loc.start, 
                          a.sc.start = a.sc.start, A.sc.start = A.sc.start, 
                          bUp = bUp, delta = delta, itmax = itmax, computeIC = FALSE)
        A1 <- Ab$A[1,1]
        A2 <- Ab$A[2,2]
        b <- Ab$b
    }else{
        A1 <- .getA1.locsc(r)
        A2 <- .getA2.locsc(r)
        b <- .getb.locsc(r)
    }

    if(rlo == 0){
        efflo <- (A1 + A2 - b^2*r^2)/(1.5*sd^2)
    }else{
        if(rlo > 80){
            Ablo <- rlsOptIC.AL(r = rlo, mean = mean, sd = sd, A.loc.start = A.loc.start, 
                                a.sc.start = a.sc.start, A.sc.start = A.sc.start, 
                                bUp = bUp, delta = delta, itmax = itmax, computeIC = FALSE)
            efflo <- (A1 + A2 - b^2*(r^2 - rlo^2))/sum(diag(Ablo$A))
        }else{
            A1lo <- .getA1.locsc(rlo)
            A2lo <- .getA2.locsc(rlo)
            efflo <- (A1 + A2 - b^2*(r^2 - rlo^2))/(A1lo + A2lo)
        }
    }

    if(rup > 80){
        Abup <- rlsOptIC.AL(r = rup, mean = mean, sd = sd, A.loc.start = A.loc.start, 
                            a.sc.start = a.sc.start, A.sc.start = A.sc.start, 
                            bUp = bUp, delta = delta, itmax = itmax, computeIC = FALSE)
        effup <- (A1 + A2 - b^2*(r^2 - rup^2))/sum(diag(Abup$A))
    }else{
        A1up <- .getA1.locsc(rup)
        A2up <- .getA2.locsc(rup)
        effup <- (A1 + A2 - b^2*(r^2 - rlo^2))/(A1up + A2up)
    }

    return(effup-efflo)
}
.getlInterval <- function(r, rlo, rup, mean, sd, bUp){
    if(r > 80){
        Ab <- rlOptIC(r = r, mean = mean, sd = sd, bUp = bUp, computeIC = FALSE)
        A <- Ab$A
        b <- Ab$b
    }else{
        A <- .getA.loc(r)
        b <- .getb.loc(r)
    }

    if(rlo == 0){
        efflo <- (A - b^2*r^2)/sd^2
    }else{
        if(rlo > 80){
            Ablo <- rlOptIC(r = rlo, mean = mean, sd = sd, bUp = bUp, computeIC = FALSE)
            efflo <- (Ab$A - Ab$b^2*(r^2 - rlo^2))/Ablo$A
        }else{
            Alo <- .getA.loc(rlo)
            efflo <- (A - b^2*(r^2 - rlo^2))/Alo
        }
    }

    if(rup > 80){
        Abup <- rlOptIC(r = rup, mean = mean, sd = sd, bUp = bUp, computeIC = FALSE)
        effup <- (Ab$A - Ab$b^2*(r^2 - rup^2))/Abup$A
    }else{
        Aup <- .getA.loc(rup)
        effup <- (A - b^2*(r^2 - rup^2))/Aup
    }

    return(effup-efflo)
}
.getsInterval <- function(r, rlo, rup, mean, sd, delta, bUp, itmax){
    if(r > 80){
        Ab <- rsOptIC(r = r, mean = mean, sd = sd, bUp = bUp, delta = delta, 
                      itmax = itmax, computeIC = FALSE)
        A <- Ab$A
        b <- Ab$b
    }else{
        A <- .getA.sc(r)
        b <- .getb.sc(r)
    }

    if(rlo == 0){
        efflo <- (A - b^2*r^2)/(0.5*sd^2)
    }else{
        if(rlo > 80){
            Ablo <- rsOptIC(r = rlo, mean = mean, sd = sd, bUp = bUp, delta = delta, 
                            itmax = itmax, computeIC = FALSE)
            efflo <- (A - b^2*(r^2 - rlo^2))/Ablo$A
        }else{
            Alo <- .getA.sc(rlo)
            efflo <- (A - b^2*(r^2 - rlo^2))/Alo
        }
    }

    if(rup > 80){
        Abup <- rsOptIC(r = rup, mean = mean, sd = sd, bUp = bUp, delta = delta, 
                        itmax = itmax, computeIC = FALSE)
        effup <- (A - b^2*(r^2 - rup^2))/Abup$A
    }else{
        Aup <- .getA.sc(rup)
        effup <- (A - b^2*(r^2 - rup^2))/Aup
    }

    return(effup-efflo)
}
###############################################################################
## optimally robust estimator for normal location and/or scale
###############################################################################
roblox <- function(x, mean, sd, eps, eps.lower, eps.upper, initial.est = "med", 
                   tol = 1e-6, A.loc.start = 1, a.sc.start = 0, A.sc.start = 0.5, 
                   bUp = 1000, itmax = 100, returnIC = FALSE){

    if(missing(x))
        stop("'x' is missing with no default")
    if(missing(eps) & missing(eps.lower) & missing(eps.upper)){
        eps.lower <- 0
        eps.upper <- 0.5
    }
    if(missing(eps)){
        if(!missing(eps.lower) & missing(eps.upper))
            eps.upper <- 0.5
        if(missing(eps.lower) & !missing(eps.upper))
            eps.lower <- 0
        if(!is.numeric(eps.lower) || !is.numeric(eps.upper) || eps.lower >= eps.upper) 
            stop("'eps.lower' < 'eps.upper' is not fulfilled")
        if((eps.lower < 0) || (eps.upper > 0.5))
            stop("'eps.lower' and 'eps.upper' have to be in [0, 0.5]")
    }else{
        if(eps == 0)
            stop("'eps = 0'! => use functions 'mean' and 'sd' for estimation")
        if((eps < 0) || (eps > 0.5))
            stop("'eps' has to be in (0, 0.5]")
    }
    if((initial.est != "ksMD") && (initial.est != "med"))
        stop("invalid 'initial.est'")

    if(missing(mean) & missing(sd)){
        if(!is.numeric(A.loc.start) || !is.numeric(a.sc.start) || !is.numeric(A.sc.start))
            stop("Starting values 'A.loc.start', 'a.sc.start' and 'A.sc.start' have to be numeric")

        if(initial.est == "ksMD"){
            KSdist <- function(param, x){
                if(param[2] <= 0) return(Inf)
                return(ks.test(x, "pnorm", mean = param[1], sd = param[2])$statistic)
            }
            res <- optim(c(0, 1), f = KSdist, method = "Nelder-Mead", 
                         control=list(reltol = tol), x = x)$par
            mean <- res[1]
            sd <- res[2]
        }
        if(initial.est == "med"){
            mean <- median(x, na.rm = TRUE)
            sd <- mad(x, na.rm = TRUE)
        }

        if(!missing(eps)){
            r <- sqrt(length(x))*eps
            if(r > 80){
                IC1 <- rlsOptIC.AL(r = r, mean = mean, sd = sd, 
                                   A.loc.start = A.loc.start, a.sc.start = a.sc.start, 
                                   A.sc.start = A.sc.start, bUp = bUp, delta = tol, 
                                   itmax = itmax)
                Infos(IC1) <- matrix(c("roblox", 
                                       "optimally robust IC for AL estimators and 'asMSE'"), 
                                     ncol = 2, dimnames = list(NULL, c("method", "message")))
            }else{
                A <- sd^2*diag(c(.getA1.locsc(r), .getA2.locsc(r)))
                a <- sd*c(0, .geta.locsc(r))
                b <- sd*.getb.locsc(r)
                mse <- sd^2*(a1 + a3)
                IC1 <- generateIC(neighbor = ContNeighborhood(radius = r), 
                                  L2Fam = NormLocationScaleFamily(mean = mean, sd = sd), 
                                  res = list(A = A, a = a, b = b, d = NULL, 
                                      risk = list(asMSE = mse, asBias = b, asCov = mse - r^2*b^2), 
                                      info = c("roblox", "optimally robust IC for AL estimators and 'asMSE'")))
            }
            robEst <- oneStepEstimator(x, IC1, c(mean, sd))
            if(returnIC)
                return(list(optIC = IC1, mean = robEst[1], sd = robEst[2]))
            else
                return(list(mean = robEst[1], sd = robEst[2]))
        }else{
            sqrtn <- sqrt(length(x))
            rlo <- sqrtn*eps.lower
            rup <- sqrtn*eps.upper
            r <- uniroot(.getlsInterval, lower = rlo+1e-5, upper = rup, 
                         tol = .Machine$double.eps^0.25, rlo = rlo, rup = rup,
                         mean = mean, sd = sd, delta = tol, A.loc.start = A.loc.start, 
                         a.sc.start = a.sc.start, A.sc.start = A.sc.start,
                         bUp = bUp, itmax = itmax)$root
            if(r > 80){
                IC1 <- rlsOptIC.AL(r = r, mean = mean, sd = sd, 
                                   A.loc.start = A.loc.start, a.sc.start = a.sc.start, 
                                   A.sc.start = A.sc.start, bUp = bUp, delta = tol, 
                                   itmax = itmax)
            }else{
                A <- sd^2*diag(c(.getA1.locsc(r), .getA2.locsc(r)))
                a <- sd*c(0, .geta.locsc(r))
                b <- sd*.getb.locsc(r)
                mse <- sd^2*(a1 + a3)
                IC1 <- generateIC(neighbor = ContNeighborhood(radius = r), 
                                  L2Fam = NormLocationScaleFamily(mean = mean, sd = sd), 
                                  res = list(A = A, a = a, b = b, d = NULL, 
                                      risk = list(asMSE = mse, asBias = b, asCov = mse - r^2*b^2)))
            }
            if(rlo == 0){
                ineff <- (sum(diag(stand(IC1))) - clip(IC1)^2*r^2)/(1.5*sd^2)
            }else{
                if(rlo > 80){
                    Ablo <- rlsOptIC.AL(r = rlo, mean = mean, sd = sd, A.loc.start = A.loc.start, 
                                        a.sc.start = a.sc.start, A.sc.start = A.sc.start, 
                                        bUp = bUp, delta = tol, itmax = itmax, computeIC = FALSE)
                    ineff <- (sum(diag(stand(IC1))) - clip(IC1)^2*(r^2 - rlo^2))/sum(diag(Ablo$A))
                }else{
                    A1lo <- .getA1.locsc(rlo)
                    A2lo <- .getA2.locsc(rlo)
                    ineff <- (sum(diag(stand(IC1))) - clip(IC1)^2*(r^2 - rlo^2))/(A1lo + A2lo)
                }
            }
            Infos(IC1) <- matrix(c(rep("roblox", 3), 
                             paste("radius-minimax IC for radius interval [", 
                               round(rlo, 3), ", ", round(rup, 3), "]", sep = ""),
                             paste("least favorable radius: ", round(r, 3), sep = ""),
                             paste("maximum MSE-inefficiency: ", round(ineff, 3), sep = "")), 
                             ncol = 2, dimnames = list(NULL, c("method", "message")))
            robEst <- oneStepEstimator(x, IC1, c(mean, sd))
            if(returnIC)
                return(list(optIC = IC1, mean = robEst[1], sd = robEst[2]))
            else
                return(list(mean = robEst[1], sd = robEst[2]))
        }
    }else{
        if(missing(mean)){
            if(initial.est == "ksMD"){
                KSdist.mean <- function(mean, x, sd){
                    return(ks.test(x, "pnorm", mean = mean, sd = sd)$statistic)
                }
                mean <- optimize(f = KSdist.mean, interval = c(min(x), max(x)), 
                            tol = tol, x = x, sd = sd)$minimum
            }
            if(initial.est == "med") mean <- median(x, na.rm = TRUE)

            if(!missing(eps)){
                r <- sqrt(length(x))*eps
                if(r > 80){
                    IC1 <- rlOptIC(r = r, mean = mean, sd = sd, bUp = bUp)
                    Infos(IC1) <- matrix(c("roblox", 
                                           "optimally robust IC for AL estimators and 'asMSE'"), 
                                         ncol = 2, dimnames = list(NULL, c("method", "message")))
                }else{
                    A <- .getA.loc(r)
                    b <- .getb.loc(r)
                    IC1 <- generateIC(neighbor = ContNeighborhood(radius = r), 
                                      L2Fam = NormLocationFamily(mean = mean, sd = sd), 
                                      res = list(A = as.matrix(A), a = 0, b = b, d = NULL, 
                                          risk = list(asMSE = A, asBias = b, asCov = A - r^2*b^2), 
                                          info = c("roblox", "optimally robust IC for AL estimators and 'asMSE'")))
                }
                robEst <- oneStepEstimator(x, IC1, mean)
                if(returnIC)
                    return(list(optIC = IC1, mean = robEst, sd = sd))
                else
                    return(list(mean = robEst, sd = sd))
            }else{
                sqrtn <- sqrt(length(x))
                rlo <- sqrtn*eps.lower
                rup <- sqrtn*eps.upper
                r <- uniroot(.getlInterval, lower = rlo+1e-5, upper = rup, 
                         tol = .Machine$double.eps^0.25, rlo = rlo, rup = rup,
                         mean = mean, sd = sd, bUp = bUp)$root
                if(r > 80){
                    IC1 <- rlOptIC(r = r, mean = mean, sd = sd, bUp = bUp)
                }else{
                    A <- .getA.loc(r)
                    b <- .getb.loc(r)
                    IC1 <- generateIC(neighbor = ContNeighborhood(radius = r), 
                                      L2Fam = NormLocationFamily(mean = mean, sd = sd), 
                                      res = list(A = as.matrix(A), a = 0, b = b, d = NULL, 
                                          risk = list(asMSE = A, asBias = b, asCov = A - r^2*b^2)))
                }
                if(rlo == 0){
                    ineff <- (as.vector(stand(IC1)) - clip(IC1)^2*r^2)/sd^2
                }else{
                    if(rlo > 80){
                        Ablo <- rlOptIC(r = rlo, mean = mean, sd = sd, bUp = bUp, computeIC = FALSE)
                        ineff <- (as.vector(stand(IC1)) - clip(IC1)^2*(r^2 - rlo^2))/Ablo$A
                    }else{
                        Alo <- .getA.loc(rlo)
                        ineff <- (as.vector(stand(IC1)) - clip(IC1)^2*(r^2 - rlo^2))/Alo
                    }
                }
                Infos(IC1) <- matrix(c(rep("roblox", 3), 
                             paste("radius-minimax IC for radius interval [", 
                               round(rlo, 3), ", ", round(rup, 3), "]", sep = ""),
                             paste("least favorable radius: ", round(r, 3), sep = ""),
                             paste("maximum MSE-inefficiency: ", round(ineff, 3), sep = "")), 
                             ncol = 2, dimnames = list(NULL, c("method", "message")))
                robEst <- oneStepEstimator(x, IC1, mean)
                if(returnIC)
                    return(list(optIC = IC1, mean = robEst, sd = sd))
                else
                    return(list(mean = robEst, sd = sd))
            }
        }
        if(missing(sd)){
            if(initial.est == "ksMD"){
                KSdist.sd <- function(sd, x, mean){
                    return(ks.test(x, "pnorm", mean = mean, sd = sd)$statistic)
                }
                sd <- optimize(f = KSdist.sd, 
                            interval = c(.Machine$double.eps^0.5, max(x)-min(x)), 
                            tol = tol, x = x, mean = mean)$minimum
            }
            if(initial.est == "med") sd <- mad(x, na.rm = TRUE)

            if(!missing(eps)){
                r <- sqrt(length(x))*eps
                if(r > 80){
                    IC1 <- rsOptIC(r = r, mean = mean, sd = sd, 
                                   bUp = bUp, delta = tol, itmax = itmax)
                    Infos(IC1) <- matrix(c("roblox", 
                                           "optimally robust IC for AL estimators and 'asMSE'"), 
                                         ncol = 2, dimnames = list(NULL, c("method", "message")))
                }else{
                    A <- .getA.sc(r)
                    a <- .geta.sc(r)
                    b <- .getb.sc(r)
                    IC1 <- generateIC(neighbor = ContNeighborhood(radius = r), 
                                      L2Fam = NormScaleFamily(mean = mean, sd = sd), 
                                      res = list(A = as.matrix(A), a = a, b = b, d = NULL, 
                                                risk = list(asMSE = A, asBias = b, asCov = A - r^2*b^2), 
                                                info = c("roblox", "optimally robust IC for AL estimators and 'asMSE'")))
                }
                robEst <- oneStepEstimator(x, IC1, sd)
                if(returnIC)
                    return(list(optIC = IC1, mean = mean, sd = robEst))
                else
                    return(list(mean = mean, sd = robEst))
            }else{
                sqrtn <- sqrt(length(x))
                rlo <- sqrtn*eps.lower
                rup <- sqrtn*eps.upper
                r <- uniroot(.getsInterval, lower = rlo+1e-5, upper = rup, 
                         tol = .Machine$double.eps^0.25, rlo = rlo, rup = rup,
                         mean = mean, sd = sd, delta = tol, bUp = bUp, 
                         itmax = itmax)$root
                if(r > 80){
                    IC1 <- rsOptIC(r = r, mean = mean, sd = sd, bUp = bUp, 
                                   delta = tol, itmax = itmax)
                }else{
                    A <- .getA.sc(r)
                    a <- .geta.sc(r)
                    b <- .getb.sc(r)
                    IC1 <- generateIC(neighbor = ContNeighborhood(radius = r), 
                                      L2Fam = NormScaleFamily(mean = mean, sd = sd), 
                                      res = list(A = as.matrix(A), a = a, b = b, d = NULL, 
                                                risk = list(asMSE = A, asBias = b, asCov = A - r^2*b^2)))
                }
                if(rlo == 0){
                    ineff <- (as.vector(stand(IC1)) - clip(IC1)^2*r^2)/(0.5*sd^2)
                }else{
                    if(rlo > 80){
                        Ablo <- rsOptIC(r = rlo, mean = mean, sd = sd, bUp = bUp, delta = tol, 
                                        itmax = itmax, computeIC = FALSE)
                        ineff <- (as.vector(stand(IC1)) - clip(IC1)^2*(r^2 - rlo^2))/Ablo$A
                    }else{
                        Alo <- .getA.sc(rlo)
                        ineff <- (as.vector(stand(IC1)) - clip(IC1)^2*(r^2 - rlo^2))/Alo
                    }
                }
                Infos(IC1) <- matrix(c(rep("roblox", 3), 
                             paste("radius-minimax IC for radius interval [", 
                               round(rlo, 3), ", ", round(rup, 3), "]", sep = ""),
                             paste("least favorable radius: ", round(r, 3), sep = ""),
                             paste("maximum MSE-inefficiency: ", round(ineff, 3), sep = "")), 
                             ncol = 2, dimnames = list(NULL, c("method", "message")))
                robEst <- oneStepEstimator(x, IC1, sd)
                if(returnIC)
                    return(list(optIC = IC1, mean = mean, sd = robEst))
                else
                    return(list(mean = mean, sd = robEst))
            }
        }
    }
}
