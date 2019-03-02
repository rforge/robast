## new helper function for make and check IC to speed up things

.preparedirectCheckMakeIC <- function(L2Fam, IC, ..., diagnostic = FALSE){

        dims <- length(L2Fam@param)
        trafo <- trafo(L2Fam@param)
        nrvalues <- nrow(trafo)
        Distr <- L2Fam@distribution

        dotsI <- .filterEargsWEargList(list(...))
        if(is.null(dotsI$useApply)) dotsI$useApply <- FALSE


        IC.v <- as(diag(nrvalues) %*% IC@Curve, "EuclRandVariable")
        L2deriv <- as(diag(dims) %*% L2Fam@L2deriv, "EuclRandVariable")

        diagn <- if(diagnostic) vector("list",(nrvalues+3)*nrvalues/2) else NULL
        if(diagnostic) dotsI$diagnostic <- TRUE
        k <- 0

        res <- numeric(nrvalues)
        for(i in 1:nrvalues){
            Eargs <- c(list(object = Distr, fun = IC.v@Map[[i]]), dotsI)
            res[i] <- buf <- do.call(E, Eargs)
            if(diagnostic){ k <- k + 1; diagn[[k]] <- attr(buf,"diagnostic") }
        }
        if(diagnostic){
           attr(res, "diagnostic") <- diagn[1:nrvalues]
           if(!is.null(diagn)) class(attr(res,"diagnostic")) <- "DiagnosticClass"
        }
        erg <- matrix(0, ncol = dims, nrow = nrvalues)

        for(i in 1:nrvalues)
            for(j in 1:dims){
                integrandA <- function(x)IC.v@Map[[i]](x)*L2deriv@Map[[j]](x)
                Eargs <- c(list(object = Distr, fun = integrandA),dotsI)
                erg[i, j] <- buf <- do.call(E, Eargs)
                if(diagnostic){ k <- k + 1; diagn[[k]] <- attr(buf,"diagnostic") }
            }
        if(diagnostic){
           attr(erg, "diagnostic") <- diagn[-(1:nrvalues)]
           if(!is.null(diagn)) class(attr(erg,"diagnostic")) <- "DiagnosticClass"
        }
        return(list(E.IC=res,E.IC.L=erg))
}



## check centering and Fisher consistency
setMethod("checkIC", signature(IC = "IC", L2Fam = "missing"),
    function(IC, out = TRUE, ..., diagnostic = FALSE){
        diagn0stic <- diagnostic
        L2Fam <- eval(IC@CallL2Fam)
        getMethod("checkIC", signature(IC = "IC", L2Fam = "L2ParamFamily"))(
              IC = IC, L2Fam = L2Fam, out = out, ..., diagnostic = diagn0stic)
    })

## check centering and Fisher consistency
setMethod("checkIC", signature(IC = "IC", L2Fam = "L2ParamFamily"),
    function(IC, L2Fam, out = TRUE, ..., diagnostic = FALSE){

        diagn0stic <- diagnostic

        D1 <- L2Fam@distribution
        if(dimension(Domain(IC@Curve[[1]])) != dimension(img(D1)))
            stop("dimension of 'Domain' of 'Curve' != dimension of 'img' of 'distribution' of 'L2Fam'")

        trafo <- trafo(L2Fam@param)

        res <- .preparedirectCheckMakeIC(L2Fam, IC, ..., diagnostic = diagn0stic)

        cent <- res$E.IC
        attr(cent,"diagnostic") <- NULL
        if(out)
            cat("precision of centering:\t", cent, "\n")


        consist <- res$E.IC.L - trafo
        attr(consist,"diagnostic") <- NULL

        if(out){
            cat("precision of Fisher consistency:\n")
            print(consist)
            cat("precision of Fisher consistency - relative error [%]:\n")
            print(100*consist/trafo)
        }

        prec <- max(abs(cent), abs(consist))

        ## PR 20190222:
		## deleting all digits beyond 1e-12 (as numeric fuzz) -- 
		## but check for relative accuracy by means of the "size" of the Fisher information 
		## measured in by the max(trafo)

        names(prec) <- "maximum deviation"
		relPrec <- 12-round(log(max(trafo),10))
		prec <- round(prec*10^relPrec)/10^relPrec

        if(diagnostic && out){
           print(attr(res$E.IC,"diagnostic"),xname="E.IC")
           print(attr(res$E.IC.L,"diagnostic"),xname="E.IC.L")
        }

        if(diagnostic){
           attr(prec,"diagnostic") <- c(attr(res$E.IC,"diagnostic"),
                                        attr(res$E.IC.L,"diagnostic"))
           if(!is.null(attr(prec,"diagnostic")))
              class(attr(prec,"diagnostic")) <- "DiagnosticClass"
        }
        return(prec)
    })


## make some L2function a pIC at a model
setMethod("makeIC", signature(IC = "IC", L2Fam = "L2ParamFamily"),
    function(IC, L2Fam, ..., diagnostic = FALSE){

        diagn0stic <- diagnostic

        dims <- length(L2Fam@param)
        if(dimension(IC@Curve) != dims)
           stop("Dimension of IC and parameter must be equal")

        D1 <- L2Fam@distribution
        if(dimension(Domain(IC@Curve[[1]])) != dimension(img(D1)))
            stop("dimension of 'Domain' of 'Curve' != dimension of 'img' of 'distribution' of 'L2Fam'")

        trafo <- trafo(L2Fam@param)

        res <- .preparedirectCheckMakeIC(L2Fam, IC, ..., diagnostic = diagn0stic)

        if(diagnostic){
           print(attr(res$E.IC,"diagnostic"), xname="E.IC")
           print(attr(res$E.IC.L,"diagnostic"), xname="E.IC.L")
        }

        IC1 <- as(diag(dimension(IC@Curve)) %*% IC@Curve, "EuclRandVariable")

        cent <- res$E.IC
        stand <- trafo %*% distr::solve(res$E.IC.L, generalized = TRUE)

        IC1.0 <- IC1 - cent
        Y <- as(stand %*% IC1.0, "EuclRandVariable")

        modifyIC <- IC@modifyIC

        if(!is.function(IC@modifyIC))
            modifyIC <- function(L2Fam, IC, withMakeIC = FALSE, ...)
                                 return(makeIC(IC,L2Fam, ...))

        CallL2Fam <- L2Fam@fam.call

        IC.0 <- IC(name = name(IC),
                  Curve = EuclRandVarList(Y),
                  Risks = list(),
                  Infos=matrix(c("IC<-",
                                 "generated by affine linear trafo to enforce consistency"),
                               ncol=2, dimnames=list(character(0), c("method", "message"))),
                  CallL2Fam = CallL2Fam,
                  modifyIC = modifyIC)

        if(diagnostic){
           attr(IC.0,"diagnostic") <- c(attr(res$E.IC,"diagnostic"),
                                        attr(res$E.IC.L,"diagnostic"))
           if(!is.null(attr(IC.0,"diagnostic")))
              class(attr(IC.0,"diagnostic")) <- "DiagnosticClass"
        }
        return(IC.0)
    })

## make some L2function a pIC at a model
setMethod("makeIC", signature(IC = "IC", L2Fam = "missing"),
    function(IC, ..., diagnostic = FALSE){
        diagn0stic <- diagnostic
        L2Fam <- eval(IC@CallL2Fam)
        getMethod("makeIC", signature(IC = "IC", L2Fam = "L2ParamFamily"))(
              IC = IC, L2Fam = L2Fam, ..., diagnostic = diagn0stic)
    })

setMethod("makeIC", signature(IC = "list", L2Fam = "L2ParamFamily"),
    function(IC, L2Fam, forceIC = TRUE, name, Risks, Infos, modifyIC = NULL,..., diagnostic = FALSE){
        mc <- match.call(call = sys.call(sys.parent(1)), expand.dots = FALSE)[-1]
        mc0 <- as.list(mc)
        mc0$IC <- NULL
        mc0$L2Fam <- NULL
        mc0$forceIC <- NULL
        mc0$diagnostic <- NULL

        diagn0stic <- diagnostic

        if(!all(as.logical(c(lapply(IC,is.function)))))
           stop("First argument must be a list of functions")

        IC.1 <- lapply(IC, function(IC.2)
                  if(length(formals(IC.2))==0) function(x) IC.2(x) else IC.2)

        mc0$Curve <- EuclRandVarList(RealRandVariable(Map = IC.1, Domain = Reals()))
        mc0$CallL2Fam <- substitute(L2Fam@fam.call)


        IC.0 <- do.call(.IC,mc0)
        if(forceIC) IC.0 <- makeIC(IC.0, L2Fam,..., diagnostic = diagn0stic)
        return(IC.0)
    })



setMethod("makeIC", signature(IC = "function", L2Fam = "L2ParamFamily"),
    function(IC, L2Fam, forceIC = TRUE, name, Risks, Infos,
             modifyIC = NULL,..., diagnostic = FALSE){
        mc <- match.call(call = sys.call(sys.parent(1)), expand.dots = FALSE)[-1]
        mc0 <- as.list(mc)
        mc0$IC <- NULL
        mc0$L2Fam <- NULL
        mc0$forceIC <- NULL
        mc0$diagnostic <- NULL
        diagn0stic <- diagnostic

        IC.1 <- if(length(formals(IC))==0) function(x) IC(x) else IC
        mc0$Curve <- EuclRandVarList(RealRandVariable(Map = list(IC.1),
                         Domain = Reals()))
        mc0$CallL2Fam <- substitute(L2Fam@fam.call)
#        print(mc0)

        IC.0 <- do.call(.IC,mc0)
#        print(IC.0)
        if(forceIC) IC.0 <- makeIC(IC.0, L2Fam,...,diagnostic=diagn0stic)
        return(IC.0)
    })
## comment 20180809: reverted changes in rev 1110

.filterEargsWEargList <- function(dots){
        dotsI <- .filterEargs(dots)
        if(!is.null(dots[["E.argList"]])){
           E.argList <- dots[["E.argList"]]
           if(is.call(E.argList)) eval(E.argList)
           if(is.list(E.argList) && length(E.argList)>0){
              nms.E.argList <- names(E.argList)
              for( item in nms.E.argList) dotsI[[item]] <- E.argList[[item]]
           }
        }

        return(dotsI)
}
