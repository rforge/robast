#if(FALSE){
## faster check for ContICs

setMethod("checkIC", signature(IC = "ContIC", L2Fam = "L2ParamFamily"),
    function(IC, L2Fam, out = TRUE, forceContICMethod = FALSE, ..., diagnostic = FALSE){

        D1 <- L2Fam@distribution
        if( dimension(Domain(IC@Curve[[1]])) != dimension(img(D1)))
            stop("dimension of 'Domain' of 'Curve' != dimension of 'img' of 'distribution' of 'L2Fam'")

         res <- .prepareCheckMakeIC(L2Fam, w = IC@weight, forceContICMethod, ..., diagnostic = diagnostic)
         ## if it pays off to use symmetry/ to compute integrals in L2deriv space
        ## we compute the following integrals:
        ## G1 = E w, G2 = E Lambda w, G3 = E Lambda Lambda' w
        ## we want to compute:
        ## Delta1 = E (A Lambda-a) w, Delta2 = E (A Lambda-a) Lambda' w
        ## where A = stand(IC), a=cent(IC)
        ## hence Delta1 = A G2 - a G1, Delta2 = A G3 - a G2'
        ### otherwise the return value is NULL and we use the standard method

        if(is.null(res))
           return(getMethod("checkIC", signature(IC = "IC",
                              L2Fam = "L2ParamFamily"))(IC,L2Fam, out = out, ..., diagnostic = diagnostic))


        A <- stand(IC);  a <- cent(IC)
        G1 <- res$G1;  G2 <- res$G2;  G3 <- res$G3
        Delta1 <- A%*%G2- a*G1
        Delta2 <- A%*%G3 - a%*%t(G2)
        trafoL <- trafo(L2Fam)
        Delta2 <- Delta2 - trafoL

        Prec <- ceiling(12-round(max(log(abs(trafoL)+1e-14,10)))/2)

    ## PR 20190407: in output in if(out)
		## deleting all digits beyond 1e-12 (as numeric fuzz) --
		## but check for relative accuracy by means of the "size" of the Fisher information
		## measured in by the max(trafo)
        if(out){
            cent.out <- round(Delta1*10^Prec)/10^Prec
            cat("precision of centering:\t", cent.out, "\n")

            oldOps <- options()
            on.exit(do.call(options,oldOps))
            consist.out <- round(Delta2*10^Prec)/10^Prec
            options(digits=5,scipen=-2)
            cat("precision of Fisher information:\n")
            print(consist.out)

            cat("precision of Fisher information - relativ error [%]:\n")
            relconsist.out <- round(Delta2/trafoL*10^(Prec+2))/10^Prec
            class(relconsist.out) <- c("relMatrix",class(consist.out))
            print(relconsist.out)

          if(diagnostic){
               print(attr(res$G1, "diagnostic"))
               print(attr(res$G2, "diagnostic"))
               print(attr(res$G3, "diagnostic"))
          }
        }

        prec <- max(abs(Delta1), abs(Delta2))
        names(prec) <- "maximum deviation"

        if(diagnostic) attr(prec, "diagnostic") <- c(attr(res$G1, "diagnostic"),
           attr(res$G2, "diagnostic"), attr(res$G3, "diagnostic"))
        return(prec)
    })

## make some L2function a pIC at a model
setMethod("makeIC", signature(IC = "ContIC", L2Fam = "L2ParamFamily"),
    function(IC, L2Fam, forceContICMethod = FALSE, ..., diagnostic = FALSE){

        D1 <- L2Fam@distribution
        if( dimension(Domain(IC@Curve[[1]])) != dimension(img(D1)))
            stop("dimension of 'Domain' of 'Curve' != dimension of 'img' of 'distribution' of 'L2Fam'")

        dims <- nrow(trafo(L2Fam))
        if(dimension(IC@Curve) != dims)
           stop("Dimension of IC and parameter must be equal")

        res <- .prepareCheckMakeIC(L2Fam, w = IC@weight, forceContICMethod, ..., diagnostic = diagnostic)

        if(diagnostic &&!is.null(res)){
               print(attr(res$G1, "diagnostic"))
               print(attr(res$G2, "diagnostic"))
               print(attr(res$G3, "diagnostic"))
        }

        ## if it pays off to use symmetry/ to compute integrals in L2deriv space
        ## we compute the following integrals:
        ## G1 = E w, G2 = E Lambda w, G3 = E Lambda Lambda' w
        ## we want to compute:
        ## Delta1 = E (A Lambda-a) w, Delta2 = E (A Lambda-a) Lambda' w
        ## where A = stand(IC), a=cent(IC)
        ## hence Delta1 = A G2 - a G1, Delta2 = A G3 - a G2'
        ### otherwise the return value is NULL and we use the standard method

        if(is.null(res))
           return(getMethod("makeIC", signature(IC = "IC",
                              L2Fam = "L2ParamFamily"))(IC,L2Fam,..., diagnostic = diagnostic))

        G1 <- res$G1;  G2 <- res$G2;  G3 <- res$G3
        trafO <- trafo(L2Fam@param)
        nrvalues <- nrow(trafO)
        dims <- ncol(trafO)

        cent0 <- c(G2/G1)
        stand1 <- trafO%*%distr::solve(G3-cent0%*%t(G2))
        cent1 <- c(stand1%*%cent0)
#        print(list(stand1,stand(IC),cent1,cent(IC)))
        L2.f <- as(diag(nrvalues) %*% L2Fam@L2deriv , "EuclRandVariable")
        D1 <- L2Fam@distribution

        IC1.f <- function(x){ indS <- liesInSupport(D1,x,checkFin=TRUE)
                              Lx <- sapply(x, function(y) evalRandVar(L2.f,y))
                              indS* (stand1%*%Lx-cent1) * weight(IC@weight)(Lx)}

        IC1.l <- vector("list",nrvalues)
        for(i in 1:nrvalues){
            IC1.l[[i]] <- function(x){}
            body(IC1.l[[i]]) <- substitute( c((IC1.s(x))[i,]), list(IC1.s=IC1.f, i=i))
        }
        IC1.c <- EuclRandVariable(Map = IC1.l, Domain = Domain(IC@Curve[[1]]),
                                Range = Reals())

        cIC1 <- new("ContIC")
        cIC1@name <- IC@name
        cIC1@Curve <- EuclRandVarList(IC1.c)
        cIC1@Risks <- IC@Risks
        cIC1@Infos <- IC@Infos
        cIC1@CallL2Fam <- L2Fam@fam.call
        cIC1@clip <- IC@clip
        cIC1@cent <- cent1
        cIC1@stand <- stand1
        cIC1@lowerCase <- IC@lowerCase
        cIC1@neighborRadius <- IC@neighborRadius
        cIC1@weight <- IC@weight
        cIC1@biastype <- IC@biastype
        cIC1@normtype <- IC@normtype
        cIC1@modifyIC <- IC@modifyIC
        addInfo(cIC1) <- c("IC<-",
                           "generated by affine linear trafo to enforce consistency")

        if(diagnostic) attr(cIC1, "diagnostic") <- c(attr(res$G1, "diagnostic"),
           attr(res$G2, "diagnostic"), attr(res$G3, "diagnostic"))

        return(cIC1)
    })

.prepareCheckMakeIC <- function(L2Fam, w, forceContICMethod, ..., diagnostic = FALSE){

        dims <- length(L2Fam@param)
        trafo <- trafo(L2Fam@param)
        nrvalues <- nrow(trafo)

        z.comp <- rep(TRUE,dims)
        A.comp <- matrix(rep(TRUE,dims^2),nrow=dims)
        to.comp.i <- (dims+1)*(dims+2)/2
        to.comp.a <- (dims+1)*nrvalues

        L2deriv <- as(diag(dims) %*% L2Fam@L2deriv, "EuclRandVariable")

        z.comp <- rep(TRUE,dims)
        A.comp <- matrix(TRUE, dims, dims)
#        print(list(z.comp,A.comp))
        # otherwise if  nrvalues > 1 # formerly: trafo == unitMatrix #
        #           may use symmetry info
        if(dims>1){
            comp <- .getComp(L2deriv, L2Fam@distrSymm, L2Fam@L2derivSymm, L2Fam@L2derivDistrSymm)
            z.comp <- comp$"z.comp"
            A.comp <- comp$"A.comp"
            t.comp.i <- sum(z.comp)+sum(A.comp)+1
        }

        if(to.comp.a < to.comp.i && !forceContICMethod) return(NULL)


        res <- .getG1G2G3Stand(L2deriv = L2deriv, Distr = L2Fam@distribution,
                               A.comp = A.comp, z.comp = z.comp, w = w, ...,
                               diagnostic = diagnostic)
        return(res)
}



.getG1G2G3Stand <- function(L2deriv, Distr, A.comp, z.comp, w, ..., diagnostic = FALSE){

        dotsI <- .filterEargs(list(...))
        if(is.null(dotsI$useApply)) dotsI$useApply <- FALSE

        w.fct <- function(x){
            weight(w)(evalRandVar(L2deriv, as.matrix(x)) [,,1])
        }


        integrand2 <- function(x, L2.i){
            return(L2.i(x)*w.fct(x))
        }

        diagn <- if(diagnostic) vector("list", sum(z.comp)+sum(A.comp))
        if(diagnostic) dotsI$diagnostic <- TRUE
        Eargs <- c(list(object = Distr, fun = w.fct), dotsI)
        res1 <- do.call(E,Eargs)

        k <- 0
        nrvalues <- length(L2deriv)
        res2 <- numeric(nrvalues)
        for(i in 1:nrvalues){
            if(z.comp[i]){
                 Eargs <- c(list(object = Distr, fun = integrand2,
                                 L2.i = L2deriv@Map[[i]]), dotsI)
                 res2[i] <- buf <- do.call(E,Eargs)
                 if(diagnostic){k <- k + 1; diagn[[k]] <- attr(buf,"diagnostic")}
            }else{
                res2[i] <- 0
            }
        }
        if(diagnostic) {k1 <- k; attr(res2, "diagnostic") <- diagn[(1:k1)]}
        cent <- res2/res1

        integrandA <- function(x, L2.i, L2.j, i, j){
            return((L2.i(x) - cent[i])*(L2.j(x) - cent[j])*w.fct(x = x))
        }

        nrvalues <- length(L2deriv)
        erg <- matrix(0, ncol = nrvalues, nrow = nrvalues)

        for(i in 1:nrvalues){
            for(j in i:nrvalues){
                if(A.comp[i,j]){
                    Eargs <- c(list(object = Distr, fun = integrandA,
                                   L2.i = L2deriv@Map[[i]],
                                   L2.j = L2deriv@Map[[j]], i = i, j = j), dotsI)
                    erg[i, j] <- buf <- do.call(E,Eargs)
                    if(diagnostic){k <- k + 1; diagn[[k]] <- attr(buf,"diagnostic")}
                }
            }
        }
        erg[col(erg) < row(erg)] <- t(erg)[col(erg) < row(erg)]
        if(diagnostic) {k1 <- k; attr(erg, "diagnostic") <- diagn[-(1:k1)]}

        return(list(G1=res1,G2=res2, G3=erg))
    }





#}
