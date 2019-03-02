###############################################################################
## k-step estimator
###############################################################################

.fix.scalename <- function(obj, scalename, estname){
        hasdim <- !is.null(dim(obj))
        n.obj <- if(hasdim) rownames(obj) else names(obj)
        if(!is.na(scalename)) if(scalename!="") {
           if((! (scalename %in% estname)) && "scale" %in% estname)
                estname[estname=="scale"] <- scalename

           if((! (scalename%in% n.obj)) && "scale" %in% n.obj){
              n.obj[n.obj=="scale"] <- scalename
              if(hasdim) rownames(obj) <- n.obj else names(obj) <- n.obj
           }else{
              if(length(n.obj)==0) n.obj <- rep("", length(estname))
              if(all(n.obj=="")) {
              if(hasdim) rownames(obj) <- estname else names(obj) <- estname
              }
           }
        }
        return(obj)
}

setMethod("neighborRadius","ANY",function(object)NA)

.addTime <- function(timold,namenew){
   tim <- rbind(timold,proc.time())
   rownames(tim) <- c(rownames(timold),namenew)
   return(tim)
}

.ensureDim2 <- function(x){
    d <- dim(x)
    if(length(d)==3L && d[3]==1L) dim(x) <- d[1:2]
    if(length(d)==4L && d[2]==1L && d[4] == 1L) dim(x) <- d[c(1,3)]
    x }

### taken from: base::system.time ::
ppt <- function(y) {
        if (!is.na(y[4L]))
            y[1L] <- y[1L] + y[4L]
        if (!is.na(y[5L]))
            y[2L] <- y[2L] + y[5L]
        paste(formatC(y[1L:3L]), collapse = " ")
}


### no dispatch on top layer -> keep product structure of dependence
kStepEstimator <- function(x, IC, start = NULL, steps = 1L,
                           useLast = getRobAStBaseOption("kStepUseLast"),
                           withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
                           IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
                           withICList = getRobAStBaseOption("withICList"),
                           withPICList = getRobAStBaseOption("withPICList"),
                           na.rm = TRUE, startArgList = NULL, ...,
                           withLogScale = TRUE, withEvalAsVar = TRUE,
                           withMakeIC = FALSE, E.argList = NULL,
                           diagnostic = FALSE){

        time <- proc.time()
        on.exit(message("Timing stopped at: ", ppt(proc.time() - time)))
## save call
        es.call <- match.call()
        es.call[[1]] <- as.name("kStepEstimator")

        if(is.null(E.argList)) E.argList <- list()
        if(is.null(E.argList$useApply)) E.argList$useApply <- FALSE
        diagn <- NULL
        if(diagnostic){
           E.argList$diagnostic <- TRUE
           diagn <- list()
        }
        if(missing(IC.UpdateInKer)) IC.UpdateInKer <- NULL

## get some dimensions
        L2Fam <- eval(CallL2Fam(IC))
        sytm <- rbind(time,"eval(CallL2Fam(IC))"=proc.time())
        colnames(sytm) <- names(time)
        Param <- param(L2Fam)

        tf <- trafo(L2Fam,Param)
        Dtau <- tf$mat
        trafoF <- tf$fct

        hasnodim.main <- is.null(dim(main(L2Fam)))
        hasnodim.nuis <- is.null(dim(nuisance(L2Fam)))

        p <- nrow(Dtau)
        k <- ncol(Dtau)

        lmx <- length(main(L2Fam))
        lnx <- length(nuisance(L2Fam))
        idx <- 1:lmx
        nuis.idx <- if(lnx) lmx + 1:lnx else NULL

        var.to.be.c <- ("asCov" %in% names(Risks(IC))) | (lnx == 0)

        fixed <- fixed(L2Fam)

## names of the estimator components
        par.names  <- names(main(L2Fam))
        if(lnx)
           par.names  <- c(par.names, names(nuisance(L2Fam)) )
        est.names   <- if(.isUnitMatrix(Dtau)) par.names else rownames(Dtau)
        u.est.names <- par.names

## check input
        if(!is.integer(steps))
          steps <- as.integer(steps)
        if(steps < 1)
            stop("steps needs to be a positive integer")
        if(! is(IC, "IC"))
           stop("Argument 'IC' must be of class 'IC'")

### transform if necessary
        x0 <- x
        x0 <- if(is.numeric(x) && ! is.matrix(x)) {
                x0 <- as.matrix(x)
                }
        completecases <- complete.cases(x0)
        if(na.rm) x0 <- na.omit(x0)

        if(missing(start)||is.null(start))
           start <- L2Fam@startPar

### use dispatch here  (dispatch only on start)
        #a.var <- if( is(start, "Estimate")) asvar(start) else NULL

        IC.UpdateInKer.0 <- if(is(start,"ALEstimate")) pIC(start) else NULL
        sytm <- .addTime(sytm,"pIC(start)")
        ## pIC(start) instead of start@pIC to potentially eval a call

        force(startArgList)

        start.val <- kStepEstimator.start(start, x=x0, nrvalues = k,
                         na.rm = na.rm, L2Fam = L2Fam,
                         startList = startArgList)
        sytm <- .addTime(sytm,"kStepEstimator.start")

### use Logtransform here in scale models
        sclname <- ""
        if(is(L2Fam, "L2ScaleUnion")) sclname <- scalename(L2Fam)
        logtrf <- is(L2Fam, "L2ScaleUnion") &
                     withLogScale & sclname %in% names(start.val)
### a starting value in k-space
#        print(start.val)
        u.theta <- start.val
        theta <- if(is(start.val,"Estimate")) estimate(start.val)
                 else trafoF(u.theta[idx])$fval
        u.start.val <- matrix(start.val,ncol=1)
        start.val <- matrix(theta,ncol=1)
        rownames(u.start.val) <- u.est.names
        rownames(start.val) <- est.names
#        print(theta)
        theta <- .fix.scalename(theta, sclname, est.names)
#        print(theta)
#        print(u.theta)
        u.theta <- .fix.scalename(u.theta, sclname, u.est.names)
#        print(u.theta)

### shall intermediate IC's / pIC's be stored?
        pICList <- if(withPICList) vector("list", steps) else NULL
        ICList  <- if(withICList)  vector("list", steps) else NULL

        cvar.fct <- function(L2, IC, dim, dimn =NULL){
                Eres <- matrix(NA,dim,dim)
                if(!is.null(dimn)) dimnames(Eres) <- dimn
                ICM <- as(diag(k)%*%IC, "EuclRandVariable")@Map
                for(i in 1: dim)
                      Eres[i,i] <- E(L2@distribution,
                           fun = function(x) ICM[[i]](x)^2,
                           useApply = FALSE)
                if(dim>1){
                  for(i in 1: (dim-1)){
                    for(j in (i+1): dim)
                        Eres[j,i] <- Eres[i,j] <- E(L2@distribution,
                           fun = function(x) ICM[[i]](x)*ICM[[j]](x),
                           useApply = FALSE)
                  }
                }
                return(Eres)}


        updStp <- 0
        ### update - function
        updateStep <- function(u.theta, theta, IC, L2Fam, Param,
                               withPreModif = FALSE,
                               withPostModif = TRUE, with.u.var = FALSE,
                               withEvalAsVar.0 = FALSE
                               ){

                updStp <<- updStp + 1
                if(withPreModif){
                   main(Param)[] <- .deleteDim(u.theta[idx])
#                   print(Param)
                   if (lnx) nuisance(Param)[] <- .deleteDim(u.theta[nuis.idx])
#                   print(Param)
#                   print(L2Fam)
                   L2Fam <- modifyModel(L2Fam, Param,
                               .withL2derivDistr = L2Fam@.withEvalL2derivDistr)
                   mmPreNm <- paste("modifyModel-PreModif-",updStp)
                   sytm <<- .addTime(sytm,mmPreNm)
                   if(diagnostic) diagn[[mmPreNm]] <<- attr(L2Fam,"diagnostic")
# print(L2Fam)

                   modifyICargs <- c(list(L2Fam, IC, withMakeIC = FALSE), E.argList)
                   IC <- do.call(modifyIC(IC),modifyICargs)
                   mmPreICNm <- paste("modifyIC-PreModif-",updStp)
                   sytm <<- .addTime(sytm,mmPreICNm)
                   if(diagnostic) diagn[[mmPreICNm]] <<- attr(IC,"diagnostic")
                   if(steps==1L && withMakeIC){
                      makeICargs <- list(IC, L2Fam, diagnostic=diagnostic, E.argList=E.argList)
                      IC <- do.call(makeIC, makeICargs)
                      mmPreMkICNm <- paste("modifyIC-makeIC-",updStp)
                      sytm <<- .addTime(sytm,mmPreMkICNm)
                      if(diagnostic) diagn[[mmPreMkICNm]] <<- attr(IC,"diagnostic")
                    }
                }

                IC.c <- as(diag(p) %*% IC@Curve, "EuclRandVariable")
                sytm <<- .addTime(sytm,paste("IC.c <- as(diag(p) %*%-",updStp))

#                print(theta)
                tf <- trafo(L2Fam, Param)
                Dtau <- tf$mat
                IC.tot.0 <- NULL
#                print(Dtau)
                if(!.isUnitMatrix(Dtau)){
                     Dminus <- distr::solve(Dtau, generalized = TRUE)
                     projker <- diag(k) - Dminus %*% Dtau

                     IC.tot1 <- Dminus %*% IC.c
#                     IC.tot2 <- 0 * IC.tot1
                     IC.tot2.isnull <- TRUE

                     if(sum(diag(projker))>0.5 && ### is EM-D^-D != 0 (i.e. rk D<p)
                        withUpdateInKer){
                            if(!is.null(IC.UpdateInKer)&&!is(IC.UpdateInKer,"IC"))
                               warning("'IC.UpdateInKer' is not of class 'IC'; we use default instead.")
                            if(is.null(IC.UpdateInKer)){
                                 getBoundedICargs <- list(L2Fam, D = projker, diagnostic=diagnostic,E.argList=E.argList)
                                 IC.tot2 <- do.call(getBoundedIC, getBoundedICargs)
                                 mmgtBDICNm <- paste("getBoundedIC-",updStp)
                                 sytm <<- .addTime(sytm,mmgtBDICNm)
                                 if(diagnostic) diagn[[mmgtBDICNm]] <<- attr(IC.tot2,"diagnostic")
                            }else{
                                 IC.tot2 <- as(projker %*% IC.UpdateInKer@Curve, "EuclRandVariable")
                                 mmgtAsPrICNm <- paste("IC.tot2<-as(projker...-",updStp)
                                 sytm <<- .addTime(sytm,mmgtAsPrICNm)
                                 if(diagnostic) diagn[[mmgtAsPrICNm]] <<- attr(IC.tot2,"diagnostic")
                            }
                            IC.tot2.isnull <- FALSE
                            IC.tot.0 <- IC.tot1 + IC.tot2
                     }else{ if(is.null(IC.UpdateInKer.0)){
                               IC.tot.0 <- NULL
                            }else{
                                if(is.call(IC.UpdateInKer.0))
                                   IC.UpdateInKer.0 <- eval(IC.UpdateInKer.0)
                                sytm <<- .addTime(sytm,paste("eval(IC.UpdateInKer.0)-",updStp))
                                IC.tot.0 <- IC.tot1 + as(projker %*%
                                         IC.UpdateInKer.0@Curve,
                                                "EuclRandVariable")
                                sytm <<- .addTime(sytm,paste("IC.tot.0 <- IC.tot1 + as(proj-",updStp))
                            }
                     }
                     IC.tot <- IC.tot1
                     if(!IC.tot2.isnull) IC.tot <- IC.tot1 + IC.tot2
                     indS <- liesInSupport(distribution(L2Fam),x0,checkFin=TRUE)
                     correct <- rowMeans(t(t(.ensureDim2(evalRandVar(IC.tot, x0)))*indS), na.rm = na.rm)
                     sytm <<- .addTime(sytm,paste("Dtau-not-Unit:correct <- rowMeans-",updStp))
                     iM <- is.matrix(u.theta)
                     names(correct) <- if(iM) rownames(u.theta) else names(u.theta)
                     if(logtrf){
                        scl <- if(iM) u.theta[sclname,1] else u.theta[sclname]
                        u.theta <- u.theta + correct
                        if(iM) u.theta[sclname,1] <- scl * exp(correct[sclname]/scl) else
                               u.theta[sclname] <- scl * exp(correct[sclname]/scl)
                     }else u.theta <- u.theta + correct

                     theta <- (tf$fct(u.theta[idx]))$fval
                }else{
                     indS <- liesInSupport(distribution(L2Fam),x0,checkFin=TRUE)
                     correct <- rowMeans(t(t(.ensureDim2(evalRandVar(IC.c, x0)))*indS), na.rm = na.rm)
                     sytm <<- .addTime(sytm,paste("Dtau=Unit:correct <- rowMeans-",updStp))
                     iM <- is.matrix(theta)
#                     print(sclname)
#                     print(names(theta))
#                     print(str(theta))
                     names(correct) <- if(iM) rownames(theta) else names(theta)
                     if(logtrf){
                        scl <- if(iM) theta[sclname,1] else theta[sclname]
                        theta <- theta + correct
                        if(iM) theta[sclname,1] <- scl * exp(correct[sclname]/scl) else
                               theta[sclname] <- scl * exp(correct[sclname]/scl)
                     }else{
                        theta <- theta + correct
                     }
                     IC.tot <- IC.c
                     u.theta <- theta
                }

                var0 <- u.var <- NULL
                if(with.u.var){
                   cnms <-  if(is.null(names(u.theta))) colnames(Dtau) else names(u.theta)
                   if(!is.null(IC.tot.0)){
                      u.var <- substitute(do.call(cfct, args = list(L2F0, IC0,
                                   dim0, dimn0)), list(cfct = cvar.fct,
                                   L2F0 = L2Fam, IC0 = IC.tot.0, dim0 = k,
                                   dimn0 = list(cnms,cnms)))
                      sytm <<- .addTime(sytm,paste("u.var-",updStp))
                      if(withEvalAsVar.0){
                         u.var <- eval(u.var)
                         uvEvnm <- paste("u.var-eval-",updStp)
                         sytm <<- .addTime(sytm,uvEvnm)
                         if(diagnostic) diagn[[uvEvnm]] <<- attr(u.var,"diagnostic")
                      }
                   }
                   if(!var.to.be.c){
                      var0 <- substitute(do.call(cfct, args = list(L2F0, IC0,
                                   dim0, dimn0)), list(cfct = cvar.fct,
                                   L2F0 = L2Fam, IC0 = IC.c, dim0 = p))
                      sytm <<- .addTime(sytm,paste("var0-",updStp))
                      if(withEvalAsVar.0) {
                         var0 <- eval(var0)
                         vEvnm <- paste("var0-eval-",updStp)
                         sytm <<- .addTime(sytm,paste("var0-eval-",updStp))
                         if(diagnostic) diagn[[vEvnm]] <<- attr(var0,"diagnostic")
                      }
                   }
                }
                if(withPostModif){
                   main(Param)[] <- .deleteDim(u.theta[idx])
                   if (lnx) nuisance(Param)[] <- .deleteDim(u.theta[nuis.idx])
                   L2Fam <- modifyModel(L2Fam, Param,
                               .withL2derivDistr = L2Fam@.withEvalL2derivDistr)
                   mmPostNm <- paste("modifyModel-PostModif-",updStp)
                   sytm <<- .addTime(sytm,mmPostNm)
                   if(diagnostic) diagn[[mmPostNm]] <<- attr(L2Fam,"diagnostic")

                   modifyICargs <- c(list(L2Fam, IC, withMakeIC = withMakeIC), E.argList)
                   IC <- do.call(modifyIC(IC),modifyICargs)
                   mmPostICNm <- paste("modifyIC-PostModif-",updStp)
                   sytm <<- .addTime(sytm,mmPostICNm)
                   if(diagnostic) diagn[[mmPostICNm]] <<- attr(IC,"diagnostic")
                }

                li <- list(IC = IC, Param = Param, L2Fam = L2Fam,
                            theta = theta, u.theta = u.theta, u.var = u.var,
                            var = var0, IC.tot = IC.tot, IC.c = IC)
                sytm <<- .addTime(sytm,paste("li <- list(IC = IC,...-",updStp))
                return(li)
        }

        Infos <- matrix(c("kStepEstimator",
                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
                        ncol = 2)
        colnames(Infos) <- c("method", "message")
        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE

        ### iteration

        ksteps  <- matrix(0,ncol=steps, nrow = p)
        uksteps <- matrix(0,ncol=steps, nrow = k)
        rownames(ksteps) <- est.names
        rownames(uksteps) <- u.est.names
        if(!is(modifyIC(IC), "NULL") ){
           for(i in 1:steps){
               if(i>1){
                  IC <- upd$IC
                  L2Fam <- upd$L2Fam
                  if((i==steps)&&withMakeIC){
                      makeICargs <- list(IC, L2Fam, diagnostic=diagnostic, E.argList=E.argList)
                      IC <- do.call(makeIC, makeICargs)
                      mkICnm <- paste("makeIC-",i)
                      sytm <- .addTime(sytm,mkICnm)
                      if(diagnostic) diagn[[mkICnm]] <- attr(IC,"diagnostic")
                  }

                  Param <- upd$Param
                  tf <- trafo(L2Fam, Param)
                  withPre <- FALSE
               }else withPre <- TRUE
               upd <- updateStep(u.theta,theta,IC, L2Fam, Param,
                                 withPreModif = withPre,
                                 withPostModif = (steps>i) | useLast,
                                 with.u.var = (i==steps),
                                 withEvalAsVar.0 = (i==steps))
#               print(upd$u.theta); print(upd$theta)
               uksteps[,i] <- u.theta <- upd$u.theta
#               print(str(upd$theta))
#               print(nrow(ksteps))
               ksteps[,i] <- theta <- upd$theta
               if(withICList)
                  ICList[[i]] <- .fixInLiesInSupport(
                                  new("InfluenceCurve",
                                      name = paste(gettext("(total) IC in step"),i),
                                      Risks = list(),
                                      Infos = matrix(c("",""),ncol=2),
                                      Curve =  EuclRandVarList(upd$IC.tot)),
                                  distr = distribution(upd$L2Fam))
               sytm <- .addTime(sytm,paste("ICList-",i))
               if(withPICList)
                  pICList[[i]] <- .fixInLiesInSupport(upd$IC.c,distribution(upd$L2Fam))
               u.var <- upd$u.var
               var0 <- upd$var
           }
           if(withICList) ICList <- new("pICList",ICList)
           if(withPICList) pICList <- new("pICList",pICList)
           if(useLast){
              IC <- upd$IC
              L2Fam <- upd$L2Fam
              Param <- upd$Param
              tf <- trafo(L2Fam, Param)
              Infos <- rbind(Infos, c("kStepEstimator",
               "computation of IC, trafo, asvar and asbias via useLast = TRUE"))
              if(withMakeIC){
                  makeICargs <- list(IC, L2Fam, diagnostic=diagnostic, E.argList=E.argList)
                  IC <- do.call(makeIC, makeICargs)
                  mkICULnm <- paste("makeIC-useLast")
                  sytm <- .addTime(sytm,mkICULnm)
                  if(diagnostic) diagn[[mkICULnm]] <- attr(IC,"diagnostic")
              }
           }else{
              Infos <- rbind(Infos, c("kStepEstimator",
               "computation of IC, trafo, asvar and asbias via useLast = FALSE"))
           }
        }else{
           if(steps > 1)
              stop("slot 'modifyIC' of 'IC' is 'NULL'!")
           upd <- updateStep(u.theta,theta,IC, L2Fam, Param,withPreModif = FALSE,
                               withPostModif = TRUE)
           theta <- upd$theta
           u.theta <- upd$u.theta
           var0 <- upd$var
           u.var <- upd$u.var
           ksteps <- NULL
           uksteps <- NULL
           if(useLast){
              warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
                     is filled with some function!")
              Infos <- rbind(Infos, c("kStepEstimator",
                          "slot 'modifyIC' of 'IC' was not filled!"))
           }else
            Infos <- rbind(Infos, c("kStepEstimator",
            "computation of IC, asvar and asbias via useLast = FALSE"))
        }

        ## if non-trivial trafo: info on how update was done
#        print(IC@Risks$asCov)
#        print(Risks(IC)$asCov)

        if(! .isUnitMatrix(trafo(L2Fam)))
             Infos <- rbind(Infos, c("kStepEstimator",
                            paste("computation of IC",
                                   ifelse(withUpdateInKer,"with","without") ,
                                   "modification in ker(trafo)")))

        ## some risks
#        print(list(u.theta=u.theta,theta=theta,u.var=u.var,var=var0))
        if(var.to.be.c){
           if("asCov" %in% names(Risks(IC))){
              asVar <- if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
                       Risks(IC)$asCov else Risks(IC)$asCov$value
           }else{
                getRiskICasVarArgs <- list(IC, risk = asCov(), withCheck = FALSE,
                                        diagnostic=diagnostic, E.argList = E.argList)
                riskAsVar <- do.call(getRiskIC, getRiskICasVarArgs)
                asVar <- riskAsVar$asCov$value
                sytm <- .addTime(sytm,"getRiskIC-Var")
                if(diagnostic) diagn[["getRiskICVar"]] <- attr(asVar,"diagnostic")
           }

        }else asVar <- var0
#        print(asVar)
        if("asBias" %in% names(Risks(IC))){
                if(length(Risks(IC)$asBias) == 1)
                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
                else
                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
                if(is.na(asBias)) asBias <- NULL
        }else{
                if(is(IC, "HampIC")){
                    r <- neighborRadius(IC)
                    asBias <- r*getRiskIC(IC, risk = asBias(),
                                          neighbor = neighbor(IC), withCheck = FALSE)$asBias$value
                    sytm <- .addTime(sytm,"getRiskIC-Bias")
                }else{
                    asBias <- NULL
                }
        }

        if(hasnodim.main) theta <- .deleteDim(theta)
        if(hasnodim.nuis) u.theta <- .deleteDim(u.theta)
        names(theta) <- est.names
        names(u.theta) <- u.est.names

        if(lnx){
          nms.theta.idx <- est.names[idx]

          theta <- theta[idx]
#          print(asVar);print(idx)
          asVar <- asVar[idx,idx,drop=FALSE]
#          print(asVar)
          names(theta) <- nms.theta.idx
          dimnames(asVar) <- list(nms.theta.idx, nms.theta.idx)
        }

        IC <- .fixInLiesInSupport(IC, distribution(L2Fam))


        estres <- new("kStepEstimate", estimate.call = es.call,
                name = paste(steps, "-step estimate", sep = ""),
                estimate = theta, samplesize = nrow(x0), asvar = asVar,
                trafo = tf, fixed = fixed, nuis.idx = nuis.idx,
                untransformed.estimate = u.theta, completecases = completecases,
                untransformed.asvar = u.var, asbias = asBias, pIC = IC,
                steps = steps, Infos = Infos, start = start,
                startval = start.val, ustartval = u.start.val, ksteps = ksteps,
                uksteps = uksteps, pICList = pICList, ICList = ICList)
        sytm <- .addTime(sytm,"new('kStepEstimate'...")
        estres <- .checkEstClassForParamFamily(L2Fam,estres)

        attr(estres,"timings") <- apply(sytm,2,diff)
        if(diagnostic){
           attr(estres,"diagnostic") <- diagn
           if(!is.null(diagn))
              class(attr(estres,"diagnostic")) <- "DiagnosticClass"
        }
        on.exit()
        return(estres)

}
#  (est1.NS <- kStepEstimator(x, IC2.NS, est0, steps = 1))

#### old method:

# setMethod("kStepEstimator", signature(x = "numeric",
#                                      IC = "IC",
#                                      start = "numeric"),
#    function(x, IC, start, steps = 1L, useLast = getRobAStBaseOption("kStepUseLast")){
#        es.call <- match.call()
#        es.call[[1]] <- as.name("kStepEstimator")
#        if(!is.integer(steps))
#          steps <- as.integer(steps)
#        if(steps < 1)
#            stop("steps needs to be a positive integer")
#
#        nrvalues <- dimension(IC@Curve)
#        if(is.list(start)) start <- unlist(start)
#        if(nrvalues != length(start))
#            stop("dimension of 'start' != dimension of 'Curve'")
#
#        res <- start + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)
#
#        L2Fam <- eval(CallL2Fam(IC))
#        Infos <- matrix(c("kStepEstimator",
#                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
#                        ncol = 2)
#        colnames(Infos) <- c("method", "message")
#        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
#
#        if(steps == 1){
#            if(useLast && !is(modifyIC(IC), "NULL") ){
#                newParam <- param(L2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(L2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                if(useLast && is(modifyIC(IC), "NULL")){
#                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
#                             is filled with some function!")
#                    Infos <- rbind(Infos, c("kStepEstimator",
#                                            "slot 'modifyIC' of 'IC' was not filled!"))
#                }
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = length(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = 1L, Infos = Infos))
#        }else{
#            if(is(modifyIC(IC), "NULL"))
#                stop("slot 'modifyIC' of 'IC' is 'NULL'!")
#            for(i in 2:steps){
#                start <- res
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- start
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                res <- start + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)
#            }
#            if(useLast){
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = length(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = steps, Infos = Infos))
#        }
#    })
#
#setMethod("kStepEstimator", signature(x = "matrix",
#                                      IC = "IC",
#                                      start = "numeric"),
#    function(x, IC, start, steps = 1, useLast = getRobAStBaseOption("kStepUseLast")){
#        es.call <- match.call()
#        es.call[[1]] <- as.name("kStepEstimator")
#        if(!is.integer(steps))
#          steps <- as.integer(steps)
#        if(steps < 1)
#            stop("steps needs to be a positive integer")
#
#        nrvalues <- dimension(IC@Curve)
#        if(is.list(start)) start <- unlist(start)
#        if(nrvalues != length(start))
#            stop("dimension of 'start' != dimension of 'Curve'")
#        if(ncol(x) != IC@Curve[[1]]@Domain@dimension)
#            stop("'x' has wrong dimension")
#
#        res <- start + rowMeans(evalIC(IC, x), na.rm = TRUE)
#
#        L2Fam <- eval(CallL2Fam(IC))
#        Infos <- matrix(c("kStepEstimator",
#                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
#                        ncol = 2)
#        colnames(Infos) <- c("method", "message")
#        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
#
#        if(steps == 1){
#            if(useLast && !is(modifyIC(IC), "NULL") ){
#                newParam <- param(L2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(L2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                if(useLast && is(modifyIC(IC), "NULL")){
#                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
#                             is filled with some function!")
#                    Infos <- rbind(Infos, c("kStepEstimator",
#                                            "slot 'modifyIC' of 'IC' was not filled!"))
#                }
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#             }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = ncol(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = 1L, Infos = Infos))
#        }else{
#            if(is(modifyIC(IC), "NULL"))
#                stop("slot 'modifyIC' of 'IC' is 'NULL'!")
#            for(i in 2:steps){
#                start <- res
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- start
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                res <- start + rowMeans(evalIC(IC, x), na.rm = TRUE)
#            }
#            if(useLast){
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = ncol(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = steps, Infos = Infos))
#        }
#    })
#setMethod("kStepEstimator", signature(x = "numeric",
#                                      IC = "IC",
#                                      start = "Estimate"),
#    function(x, IC, start, steps = 1, useLast = getRobAStBaseOption("kStepUseLast")){
#        es.call <- match.call()
#        es.call[[1]] <- as.name("kStepEstimator")
#        if(!is.integer(steps))
#          steps <- as.integer(steps)
#        if(steps < 1)
#            stop("steps needs to be a positive integer")
#
#        nrvalues <- dimension(IC@Curve)
#        start0 <- estimate(start)
#        if(is.list(start0)) start0 <- unlist(start0)
#        if(nrvalues != length(start0))
#            stop("dimension of slot 'estimate' of 'start' != dimension of 'Curve'")
#
#        res <- start0 + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)
#
#        L2Fam <- eval(CallL2Fam(IC))
#        Infos <- matrix(c("kStepEstimator",
#                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
#                        ncol = 2)
#        colnames(Infos) <- c("method", "message")
#        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
#
#        if(steps == 1){
#            if(useLast && !is(modifyIC(IC), "NULL") ){
#                newParam <- param(L2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(L2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                if(useLast && is(modifyIC(IC), "NULL")){
#                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
#                             is filled with some function!")
#                    Infos <- rbind(Infos, c("kStepEstimator",
#                                            "slot 'modifyIC' of 'IC' was not filled!"))
#                }
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = length(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = 1L, Infos = Infos))
#        }else{
#            if(is(modifyIC(IC), "NULL"))
#                stop("slot 'modifyIC' of 'IC' is 'NULL'!")
#            for(i in 2:steps){
#                start0 <- res
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- start0
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                res <- start0 + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)
#            }
#            if(useLast){
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                 Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = length(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = steps, Infos = Infos))
#        }
#    })
#setMethod("kStepEstimator", signature(x = "matrix",
#                                      IC = "IC",
#                                      start = "Estimate"),
#    function(x, IC, start, steps = 1, useLast = getRobAStBaseOption("kStepUseLast")){
#        es.call <- match.call()
#        es.call[[1]] <- as.name("kStepEstimator")
#        if(!is.integer(steps))
#          steps <- as.integer(steps)
#        if(steps < 1)
#            stop("steps needs to be a positive integer")
#
#        nrvalues <- dimension(IC@Curve)
#        start0 <- estimate(start)
#        if(is.list(start0)) start0 <- unlist(start0)
#        if(nrvalues != length(start0))
#            stop("dimension of slot 'estimate' of 'start' != dimension of 'Curve'")
#        if(ncol(x) != IC@Curve[[1]]@Domain@dimension)
#            stop("'x' has wrong dimension")
#
#        res <- start0 + rowMeans(evalIC(IC, x), na.rm = TRUE)
#
#        L2Fam <- eval(CallL2Fam(IC))
#        Infos <- matrix(c("kStepEstimator",
#                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
#                        ncol = 2)
#        colnames(Infos) <- c("method", "message")
#        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
#
#        if(steps == 1){
#            if(useLast && !is(modifyIC(IC), "NULL") ){
#                newParam <- param(L2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(L2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                if(useLast && is(modifyIC(IC), "NULL")){
#                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
#                             is filled with some function!")
#                    Infos <- rbind(Infos, c("kStepEstimator",
#                                            "slot 'modifyIC' of 'IC' was not filled!"))
#                }
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = ncol(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = 1L, Infos = Infos))
#        }else{
#            if(is(modifyIC(IC), "NULL"))
#                stop("slot 'modifyIC' of 'IC' is 'NULL'!")
#            for(i in 2:steps){
#                start0 <- res
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- start0
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                res <- start0 + rowMeans(evalIC(IC, x), na.rm = TRUE)
#            }
#            if(useLast){
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = ncol(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = steps, Infos = Infos))
#        }
#    })

