
setMethod("getStartIC",signature(model = "L2LocationFamily", risk = "interpolRisk"),
           function(model, risk, ...) .getStIC(model, risk, ..., intfct=.getPsi.loc, pkg="ROptEst"))

setMethod("getStartIC",signature(model = "L2ScaleFamily", risk = "interpolRisk"),
           function(model, risk, ...) .getStIC(model, risk, ..., intfct=.getPsi.sca, pkg="ROptEst"))

setMethod("getStartIC",signature(model = "L2LocationScaleFamily", risk = "interpolRisk"),
           function(model, risk, ...) .getStIC(model, risk, ..., intfct=.getPsi.lsc, pkg="ROptEst"))

.getStIC <- function(model,risk, ..., intfct, pkg="ROptEst"){

    mc <- match.call(call = sys.call(sys.parent(1)))
    dots <- match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)$"..."

    gridn <- gsub("\\.","",type(risk))

    nam <- paste(".",gsub("[F,f]amily","",gsub(" ","",name(model))),sep="")
    param1 <- param(model)

    nsng <- character(0)
    famg <- try(getFromNamespace(nam, ns = pkg), silent=TRUE)
    if(!is(famg,"try-error")) nsng <- names(famg)
    if(length(nsng)){
       if(gridn %in% nsng){
          LMref <- famg[[gridn]]
          .modifyIC0 <- function(L2Fam, IC){
                    para <- param(L2Fam)
                    return(intfct(para, LMref, L2Fam, type(risk)))
          }

          .modifyIC <- function(L2Fam,IC, withMakeIC = FALSE, ...){
               psi.0 <- .modifyIC0(L2Fam,IC)
               psi.0@modifyIC <- .modifyIC
               return(psi.0)
          }

          IC0 <- intfct(param1, LMref, model, type(risk))
          IC0@modifyIC <- .modifyIC
          return(IC0)
       }
    }
    mc1 <- as.list(mc)[-1]
    mc1[["risk"]] <- if(type(risk)==".MBRE") asBias() else asMSE()
    mc1[["neighbor"]] <- ContNeighborhood(radius=0.5)
    mc1[["verbose"]] <- FALSE
    if(type(risk)==".MBRE") mc1[["eps"]] <- list(e=40)
    if(type(risk)==".OMSE"){
        n <- length(get("x", envir=parent.frame(2)))
        eps <- list("e" =0.5/sqrt(n), "sqn"= sqrt(n))
        mc1[["eps"]] <- eps
    }
    if(type(risk)==".RMXE"){
        n <- length(get("x", envir=parent.frame(2)))
        eps <- list("eps.lower"=0, "eps.upper"=20, "sqn"= sqrt(n))
        mc1[["eps"]] <- eps
    }
    IC <- do.call(getStartIC, mc1, envir=parent.frame(2))
    return(IC)
}

.getPsi.loc <- function(param, LMref, L2Fam , type){
   LM0 <- LMref
   if(type==".MBRE") LM0 <- .symmIA.gen(LM0)

   return(.mkPsi.gen(LM0,L2Fam,type))
}


.getPsi.sca <- function(param, LMref, L2Fam , type){

   scnam <- scalename(L2Fam)
   xi <- main(param)[scnam] #[["shape"]]

   LM0 <- LMref
   if(type==".MBRE") LM0 <- .symmIA.gen(LM0)

   return(.mkPsi.gen(.xiMkLM(LM0,xi),L2Fam,type))
}


.getPsi.lsc <- function(param, LMref, L2Fam , type){

   scnam <- locscalename(L2Fam)["scale"]
   xi <- if(scnam %in% names(main(param)))
            main(param)[scnam] else nuisance(param)[scnam]

   LM0 <- LMref
   if(type==".MBRE") LM0 <- .symmIA.gen(LM0)

   return(.mkPsi.gen(.xiMkLM(LM0,xi),L2Fam,type))
}

.symmIA.gen <- function(LMset){
         LM0 <- LMset
         aa <- LMset[["a"]]
         zi <- LMset[["aw"]]
         Aa <- LMset[["A"]]
         Ai <- LMset[["Aw"]]
         LM0$Aw <- LM0$A <- (Aa+Ai+t(Ai)+t(Aa))/4
         ai <- Ai %*% zi
         LM0$a <- (ai+aa)/2
         LM0$aw <- distr::solve(LM0$A, LM0$a)
}

.xiMkLM <- function(LMset,xi){
   LMset$b  <- xi*LMset[["b"]]
   LMset$a  <- xi*LMset[["a"]]
   LMset$aw <- LMset[["aw"]]/xi
   LMset$A  <- xi^2*LMset[["A"]]
   LMset$Aw <- xi^2*LMset[["Aw"]]
   return(LMset)
}

.mkPsi.gen <- function(res0,L2Fam, type){
   normt <- NormType()
   biast <- symmetricBias()
   nb <- ContNeighborhood(radius=0.5)
   ICT <- paste("optimally robust IC for", switch(type,
                      ".OMSE"="maxMSE",".RMXE"="RMX", ".MBRE"="maxBias"))
   riskT <- if(type!=".MBRE") "asGRisk" else "asBias"

   w <- new("HampelWeight")
      stand(w) <- res0$Aw
      cent(w) <- res0$aw
      clip(w) <- res0$b

   if(type!=".MBRE"){
        weight(w) <- getweight(w, neighbor = nb, biastype = biast,
                          normW = normt)
   }else weight(w) <- minbiasweight(w, neighbor = nb, biastype = biast,
                          normW = normt)

   res <- list(a = res0$a, A = res0$A, b = res0$b, d = 0*res0$a,
               normtype = normt, biastype = biast, w = w,
               info = c("optIC", ICT), risk = list(),
               modifyIC = NULL)

   IC <- generateIC(nb, L2Fam, res)
   return(IC)
}
