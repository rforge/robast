.is.na.Psi <- function(param, fct, nam = "shape"){
   xi <- main(param)[nam]
   return(is.na(fct[[1]](xi)))
}
.getPsi <- function(param, fct, L2Fam , type){

   scshnm <- scaleshapename(L2Fam)
   shnam <- scshnm["shape"]
   scnam <- scshnm["scale"]
   xi <- main(param)[shnam] #[["shape"]]
   beta <- main(param)[scnam] #[scaleshapename(model)["scale"]]

   #print(param)
   #L2deriv <- L2Fam@L2deriv # .fct(param)
   #print(get("tr",environment(get("Lambda1", environment(L2deriv[[1]]@Map[[1]])))))
   #print(get("k",environment(get("Lambda1", environment(L2deriv[[1]]@Map[[1]])))))
   #print(get("sc",environment(get("Lambda1", environment(L2deriv[[1]]@Map[[1]])))))

   .dbeta <- diag(c(beta,1)); .dbeta1 <- diag(c(1/beta,1))
   b <- fct[[1]](xi)

   aa <-  c(fct[[2]](xi),fct[[3]](xi))
   zi <-  c(fct[[4]](xi),fct[[5]](xi))
   am <- mean(c(fct[[7]](xi),fct[[8]](xi)))
   Aa <- matrix(c(fct[[6]](xi),am,am,fct[[9]](xi)),2,2)
   am <- mean(c(fct[[11]](xi),fct[[12]](xi)))
   Ai <- matrix(c(fct[[10]](xi),am,am,fct[[13]](xi)),2,2)
   if(type==".MBRE"){
      ai <- Ai %*% zi
      Am <- (Ai+Aa)/2; Ai <- Aa <- Am
      am <- (ai+aa)/2; ai <- aa <- am
      zi <- distr::solve(Ai,ai)
   }
   a <-  c(.dbeta%*%aa)
   aw <- c(.dbeta1%*%zi)
   A <-  .dbeta%*%Aa%*%.dbeta
   Aw <- Ai%*%.dbeta

   normt <- NormType()
   biast <- symmetricBias()
   nb <- ContNeighborhood(radius=0.5)
   ICT <- paste("optimally robust IC for", switch(type,
                      ".OMSE"="maxMSE",".RMXE"="RMX", ".MBRE"="maxBias"))
   riskT <- if(type!=".MBRE") "asGRisk" else "asBias"

   w <- new("HampelWeight")
      stand(w) <- Aw
      cent(w) <- aw
      clip(w) <- b

   if(type!=".MBRE"){
        weight(w) <- getweight(w, neighbor = nb, biastype = biast,
                          normW = normt)
   }else weight(w) <- minbiasweight(w, neighbor = nb, biastype = biast,
                          normW = normt)
   Risk <- list(asBias = list(value = b, biastype = biast,
                                       normtype = normt,
                                       neighbortype = class(nb)))

   res <- list(a = a, A = A, b = b, d = 0*a,
               normtype = normt, biastype = biast, w = w,
               info = c("optIC", ICT), risk = Risk,
               modifyIC = NULL)


   IC <- generateIC(nb, L2Fam, res)
   return(IC)
}


.getPsi.wL <- function(param, fct, L2Fam , type){

   scshnm <- scaleshapename(L2Fam)
   shnam <- scshnm["shape"]
   scnam <- scshnm["scale"]
   xi <- main(param)[shnam] #[["shape"]]
   beta <- main(param)[scnam] #[scaleshapename(model)["scale"]]

   #print(param)
   #L2deriv <- L2Fam@L2deriv # .fct(param)
   #print(get("tr",environment(get("Lambda1", environment(L2deriv[[1]]@Map[[1]])))))
   #print(get("k",environment(get("Lambda1", environment(L2deriv[[1]]@Map[[1]])))))
   #print(get("sc",environment(get("Lambda1", environment(L2deriv[[1]]@Map[[1]])))))

   .dbeta <- diag(c(beta,beta,1)); .dbeta1 <- diag(c(1/beta,1/beta,1))
   b <- fct[[1]](xi)
   aa <-  c(fct[[2]](xi),fct[[3]](xi),fct[[4]](xi))
   zi <-  c(fct[[5]](xi),fct[[6]](xi),fct[[7]](xi))
   am1 <- mean(c(fct[[9]](xi),fct[[11]](xi)))
   am2 <- mean(c(fct[[10]](xi),fct[[14]](xi)))
   am3 <- mean(c(fct[[13]](xi),fct[[15]](xi)))
   Aa <- matrix(c(fct[[8]](xi),am1,am2,am1,fct[[12]](xi),am3,am2,am3,fct[[16]](xi)),3,3)
   am1 <- mean(c(fct[[18]](xi),fct[[20]](xi)))
   am2 <- mean(c(fct[[19]](xi),fct[[23]](xi)))
   am3 <- mean(c(fct[[22]](xi),fct[[24]](xi)))
   Ai <- matrix(c(fct[[8]](xi),am1,am2,am1,fct[[17]](xi),am3,am2,am3,fct[[25]](xi)),3,3)
   if(type==".MBRE"){
      ai <- Ai %*% zi
      Am <- (Ai+Aa)/2; Ai <- Aa <- Am
      am <- (ai+aa)/2; ai <- aa <- am
      zi <- distr::solve(Ai,ai)
   }
   a <-  c(.dbeta%*%aa)
   aw <- c(.dbeta1%*%zi)
   A <-  .dbeta%*%Aa%*%.dbeta
   Aw <- Ai%*%.dbeta

   normt <- NormType()
   biast <- symmetricBias()
   nb <- ContNeighborhood(radius=0.5)
   ICT <- paste("optimally robust IC for", switch(type,
                      ".OMSE"="maxMSE",".RMXE"="RMX", ".MBRE"="maxBias"))
   riskT <- if(type!=".MBRE") "asGRisk" else "asBias"

   w <- new("HampelWeight")
      stand(w) <- Aw
      cent(w) <- aw
      clip(w) <- b

   if(type!=".MBRE"){
        weight(w) <- getweight(w, neighbor = nb, biastype = biast,
                          normW = normt)
   }else weight(w) <- minbiasweight(w, neighbor = nb, biastype = biast,
                          normW = normt)

   Risk <- list(asBias = list(value = b, biastype = biast,
                                       normtype = normt,
                                       neighbortype = class(nb)))

   res <- list(a = a, A = A, b = b, d = 0*a,
               normtype = normt, biastype = biast, w = w,
               info = c("optIC", ICT), risk = Risk,
               modifyIC = NULL)


   IC <- generateIC(nb, L2Fam, res)
   return(IC)
}

if(FALSE){
   res <- list(a = a, A = A, b = b, d = 0*a, w = w)

   IC <- ContIC(name = "interpolated IC of contamination type",
                CallL2Fam = L2Fam@fam.call,
                Curve = generateIC.fct(nb, L2Fam, res),
                clip = b,
                cent = a,
                stand = A,
                lowerCase =0*a,
                w = w,
                neighborRadius = nb@radius,
                modifyIC = NULL,
                normtype = normt,
                biastype = biast,
                Risks = list(),
                Infos = matrix(c("optIC", ICT), ncol = 2,
                            dimnames = list(character(0), c("method", "message"))
   ))
}
