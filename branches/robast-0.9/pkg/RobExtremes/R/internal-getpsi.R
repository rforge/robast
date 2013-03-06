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

   .dbeta <- diag(c(beta,1))
   b <- fct[[1]](xi)
   a <-  c(.dbeta%*%c(fct[[2]](xi),fct[[3]](xi)))
   aw <- c(.dbeta%*%c(fct[[4]](xi),fct[[5]](xi)))
   am <- mean(c(fct[[7]](xi),fct[[8]](xi)))
   A <-  .dbeta%*%matrix(c(fct[[6]](xi),am,am,fct[[9]](xi)),2,2)%*%.dbeta
   am <- mean(c(fct[[11]](xi),fct[[12]](xi)))
   Aw <- .dbeta%*%matrix(c(fct[[10]](xi),am,am,fct[[13]](xi)),2,2)%*%.dbeta



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

   res <- list(a = a, A = A, b = b, d = 0*a,
               normtype = normt, biastype = biast, w = w,
               info = c("optIC", ICT), risk = list(),
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