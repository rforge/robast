.getpsi <- function(xi, fct, L2Fam , type){

   L2deriv <- L2deriv(L2Fam)
   b <- fct(xi,1)
   a <- c(fct(xi,2),fct(xi,3))
   aw <- c(fct(xi,4),fct(xi,5))
   am <- mean(c(fct(xi,7),fct(xi,8)))
   A <-  matrix(c(fct(xi,6),am,am,fct(xi,9)),2,2)}
   am <- mean(c(fct(xi,11),fct(xi,12)))
   Aw <- matrix(c(fct(xi,10),am,am,fct(xi,13)),2,2)

   normt <- NormType()
   biast <- symmetricBias()
   ICT <- paste("optimally robust IC for", switch(type,
                      c(".OMSE"="maxMSE",".RMXE"="RMX", ".MBRE"="maxBias")))
   riskT <- if(nameInSysdata!=".MBRE") "asGRisk" else "asBias"

   w <- new("HampelWeight")
      stand(w) <- Aw(xi0)
   cent(w) <- aw(xi0)
   clip(w) <- b(xi0)
   if(type!=".MBRE")
        weight(w) <- getweight(w, neighbor = neighbor, biastype = biast,
                          normW = normt)
   else weight(w) <- minbiasweight(w, neighbor = neighbor, biastype = biast,
                          normW = normt)

   res <- list(a = a(xi0), A = A(xi0), b = b(xi0), d = 0,
               normtype = normt, biastype = biast, w = w,
               info = c("optIC", ICT), risk = riskT,
               modifyIC = function(L2Fam, IC){
                   para <- param(L2Fam)
                   xi0 <- main(para)[scaleshapename(L2Fam)["scale"]]
                   L2deriv0 <-  EuclRandVarList(RealRandVariable(
                               L2Fam@L2deriv.fct(para),
                               Domain = Reals()))
                   .getpsi(xi0,fct, L2deriv0, type)
               }
            )
   return(generateIC(ContNeighborhood(r=0.5), L2Fam, res))
}
