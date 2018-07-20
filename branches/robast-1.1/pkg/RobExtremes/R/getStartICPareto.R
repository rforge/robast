setMethod("getStartIC",signature(model = "ParetoFamily", risk = "interpolRisk"),
           function(model, risk, ...){

    param1 <- param(model)
    xi <- main(param1)
    .modifyIC0 <- function(L2Fam, IC){
              xi0 <- main(param(L2Fam))
              return(.getPsi.P(xi0, type(risk)))
    }
    .modifyIC <- function(L2Fam,IC){
         psi.0 <- .modifyIC0(L2Fam,IC)
         psi.0@modifyIC <- .modifyIC
         return(psi.0)
    }
    IC0 <- .getPsi.P(xi, type(risk))
    IC0@modifyIC <- .modifyIC
    return(IC0)
    })

.getPsi.P <- function(xi, type){
   ## the respective LMs have been computed ahead of time
   ## and stored in sysdata.rda of this package
   ## the code for this computation is in AddMaterial/getLMPareto.R
   if(type==".MBRE"){
         b  <- xi*.ParetoLM$MBR["b"]
         a  <- xi*.ParetoLM$MBR["a"]
         aw <-    .ParetoLM$MBR["aw"]
         A  <- xi*.ParetoLM$MBR["A"]
         Aw <- xi*.ParetoLM$MBR["Aw"]
   }else{if(type==".RMXE"){
         b  <- xi*.ParetoLM$RMX["b"]
         a  <- xi*.ParetoLM$RMX["a"]
         aw <-    .ParetoLM$RMX["aw"]
         A  <- xi*.ParetoLM$RMX["A"]
         Aw <- xi*.ParetoLM$RMX["Aw"]
      }else{if(type==".OMSE"){
         b  <- xi*.ParetoLM$OMS["b"]
         a  <- xi*.ParetoLM$OMS["a"]
         aw <-    .ParetoLM$OMS["aw"]
         A  <- xi*.ParetoLM$OMS["A"]
         Aw <- xi*.ParetoLM$OMS["Aw"]
         }
      }
   }
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
