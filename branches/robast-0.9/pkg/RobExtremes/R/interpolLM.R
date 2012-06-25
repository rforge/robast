.modify.xi.PFam.call <- function(xi, PFam){
      Param <- param(PFam)
      param <- main(Param)
      param["shape"] <- xi
      main(Param) <- param
      nModel <- modifyModel(PFam, Param)
#GParetoFamily(shape=xi,scale=1)
}

.RMXE.xi <- function(xi, PFam){
      PFam <- .modify.xi.PFam.call(xi,PFam)
      IC <- radiusMinimaxIC(L2Fam=PFam, neighbor= ContNeighborhood(),
                            risk = asMSE(), verbose = FALSE)
      return(c(b=clip(IC), a=cent(IC), a.w = cent(weight(IC)),
                           A=stand(IC),  A.w = stand(weight(IC))))
}

.MBRE.xi <- function(xi, PFam){
      PFam <- .modify.xi.PFam.call(xi,PFam)
      RobM <- InfRobModel(center = PFam, neighbor = ContNeighborhood(radius = 15))
      IC <- optIC(model = RobM, risk = asBias(), verbose = FALSE)
      mA <- max(stand(IC))
      mAw <- max(stand(weight(IC)))
      return(c(b=clip(IC), a=cent(IC), aw=cent(weight(IC)),
               A=stand(IC)/mA, Aw=stand(weight(IC))/mAw))
}

.OMSE.xi <- function(xi, PFam){
      PFam <- .modify.xi.PFam.call(xi,PFam)
      RobM <- InfRobModel(center = PFam, neighbor = ContNeighborhood(radius = .5))
      IC <- optIC(model = RobM, risk = asMSE(), verbose = FALSE)
      res=c(b=clip(IC), a=cent(IC), a.w = cent(weight(IC)),
                A=stand(IC), A.w = stand(weight(IC)))
      return(res)
}


.getLMGrid <- function(xiGrid = getShapeGrid(),
                      PFam = GParetoFamily(scale=1,shape=2),
                      optFct = .RMXE.xi,
                      withSmooth = TRUE,
                      withPrint = FALSE){
   print(match.call())
   call <- match.call()
   itLM <- 0
   getLM <- function(xi){
               itLM <<- itLM + 1
               if(withPrint) cat("Evaluation Nr.", itLM," at xi = ",xi,"\n")
               a <- try(optFct(xi,PFam), silent=TRUE)
               if(is(a,"try-error")) a <- rep(NA,13)
               return(a)
               }

   distroptions.old <- distroptions()
   distrExOptions.old <- distrExOptions()
   distroptions("withgaps"=FALSE)
   distrExOptions( MCIterations=1e6,
                   GLIntegrateTruncQuantile=.Machine$double.eps,
                   GLIntegrateOrder=1000,
                   ElowerTruncQuantile=1e-7,
                   EupperTruncQuantile=1e-7,
                   ErelativeTolerance = .Machine$double.eps^0.4,
                   m1dfRelativeTolerance = .Machine$double.eps^0.4,
                   m2dfRelativeTolerance = .Machine$double.eps^0.4,
                   nDiscretize = 300, IQR.fac = 20)
   on.exit({do.call(distrExOptions,args=distrExOptions.old)
            do.call(distroptions,args=distroptions.old)

            })
   LMGrid <- sapply(xiGrid,getLM)
   save(LMGrid, file="LMGrid.Rdata")
   res <- .MakeGridList(xiGrid, Y=t(LMGrid), withSmooth = withSmooth)
   print(res)
   return(list(grid = res$grid,
               fct = res$fct, call = call))
}

.MakeGridList <- function(xiGrid, Y, withSmooth = TRUE){
  if(length(dim(Y))==3)
     LMGrid <- Y[,1,,drop=TRUE]
  else LMGrid <- Y[,drop=FALSE]

   iNA <- apply(LMGrid,1, function(u) any(is.na(u)))
   LMGrid <- LMGrid[!iNA,,drop=FALSE]
   xiGrid <- xiGrid[!iNA]
   if(withSmooth)
      LMGrid2 <- apply(LMGrid,2,function(u) smooth.spline(xiGrid,u)$y)

   fctL <- vector("list",ncol(LMGrid))
   for(i in 1:ncol(LMGrid)){
       LMG <- LMGrid[,i]
       fct <- splinefun(x=xiGrid,y=LMG)
       xm <- xiGrid[1]
       ym <- LMG[1]
       dym <- (LMG[2]-LMG[1])/(xiGrid[2]-xiGrid[1])
       xM <- (rev(xiGrid))[1]
       yM <- ym
       dyM <- dym
       yM <- (rev(LMG))[1]
       dyM <- ((rev(LMG))[2]-(rev(LMG))[1])/((rev(xiGrid))[2]-(rev(xiGrid))[1])
       fctL[[i]] <- function(x){
            y0 <- fct(x)
            y1 <- y0
            y1[x<xm] <- ym+dym*(x[x<xm]-xm)
            y1[x>xM] <- yM+dyM*(x[x>xM]-xM)
            if(any(is.na(y0)))
               warning("There have been xi-values out of range of the interpolation grid.")
            return(y1)
       }
       environment(fctL) <- new.env()
       assign("fct", fct, envir=environment(fctL))
   }
   if(ncol(LMGrid)==1) fctL <- fctL[[1]]

   return(list(grid = cbind(xi=xiGrid,LM=LMGrid),
               fct = fctL))
}

.myFolder <- "C:/rtest/RobASt/branches/robast-0.9/pkg/ROptEst/R"
.svInt <- function(optF = .RMXE.xi, nam = ".RMXE")
             .saveInterpGrid(xiGrid = getShapeGrid(200,
                  cutoff.at.0=0.01),
                  PFam = GParetoFamily(shape=1,scale=2),
                  sysRdaFolder=.myFolder, optFct = optF,
                  nameInSysdata = nam, getFun = .getLMGrid,
                  withSmooth = TRUE, withPrint = TRUE)

