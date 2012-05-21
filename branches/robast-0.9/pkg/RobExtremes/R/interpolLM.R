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
               if(is(a,"try-error")) a <- rep(NA,14)
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
   res <- .MakeGridList(xiGrid, Y=t(LMGrid), withSmooth = withSmooth)
   print(res)
   return(list(grid = res$grid,
               fct = res$fct, call = call))
}

.MakeGridList <- function(xiGrid, Y, withSmooth = TRUE){
  if(length(dim(Y))==3)
     LMGrid <- Y[,1,,drop=TRUE]
  else LMGrid <- Y

   iNA <- apply(LMGrid,1, function(u) any(is.na(u)))
   LMGrid <- LMGrid[!iNA,]
   xiGrid <- xiGrid[!iNA]
   if(withSmooth)
      LMGrid2 <- apply(LMGrid,2,function(u) smooth.spline(xiGrid,u)$y)

   print(LMGrid2)
   fct0 <- function(x,i) (splinefun(x=xiGrid,y=LMGrid[,i]))(x)

   xm <- xiGrid[1]
   ym <- LMGrid[1,]
   dym <- (LMGrid[2,]-LMGrid[1,])/(xiGrid[2]-xiGrid[1])
   xM <- (rev(xiGrid))[1]
   yM <- ym
   dyM <- dym
   for(i in 1:ncol(LMGrid)){
       yM[i] <- (rev(LMGrid[,i]))[1]
       dyM[i] <- ((rev(LMGrid[,i]))[2]-(rev(LMGrid[,i]))[1])/
                 ((rev(xiGrid))[2]-(rev(xiGrid))[1])
   }
   fct <- function(x,i){
       y0 <- fct0(x,i)
       y1 <- y0
       y1[x<xm] <- ym[i]+dym[i]*(x[x<xm]-xm)
       y1[x>xM] <- yM[i]+dyM[i]*(x[x>xM]-xM)
       if(any(is.na(y0)))
          warning("There have been xi-values out of range of the interpolation grid.")
       return(y1)
   }

   return(list(grid = cbind(xi=xiGrid,LM=LMGrid),
               fct = fct))
}

.svInt <- function(optF = .RMXE.xi, nam = ".RMXE")
                  .saveInterpGrid(xiGrid = getShapeGrid(400,
                  cutoff.at.0=0.005),
                  PFam = GParetoFamily(shape=1,scale=2),
                  sysRdaFolder = .myFolder, optFct = optF,
                  nameInSysdata = nam, getFun = .getLMGrid,
                  withSmooth = TRUE, withPrint = TRUE)

if(FALSE){
.myFolder <- "C:/rtest/RobASt/branches/robast-0.9/pkg/ROptEst/R"
svInt <- RobExtremes:::.svInt;
.OMSE.xi <- RobExtremes:::.OMSE.xi
.MBRE.xi <- RobExtremes:::.MBRE.xi
.RMXE.xi <- RobExtremes:::.RMXE.xi
.svInt(.OMSE.xi, ".OMSE")
.svInt(.MBRE.xi, ".MBRE")
.svInt(.RMXE.xi, ".RMXE")
}
