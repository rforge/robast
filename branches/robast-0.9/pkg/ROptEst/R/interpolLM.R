.RMXE.th <- function(th, PFam, modifyfct){
      PFam <- modifyfct(th,PFam)
      IC <- radiusMinimaxIC(L2Fam=PFam, neighbor= ContNeighborhood(),
                            risk = asMSE(), verbose = FALSE)
      return(c(b=clip(IC), a=cent(IC), a.w = cent(weight(IC)),
                           A=stand(IC),  A.w = stand(weight(IC))))
}

.MBRE.th <- function(th, PFam, modifyfct){
      PFam <- modifyfct(th,PFam)
      RobM <- InfRobModel(center = PFam, neighbor = ContNeighborhood(radius = 15))
      IC <- optIC(model = RobM, risk = asBias(), verbose = FALSE)
      mA <- max(stand(IC))
      mAw <- max(stand(weight(IC)))
      return(c(b=clip(IC), a=cent(IC), aw=cent(weight(IC)),
               A=stand(IC)/mA, Aw=stand(weight(IC))/mAw))
}

.OMSE.th <- function(th, PFam, modifyfct){
      PFam <- modifyfct(th,PFam)
      RobM <- InfRobModel(center = PFam, neighbor = ContNeighborhood(radius = .5))
      IC <- optIC(model = RobM, risk = asMSE(), verbose = FALSE)
      res=c(b=clip(IC), a=cent(IC), a.w = cent(weight(IC)),
                A=stand(IC), A.w = stand(weight(IC)))
      return(res)
}

.getLMGrid <- function(thGrid, PFam, optFct = .RMXE.th, modifyfct,
                       GridFileName="LMGrid.Rdata",
                       withSmooth = TRUE, withPrint = FALSE, withCall = FALSE){
   print(match.call())
   call <- match.call()
   thGrid <- unique(sort(thGrid))
   itLM <- 0
   getLM <- function(th){
               itLM <<- itLM + 1
               if(withPrint) cat("Evaluation Nr.", itLM," at th = ",th,"\n")
               a <- try(optFct(th=th,PFam=PFam,modifyfct=modifyfct), silent=TRUE)
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
   LMGrid <- sapply(thGrid,getLM)
   if(GridFileName!="") save(LMGrid, file=GridFileName)
   res <- .MakeGridList(thGrid, Y=t(LMGrid), withSmooth = withSmooth)
   print(res)
   rm(itLM,getLM)
   if(withCall) rm(call)
   return(list(grid = res$grid,
               fct = res$fct, call = if(withCall) call else NULL))
}

.MakeGridList <- function(thGrid, Y, withSmooth = TRUE){
  if(length(dim(Y))==3)
     LMGrid <- Y[,1,,drop=TRUE]
  else LMGrid <- Y[,drop=FALSE]

   iNA <- apply(LMGrid,1, function(u) any(is.na(u)))
   LMGrid <- LMGrid[!iNA,,drop=FALSE]
   thGrid <- thGrid[!iNA]
   oG <- order(thGrid)
   thGrid <- thGrid[oG]
   LMGrid <- LMGrid[oG,,drop=FALSE]
   if(withSmooth)
      LMGrid2 <- apply(LMGrid,2,function(u) smooth.spline(thGrid,u)$y)

   fctL <- vector("list",ncol(LMGrid))
   xm <- thGrid[1]
   xM <- (rev(thGrid))[1]
   for(i in 1:ncol(LMGrid)){
       LMG <- LMGrid[,i]
       fct <- splinefun(x=thGrid,y=LMG)
       ym <- LMG[1]
       dym <- (LMG[2]-LMG[1])/(thGrid[2]-thGrid[1])
       yM <- (rev(LMG))[1]
       dyM <- ((rev(LMG))[2]-(rev(LMG))[1])/((rev(thGrid))[2]-(rev(thGrid))[1])
       fctX <- function(x){
            y0 <- fct(x)
            y1 <- y0
            y1[x<xm] <- ym+dym*(x[x<xm]-xm)
            y1[x>xM] <- yM+dyM*(x[x>xM]-xM)
            if(any(is.na(y0)))
               warning("There have been xi-values out of range of the interpolation grid.")
            return(y1)
       }
       environment(fctX) <- nE <- new.env()
       assign("fct",fct, envir=nE)
       assign("yM",yM, envir=nE)
       assign("ym",ym, envir=nE)
       assign("dyM",dyM, envir=nE)
       assign("dym",dym, envir=nE)
       fctL[[i]] <- fctX
   }
   if(ncol(LMGrid)==1) fctL <- fctL[[1]]
   rm(LMG,fct,fctX,iNA,ym,yM,dym,dyM)
   return(list(grid = cbind(xi=thGrid,LM=LMGrid),
               fct = fctL))
}


.saveInterpGrid <- function(thGrid, PFam, sysRdaFolder,
            sysdataWriteFile = "sysdata.rda", getFun = .getLMGrid, ...,
            modifyfct, nameInSysdata, GridFileName, withSmooth = TRUE,
            withPrint = TRUE, withCall = FALSE, Y = NULL, elseFun = NULL){
  if(missing(sysRdaFolder)) stop("You must specify argument 'sysRdaFolder'.")

  if(missing(GridFileName))
     GridFileName <- paste(sub("^\\.(.+)","\\1",nameInSysdata),".Rdata",sep="")
  newEnv <- new.env()
  sysdataFile <- file.path(sysRdaFolder, sysdataWriteFile)
  cat("sysdataFile = ", sysdataFile, "\n")

  if(file.exists(sysdataFile)){
     load(file=sysdataFile,envir=newEnv)
     whatIsThereAlready <- ls(envir=newEnv, all.names=TRUE)
  }else whatIsThereAlready <- character(0)

  cat("whatIsThereAlready = ", head(whatIsThereAlready), "\n")

  if(exists(.versionSuff(nameInSysdata),envir=newEnv,inherits=FALSE)){
    InterpGrids <- get(.versionSuff(nameInSysdata), envir=newEnv)
    namesInterpGrids <- names(InterpGrids)
    cat(gettext("Names of existing grids:\n"))
    cat(paste("   ", namesInterpGrids , "\n"))

    if(name(PFam)%in% namesInterpGrids){
       cat(gettext("There already is a grid for family "), name(PFam),".\n",sep="")
       if(!is.null(InterpGrids[[name(PFam)]]$call)){
           cat(gettextf("It was generated by\n"),sep="")
           print(InterpGrids[[name(PFam)]]$call)
       }
       cat("\n",
           gettext("Do you really want to overwrite it (yes=1/no else)?"),"\n",
           sep="")
       ans <- try(scan(what=integer(1)), silent = TRUE)
       if(is(ans,"try-error")) ans <- 0
       if(ans==1){
          if(is.null(Y)) {
              InterpGrids[[name(PFam)]] <- getFun(thGrid = thGrid, PFam = PFam,
                         ..., modifyfct = modifyfct, withSmooth = withSmooth,
                         withPrint = withPrint, withCall = withCall,
                         GridFileName = GridFileName)
          }else{ if(!is.null(elseFun)){
                   InterpGrids[[name(PFam)]] <- elseFun(thGrid, Y,
                                                      withSmooth = withSmooth)
                 }else return(NULL)
          }

          l.ng <- -1
          cat(gettext("SnGrid successfully produced.\n"))
       }else l.ng <- -2
    }else l.ng <- length(InterpGrids)+1
  }else{
    l.ng <- 1
    InterpGrids <- vector("list",1)
    whatIsThereAlready <- c(whatIsThereAlready,.versionSuff(nameInSysdata))
  }

  if(l.ng>0){
     if(is.null(Y)) {
           InterpGrids[[l.ng]] <- getFun(thGrid = thGrid, PFam = PFam,
                         ..., modifyfct = modifyfct, withSmooth = withSmooth,
                         withPrint = withPrint, withCall = withCall,
                         GridFileName = GridFileName)
     }else{ if(!is.null(elseFun)){
               InterpGrids[[l.ng]] <- elseFun(thGrid, Y, withSmooth = withSmooth)
            }else return(NULL)
     }
     cat(gettext("Grid successfully produced.\n"))
     names(InterpGrids)[l.ng] <- name(PFam)
  }

  if(l.ng> -2){
     assign(.versionSuff(nameInSysdata), InterpGrids, envir=newEnv)
     save(list=whatIsThereAlready, file=sysdataFile, envir=newEnv)
     tools::resaveRdaFiles(sysdataFile)
     cat(gettextf("%s successfully written to sysdata.rda file.\n",
            nameInSysdata))
  }
  rm(list=whatIsThereAlready,envir=newEnv)
  gc()
  return(invisible(NULL))
}

