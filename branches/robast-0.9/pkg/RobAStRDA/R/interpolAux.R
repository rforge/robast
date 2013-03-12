.versionSuff <- function(name){
    paste(sep="", name, if(getRversion()<"2.16") ".O" else ".N")
}

.MakeSmoothGridList <- function(thGrid, Y, df=NULL){
   if(length(dim(Y))==3)
      LMGrid <- Y[,1,,drop=TRUE]
   else LMGrid <- Y[,drop=FALSE]

   iNA <- apply(LMGrid,1, function(u) any(is.na(u)))
   LMGrid <- LMGrid[!iNA,,drop=FALSE]
   thGrid <- thGrid[!iNA]
   oG <- order(thGrid)
   thGrid <- thGrid[oG]
   LMGrid <- LMGrid[oG,,drop=FALSE]

   LMGrid <- apply(LMGrid,2,function(u) if(is.null(df))
                  smooth.spline(thGrid,u)$y else smooth.spline(thGrid,u,df=df)$y
                  )
   return(cbind(xi=thGrid,LM=LMGrid))
}

.readGridFromCSV <- function(fromFileCSV){
  Grid <- as.matrix(read.csv(fromFileCSV)); dimnames(Grid) <- NULL
  fromFileTXT <- gsub("(.+\\.)csv$","\\1txt",fromFileCSV)
  res2 <- scan(file=fromFileTXT, what=c("character","character"))
  return(list(Grid=Grid, namPFam=res2[1], namInSysdata=res2[2]))
}

############################################################################
# .generateInterpolators generates the interpolators to a given grid
#     and returns a list of the given grid and the function list
############################################################################
.generateInterpolators <- function(Grid, approxOrspline = "spline"){
  thGrid <- Grid[,1]
  LMGrid <- Grid[,-1,drop=FALSE]
  fctL <- vector("list",ncol(LMGrid))
  xm <- thGrid[1]
  xM <- (rev(thGrid))[1]
  for(i in 1:ncol(LMGrid)){
       LMG <- LMGrid[,i]
       fct <- if(approxOrspline=="spline")
                  splinefun(x=thGrid,y=LMG) else approxfun(x=thGrid,y=LMG)
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
               warning(paste("There have been xi-values out of range ",
                             "of the interpolation grid.", sep = ""))
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
  rm(LMG,fct,fctX,ym,yM,dym,dyM)
  return(list(grid = Grid, fct = fctL))
}


############################################################################
# .saveGridToRda loads in one or more grids from one ore more csv file(s)
#   (argument fromFileCSV) and writes the respective merged grid to an
#    rda-file generated from toFileRDA, sysRdaFolder
#    if withMerge == FALSE corresponding entries are not merged but overwritten
############################################################################
.saveGridToRda <- function(fromFileCSV, toFileRDA = "sysdata.rda",
                           withMerge =FALSE, withPrint = TRUE,
                           withSmooth = TRUE, df = NULL){

  ### check whether input is complete
  if(missing(fromFileCSV)) stop("You must specify argument 'fromFileCSV'.")
  if(missing(toFileRDA)) stop("You must specify argument 'toFileRDA'.")

  ## new environment to store all merged information from sysdata.rda-type file
  ## and new grids
  newEnv <- new.env()


  ### determine what objects already exist in sysdata.rda - type file
  if(file.exists(toFileRDA)){
     load(file=toFileRDA,envir=newEnv)
     whatIsThereAlready <- ls(envir=newEnv, all.names=TRUE)
  }else whatIsThereAlready <- character(0)

  ### load precomputed grids from file
  le <- length(fromFileCSV)
  CSVlist <- vector("list",le)
  if(le>0) for(i in 1:le){
      CSVlist[[i]] <- .readGridFromCSV(fromFileCSV[i])
      nameInSysdata <- CSVlist[[i]]$namInSysdata
      namPFam <- CSVlist[[i]]$namPFam
      Grid <- CSVlist[[i]]$Grid
      GridFileName <- paste(sub("^\\.(.+)","\\1",nameInSysdata),".Rdata",sep="")

      ### check whether object nameInSysdata already exists (ie. some
      ##   grids for this family already exist) or not
      if(!exists(nameInSysdata,envir=newEnv,inherits=FALSE)){
          l.ng <- 1
          InterpGrids <- vector("list",1)
          whatIsThereAlready <- c(whatIsThereAlready,nameInSysdata)
      }else{ ## already exists -> some merging necessary
          InterpGrids <- get(nameInSysdata, envir=newEnv)
          namesInterpGrids <- names(InterpGrids)
          cat(gettext("Names of existing grids:\n"))
          cat(paste("   ", namesInterpGrids , "\n"))
          if(namPFam %in% namesInterpGrids){
             cat(gettext("There already is a grid for family "),
                 namPFam,".\n",sep="")
             cat("\n",
                 gettext("Do you really want to overwrite/merge it (yes=1/no else)?"),"\n",
                 sep="")
             ans <- try(scan(what=integer(1)), silent = TRUE)
             if(is(ans,"try-error")) ans <- 0
             if(ans==1){
                if(withMerge){
                   gr0 <- .mergeGrid(InterpGrids[[namPFam]]$grid, Grid)
                   InterpGrids[[namPFam]]$grid <- gr0
                   if(withSmooth)
                      InterpGrids[[namPFam]]$gridS <-
                        .MakeSmoothGridList(gr0[,1],gr0[,-1,drop=FALSE],
                                            df = df)
                   cat(gettext("Grid successfully merged.\n"))
                }else{
                   InterpGrids[[namPFam]]$grid <- Grid
                   InterpGrids[[namPFam]]$gridS <-
                        .MakeSmoothGridList(Grid[,1],Grid[,-1,drop=FALSE],
                                            df = df)
                   cat(gettext("Grid successfully overwritten.\n"))
                }
                l.ng <- -1
             }else l.ng <- -2
          }else l.ng <- length(InterpGrids)+1
      }
      if(l.ng>0){ ## a new family is entered
         InterpGrids[[l.ng]]$grid <- Grid
         InterpGrids[[l.ng]]$gridS <-
              .MakeSmoothGridList(Grid[,1], Grid[,-1,drop = FALSE], df = df)
         cat(gettext("New Grid successfully produced.\n"))
         names(InterpGrids)[l.ng] <- namPFam
      }

      if(l.ng> -2){
         assign(nameInSysdata, InterpGrids, envir=newEnv)
      }
    }
  save(list=whatIsThereAlready, file=toFileRDA, envir=newEnv)
  tools::resaveRdaFiles(toFileRDA)
  rm(list=whatIsThereAlready,envir=newEnv)
  gc()
  return(invisible(NULL))
}

############################################################################
# .mergeGrid merges two grids according to the lines
############################################################################
.mergeGrid <- function(Grid1, Grid2){
   if(ncol(Grid1)==ncol(Grid2))
      stop("Grids must have the same number of columns")
   Grid <- rbind(Grid1, Grid2)
   on <- order(Grid[,1])
   Grid <- Grid[on,,drop=FALSE]
   dup <- duplicated(Grid[,1])
   Grid <- Grid[!dup,,drop=FALSE]
   return(Grid)
}

############################################################################
# .computeInterpolators takes one ore more given  sysdatafiles and produces
#    the respective interpolators writing them to file
############################################################################
.computeInterpolators <- function(sysdataFiles,  toFileRDA = "sysdata.rda",
                                   includeGrids = NULL, includeNams = NULL,
                                   excludeGrids = NULL, excludeNams = NULL,
                                   withPrint = TRUE, withSmoothFct = FALSE,
                                   approxOrspline = "spline"){

  wprint <- function(...){ if (withPrint) print(...)}

  samEnv <- new.env()
  toEnv  <- new.env()
  for(File in sysdataFiles) .mergeF(File, envir = samEnv,
            includeGrids = includeGrids , includeNams = includeNams,
            excludeGrids = excludeGrids , excludeNams = excludeNams)

  funN <- paste("fun.", if(getRversion()>="2.16") "N" else "O", sep = "")

  whatIsThereAlready <-  ls(all.names=TRUE, envir=samEnv)
  wprint(whatIsThereAlready)

  for(what in whatIsThereAlready){
      whatG <- get(what, envir=samEnv)
      for(Fam in names(whatG)){
          Grid <- whatG[[Fam]]$grid
          if(withSmoothFct && !is.null(whatG[[Fam]]$gridS))
             Grid <- whatG[[Fam]]$gridS
          whatG[[Fam]][[funN]] <- .generateInterpolators(Grid, approxOrspline)$fct
      }
      assign(what,whatG,toEnv)
  }

  save(list=whatIsThereAlready, envir=toEnv, file=toFileRDA)
  tools::resaveRdaFiles(toFileRDA)
  return(invisible(NULL))
}


.mergeF <- function(file,envir, includeGrids = NULL, includeNams = NULL,
                    excludeGrids = NULL, excludeNams = NULL){
  envir2 <- new.env()
  load(file,envir=envir2)
  what0 <- ls(all.names=TRUE,envir=envir2)
  rm(list=what0[! what0 %in% includeGrids], envir=envir2)
  rm(list=excludeGrids, envir=envir2)
  what1 <- ls(all.names=TRUE,envir=envir)
  what2 <- ls(all.names=TRUE,envir=envir2)
  for(w2 in what2){
      wG2 <- get(w2, envir=envir2)
      if(w2 %in% what1){
         wG1 <- get(w2, envir=envir)
         nwG1 <- names(wG1)
         if(!is.null(includeNams)) nwG1 <- nwG1[nwG1 %in% includeNams]
         if(length(nwG1))
         for(Fam1 in nwG1){
             if( Fam1 %in% excludeNams)   wG2[[Fam1]] <- NULL
             if( ! Fam1 %in% names(wG2))  wG2[[Fam1]] <- wG1[[Fam1]]
         }
      }
      assign(w2,wG2,envir=envir)
  }
  return(invisible(NULL))
}

.copy_smoothGrid <- function(gridEntry = NULL, rdafileOld, gridnamOld, FamnamOld,
                 rdafileNew, gridnamNew, FamnamNew, withSmooth = FALSE, df = NULL){

  if(missing(rdafileOld)) stop("Argument 'rdafileOld' must not be missing.")
  if(missing(gridnamOld)) stop("Argument 'gridnamOld' must not be missing.")
  if(missing(FamnamOld)) stop("Argument 'FamnamOld' must not be missing.")
  if(missing(rdafileNew)) rdafileNew <- rdafileOld
  if(missing(gridnamNew)) gridnamNew <- gridnamOld
  if(missing(FamnamNew)) FamnamNew <- FamnamOld
  nE <- new.env()
  load(rdafileOld,envir=nE)
  gr <- get(gridnamOld,envir=nE)
  gr[[FamnamNew]] <- gr[[FamnamOld]]
  if(is.null(gridEntry))
     gridEntry <- gr[[FamnamOld]]$grid  else  gr[[FamnamNew]]$grid <- gridEntry
  if(withSmooth){
        gr[[FamnamNew]]$gridS <- .MakeSmoothGridList(gridEntry[,1],
                                    gridEntry[,-1, drop = FALSE],  df = df)
  }else gr[[FamnamNew]]$gridS <- NULL
  assign(gridnamNew,gr,envir=nE)
  what <- ls(envir=nE, all.names = TRUE)
  save(list=what, file= rdafileNew, envir=nE)
  tools::resaveRdaFiles(rdafileNew)
  return(invisible(NULL))
}

.renameGridName <- function(rdafileOld, gridnamOld, FamnamOld,
                            rdafileNew, gridnamNew, FamnamNew){
  if(missing(rdafileOld)) stop("Argument 'rdafileOld' must not be missing.")
  if(missing(gridnamOld)) stop("Argument 'gridnamOld' must not be missing.")
  if(missing(FamnamOld)) stop("Argument 'FamnamOld' must not be missing.")
  if(missing(rdafileNew)) rdafileNew <- rdafileOld
  if(missing(gridnamNew)) gridnamNew <- gridnamOld
  if(missing(FamnamNew)) FamnamNew <- FamnamOld
   nE <- new.env()
   load(rdafileOld,envir=nE)
   what <- ls(all.names=TRUE,envir=nE)
   a <- get(gridnamOld, envir=nE)
   na <- names(a)
   wi <- which(FamnamOld==na)
   na[wi] <- FamnamNew
   names(a) <- na
   assign(gridnamNew,a,envir=nE)
   save(list=what, file=rdafileNew, envir=nE)
   tools::resaveRdaFiles(rdafileNew)
   return(invisible(NULL))
}

