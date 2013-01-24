.recomputeInterpolators <- function(sysdataFiles, sysRdaFolder = ".",
                                   others = FALSE, onlyothers = FALSE,
                                   overwrite = TRUE, integrateto = FALSE,
                                   onlyCurrent = FALSE, withPrint =TRUE,
                                   withSmooth = TRUE,
                                   debug = FALSE){

  wprint <- function(...){ if (withPrint) print(...)}

  sam <- new.env()
  for(File in sysdataFiles) load(File, envir = sam)

  keep <- if(getRversion()>="2.16") "N" else "O"
  todo <- if(getRversion()>="2.16") "O" else "N"

  whatIsThereAlready <-  ls(all.names=TRUE, envir=sam)
  whatIsThereAlready.N <- grep(paste("^\\..+\\.",keep,"$",sep=""),
                              whatIsThereAlready,value=T)

  whatIsThereAlready.O <- grep(paste("^\\..+\\.",todo,"$",sep=""),
                             whatIsThereAlready,value=T)
  whatIsThereAlready.E <- setdiff(setdiff(whatIsThereAlready,
                                 whatIsThereAlready.N),whatIsThereAlready.O)
  
  wprint(whatIsThereAlready.N)

  only.grid <- new.env()

  if(others){
     wprint("recomputed anew from neither .O nor .N architecture")

     for(what in whatIsThereAlready.E){
       wprint(what)
       what.to <- paste(what,".",keep,sep="")
       vec <- get(what, envir=sam)
       for(Fam in names(vec)){
             wprint(Fam)
             grid <- vec[[Fam]]$grid
             wprint(head(grid))
             a0 <- .MakeGridList(grid[,1], Y=grid[,-1,drop=FALSE],
                                 withSmooth = withSmooth)
             vec[[Fam]] <- a0
       }
       assign(what.to, vec, envir=only.grid)
     }
     lsA <- ls(all.names=T,envir=only.grid)
     wprint(lsA)
  }

  if(!onlyothers){

    wprint("copied/recomputed anew from from current architecture")

    for(what in whatIsThereAlready.N){
      wprint(what)
      vec <- get(what, envir=sam)
      if(overwrite){
         for(Fam in names(vec)){
            wprint(Fam)
            grid <- vec[[Fam]]$grid
            wprint(head(grid))
            a0 <- .MakeGridList(grid[,1], Y=grid[,-1,drop=FALSE],
                                withSmooth = withSmooth)
            vec[[Fam]] <- a0
         }
      }
      if(integrateto){
         vec.E <- get(what, envir = only.grid)
         for(Fam in names(vec)){
            wprint(Fam)
            grid.E <- vec.E[[Fam]]$grid
            grid <- vec[[Fam]]$grid
            grid.0 <- rbind(grid.E, grid)
            oI <- order(grid.0[,1])
            wI <- !duplicated(grid.0[oI,1])
            grid <- grid.0[wI,]
            wprint(head(grid))
            a0 <- .MakeGridList(grid[,1], Y=grid[,-1,drop=FALSE],
                                withSmooth = withSmooth)
            vec[[Fam]] <- a0
         }
      }
      assign(what, vec, envir=only.grid)
    }
    lsA <- ls(all.names=T,envir=only.grid)
    wprint(lsA)

    if(!onlyCurrent){
       wprint("copy foreign architecture")
       for(what in whatIsThereAlready.O){
           wprint(what)
           assign(what, get(what, envir=sam), envir=only.grid)
           }
    }
    lsA <- ls(all.names=T,envir=only.grid)
    wprint(lsA)

    for(what in whatIsThereAlready.O){
        wprint("translating foreign to current architecture")
        what.N <- sub(paste("\\.", todo, "$", sep=""),
                      paste(".", keep, sep=""),what)
        wprint(c(from=what, to=what.N))

        wG <- get(what, envir=sam)
        anyFam <- FALSE
        vec <- NULL
        if(onlyCurrent) if(what.N %in% whatIsThereAlready.N)
                           vec <- get(what.N,envir=sam)
        for(Fam in names(wG)){
            wprint(Fam)
            if(! Fam %in% names(vec)){
               anyFam <- TRUE
               grid <- wG[[Fam]]$grid
               wprint(head(grid))
               a0 <- .MakeGridList(grid[,1], Y=grid[,-1,drop=FALSE],
                                withSmooth = withSmooth)
               vec[[Fam]] <- a0
            }
        }
        if(integrateto){
            vec.E <- get(what, envir = only.grid)
            anyFam <- TRUE
            for(Fam in names(wG)){
               wprint(Fam)
               grid.E <- vec.E[[Fam]]$grid
               grid <- wG[[Fam]]$grid
               grid.0 <- rbind(grid.E, grid)
               oI <- order(grid.0[,1])
               wI <- !duplicated(grid.0[oI,1])
               grid <- grid.0[wI,]
               wprint(head(grid))
               a0 <- .MakeGridList(grid[,1], Y=grid[,-1,drop=FALSE],
                                withSmooth = withSmooth)
               vec[[Fam]] <- a0
            }
        }
        if(anyFam) assign(what.N, vec, envir=only.grid)
    }
    lsA <- ls(all.names=T,envir=only.grid)
    wprint(lsA)
  }

  sysFile <- file.path(sysRdaFolder,"sysdata.rda")

  if(!debug){
     save(list=lsA, envir=only.grid, file=sysFile)
     tools::resaveRdaFiles(sysRdaFolder)
  }else{
     print(paste("save(list=lsA, envir=only.grid, file=", sysFile,")", sep=""))
  }
}

if(FALSE){
  source("makegridlist.R")
 .myFolder <- "C:/rtest/RobASt/branches/robast-0.9/pkg"
  source(file.path(.myFolder,"RobExtremes/R","recomputeinterpolators.R"))
 .myfiles1 <- file.path(.myFolder,
               c("ROptEst/R", "RobExtremes/R", "RobExtremesBuffer"),
               "sysdata.rda")

 .myfiles <- file.path(.myFolder, "RobExtremes/R/sysdata.rda")
 
 .recomputeInterpolators(file.path(.myFolder,"RobExtremes/R/sysdata.rda"),
        sysRdaFolder = file.path(.myFolder,"RobExtremes/R"), debug = TRUE)

  wha <- c(".OMSE",".RMXE",".MBRE",".SnGrids")
  for(w in wha) assign(w,get(w, envir=asNamespace("RobExtremes")),envir=nE)
  save(list=wha,envir=nE, file="sysdata-oold.rda")

 .recomputeInterpolators(c("sysdata-ooold.rda","sysdata-0.rda"), others = TRUE, onlyothers = FALSE,    sysRdaFolder = ".", integrateto = TRUE)

}