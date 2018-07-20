################################################################################
###
### merge and thin out results on rdafile
###
require(RobAStRDA)
.saveGridToRda <- RobAStRDA:::.saveGridToRda
.readGridFromCSV <- RobAStRDA:::.readGridFromCSV
.mergeGrid <- RobAStRDA:::.mergeGrid
.MakeSmoothGridList <- RobAStRDA:::.MakeSmoothGridList
.computeInterpolators <- RobAStRDA:::.computeInterpolators
.mergeF <- RobAStRDA:::.computeInterpolators
.generateInterpolators <- RobAStRDA:::.generateInterpolators
#####

oldwd <- getwd()
.basepath <- "C:/rtest/RobASt/branches/robast-1.1/pkg"
.myFolderFrom <- file.path(.basepath,"RobExtremesBuffer")
myRDA0 <- file.path(.basepath,"RobExtremesBuffer/sysdata.rda")
#myRDA <- file.path(.basepath,"RobExtremesBuffer/sysdata.rda")
#myRDA0 <- file.path(.basepath,"RobAStRDA/R/sysdata0.rda")
myRDA <- file.path(.basepath,"RobAStRDA/R/sysdata.rda")
CSVFiles <- grep("\\.csv$", dir(.myFolderFrom), value=TRUE)
CSVFiles <- paste(.myFolderFrom, CSVFiles, sep="/")

.saveGridToRda(CSVFiles, toFileRDA = myRDA0, withMerge = FALSE,
               withPrint = TRUE, withSmooth = TRUE, df = NULL)
##
.computeInterpolators(myRDA0, myRDA,withSmoothFct = TRUE)
###

if(FALSE){
#---------------------------------------------------------
# (1) load package in R>3.0
#---------------------------------------------------------
  if(getRversion()>"3.0"){
    require(RobAStRDA)
    .basepath <- "C:/rtest/RobASt/branches/robast-1.1/pkg"
    .myFolderFrom <- file.path(.basepath,"RobExtremesBuffer")
    (myRDAg30 <- file.path(.basepath,"RobExtremesBuffer/sysdataOnlyGridsOnlyR-3.5.1rc.rda"))
    (myRDAg3 <- file.path(.basepath,"RobExtremesBuffer/sysdataWithInterpOnlyR-3.5.1rc.rda"))
    file.remove(myRDAg3)
    file.remove(myRDAg30)
    CSVFiles <- grep("\\.csv$", dir(.myFolderFrom), value=TRUE)
    (CSVFiles <- paste(.myFolderFrom, CSVFiles, sep="/"))
    RobAStRDA:::.saveGridToRda(CSVFiles, toFileRDA = myRDAg30, withMerge = FALSE,
                               withPrint = TRUE, withSmooth = TRUE, df = NULL)
    RobAStRDA:::.computeInterpolators(myRDAg30, myRDAg3,withSmoothFct = TRUE)
    ###
    nEg3 <- new.env()
    load(myRDAg3,env=nEg3)
    nEg3L <- ls(all=TRUE,env=nEg3)
    for(nam in nEg3L){
       loc <- get(nam,env=nEg3)
       namU <- names(loc)
       for(nams in namU){
           print(c(nam,nams,names(loc[[nams]])))
       }
    }
  }
#---------------------------------------------------------
# (2) load package in R<2.15
#---------------------------------------------------------
  if(getRversion()<"2.16"){
    require(RobAStRDA)
    .basepath <- "C:/rtest/RobASt/branches/robast-1.1/pkg"
    .myFolderFrom <- file.path(.basepath,"RobExtremesBuffer")
    myRDAs30 <- file.path(.basepath,"RobExtremesBuffer/sysdataOnlyGridsOnlyR-2.15.1.rda")
    myRDAs3 <- file.path(.basepath,"RobExtremesBuffer/sysdataWithInterpOnlyR-2.15.1.rda")
    file.remove(myRDAs3)
    file.remove(myRDAs30)
    CSVFiles <- grep("\\.csv$", dir(.myFolderFrom), value=TRUE)
    CSVFiles <- paste(.myFolderFrom, CSVFiles, sep="/")
    RobAStRDA:::.saveGridToRda(CSVFiles, toFileRDA = myRDAs30, withMerge = FALSE,
                               withPrint = TRUE, withSmooth = TRUE, df = NULL)
    RobAStRDA:::.computeInterpolators(myRDAs30, myRDAs3,withSmoothFct = TRUE)
    ###
    nEs3 <- new.env()
    load(myRDAs3,env=nEs3)
    nEs3L <- ls(all=TRUE,env=nEs3)
    for(nam in nEs3L){
       loc <- get(nam,env=nEs3)
       namU <- names(loc)
       for(nams in namU){
           print(c(nam,nams,names(loc[[nams]])))
       }
    }
  }
#---------------------------------------------------------
# (3) back in R>3.0 merge grids
#---------------------------------------------------------
  if(getRversion()>"3.0"){
    myRDAs30 <- file.path(.basepath,"RobExtremesBuffer/sysdataOnlyGridsOnlyR-2.15.1.rda")
    myRDAs3 <- file.path(.basepath,"RobExtremesBuffer/sysdataWithInterpOnlyR-2.15.1.rda")
    nEs3 <- new.env()
    mergeE <- new.env()
    load(myRDAs3,env=nEs3)
    (nEs3L <- ls(all=TRUE,env=nEs3))
    names(get(".Gamma", env = nEs3)$MBRE)

    for(nam in nEg3L){
       loc <- get(nam,env=nEg3)
       namU <- names(loc)
       for(nams in namU){
           print(c(nam,nams,names(loc[[nams]])))
       loc[[nams]]$fun.O <- get(nam,env=nEs3)[[nams]][["fun.O"]]
       }
       assign(nam,loc,env=mergeE)
    }
    (mergeEL <- ls(all=TRUE,env=mergeE))

    for(nam in mergeEL){
       loc <- get(nam,env=mergeE)
       namU <- names(loc)
       for(nams in namU){
           print(c(nam,nams,names(loc[[nams]])))
       }
    }
  }
#---------------------------------------------------------
# (4) save merged files and zip them
#---------------------------------------------------------
  if(getRversion()>"3.0"){

    myRDAmerge <- file.path(.basepath,"RobExtremesBuffer/sysdataWithInterpMerge.rda")
    myRDAmergeZip <- file.path(.basepath,"RobExtremesBuffer/sysdataWithInterpMergeZip.rda")
    myRDA <- file.path(.basepath,"RobAStRDA/R/sysdata.rda")


    mergeEL <- ls(all=TRUE,env=mergeE)
    save(list=mergeEL,envir=mergeE,file=myRDAmerge)
    file.copy(from = myRDAmerge, to = myRDAmergeZip, overwrite = TRUE)
    tools::resaveRdaFiles(myRDAmergeZip)
    file.copy(from = myRDAmergeZip, to = myRDA, overwrite = TRUE)
    }

#---------------------------------------------------------
# end of if(FALSE)
#---------------------------------------------------------
}