#############################################################################
### GPD
#############################################################################
##---------------------------------------------------------------------
## on R-3.0.0
##---------------------------------------------------------------------
## install new versions of distr-family and robast-family of pkgs
### open R session
require(RobExtremes)
### -> change this according to where you checked out the svn repo:
.basepath <- "C:/rtest/RobASt/branches/robast-0.9/pkg"
## <-
.saveInterpGrid <- RobExtremes:::.saveInterpGrid
.OMSE.xi <- RobExtremes:::.OMSE.xi
.MBRE.xi <- RobExtremes:::.MBRE.xi
.RMXE.xi <- RobExtremes:::.RMXE.xi
.modify.xi.PFam.call <- RobExtremes:::.modify.xi.PFam.call
.myFolder <- file.path(.basepath,"RobAStRDA/R")
.myFolderTo <- file.path(.basepath,"RobExtremesBuffer")
#chkExist <- function(fN) if(!file.exists(fN)) dir.create(fN, recursive = TRUE)
#sapply(c(.myFolder1,.myFolder2,.myFolder3), chkExist)
PF <- GParetoFamily(); n0 <- "GPD"
###
toS <- gsub("XXXX",n0,"sysdataSnXXXX.rda")
.svInt0 <- function(what){
    RobExtremes:::.svInt(optF = what, PFam = PF, sysRdafolder = .myFolderTo,
                         nam = n0)}
## to make this parallel, start this on several processors
.saveInterpGrid(getShapeGrid(gridsize=500, cutoff.at.0=0.005),
                sysRdaFolder = .myFolderTo, sysdataWriteFile=toS,
                accuracy = 5000,upp=10, PFam = PF)
.svInt0(.OMSE.xi)
.svInt0(.MBRE.xi)
.svInt0(.RMXE.xi)

#############################################################################
### GEVD
#############################################################################
##---------------------------------------------------------------------
## on R-3.0.0
##---------------------------------------------------------------------
## install new versions of distr-family and robast-family of pkgs
### open R session
## install new versions of distr-family and robast-family of pkgs
### open R session
require(RobExtremes)
### -> change this according to where you checked out the svn repo:
.basepath <- "C:/rtest/RobASt/branches/robast-0.9/pkg"
## <-
.saveInterpGrid <- RobExtremes:::.saveInterpGrid
.OMSE.xi <- RobExtremes:::.OMSE.xi
.MBRE.xi <- RobExtremes:::.MBRE.xi
.RMXE.xi <- RobExtremes:::.RMXE.xi
.modify.xi.PFam.call <- RobExtremes:::.modify.xi.PFam.call
.myFolder <- file.path(.basepath,"RobAStRDA/R")
.myFolderTo <- file.path(.basepath,"RobExtremesBuffer")
#chkExist <- function(fN) if(!file.exists(fN)) dir.create(fN, recursive = TRUE)
#sapply(c(.myFolder1,.myFolder2,.myFolder3), chkExist)
PF <- GEVFamily(); n0 <- "GEV"
###
toS <- gsub("XXXX",n0,"sysdataSnXXXX.rda")
.svInt0 <- function(what){
    RobExtremes:::.svInt(optF = what, PFam = PF, sysRdafolder = .myFolderTo,
                         nam = n0)}
## to make this parallel, start this on several processors
.saveInterpGrid(getShapeGrid(gridsize=500, cutoff.at.0=0.005),
                sysRdaFolder = .myFolderTo, sysdataWriteFile=toS,
                accuracy = 5000,upp=10, PFam = PF)
.svInt0(.OMSE.xi)
.svInt0(.MBRE.xi)
.svInt0(.RMXE.xi)

#############################################################################
### Weibull
#############################################################################
##---------------------------------------------------------------------
## on R-3.0.0
##---------------------------------------------------------------------
## install new versions of distr-family and robast-family of pkgs
### open R session
require(RobExtremes)
### -> change this according to where you checked out the svn repo:
.basepath <- "C:/rtest/RobASt/branches/robast-0.9/pkg"
## <-
.saveInterpGrid <- RobExtremes:::.saveInterpGrid
.OMSE.xi <- RobExtremes:::.OMSE.xi
.MBRE.xi <- RobExtremes:::.MBRE.xi
.RMXE.xi <- RobExtremes:::.RMXE.xi
.modify.xi.PFam.call <- RobExtremes:::.modify.xi.PFam.call
.myFolder <- file.path(.basepath,"RobAStRDA/R")
.myFolderTo <- file.path(.basepath,"RobExtremesBuffer")
#chkExist <- function(fN) if(!file.exists(fN)) dir.create(fN, recursive = TRUE)
#sapply(c(.myFolder1,.myFolder2,.myFolder3), chkExist)
PF <- WeibullFamily(); n0 <- "Weib"
###
toS <- gsub("XXXX",n0,"sysdataSnXXXX.rda")
.svInt0 <- function(what){
    RobExtremes:::.svInt(optF = what, PFam = PF, sysRdafolder = .myFolderTo,
                         nam = n0)}
## to make this parallel, start this on several processors
.saveInterpGrid(getShapeGrid(gridsize=500, cutoff.at.0=0.005),
                sysRdaFolder = .myFolderTo, sysdataWriteFile=toS,
                accuracy = 5000,upp=10, PFam = PF)
.svInt0(.OMSE.xi)
.svInt0(.MBRE.xi)
.svInt0(.RMXE.xi)

#############################################################################
### Gamma
#############################################################################
##---------------------------------------------------------------------
## on R-3.0.0
##---------------------------------------------------------------------
## install new versions of distr-family and robast-family of pkgs
### open R session
require(RobExtremes)
### -> change this according to where you checked out the svn repo:
.basepath <- "C:/rtest/RobASt/branches/robast-0.9/pkg"
## <-
.saveInterpGrid <- RobExtremes:::.saveInterpGrid
.OMSE.xi <- RobExtremes:::.OMSE.xi
.MBRE.xi <- RobExtremes:::.MBRE.xi
.RMXE.xi <- RobExtremes:::.RMXE.xi
.modify.xi.PFam.call <- RobExtremes:::.modify.xi.PFam.call
.myFolder <- file.path(.basepath,"RobAStRDA/R")
.myFolderTo <- file.path(.basepath,"RobExtremesBuffer")
#chkExist <- function(fN) if(!file.exists(fN)) dir.create(fN, recursive = TRUE)
#sapply(c(.myFolder1,.myFolder2,.myFolder3), chkExist)
PF <- GammaFamily(); n0 <- "Gam"
###
toS <- gsub("XXXX",n0,"sysdataSnXXXX.rda")
.svInt0 <- function(what){
    RobExtremes:::.svInt(optF = what, PFam = PF, sysRdafolder = .myFolderTo,
                         nam = n0)}
## to make this parallel, start this on several processors
.saveInterpGrid(getShapeGrid(gridsize=500, cutoff.at.0=0.005),
                sysRdaFolder = .myFolderTo, sysdataWriteFile=toS,
                accuracy = 5000,upp=10, PFam = PF)
.svInt0(.OMSE.xi)
.svInt0(.MBRE.xi)
.svInt0(.RMXE.xi)
###
################################################################################
###
### merge and thin out results on R-3.0.0
###
rdafiles <- file.path(c(.myFolder,.myFolder1,.myFolder2,.myFolder3),"sysdata.rda")
.recomputeInterpolators(rdafiles, sysRdaFolder = myFolder)
### close R session;
##  R CMD build RobExtremes
##---------------------------------------------------------------------
## on R-2.15.2
##---------------------------------------------------------------------
## install new versions of distr-family and robast-family of pkgs
##  install RobExtremes from source on R-2.15.2
require(RobExtremes); RobExtremes:::.recomputeInterpolators(rdafiles[1], sysRdaFolder = .myFolder)
### close R session;
##  R CMD build RobExtremes
##  R CMD install RobExtremes from source
##---------------------------------------------------------------------
## on R-3.0.0
##---------------------------------------------------------------------
##  R CMD install RobExtremes from source
###

#### Fix Sn for GEV (which was wrong ...)
require(RobExtremes)
### -> change this according to where you checked out the svn repo:
.basepath <- "C:/rtest/RobASt/branches/robast-0.9/pkg"
## <-
.saveInterpGrid <- RobExtremes:::.saveInterpGrid
.svInt <- RobExtremes:::.svInt
.OMSE.xi <- RobExtremes:::.OMSE.xi
.MBRE.xi <- RobExtremes:::.MBRE.xi
.RMXE.xi <- RobExtremes:::.RMXE.xi
.modify.xi.PFam.call <- RobExtremes:::.modify.xi.PFam.call
.myFolder <- file.path(.basepath,"RobExtremes/R")
.myFolder1 <- file.path(.basepath,"RobExtremesBuffer/tmp1")
PF <- GEVFamily()
.saveInterpGrid(getShapeGrid(gridsize=500, cutoff.at.0=0.005),
                sysRdaFolder = .myFolder1, accuracy = 5000,upp=10,
                PFam = PF)


###############################################################################
### to merge files do
#
### on R-3.0.0
require(RobExtremes)
### -> change this according to where you checked out the svn repo:
.basepath <- "C:/rtest/RobASt/branches/robast-0.9/pkg"
## <-
.myFolderFrom <- file.path(.basepath,"RobExtremesBuffer")
rda0 <- dir(.myFolderFrom)
rdafiles <- rda0[grepl("\\.rda$",rda0)]
.recomputeInterpolators(rdafiles, sysdataWriteFile = "sysdata30.rda",
                        sysRdaFolder = myFolderFrom, overwrite=TRUE,
                        translate=FALSE)
#
### on R-2.15.0
require(RobExtremes)
### -> change this according to where you checked out the svn repo:
.basepath <- "C:/rtest/RobASt/branches/robast-0.9/pkg"
## <-
.myFolderFrom <- file.path(.basepath,"RobExtremesBuffer")
rdafile <- file.path(.basepath,"RobExtremesBuffer","sysdata30.rda")
.recomputeInterpolators(rdafile, sysdataWriteFile = "sysdataAll.rda",
                        sysRdaFolder = myFolderFrom)
### after check
### -> change this according to where you checked out the svn repo:
.basepath <- "C:/rtest/RobASt/branches/robast-0.9/pkg"
## <-
.myFolderFrom <- file.path(.basepath,"RobExtremesBuffer")
.myFolderTo <- file.path(.basepath,"RobAStRDA/R")
file.copy(from=file.path(.myFolderFrom, "sysdataAll.rda"),
          to=file.path(.myFolderTo, "sysdata.rda"),