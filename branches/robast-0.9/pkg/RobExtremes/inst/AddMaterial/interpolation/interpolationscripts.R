require(RobExtremes)

.myFolder <- "C:/rtest/RobASt/branches/robast-0.9/pkg/RobExtremesBuffer"
.myFolder1 <- "C:/rtest/RobASt/branches/robast-0.9/pkg/RobExtremes/R"
.myFolder2 <- "C:/rtest/RobASt/branches/robast-0.9/pkg/RobExtremes"

### produce Sn grid
.saveInterpGrid <- RobExtremes:::.saveInterpGrid
.saveInterpGrid(getShapeGrid(gridsize=500, cutoff.at.0=0.005),
                sysRdaFolder = .myFolder, accuracy = 5000,upp=10)

### produce LM grids
.saveInterpGrid <- RobExtremes:::.saveInterpGrid
.svInt <- RobExtremes:::.svInt
.OMSE.xi <- RobExtremes:::.OMSE.xi
.MBRE.xi <- RobExtremes:::.MBRE.xi
.RMXE.xi <- RobExtremes:::.RMXE.xi
.modify.xi.PFam.call <- RobExtremes:::.modify.xi.PFam.call
.svInt(.OMSE.xi, ".OMSE", sysRdaFolder = .myFolder)
.svInt(.MBRE.xi, ".MBRE", sysRdaFolder = .myFolder)
.svInt(.RMXE.xi, ".RMXE", sysRdaFolder = .myFolder)

#### now this is still much too large.
#### therefore, in a first step, we only use the saved grids
### and by .MakeGridList, separately for R<2.16 and R>2.16, we thin out
### the references

 .myfiles <- file.path(.myFolder, "sysdata.rda")

 .myfiles1 <- file.path(.myFolder1, "RobExtremes/R/sysdata.rda")

 ## on R-3.0.0
RobExtremes:::.recomputeInterpolators(.myfiles, sysRdaFolder = .myFolder)
 ## on R-2.15.2
RobExtremes:::.recomputeInterpolators(.myfiles1, sysRdaFolder = .myfolder2)

## some check (R-3.0.0) : fct[[1]] and fct[[2]] should be different...
require(RobExtremes); RobExtremes:::.recomputeInterpolators("sysdata.rda", sysRdaFolder = ".")
fct <- getFromNamespace(".OMSE.N", "RobExtremes")[[1]]$fct
fct[[1]](2);fct[[2]](2)


#############################################################################
### GEVD
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
.svInt <- RobExtremes:::.svInt
.OMSE.xi <- RobExtremes:::.OMSE.xi
.MBRE.xi <- RobExtremes:::.MBRE.xi
.RMXE.xi <- RobExtremes:::.RMXE.xi
.modify.xi.PFam.call <- RobExtremes:::.modify.xi.PFam.call
.myFolder <- file.path(.basepath,"RobExtremes/R")
.myFolder1 <- file.path(.basepath,"RobExtremesBuffer/tmp1")
.myFolder2 <- file.path(.basepath,"RobExtremesBuffer/tmp2")
.myFolder3 <- file.path(.basepath,"RobExtremesBuffer/tmp3")
chkExist <- function(fN) if(!file.exists(fN)) dir.create(fN, recursive = TRUE)
sapply(c(.myFolder1,.myFolder2,.myFolder3), chkExist)
PF <- GEVFamily()
.saveInterpGrid(getShapeGrid(gridsize=500, cutoff.at.0=0.005),
                sysRdaFolder = .myFolder, accuracy = 5000,upp=10,
                PFam = PF)
## to make this parallel, we write the results to different folders:
.svInt(.OMSE.xi, ".OMSE", PFam = PF, sysRdafolder = .myFolder1)
.svInt(.MBRE.xi, ".MBRE", PFam = PF, sysRdafolder = .myFolder2)
.svInt(.RMXE.xi, ".RMXE", PFam = PF, sysRdafolder = .myFolder3)

### merge and thin out results on R-3.0.0
rdafiles <- file.path(c(myFolder,.myFolder1,myFolder2,.myFolder3),"sysdata.rda")
.recomputeInterpolators("sysdata.rda", sysRdaFolder = myFolder)
### close R session;
##  R CMD build RobExtremes
##---------------------------------------------------------------------
## on R-2.15.2
##---------------------------------------------------------------------
## install new versions of distr-family and robast-family of pkgs
##  install RobExtremes from source on R-2.15.2
require(RobExtremes); RobExtremes:::.recomputeInterpolators(rdafiles, sysRdaFolder = .myFolder)
### close R session;
##  R CMD build RobExtremes
##  R CMD install RobExtremes from source
##---------------------------------------------------------------------
## on R-3.0.0
##---------------------------------------------------------------------
##  R CMD install RobExtremes from source
###

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
.svInt <- RobExtremes:::.svInt
.OMSE.xi <- RobExtremes:::.OMSE.xi
.MBRE.xi <- RobExtremes:::.MBRE.xi
.RMXE.xi <- RobExtremes:::.RMXE.xi
.modify.xi.PFam.call <- RobExtremes:::.modify.xi.PFam.call
.myFolder <- file.path(.basepath,"RobExtremes/R")
.myFolder0 <- file.path(.basepath,"RobExtremesBuffer/tmp0")
.myFolder1 <- file.path(.basepath,"RobExtremesBuffer/tmp1")
.myFolder2 <- file.path(.basepath,"RobExtremesBuffer/tmp2")
.myFolder3 <- file.path(.basepath,"RobExtremesBuffer/tmp3")
chkExist <- function(fN) if(!file.exists(fN)) dir.create(fN, recursive = TRUE)
sapply(c(.myFolder0,.myFolder1,.myFolder2,.myFolder3), chkExist)
PF <- WeibullFamily()
.saveInterpGrid(getShapeGrid(gridsize=500, cutoff.at.0=0.005),
                sysRdaFolder = .myFolder0, accuracy = 5000,upp=10,
                PFam = PF)
## to make this parallel, we write the results to different folders:
.svInt(.OMSE.xi, ".OMSE", PFam = PF, sysRdafolder = .myFolder1)
.svInt(.MBRE.xi, ".MBRE", PFam = PF, sysRdafolder = .myFolder2)
.svInt(.RMXE.xi, ".RMXE", PFam = PF, sysRdafolder = .myFolder3)

### merge and thin out results on R-3.0.0
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

#############################################################################
### Gamma
#############################################################################
##---------------------------------------------------------------------
## on R-3.0.0
##---------------------------------------------------------------------
## install new versions of distr-family and robast-family of pkgs
### open R session
require(RobExtremes)
.saveInterpGrid <- RobExtremes:::.saveInterpGrid
.svInt <- RobExtremes:::.svInt
.OMSE.xi <- RobExtremes:::.OMSE.xi
.MBRE.xi <- RobExtremes:::.MBRE.xi
.RMXE.xi <- RobExtremes:::.RMXE.xi
.modify.xi.PFam.call <- RobExtremes:::.modify.xi.PFam.call
### -> change this according to where you checked out the svn repo:
.basepath <- "C:/rtest/RobASt/branches/robast-0.9/pkg"
## <-
.myFolder0 <- file.path(.basepath,"RobExtremesBuffer/tmp0")
.myFolder1 <- file.path(.basepath,"RobExtremesBuffer/tmp1")
.myFolder2 <- file.path(.basepath,"RobExtremesBuffer/tmp2")
.myFolder3 <- file.path(.basepath,"RobExtremesBuffer/tmp3")
chkExist <- function(fN) if(!file.exists(fN)) dir.create(fN, recursive = TRUE)
sapply(c(.myFolder0,.myFolder1,.myFolder2,.myFolder3), chkExist)
PF <- GammaFamily()
.saveInterpGrid(getShapeGrid(gridsize=500, cutoff.at.0=0.005),
                sysRdaFolder = .myFolder0, accuracy = 5000,upp=10,
                PFam = PF)
## to make this parallel, we write the results to different folders:
.svInt(.OMSE.xi, ".OMSE", PFam = PF, sysRdafolder = .myFolder1)
.svInt(.MBRE.xi, ".MBRE", PFam = PF, sysRdafolder = .myFolder2)
.svInt(.RMXE.xi, ".RMXE", PFam = PF, sysRdafolder = .myFolder3)

### merge and thin out results on R-3.0.0
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
