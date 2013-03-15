#############################################################################
# just creating csv files
#############################################################################
## install new versions of distr-family and robast-family of pkgs
### open R session
require(RobExtremes)
### -> change this according to where you checked out the svn repo:
.basepath <- "C:/rtest/RobASt/branches/robast-0.9/pkg"
## <-
oldwd <- getwd()
.myFolderTo <- file.path(.basepath,"RobExtremesBuffer")
setwd(.myFolderTo)
.OMSE.th <- ROptEst:::.OMSE.th
.MBRE.th <- ROptEst:::.MBRE.th
.RMXE.th <- ROptEst:::.RMXE.th
.modify.xi.PFam.call <- RobExtremes:::.modify.xi.PFam.call
#
#PF <- GParetoFamily()
#PF <- GEVFamily()
PF <- GammaFamily()
#PF <- WeibullFamily()
###
.svInt <- RobExtremes:::.svInt
#.svInt1 <- function(){
#    RobExtremes:::.generateInterpGridSn(PFam = PF)}
## to make this parallel, start this on several processors
#.svInt1()
#.svInt(.OMSE.th, PFam=PF)
#.svInt(.MBRE.th, PFam=PF)
.svInt(.RMXE.th, PFam=PF)
setwd(oldwd)
