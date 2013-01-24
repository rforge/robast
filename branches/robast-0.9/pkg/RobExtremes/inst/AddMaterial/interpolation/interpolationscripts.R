require(RobExtremes)

.myFolder <- "C:/rtest/RobASt/branches/robast-0.9/pkg/RobExtremesBuffer"
.myFolder1 <- "C:/rtest/RobASt/branches/robast-0.9/pkg/RobExtremes/R"
.myFolder2 <- "C:/rtest/RobASt/branches/robast-0.9/pkg/RobExtremes"

### produce Sn grid
.saveInterpGrid <- RobExtremes:::.saveInterpGrid
.saveInterpGrid(getShapeGrid(gridsize=500, cutoff.at.0=0.005),
                sysRdaFolder = .myFolder, accuracy = 5000,upp=10)

### produce LM grids
.svInt <- RobExtremes:::.svInt
.OMSE.xi <- RobExtremes:::.OMSE.xi
.MBRE.xi <- RobExtremes:::.MBRE.xi
.RMXE.xi <- RobExtremes:::.RMXE.xi
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
