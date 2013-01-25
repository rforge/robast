.basepath <- "C:/rtest/RobASt/branches/robast-0.9/pkg"
.myFolder <- file.path(.basepath,"RobExtremes/R")
.myFolder1 <- file.path(.basepath,"RobExtremesBuffer/tmp1")
.myFolder2 <- file.path(.basepath,"RobExtremesBuffer/tmp2")
.myFolder3 <- file.path(.basepath,"RobExtremesBuffer/tmp3")
.myFoldera <- file.path(.basepath,"RobExtremesBuffer/tmpa")
rdafiles <- file.path(c(.myFolder,.myFolder1,.myFolder2,.myFolder3),"sysdata.rda")
require(RobExtremes); RobExtremes:::.recomputeInterpolators(rdafiles[1], sysRdaFolder = .myFolder)


require(RobExtremes)
.basepath <- "C:/rtest/RobASt/branches/robast-0.9/pkg"
.myFolder <- file.path(.basepath,"RobExtremes/R")
.myFolder0 <- file.path(.basepath,"RobExtremesBuffer/tmp0")
.myFoldera <- file.path(.basepath,"RobExtremesBuffer/tmpa")
fo=file.path(.myFolder0,"sysdata.rda")
fn=file.path(.myFolder,"sysdata-1.rda")
fi=file.path(.myFolder,"sysdata.rda")
PF <- GEVFamily()

RobExtremes:::.renameGridName(".SnGrids.N","Generalized Pareto Family", name(PF),fo,fn)
RobExtremes:::.recomputeInterpolators(c(fi,fn), sysRdaFolder = .myFoldera, overwrite=TRUE, translate=FALSE)
   nE <- new.env()
   load(file.path(.myFoldera,"sysdata.rda"),envir=nE)
   ls(all.names=TRUE,envir=nE)
   str(get(".SnGrids.N",envir=nE))
   head(get(".SnGrids.N",envir=nE)[[1]]$grid)
   nE <- new.env()
   load(fn,envir=nE)
   ls(all.names=TRUE,envir=nE)
   str(get(".SnGrids.N",envir=nE))
   head(get(".SnGrids.N",envir=nE)[[1]]$grid)
   nE <- new.env()
   load(fo,envir=nE)
   ls(all.names=TRUE,envir=nE)
   str(get(".SnGrids.N",envir=nE))
   head(get(".SnGrids.N",envir=nE)[[1]]$grid)
   nE <- new.env()
   load(fi,envir=nE)
   ls(all.names=TRUE,envir=nE)
   str(get(".SnGrids.N",envir=nE))
   head(get(".SnGrids.N",envir=nE)[[1]]$grid)


assignInNamespace(".renameGridName",.renameGridName,ns="RobExtremes")

.basepath <- "C:/rtest/RobASt/branches/robast-0.9/pkg"
.myFolder <- file.path(.basepath,"RobExtremes/R")
.myFolder0 <- file.path(.basepath,"RobExtremesBuffer/tmp0")
.myFolder1 <- file.path(.basepath,"RobExtremesBuffer/tmp1")
.myFolder2 <- file.path(.basepath,"RobExtremesBuffer/tmp2")
.myFolder3 <- file.path(.basepath,"RobExtremesBuffer/tmp3")
.myFoldera <- file.path(.basepath,"RobExtremesBuffer/tmpa")
require(RobExtremes);
fo=file.path(.myFolder0,"sysdata.rda")
fn=file.path(.myFolder,"sysdata-1.rda")
fn2=file.path(.myFoldera,"sysdata-1.rda")
fi=file.path(.myFolder,"sysdata.rda")
PF <- GEVFamily()
nE <- new.env()
load("MBRE.Rdata", envir=nE)
ls(all=T, envir=nE)
grid <- get("LMGrid", envir=nE)
RobExtremes:::.copyGrid(grid, ".MBRE.N", name(PF), name(PF), fo, fn2)
nE <- new.env()
load(fn2, envir=nE)
str(get(".MBRE.N",envir=nE))

nE <- new.env()

f <- function(a=1,env=nE){
   assign("K",a,envir=env)
}

ls(envir=nE); f(); ls(envir=nE)
