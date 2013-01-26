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


fu <- function(xi){
  ext <- if(getRversion()<"2.16") "\\.O$" else "\\.N$"
  lapply(grep(ext,w,val=T), function(x) {
    print(x)
    u <- get(x,envir=nE);
    for(i in 1:length(u)){
        ni <- names(u)[i]
        print(ni)
        print(!grepl("Sn",x))
        if(!grepl("Sn",x)){
           len <- length(u[[i]]$fct)
           yi <- sapply(1:len, function(j) u[[i]]$fct[[j]](xi))
           if(length(xi)==1) yi <- matrix(yi,ncol=len)
           colnames(yi) <- c("b","a1.a", "a2.a", "a1.i", "a2.i", "A11.a",
                   "A12.a", "A21.a", "A22.a", "A11.i", "A12.i", "A21.i", "A22.i")
           print(cbind(xi,yi))
        }else{
           Sn <- u[[i]]$fct(xi)
           print(cbind(xi,Sn))
        }
    }
    return(invisible(NULL))
    })
  return(invisible(NULL))
}

require(RobExtremes)
.basepath <- "C:/rtest/RobASt/branches/robast-0.9/pkg"
.myFolder <- file.path(.basepath,"RobExtremesBuffer")
.myFolderA <- file.path(.basepath,"RobExtremesBuffer/all2")
.myFolderW <- file.path(.basepath,"RobExtremesBuffer/WTS2")
fn2 <- file.path(.myFolder,"tmp2/sysdata.rda")
fn3 <- file.path(.myFolder,"tmp3/sysdata.rda")
fn00=file.path(.myFolderW,"tmp0/sysdata.rda")
fn01=file.path(.myFolderW,"tmp1/sysdata.rda")
fn02=file.path(.myFolderW,"tmp2/sysdata.rda")
fn03=file.path(.myFolderW,"tmp3/sysdata.rda")
fnA <- file.path(.myFolderA,"sysdata.rda")
#fn2=file.path(.myFoldera,"sysdata-1.rda")
#RobExtremes:::.recomputeInterpolators(c(fn01,fn02, fn1), sysRdaFolder = .myFolderA, overwrite=TRUE, translate=FALSE)
file.copy(fnA,fn1, overwrite=T)
RobExtremes:::.recomputeInterpolators(c(fn3,fn1), sysRdaFolder = .myFolderA, overwrite=TRUE, translate=FALSE)
nE= new.env()
load(fnA,envir=nE)
w = ls(all=T,envir=nE)
lapply(w, function(x) {u=get(x,envir=nE); print(x);print(names(u))})
#lapply(grep("\\.N$",w,val=T), function(x) {u=get(x,envir=nE); for(i in 1:length(u)){if(length(u)<4){ print(u[[i]]$fct[[1]](0.3)); print(u[[i]]$fct[[2]](0.3))}else{print(u[[i]]$fct(0.3))}}})

.basepath <- "C:/rtest/RobASt/branches/robast-0.9/pkg"
.myFolderA <- file.path(.basepath,"RobExtremesBuffer/all2")
fnA0 <- file.path(.myFolderA,"sysdata0.rda")
fnA <- file.path(.myFolderA,"sysdata.rda")
file.copy(fnA,fnA0, overwrite=T)
#fn2=file.path(.myFoldera,"sysdata-1.rda")
require(RobExtremes); RobExtremes:::.recomputeInterpolators(fnA0, sysRdaFolder = .myFolderA)
nE= new.env()
load(fnA,envir=nE)
w = ls(all=T,envir=nE)
lapply(w, function(x) {u=get(x,envir=nE); print(x);print(names(u))})
#lapply(grep("\\.O$",w,val=T), function(x) {u=get(x,envir=nE); for(i in 1:length(u)){if(length(u)<4){ print(u[[i]]$fct[[1]](0.3)); print(u[[i]]$fct[[2]](0.3))}else{print(u[[i]]$fct(0.3))}}})
fu(c(0.05,0.2,0.7,1,1.5,2,4))
