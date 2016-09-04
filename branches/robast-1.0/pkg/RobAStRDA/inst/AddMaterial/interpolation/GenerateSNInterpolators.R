### session to supply Sn Grids on Sep 04 2016

require(RobAStRDA)
require(RobExtremes)

s = getFromNamespace(".RMXE", ns = "RobAStRDA")
nE = new.env()
oldwd <- getwd()
.basepath <- "C:/rtest/RobASt/branches/robast-1.0/pkg"
.myFolderFrom <- file.path(.basepath,"RobExtremesBuffer/neu")
myRDA0 <- file.path(.basepath,"RobExtremesBuffer/sysdata.rda")
myRDA <- file.path(.basepath,"RobAStRDA/R/sysdata.rda")

load(myRDA,nE)
ls(nE,all=T)
ga=get(".Generalized",nE)
names(ga)

nEalt = new.env()
load("C:/rtest/robast/branches/robast-0.9/pkg/RobAStRDA/R/sysdata.rda",nEalt)
ls(nEalt, all=TRUE)

SN=get(".Sn",nEalt)
names(SN)
SN[[1]]
SN[[1]][[1]]
names(SN)

listGrid <- list()
listGrid[[1]] <- SN[[1]][[1]]
listGrid[[2]] <- SN[[2]][[1]]
listGrid[[3]] <- SN[[3]][[1]]
listGrid[[4]] <- SN[[4]][[1]]
names(listGrid) <- names(SN)

xiseq <- seq(-0.5,-0.00001,length.out=150)


SnV <- numeric(150)
for(i in 1:150){
  Fam <- GPareto(shape=xiseq[i])
  print(i)
  SnV[i] <- Sn(Fam)
}

newGrid <- rbind(cbind(xiseq,SnV),listGrid[[2]])
newGrid <- cbind(newGrid,predict(smooth.spline(newGrid[,1], newGrid[,2]),newGrid[,1])$y)

matplot(newGrid[,1], newGrid[,c(2,3)],type="l", log="y")

SnV2 <- SnV
for(i in 1:150){
  Fam <- GEV(shape=xiseq[i])
  print(i)
  SnV2[i] <- Sn(Fam)
}

newGrid2 <- rbind(cbind(xiseq,SnV2),listGrid[[3]])
newGrid2 <- cbind(newGrid2,predict(smooth.spline(newGrid2[,1], newGrid2[,2]),newGrid2[,1])$y)

matplot(newGrid2[,1], newGrid2[,c(2,3)],type="l", log="y")

newGrid0 <- listGrid[[1]]
newGrid0 <- cbind(newGrid0,predict(smooth.spline(newGrid0[,1], newGrid0[,2]),newGrid0[,1])$y)
matplot(newGrid0[,1], newGrid0[,c(2,3)],type="l", log="y")


newGrid3 <- listGrid[[4]]
newGrid3 <- cbind(newGrid3,predict(smooth.spline(newGrid3[,1], newGrid3[,2]),newGrid3[,1])$y)
matplot(newGrid3[,1], newGrid3[,c(2,3)],type="l", log="y")

save(newGrid0, newGrid, newGrid2, newGrid3, listGrid, file="SnGridAll.rda")

tu <- function(Grid) {gr <- RobAStRDA:::.generateInterpolators(Grid[,c(1,3)])
					  gr0 <- list()
					  gr0[["grid"]] <- Grid[,c(1,2)]
					  gr0[["gridS"]] <- Grid[,c(1,3)]
					  gr0[["fun.N"]] <- gr$fct
					  return(gr0)}


GammaF <- tu(newGrid0)
GPDF <- tu(newGrid)
GEVF <- tu(newGrid2)
WeibullF <- tu(newGrid3)

.basepath <- "C:/rtest/RobASt/branches/robast-1.0/pkg"
myRDA <- file.path(.basepath,"RobAStRDA/R/sysdata.rda")
nE = new.env()
load(myRDA,nE)
ls(nE,all=T)

nams <- ls(nE,all=TRUE)
zw <- get(nams[1],nE)
zw[["Sn"]] <- GammaF
assign(nams[1],zw,envir=nE)
zw <- get(nams[2],nE)
zw[["Sn"]] <- GPDF
assign(nams[2],zw,envir=nE)
zw <- get(nams[3],nE)
zw[["Sn"]] <- GEVF
assign(nams[3],zw,envir=nE)
zw <- get(nams[4],nE)
zw[["Sn"]] <- GEVF
assign(nams[4],zw,envir=nE)
zw <- get(nams[5],nE)
zw[["Sn"]] <- WeibullF
assign(nams[5],zw,envir=nE)

zw <- get(nams[2],nE)
assign(".GPareto",zw,envir=nE)
rm(".Generalized",envir=nE)
(namesNew <- ls(nE,all=T))

save(list=namesNew, envir=nE, file=myRDA)
tools::resaveRdaFiles(myRDA)
