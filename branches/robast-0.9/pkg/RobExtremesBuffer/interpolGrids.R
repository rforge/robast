###generate first InterpolGrids
if(FALSE){
### if grid were dense enough....
.myFolder <- "C:/rtest/RobASt/branches/robast-0.9/pkg/ROptEst/R"
.myfolData <- "D:/Documents/Arbeit/Projekte/NataliyaDiss/Paper R-files/Extrapolation"
newEnv1 <- new.env()
load(file.path(.myfolData,"ExtrapolationMC_RMX1.RData"),envir=newEnv1)
whatIn <- ls(envir=newEnv1, all.names=TRUE)
xi <- get("xi.g",envir=newEnv1)
Y <- get("Y",envir=newEnv1)
rm(list=whatIn,envir=newEnv1)
gc()
.saveInterpGrid(xiGrid = xi,
                PFam = GParetoFamily(),
                sysRdaFolder = .myFolder,
                getFun = getLMGrid, optFct = .RMXE.xi, nameInSysdata = ".RMXE",
                withSmooth = TRUE, withPrint = TRUE,
                Y=Y, elseFun= .MakeGridList)

newEnv1 <- new.env()
.myfolData <- "D:/Documents/Arbeit/Projekte/NataliyaDiss/Paper JOP"
load(file.path(.myfolData,"ExtrapolationMC_OMSE.RData"),envir=newEnv1)
whatIn <- ls(envir=newEnv1, all.names=TRUE)
xi0 <- get("xi.g",envir=newEnv1)
Y0 <- get("Y",envir=newEnv1)
rm(list=whatIn,envir=newEnv1)
gc()
load(file.path(.myfolData,"ExtrapolationMC_OMSE_infmean.RData"),envir=newEnv1)
whatIn <- ls(envir=newEnv1, all.names=TRUE)
xi1 <- get("xi.g",envir=newEnv1)
Y1 <- get("Y",envir=newEnv1)
rm(list=whatIn,envir=newEnv1)
gc()
xi <- rbind(xi0,xi1)
on <- order(xi)
xi <- xi[on]
Y <- (rbind(Y0[,1,],Y1[,1,]))[on,]
.saveInterpGrid(xiGrid = xi,
                PFam = GParetoFamily(),
                sysRdaFolder = .myFolder,
                getFun = getLMGrid, optFct = .OMSE.xi, nameInSysdata = ".OMSE",
                withSmooth = TRUE, withPrint = TRUE,
                Y=Y, elseFun= .MakeGridList)

newEnv1 <- new.env()
load(file.path(.myfolData,"ExtrapolationMC_OBRE.RData"),envir=newEnv1)
whatIn <- ls(envir=newEnv1, all.names=TRUE)
xi0 <- get("xi.g",envir=newEnv1)
Y0 <- get("Y",envir=newEnv1)
rm(list=whatIn,envir=newEnv1)
gc()
load(file.path(.myfolData,"ExtrapolationMC_OBRE_infmean.RData"),envir=newEnv1)
whatIn <- ls(envir=newEnv1, all.names=TRUE)
xi1 <- get("xi.g",envir=newEnv1)
Y1 <- get("Y",envir=newEnv1)
rm(list=whatIn,envir=newEnv1)
gc()
xi <- rbind(xi0,xi1)
on <- order(xi)
xi <- xi[on]
Y <- (rbind(Y0[,1,],Y1[,1,]))[on,]
.saveInterpGrid(xiGrid = xi,
                PFam = GParetoFamily(),
                sysRdaFolder = .myFolder,
                getFun = getLMGrid, optFct = .MBRE.xi, nameInSysdata = ".MBRE",
                withSmooth = TRUE, withPrint = TRUE,
                Y=Y, elseFun= .MakeGridList)
}
