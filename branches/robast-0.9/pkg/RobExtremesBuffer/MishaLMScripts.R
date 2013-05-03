### preparations:
# (0) R-forge checkout von distr und robast machen; Pakete installieren
######
# Reihenfolge
#### *: von r-forge, **: von CRAN, ***: von BioConductor
# vorab:
# CRAN: **  sfsmisc, setRNG, fBasics, fGarch, mvtnorm, lattice, RColorBrewer
# BioConductor: *** Biobase, affy, beadarray
# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite(c("affy", "beadarray"))
#
# *   RobAStRDA
# *   startupmsg
# *   SweaveListingUtils
# *   distr
# *   distrEx
# *   distrTeach
# *   distrRmetrics
# *   distrSim
# *   distrEllipse
# *   distrTEst
# *   RandVar
# *   distrMod
# *   distrDoc
# *   RobAStBase
# *   ROptEst
# *   RobExtremes
# *   RobLox
# *   RobLoxBioC
# *   ROptEstOld
# *   ROptRegTS
# *   RobRex
#
## evtl naechste Zeile modifizieren
baseDir0 <- "D:/SVN repositories/robast"
interpolDir <- "branches/robast-0.9/pkg/RobExtremes/inst/AddMaterial/interpolation"
interpolFile <- "plotInterpol.R"
##
# (1) Paket laden
# sourceDir <- function(path, trace = TRUE, ...) {
#   for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
#     if(trace) cat(nm,":")
#     try(source(file.path(path, nm), ...))
#     if(trace) cat("\n")
#   }
# }
# 
# path = "D:/SVN repositories/robast/branches/robast-0.9/pkg/RandVar/R/"
# setwd(path)
# sourceDir(path)
# 
# path = "D:/SVN repositories/robast/branches/robast-0.9/pkg/RobAStBase/R/"
# setwd(path)
# sourceDir(path)
# 
# path = "D:/SVN repositories/robast/branches/robast-0.9/pkg/RobAStRDA/R/"
# setwd(path)
# sourceDir(path)
# 
# path = "D:/SVN repositories/robast/branches/robast-0.9/pkg/RobExtremes/R/"
# setwd(path)
# sourceDir(path)
# 
# path = "D:/SVN repositories/robast/branches/robast-0.9/pkg/RobExtremesBuffer/R/"
# setwd(path)
# sourceDir(path)
# 
# path = "D:/SVN repositories/robast/branches/robast-0.9/pkg/RobLox/R/"
# setwd(path)
# sourceDir(path)
# 
# path = "D:/SVN repositories/robast/branches/robast-0.9/pkg/RobLoxBioC/R/"
# setwd(path)
# sourceDir(path)
# 
# path = "D:/SVN repositories/robast/branches/robast-0.9/pkg/RobRex/R/"
# setwd(path)
# sourceDir(path)
# 
# path = "D:/SVN repositories/robast/branches/robast-0.9/pkg/ROptEst/R/"
# setwd(path)
# sourceDir(path)
# 
# path = "D:/SVN repositories/robast/branches/robast-0.9/pkg/ROptEstOld/R/"
# setwd(path)
# sourceDir(path)
# 
# path = "D:/SVN repositories/robast/branches/robast-0.9/pkg/ROptReg/R/"
# setwd(path)
# sourceDir(path)
# 
# path = "D:/SVN repositories/robast/branches/robast-0.9/pkg/ROptRegTS/R/"
# setwd(path)
# sourceDir(path)


require(RobExtremes)
##
## in \branches\robast-0.9\pkg\RobExtremes\inst\AddMaterial\interpolation
## file plotInterpol.R einsourcen
source(file.path(baseDir0,interpolDir, interpolFile))

### .saveGridToRDA und .computeInterpolators aus Namespace holen:
.saveGridToRda <- RobAStRDA:::.saveGridToRda
.computeInterpolators <- RobAStRDA:::.computeInterpolators

## Risiken auf P+M+B+G+MP+D (jeder 22)
#P OMSE.GEV, OMSE.Gamma
#MP MBRE.GEV, MBRE.Gamma,
#M RMXE.GEV, RMXE.Gamma
#G OMSE.GPD, OMSE.Weibull
#D MBRE.GPD MBRE.Weibull
#B RMXE.GPD RMXE.Weibull

## in den Plots: schwarz: ungeglättet;
##               rot: bereits im Gitter vorhandene Glättung;
##               grün: aktuelle TestGlättung

## Definition von Shortcuts
## Peter: / bei Euch entsprechend erste beide Argumente von myplot2, myplot3, zu ersetzen
myplot2 <- function(whichLM, plotGridRestriction = NULL,
               df = NULL, gridRestrForSmooth = NULL, withSmooth=TRUE, ...)
       plotLM("MBRE",Famnam="GEV",whichLM=whichLM, baseDir=baseDir0, withSmooth=withSmooth,
               plotGridRestriction=plotGridRestriction,
               smoothtry = TRUE, df = df,
               gridRestrForSmooth = gridRestrForSmooth, ...)
myplot3 <- function(whichLM, plotGridRestriction = NULL,
               df = NULL, gridRestrForSmooth = NULL, withSmooth=TRUE, ...)
       plotLM("MBRE",Famnam="Gam",whichLM=whichLM, baseDir=baseDir0, withSmooth=withSmooth,
               plotGridRestriction=plotGridRestriction,
               smoothtry = TRUE, df = df,
               gridRestrForSmooth = gridRestrForSmooth, ...)

### folder setzen
oldwd <- getwd()
.basepath <- file.path(baseDir0, "branches/robast-0.9/pkg")
.myFolderFrom <- file.path(.basepath,"RobExtremesBuffer")
### Zwischenspeichern des rda-files
myRDA1 <- file.path(.basepath,"RobExtremesBuffer/sysdata.rda")
### Endort des rda-files
myRDA <- file.path(.basepath,"RobAStRDA/R/sysdata.rda")
CSVFiles <- grep("\\.csv$", dir(.myFolderFrom), value=TRUE)
CSVFiles <- paste(.myFolderFrom, CSVFiles, sep="/")
CSVFiles2 <- file.path(.myFolderFrom,"interpolMBREGEVFamily.csv")
CSVFiles3 <- file.path(.myFolderFrom,"interpolMBREGammafamily.csv")
file.copy(from=myRDA,to=myRDA1)

### 1. Runde
### "MBRE"-"GEV"
## df und gridR Werte durch Ausprobieren gewonnen
myplot2(1, df = 10, gridR = -(1:270))
myplot2(2, df = 12, gridR = -(1:270))
myplot2(3, df = 10, gridR = -(1:270))
myplot2(4, df = 10, gridR = -(1:270))
myplot2(5, df = 10, gridR = -(1:270))
myplot2(6, df = 20, gridR = -(1:270))
myplot2(7, df = 20, gridR = -(1:270))
myplot2(8, df = 20, gridR = -(1:270))
myplot2(9, df = 20, gridR = -(1:270))
myplot2(10, df = 20, gridR = -(1:270))
myplot2(11, df = 20, gridR = -(1:270))
myplot2(12, df = 20, gridR = -(1:270))
myplot2(13, df = 20, gridR = -(1:270))

### sammeln der gridR und df Werte (ggf in listen)
gridR2 <- -(1:275)
dfR2 <- 20

### alle Plotten zur Kontrolle
myplot2("all", df=20, gridR=gridR2, withSmooth=FALSE, pre=windows())

### schreiben der geglätteten Gitter ins rda-file,
##      aber zunächst noch woanders (myRDA1) gespeichert:
.saveGridToRda(CSVFiles2, toFileRDA = myRDA1, withMerge = TRUE,
               withPrint = TRUE, withSmooth = TRUE, df = dfR2,
               gridRestrForSmooth=gridR2)

### 1. Runde
### "MBRE"-"Gamma"
## df und gridR Werte durch Ausprobieren gewonnen
myplot3(1, df = 4, gridR = -(1:260), plotG=-(1:20))
myplot3(2, df = 4, gridR = -(1:260), plotG=-(1:10))
myplot3(3, df = 4, gridR = -(1:260), plotG=-(1:20))
myplot3(4, df = 4, gridR = -(1:260))
myplot3(5, df = 4, gridR = -(1:260), plotG=-(1:20))
myplot3(6, df = 5, gridR = -(1:150), plotG=-(1:20), withSmooth=FALSE)
myplot3(7, df = 2, gridR = -(1:260))
myplot3(8, df = 2, gridR = -(1:260), plotG=-(1:20))
myplot3(9, df = 5, gridR = -(1:260), plotG=-(1:10))
myplot3(10, df = 5, gridR = -(1:150), plotG=-(1:20), withSmooth=FALSE)
myplot3(11, df = 2, gridR = -(1:260), plotG=-(1:10))
myplot3(12, df = 2, gridR = -(1:260), plotG=-(1:20))
myplot3(13, df = 4, gridR = -(1:260), plotG=-(1:10))

### sammeln der gridR, plotR  und df Werte (ggf in listen)
plotR3 <- list(-(1:20),-(1:20),-(1:10),NULL,-(1:20),
                                  -(1:20), NULL, -(1:20), -(1:10),-(1:20),
                                  -(1:10),-(1:20),-(1:20))
gridR3 <- list(-(1:260),-(1:260),-(1:260),-(1:260),-(1:260),-(1:150),
             -(1:260),-(1:260),-(1:260),-(1:150),-(1:260),-(1:260),-(1:260))
dfR3 <- c(4,4,4,4,4,5,2,2,5,5,2,2,4)
### alle Plotten zur Kontrolle
myplot3("all", df=dfR3, gridR=gridR3, plotG=plotR3, withSmooth=FALSE, pre=windows())
### schreiben der geglätteten Gitter ins rda-file,
##      aber zunächst noch woanders (myRDA1) gespeichert:
.saveGridToRda(CSVFiles3, toFileRDA = myRDA1, withMerge = TRUE,
               withPrint = TRUE, withSmooth = TRUE, df = dfR3,
               gridRestrForSmooth=gridR3)

if(getRversion()>"2.16"){
### generierung der Interpolatoren (in R>3.0)
.computeInterpolators(myRDA1, myRDA,withSmoothFct = TRUE)
}
######################################---bis hierher mit R-3.0.0 laufen lassen ##

######################################---ab hier mit R-2.15.2 laufen lassen ##
if(getRversion()<"2.16"){
### generierung der Interpolatoren (in R<=2.15)
## folgenden Code einsourcen:
### change adequately:
.baseDir.loc <- "C:/rtest/RobASt"
.basepath <- file.path(.baseDir.loc,"branches/robast-0.9/pkg")
myRDA <- file.path(.basepath,"RobAStRDA/R/sysdata.rda")
require(RobAStRDA)
RobAStRDA:::.computeInterpolators(myRDA, myRDA,withSmoothFct = TRUE)
}

