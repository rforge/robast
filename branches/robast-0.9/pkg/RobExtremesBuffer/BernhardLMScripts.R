### preparations:
# (0) R-forge checkout von distr und robast machen; Pakete installieren
######
# Reihenfolge
#### *: von r-forge, **: von CRAN, ***: von BioConductor
# vorab:
# CRAN: **  sfsmisc, setRNG, fBasics, fGarch, mvtnorm, lattice, RColorBrewer,
#           rrcov, evd, actuar
# BioConductor: *** Biobase, affy, beadarray
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
baseDir0 <- "/home/bernhard/university/svn/r-forge/robast"
interpolDir <- "branches/robast-0.9/pkg/RobExtremes/inst/AddMaterial/interpolation"
interpolFile <- "plotInterpol.R"
##
# (1) Paket laden
require(RobExtremes)
##
## in \branches\robast-0.9\pkg\RobExtremes\inst\AddMaterial\interpolation
## file plotInterpol.R einsourcen
source(file.path(baseDir0,interpolDir, interpolFile))

###############################################
###  WARNING: new .MakeSmoothGridList()!!!  ###
source("MakeSmoothGridListBe.R")
###############################################


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
       plotLM("RMXE",Famnam="Generalized",whichLM=whichLM, baseDir=baseDir0, withSmooth=withSmooth,
               plotGridRestriction=plotGridRestriction,
               smoothtry = TRUE, df = df,
               gridRestrForSmooth = gridRestrForSmooth, ...)
myplot3 <- function(whichLM, plotGridRestriction = NULL,
               df = NULL, gridRestrForSmooth = NULL, withSmooth=TRUE, ...)
       plotLM("RMXE",Famnam="Weibull",whichLM=whichLM, baseDir=baseDir0, withSmooth=withSmooth,
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
CSVFiles2 <- file.path(.myFolderFrom,"interpolRMXEGeneralizedParetoFamily.csv")
CSVFiles3 <- file.path(.myFolderFrom,"interpolRMXEWeibullFamily.csv")
file.copy(from=myRDA,to=myRDA1)

### 1. Runde
### "RMXE"-"GPD"
## df und gridR Werte durch Ausprobieren gewonnen

#myplot2(1, df = NULL, gridR = NULL)
#myplot2(2, df = NULL, gridR = NULL)
#myplot2(3, df = NULL, gridR = NULL)
#myplot2(4, df = NULL, gridR = NULL)
#myplot2(5, df = NULL, gridR = NULL)
#myplot2(6, df = NULL, gridR = NULL)
#myplot2(7, df = NULL, gridR = NULL)
#myplot2(8, df = NULL, gridR = NULL)
#myplot2(9, df = NULL, gridR = NULL)
#myplot2(10, df = NULL, gridR = NULL)
#myplot2(11, df = NULL, gridR = NULL)
#myplot2(12, df = NULL, gridR = NULL)
#myplot2(13, df = NULL, gridR = NULL)

myplot2(1)
myplot2(2, gridR=-(106:180), plotG=1:300, withSmooth=FALSE)
myplot2(3, gridR=-(106:180), plotG=1:300, withSmooth=FALSE)
myplot2(4, gridR=-(106:180), plotG=1:300, withSmooth=FALSE)
myplot2(5)
myplot2(6)
myplot2(7)
myplot2(8)
myplot2(9)
myplot2(10)
myplot2(11)
myplot2(12)
myplot2(13)

### sammeln der gridR und df Werte (ggf in listen)
gridR2 <- list(NULL, -(106:180), -(106:180), -(106:180), NULL, NULL,
               NULL, NULL, NULL, NULL, NULL, NULL,
               NULL)
dfR2 <- NULL

### alle Plotten zur Kontrolle
myplot2("all", df=dfR2, gridR=gridR2, withSmooth=FALSE, pre=X11())

### schreiben der geglätteten Gitter ins rda-file,
##      aber zunächst noch woanders (myRDA1) gespeichert:
.saveGridToRda(CSVFiles2, toFileRDA = myRDA1, withMerge = TRUE,
               withPrint = TRUE, withSmooth = TRUE, df = dfR2,
               gridRestrForSmooth=gridR2)

### 2. Runde
### "RMXE"-"Weibull"
## df und gridR Werte durch Ausprobieren gewonnen

#myplot3(1, df=25, gridR=-c(1:3, 200:670), plotG=3:200)
#myplot3(1, df=12, gridR=-c(1:3, 100:670), plotG=3:200)

myplot3(1, gridR=-(1:2), plotG=1:200, withSmooth=FALSE)
myplot3(2, df=NULL, gridR=-2, plotG=1:200, withSmooth=FALSE)
myplot3(3, gridR=-(1:2), plotG=1:100, withSmooth=FALSE)
myplot3(4)
myplot3(5, gridR=-(1:2), plotG=1:100, withSmooth=FALSE)
myplot3(6, gridR=-2, plotG=1:100, withSmooth=FALSE)
myplot3(7, gridR=-(1:2), plotG=c(1, 3:100), withSmooth=FALSE)
myplot3(8, gridR=-(1:2), plotG=c(1, 3:100), withSmooth=FALSE)
myplot3(9)
myplot3(10, gridR=-2, plotG=1:100, withSmooth=FALSE)
myplot3(11, gridR=-(1:2), plotG=1:100, withSmooth=FALSE)
myplot3(12, gridR=-(1:2), plotG=1:100, withSmooth=FALSE)
myplot3(13)

### sammeln der gridR, plotR  und df Werte (ggf in listen)
gridR3 <- list(-(1:2),    -2, -(1:2), NULL, -(1:2),     -2,
               -(1:2),-(1:2),   NULL,   -2, -(1:2), -(1:2),
               NULL)
dfR3 <- NULL
### alle Plotten zur Kontrolle
myplot3("all", df=dfR3, gridR=gridR3, withSmooth=FALSE, pre=X11())
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

