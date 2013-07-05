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
#Peter
baseDir0 <- "C:/rtest/RobASt"
#Gerald -> bitte checken
#baseDir0 <- "/home/bernhard/university/svn/r-forge/robast"
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
#setwd(file.path(baseDir0,"branches/robast-0.9/pkg/RobExtremesBuffer"))
source(file.path(baseDir0,"branches/robast-0.9/pkg/RobExtremesBuffer","MakeSmoothGridListBe.R"))
###############################################


### .saveGridToRDA und .computeInterpolators aus Namespace holen:
.saveGridToRda <- RobAStRDA:::.saveGridToRda
.computeInterpolators <- RobAStRDA:::.computeInterpolators

## Risiken auf P+M+B+G+MP+D (jeder 22)
#P OMSE.GEV, OMSE.Gamma   oK
#MP MBRE.GEV, MBRE.Gamma, oK
#M RMXE.GEV, RMXE.Gamma   oK
#G OMSE.GPD, OMSE.Weibull oK
#D MBRE.GPD MBRE.Weibull
#B RMXE.GPD RMXE.Weibull

## in den Plots: schwarz: ungeglättet;
##               rot: bereits im Gitter vorhandene Glättung;
##               grün: aktuelle TestGlättung

## Definition von Shortcuts
## Peter: / bei Euch entsprechend erste beide Argumente von myplot2, myplot3, zu ersetzen
myplot2 <- function(whichLM, plotGridRestriction = NULL,
               df = NULL, gridRestrForSmooth = NULL, withSmooth=TRUE, ...)
       plotLM("OMSE",Famnam="Generalized",whichLM=whichLM, baseDir=baseDir0, withSmooth=withSmooth,
               plotGridRestriction=plotGridRestriction,
               smoothtry = TRUE, df = df,
               gridRestrForSmooth = gridRestrForSmooth, ...)
myplot3 <- function(whichLM, plotGridRestriction = NULL,
               df = NULL, gridRestrForSmooth = NULL, withSmooth=TRUE, ...)
       plotLM("OMSE",Famnam="Weibull",whichLM=whichLM, baseDir=baseDir0, withSmooth=withSmooth,
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
CSVFiles2 <- file.path(.myFolderFrom,"interpolOMSEGeneralizedParetoFamily.csv")
CSVFiles3 <- file.path(.myFolderFrom,"interpolOMSEWeibullFamily.csv")
file.copy(from=myRDA,to=myRDA1)

### 1. Runde
### "OMSE"-"GPD"
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

myplot2(1, df=8)
myplot2(2, withSmooth=FALSE, df=8)
myplot2(3, df=8)
myplot2(4)
myplot2(5)
myplot2(6)
myplot2(7, df=8)
myplot2(8, df=8)
myplot2(9, df=14)
myplot2(10, df=6)
myplot2(11, df=8)
myplot2(12, df=8)
myplot2(13, df=14)

### sammeln der gridR und df Werte (ggf in listen)
gridR2 <- list(NULL, NULL, NULL, NULL, NULL, NULL,
               NULL, NULL, NULL, NULL, NULL, NULL,
               NULL)
dfR2 <- list(8,8,8,NULL,NULL,NULL,8,8,14,6,8,8,14)

### alle Plotten zur Kontrolle
## linux:
# myplot2("all", df=dfR2, gridR=gridR2, withSmooth=FALSE, pre=x11())
myplot2("all", df=dfR2, gridR=gridR2, withSmooth=FALSE, pre=windows())
myplot2("all", df=dfR2, gridR=gridR2, withSmooth=FALSE,  pre=pdf("GPD-OMSE-s-Gerald.pdf"), post=dev.off())

### schreiben der geglätteten Gitter ins rda-file,
##      aber zunächst noch woanders (myRDA1) gespeichert:
.saveGridToRda(CSVFiles2, toFileRDA = myRDA1, withMerge = FALSE,
               withPrint = TRUE, withSmooth = TRUE, df = dfR2,
               gridRestrForSmooth=gridR2)

### 2. Runde
### "OMSE"-"Weibull"
## df und gridR Werte durch Ausprobieren gewonnen

#myplot3(1, df=25, gridR=-c(1:3, 200:670), plotG=3:200)
#myplot3(1, df=12, gridR=-c(1:3, 100:670), plotG=3:200)

myplot3(1, plotG=-c(1:15), gridR=-c(1:18), withS=FALSE)
myplot3(2, plotG=-c(1:28), gridR=-c(1:25), withS=FALSE)
myplot3(3, df=8, withS=FALSE)
myplot3(4, withS=FALSE)
myplot3(5, plotG=-c(1:15), gridR=-c(1:18), withS=FALSE)
myplot3(6, plotG=-c(1:20), gridR=-c(1:22), withS=FALSE)
myplot3(7, plotG=-c(1:20), gridR=-c(1:10), withS=FALSE)
myplot3(8, plotG=-c(1:20), gridR=-c(1:10), withS=FALSE)
myplot3(9, plotG=-c(1:10), gridR=-c(1:3), withS=FALSE, xlim=c(0,2), ylim=c(0,3))
myplot3(10, plotG=-c(1:1), gridR=-c(1:10), withS=FALSE, xlim=c(0,6), ylim=c(0,5))
myplot3(11, plotG=-c(1:20), gridR=-c(1:10), withS=FALSE)
myplot3(12, plotG=-c(1:20), gridR=-c(1:10), withS=FALSE)
myplot3(13, plotG=-c(1:10), gridR=-c(1:3), withS=FALSE, xlim=c(0,2), ylim=c(0,3))


### sammeln der gridR, plotR  und df Werte (ggf in listen)
gridR3 <- list(-(1:18),-(1:25),NULL,NULL,-(1:18),-(1:22),-(1:10),-(1:10),-(1:3),
               -(1:10), -(1:10), -(1:10),-(1:3))
dfR3 <- list(NULL,NULL,8,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
             NULL,NULL,NULL)
### alle Plotten zur Kontrolle
#linux
#myplot3("all", df=dfR3, gridR=gridR3, withSmooth=FALSE, pre=X11())
myplot3("all", df=dfR3, gridR=gridR3, withSmooth=FALSE, pre=windows())
myplot3("all", df=dfR3, gridR=gridR3, withSmooth=FALSE, pre=windows(), plotG=-c(1:15))
myplot3("all", df=dfR3, gridR=gridR3, withSmooth=FALSE, plotG=-c(1:15), pre=pdf("Weibull-OMSE-s-Gerald.pdf"), post=dev.off())
### schreiben der geglätteten Gitter ins rda-file,
##      aber zunächst noch woanders (myRDA1) gespeichert:
.saveGridToRda(CSVFiles3, toFileRDA = myRDA1, withMerge = FALSE,
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

