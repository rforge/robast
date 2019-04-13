### preparations:
# (0) R-forge checkout von distr und robast machen; Pakete installieren
######
# Reihenfolge
#### *: von r-forge, **: von CRAN, ***: von BioConductor
# vorab:
# CRAN: **  sfsmisc, setRNG, fBasics, fGarch, mvtnorm, lattice, RColorBrewer
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
baseDir0 <- "C:/rtest/RobASt"
interpolDir <- "branches/robast-1.1/pkg/RobExtremes/inst/AddMaterial/interpolation"
interpolFile <- "plotInterpol.R"
##
# (1) Paket laden
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
       plotLM("OMSE",Famnam="GEV Fam",whichLM=whichLM, baseDir=baseDir0, withSmooth=withSmooth,
               plotGridRestriction=plotGridRestriction,
               smoothtry = TRUE, df = df,
               gridRestrForSmooth = gridRestrForSmooth, ...)
myplot3 <- function(whichLM, plotGridRestriction = NULL,
               df = NULL, gridRestrForSmooth = NULL, withSmooth=TRUE, ...)
       plotLM("OMSE",Famnam="Gam",whichLM=whichLM, baseDir=baseDir0, withSmooth=withSmooth,
               plotGridRestriction=plotGridRestriction,
               smoothtry = TRUE, df = df,
               gridRestrForSmooth = gridRestrForSmooth, ...)

### folder setzen
oldwd <- getwd()
.basepath <- file.path(baseDir0, "branches/robast-1.0/pkg")
.myFolderFrom <- file.path(.basepath,"RobExtremesBuffer")
### Zwischenspeichern des rda-files
myRDA1 <- file.path(.basepath,"RobExtremesBuffer/sysdata.rda")
### Endort des rda-files
myRDA <- file.path(.basepath,"RobAStRDA/R/sysdata.rda")
CSVFiles <- grep("\\.csv$", dir(.myFolderFrom), value=TRUE)
CSVFiles <- paste(.myFolderFrom, CSVFiles, sep="/")
CSVFiles2 <- file.path(.myFolderFrom,"interpolOMSEGEVFamily.csv")
CSVFiles3 <- file.path(.myFolderFrom,"interpolOMSEGammafamily.csv")
file.copy(from=myRDA,to=myRDA1)

### 1. Runde
### "OMSE"-"GEV"
## df und gridR Werte durch Ausprobieren gewonnen

myplot2(1, gridR=c(1:183,228:234,242:286,293), withS=FALSE)
myplot2(2, gridR=c(1:183,228:234,242:286,293), df=30,plotG=1:290, withSmooth=FALSE)
myplot2(3, gridR=c(1:183,228:234,242:286,293), df=30,plotG=1:290, withSmooth=FALSE)
myplot2(4, gridR=c(1:183,228:234,242:286,293), df=30,plotG=1:290, withSmooth=FALSE)
myplot2(5, gridR=c(1:183,228:234,242:286,293), df=30,plotG=1:290, withSmooth=FALSE)
myplot2(6, gridR=c(1:183,228:234,242:286,293), df=30,plotG=1:290, withSmooth=FALSE)
myplot2(7, gridR=c(1:183,228:234,242:286,293), df=14,plotG=1:290, withSmooth=FALSE)
myplot2(8, gridR=c(1:183,228:234,242:286,293), df=14,plotG=1:290, withSmooth=FALSE)
myplot2(9, gridR=c(1:183,228:234,242:286,293), df=18,plotG=1:290, withSmooth=FALSE)
myplot2(10, gridR=c(1:183,228:234,242:286,293), df=18,plotG=1:290, withSmooth=FALSE)
myplot2(11, gridR=c(1:183,228:234,242:286,293), df=8,plotG=1:290, withSmooth=FALSE)
myplot2(12, gridR=c(1:183,228:234,242:286,293), df=8,plotG=1:290, withSmooth=FALSE)
myplot2(13, gridR=c(1:183,228:234,242:286,293), df=18,plotG=1:290, withSmooth=FALSE)


### sammeln der gridR und df Werte (ggf in listen)
gridR20 <- c(1:183,228:234,242:286,293)
gridR2 <- list(gridR20,gridR20,gridR20,gridR20,gridR20,gridR20,gridR20,gridR20,
               gridR20,gridR20,gridR20,gridR20,gridR20)
dfR2 <- list(NULL,30,30,30,30,30,14,14,18,18,8,8,18)

### alle Plotten zur Kontrolle
myplot2("all", df=20, gridR=gridR2, withSmooth=FALSE, pre=windows())

### schreiben der geglätteten Gitter ins rda-file,
##      aber zunächst noch woanders (myRDA1) gespeichert:
.saveGridToRda(CSVFiles2, toFileRDA = myRDA1, withMerge = FALSE,
               withPrint = TRUE, withSmooth = TRUE, df = dfR2,
               gridRestrForSmooth=gridR2)

### 1. Runde
### "OMSE"-"Gamma"
## df und gridR Werte durch Ausprobieren gewonnen
gridR30=-c(1:2,176:179,183:213,268,276,284:286,295:296,300,305,310:324,339:345,351:362,368:415,419:472,508:514,520:528,523:543)
myplot3(1, gridR=gridR30,df=95, plotG=2:530, withSmooth=FALSE, pre=substitute(print(gr0[,c(1,3)])))
myplot3(2, gridR=gridR30,df=65, plotG=2:530, withSmooth=FALSE, pre=substitute(print(gr0[,c(1,3)])))
myplot3(3, gridR=gridR30,df=55, plotG=2:530, withSmooth=FALSE, pre=substitute(print(gr0[,c(1,3)])))
myplot3(4, gridR=gridR30,df=50, plotG=1:30, withSmooth=FALSE, pre=substitute(print(gr0[,c(1,3)])))
myplot3(5, gridR=gridR30,df=85, plotG=1:30, withSmooth=FALSE, pre=substitute(print(gr0[,c(1,3)])))
myplot3(6, gridR=gridR30,df=85,withSmooth=FALSE,plotG=160:300, pre=substitute(print(gr0[,c(1,7)])))
myplot3(7, gridR=gridR30,withSmooth=FALSE,plotG=160:300, df=29, pre=substitute(print(gr0[,c(1,8)])))
myplot3(8, gridR=gridR30,withSmooth=FALSE,plotG=160:300, df=29, pre=substitute(print(gr0[,c(1,8)])))
myplot3(9, gridR=gridR30,withSmooth=FALSE,plotG=160:300, df=45, pre=substitute(print(gr0[,c(1,8)])))
myplot3(10,gridR=gridR30,withSmooth=FALSE,plotG=1:30, df=80, pre=substitute(print(gr0[,c(1,11)])))
myplot3(11,gridR=gridR30,withSmooth=FALSE,plotG=1:30, df=40, pre=substitute(print(gr0[,c(1,11)])))
myplot3(12,gridR=gridR30,withSmooth=FALSE,plotG=1:30, df=40, pre=substitute(print(gr0[,c(1,11)])))
myplot3(13,gridR=gridR30,withSmooth=FALSE,plotG=1:30, df=90, pre=substitute(print(gr0[,c(1,11)])))

### sammeln der gridR, plotR  und df Werte (ggf in listen)
gridR3 <- list(gridR30,gridR30,gridR30,gridR30,gridR30,gridR30,gridR30,gridR30,
               gridR30,gridR30,gridR30,gridR30)
dfR3 <- list(95,65,55,50,85,85,29,29,45,80,40,40,90)

### sammeln der gridR, plotR  und df Werte (ggf in listen)
### alle Plotten zur Kontrolle
myplot3("all", df=dfR3, gridR=gridR3, plotG=-c(1:2), withSmooth=FALSE, pre=windows())
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
.basepath <- file.path(.baseDir.loc,"branches/robast-1.1/pkg")
myRDA <- file.path(.basepath,"RobAStRDA/R/sysdata.rda")
require(RobAStRDA)
RobAStRDA:::.computeInterpolators(myRDA, myRDA,withSmoothFct = TRUE)
}

