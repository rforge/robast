###############################################################################
## Illumina Example
###############################################################################

###############################################################################
## References:
## Dunning, M.J., Smith, M.L., Ritchie, M.E., Tavare, S.:
## beadarray: R classes and methods for Illumina bead-based data. 
## Bioinformatics 2007, 23(16):2183-4.
##
## M.J. Dunning, N.L. Barbosa-Morais, A.G. Lynch, S. Tavaré and M.E. Ritchie:
## Statistical issues in the analysis of Illumina data.
## BMC Bioinformatics 2008, 9:85.
###############################################################################

###############################################################################
## Data:
## Can be obtained via
## http://www.compbio.group.cam.ac.uk/Resources/spike/index.html
###############################################################################

## Load the required packages
library(beadarray)
library(gplots)
library(RColorBrewer)
library(RobLoxBioC)


###########################################################
## Read targets information
targets <- read.table("./SpikeInData/spike_targets.txt",header=TRUE)
arraynms <- as.character(targets$ArrayNo)

## Use sharpened, subtracted data from text files
spikeInData <- readIllumina(path = "./SpikeInData", arrayNames=arraynms[1:2], 
                            useImages=FALSE, textType=".csv")
#save(spikeInData, compress = TRUE, file = "spikeInData.RData")
#load(file = "spikeInData.RData")

## takes about 9 hours on Intel P9500 (64bit Linux, 4 GByte RAM)
system.time(minKD.Illumina <- KolmogorovMinDist(spikeInData, Norm(), imagesPerArray = 2))
save(minKD.Illumina, compress = TRUE, file = "minKD_Illumina.RData")

## takes about 9 hours on Intel P9500 (64bit Linux, 4 GByte RAM)
system.time(minKD.Illumina.log <- KolmogorovMinDist(spikeInData, Norm(), log = TRUE, imagesPerArray = 2))
save(minKD.Illumina.log, compress = TRUE, file = "minKD_Illumina_log.RData")

## load the results from R-forge ...
con <- url("http://robast.r-forge.r-project.org/data/minKD_Illumina.RData")
load(file = con)
close(con)
con <- url("http://robast.r-forge.r-project.org/data/minKD_Illumina_log.RData")
load(file = con)
close(con)

## takes more than 90 min on Intel P9500 (64bit Linux, 4 GByte RAM)
ns <- c(10:70)
M <- length(ns)
minKD.Illumina.norm <- matrix(NA, nrow = 50000, ncol = M)
colnames(minKD.Illumina.norm) <- ns
#for(i in seq_len(M)){
for(i in 27:61){
    tm <- proc.time()
    print(ns[i])
    temp <- matrix(rnorm(50000*ns[i]), ncol = ns[i])
    minKD.Illumina.norm[,i] <- KolmogorovMinDist(temp, Norm())$dist
    cat("Dauer:\t", proc.time()-tm, "\n")
    save(minKD.Illumina.norm, compress = TRUE, file = "minKD_Illumina_norm1.RData")
}

## load the results from R-forge
con <- url("http://robast.r-forge.r-project.org/data/minKD_Illumina_norm.RData")
load(file = con)
close(con)

#######################################
## Figure in Kohl and Deigner (2009)
#######################################
res1 <- split(as.vector(minKD.Illumina$dist), as.vector(minKD.Illumina$n))[20:60]
res2 <- split(as.vector(minKD.Illumina.log$dist), as.vector(minKD.Illumina.log$n))[20:60]
res3 <- lapply(as.data.frame(minKD.Illumina.norm[,11:51]), function(x) x)
uni.n <- rep(20:60, 3)

postscript(file = "minKDIllumina.eps", height = 6, width = 9, paper = "special", 
           horizontal = TRUE)
par(mar = c(4, 4, 3, 1))
plot(0, 0, type = "n", ylim = c(-0.01, 0.4), xlim = c(0.5, 125.5), 
     panel.first = abline(h = seq(0, 0.35, by = 0.05), lty = 2, col = "grey"), 
     main = "Minimum Kolmogorov distance", 
     ylab = "minimum Kolmogorov distance", 
     xlab = "sample size", axes = FALSE)
axis(1, c(1:41, 43:83, 85:125), labels = uni.n, cex.axis = 0.6)
axis(2, seq(0, 0.35, by = 0.05), labels = seq(0, 0.35, by = 0.05), las = 2,
     cex.axis = 0.8)
box()
boxplot(c(res1, res2, res3), at = c(1:41, 43:83, 85:125), add = TRUE, pch = 20, 
        names = FALSE, axes = FALSE)
abline(v = c(42, 84), lwd = 1.5)
text(c(20, 63, 105), rep(0.38, 3), labels = c("Bead Level Data", "log Bead Level Data", "Normal Samples"),
     font = 2)
lines(1:41, 1/(2*(20:60)), lwd = 2)
lines(43:83, 1/(2*(20:60)), lwd = 2)
lines(85:125, 1/(2*(20:60)), lwd = 2)
legend("bottomleft", legend = "minimal possible distance", lty = 1, 
       bg = "white", cex = 0.8)
dev.off()

## Comparison of median distances
## Figure in Kohl and Deigner (2009)
res1 <- split(as.vector(minKD.Illumina$dist), as.vector(minKD.Illumina$n))[10:70]
res2 <- split(as.vector(minKD.Illumina.log$dist), as.vector(minKD.Illumina.log$n))[10:70]
res3 <- lapply(as.data.frame(minKD.Illumina.norm), function(x) x)

postscript(file = "minKDIlluminaQuant.eps", height = 6, width = 9, paper = "special", 
           horizontal = TRUE)
par(mar = c(4, 4, 3, 1))
plot(10:70, sapply(res3, quantile, prob = 0.99), type = "l", lwd = 2, xlab = "sample size", 
     ylab = "quantile of mimimum Kolmogorov distances",
     main = "50% and 99% quantiles of minimum Kolmogorov distances", ylim = c(0.05, 0.23))
lines(10:70, sapply(res1, quantile, prob = 0.99), lwd = 2, lty = 2)
lines(10:70, sapply(res2, quantile, prob = 0.99), lwd = 2, lty = 3)
lines(10:70, sapply(res3, quantile, prob = 0.5), lwd = 2, lty = 1)
lines(10:70, sapply(res1, quantile, prob = 0.5), lwd = 2, lty = 2)
lines(10:70, sapply(res2, quantile, prob = 0.5), lwd = 2, lty = 3)
text(22, 0.18, "99% quantiles", font = 2)
text(22, 0.115, "50% quantiles", font = 2)
legend("topright", legend = c("normal samples", "bead level data", "log bead level data"),
       lty = 1:3, lwd = 2)
dev.off()

###############################################################################
## The following example is based on the R code of Mark Dunning and Matt Ritchie
## available under
## http://www.compbio.group.cam.ac.uk/Resources/spike/scripts/Analysis.R
##
## This file was slightly adapted and code for the computation of 
## rmx-estimators was added.
###############################################################################

########################################
## Analysis of Illumina in Spike data set
##
## November 2007
##
## Mark Dunning and Matt Ritchie
########################################

###############################################################################
## Extract all *.zip file to directory "SpikeInData".
## Copy spike_targets.txt to directory "SpikeInData".
##
## Code to read the bead level data from the directory "SpikeInData"
##
## NOTE: reading in the raw data for the entire experiment requires at
## least 4Gb of RAM for each processing method.  
## For this reason the data is read in sequentially saved, then removed
## before moving on to the next method 
###############################################################################


###########################################################
## Read targets information
targets <- read.table("./SpikeInData/spike_targets.txt",header=TRUE)


###########################################################
## Read text and tif files from all BeadChips from spike experiment
## sharpen images, no local background subtraction
bld.sharpen.nobgc <- readIllumina(path = "./SpikeInData", textType=".csv", 
                                  arrayNames=targets$ArrayNo, 
                                  imageManipulation="sharpen", 
                                  backgroundMethod="none")
save(bld.sharpen.nobgc, file="bld.sharpen.nobgc.rda")
rm(bld.sharpen.nobgc); gc()



## Take the values found in the bead level text files as sharpened, subtracted intensities
bld.sharpen.subtract <- readIllumina(path = "./SpikeInData", textType=".csv", 
                                     arrayNames=targets$ArrayNo, 
                                     useImages=FALSE, 
                                     backgroundMethod="none")
save(bld.sharpen.subtract, file="bld.sharpen.subtract.rda")
rm(bld.sharpen.subtract); gc()

## sharpen images, with local background subtraction and
## normal-exponential convolution model (avoids negative background
## corrected intensities)
bld.sharpen.normexp <- readIllumina(path = "./SpikeInData", textType=".csv", 
                                    arrayNames=targets$ArrayNo, 
                                    useImages=FALSE, 
                                    backgroundMethod="normexp")
save(bld.sharpen.normexp, file="bld.sharpen.normexp.rda")
rm(bld.sharpen.normexp); gc()


## no sharpening, no local background subtraction
bld.nosharpen.nobgc <- readIllumina(path = "./SpikeInData", textType=".csv", 
                                    arrayNames=targets$ArrayNo,
                                    imageManipulation="none", 
                                    backgroundMethod="none")
save(bld.nosharpen.nobgc, file="bld.nosharpen.nobgc.rda")
rm(bld.nosharpen.nobgc); gc()

## no sharpening with local background subtraction
bld.nosharpen.subtract <- readIllumina(path = "./SpikeInData", textType=".csv", 
                                       arrayNames=targets$ArrayNo,
                                       imageManipulation="none", 
                                       backgroundMethod="subtract")
save(bld.nosharpen.subtract, file="bld.nosharpen.subtract.rda")
rm(bld.nosharpen.subtract); gc()


###########################################################
## Figure 1 - plots of foreground and background intensities
###########################################################
## Read text and tif files from first array only (to save time and memory)
## sharpen images, no local background subtraction
B <- readIllumina(path = "./SpikeInData", textType=".csv", 
                  arrayNames=targets$ArrayNo[1:12],
                  imageManipulation="sharpen", 
                  backgroundMethod="none")

## sharpen images, with local background subtraction
B.bc <- backgroundCorrect(B, method="subtract")

## Setup colours
pal <- brewer.pal(6, "Set1")

ylim <- c(8.5,10.5)
names <- c("Array1 Strip1", "Array1 Strip2", "Array2 Strip1", "Array2 Strip2", 
           "Array3 Strip1", "Array3 Strip2", "Array4 Strip1", "Array4 Strip2", 
           "Array5 Strip1", "Array5 Strip2", "Array6 Strip1", "Array6 Strip2")

pdf("Figure1.pdf", width=12, height=8)
par(mfrow=c(1,2), mar=c(6,3,1.5,0.15), oma=c(0,1,0,0))
boxplotBeads(B, names=names, las=2, outline=FALSE, col=rep(pal, each=2), 
             main="Raw Foreground", what="G", ylim=ylim, ylab="", medlwd=1)
boxplotBeads(B, names=names, las=2, outline=FALSE, col=rep(pal, each=2), 
             main="Raw Background", what="Gb", ylim=ylim, ylab="", medlwd=1)
mtext("log2 intensity", side=2, outer=TRUE)
dev.off()


###########################################################
## Derive low-level properties of Illumina data quoted in the text
## Number of negative corrected intensities from first BeadChip
negs <- numeric(12)
for(i in 1:12)
    negs[i] <- sum(B.bc[[i]]$G<0)

## Percentage negative from first BeadChip - sharpened
round(negs/numBeads(B.bc)*100,2)
# [1] 1.04 0.76 0.93 0.68 0.92 0.68 0.94 0.82 0.84 0.74 0.90 0.82
# old figures [1] 7.87 7.53 7.32 6.95 6.94 6.29 6.18 5.91 5.54 5.25 4.80 4.92

## quartiles
quantile(round(negs/numBeads(B.bc)*100,2), c(0,0.25,0.5,0.75,1))
#  0%    25%    50%    75%   100% 
# 0.6800 0.7550 0.8300 0.9225 1.0400

## How many missing beads (contaminated left-hand edge - 0 intensities)
artefact = NULL
for(i in 1:12)
    artefact[i] = sum(B[[i]]$G==0)

## Percentage negative from first BeadChip - sharpened
round(artefact/numBeads(B)*100,2)
# [1] 7.32 7.21 6.88 6.71 6.49 6.08 5.70 5.56 5.16 4.96 4.42 4.67

## quartiles
quantile(round(artefact/numBeads(B)*100,2), c(0,0.25,0.5,0.75,1))
#    0%    25%    50%    75%   100% 
# 4.4200 5.1100 5.8900 6.7525 7.3200 

## Calculate median background level on both original and log2-scale
bgorig = NULL
bglog = NULL
for(i in 1:12) {
    bgorig[[i]] = getArrayData(B, wh="Gb", log=FALSE, array=i)
    bglog[[i]] = getArrayData(B, wh="Gb", log=TRUE, array=i)
}

## log2-scale
sapply(bglog, FUN="median", na.rm=TRUE)
# [1] 9.308339 9.308339 9.308339 9.306062 9.306062 9.306062 9.308339 9.306062 9.308339 9.306062 9.308339 9.308339

## original scale
sapply(bgorig, FUN="median")
# [1] 634 634 634 633 633 633 634 633 634 633 634 634


###########################################################
## Number of negative corrected intensities from full data set - sharpened
load("bld.sharpen.subtract.rda")
narrays <- length(arrayNames(bld.sharpen.subtract))

negs <- numeric(12)
for(i in 1:narrays)
    negs[i] = sum(bld.sharpen.subtract[[i]]$G<0)

## Percentage negative from full data set
round(negs/numBeads(bld.sharpen.subtract)*100,2)
# [1] 1.04 0.76 0.93 0.68 0.92 0.68 0.94 0.82 0.84 0.74 0.90 0.82 0.29 0.13 0.54
#[16] 0.15 0.41 0.22 0.29 0.15 0.38 0.14 0.41 0.13 0.61 0.19 0.46 0.13 0.52 0.18
#[31] 0.38 0.21 0.38 0.17 0.40 0.16 0.29 0.10 0.48 0.11 0.28 0.12 0.32 0.15 0.39
#[46] 0.15 0.32 0.14 0.50 0.17 0.46 0.26 0.64 0.22 0.57 0.16 0.52 0.34 0.49 0.22
#[61] 0.61 0.14 0.63 0.29 0.61 0.20 0.54 0.14 0.39 0.20 0.48 0.19 0.50 0.28 0.79
#[76] 0.32 0.75 0.30 0.87 0.21 0.63 0.32 0.81 0.30 0.22 0.18 0.41 0.22 0.57 0.14
#[91] 0.46 0.23 0.22 0.12 0.55 0.13

# quartiles
quantile(round(negs/numBeads(bld.sharpen.subtract)*100,2), c(0,0.25,0.5,0.75,1))
#   0%   25%   50%   75%  100%
#0.100 0.190 0.320 0.555 1.040


###########################################################
## Create bead summary data for each pre-processing option using
## the default Illumina method of removing outliers > 3 MADs 
## from the median before averaging the intensities

## sharpen images, no local background subtraction
load("bld.sharpen.nobgc.rda")
bsd.sharpen.nobgc <- createBeadSummaryData(bld.sharpen.nobgc, log=TRUE, 
                                           imagesPerArray=2, what="G", 
                                           method="illumina", n=3)
rm(bld.sharpen.nobgc)

## sharpen images, with local background subtraction
load("bld.sharpen.subtract.rda")
bsd.sharpen.subtract <- createBeadSummaryData(bld.sharpen.subtract, log=TRUE, 
                                              imagesPerArray=2, what="G", 
                                              method="illumina", n=3)
bsd.sharpen.subtract.raw <- createBeadSummaryData(bld.sharpen.subtract, log=FALSE, 
                                                  imagesPerArray=2, what="G", 
                                                  method="illumina", n=3)
rm(bld.sharpen.subtract)


## sharpen images, with local background subtraction and
## normal-exponential convolution model (avoids negative background
## corrected intensities)
load("bld.sharpen.normexp.rda")
bsd.sharpen.normexp <- createBeadSummaryData(bld.sharpen.normexp, log=TRUE, 
                                             imagesPerArray=2, what="G", 
                                             method="illumina", n=3)
rm(bld.sharpen.normexp)

## no sharpening, no local background subtraction
load("bld.nosharpen.nobgc.rda")
bsd.nosharpen.nobgc <- createBeadSummaryData(bld.nosharpen.nobgc, log=TRUE, imagesPerArray=2, what="G",  method="illumina", n=3)
rm(bld.nosharpen.nobgc)

## no sharpening, with local background subtraction
load("bld.nosharpen.subtract.rda")
bsd.nosharpen.subtract = createBeadSummaryData(bld.nosharpen.subtract, log=TRUE, imagesPerArray=2, what="G",  method="illumina", n=3)
rm(bld.nosharpen.subtract)

## Save summary data
save(bsd.sharpen.nobgc, bsd.sharpen.subtract, bsd.sharpen.normexp, bsd.nosharpen.nobgc, bsd.nosharpen.subtract, file="bsd.log.ill.rda")


###########################################################
## Figure 2 - plot results from simulation study (varying number of outliers 
## were introduced at the saturation level to assess how many outliers can be 
## tolerated before serious bias is introduced in to the ## expression measures)

targets <- read.table("./SpikeInData/spike_targets.txt", header=TRUE, sep=" ")
arraynms <- as.character(targets$ArrayNo)
narrays <- length(arraynms)

## introduce outliers in 2nd BeadChip, at varying %
## Use sharpened, subtracted data from text files
bld.sharpen.array2 <- readIllumina(path = "./SpikeInData", arrayNames=arraynms[13:24], 
                                   useImages=FALSE, textType=".csv")

## percentage of beads to be saturated
per.out <- c(0, 0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4)

nchips <- 12
arraynms <- arrayNames(bld.sharpen.array2)

## Saturated
ill.summ.array2 <- list(NoBeads=list(), exprs=list(), se.exprs=list()) 
ill.summ2.array2 <- med.summ.array2 <- mean.summ.array2 <- trim.summ.array2 <- win.summ.array2 <- ill.summ.array2
rmx.summ.array2 <- ill.summ.array2

set.seed(06092007)
for(l in 1:length(per.out)) {
    cat(l, "of 9\n")
    simdataarray2 <- copyBeadLevelList(bld.sharpen.array2)
    if(per.out[l]!=0) {
        for(i in 1:nchips) {
            cat(i, "of 12\n")
        ## no gross error model!
#            nsamp <- round(per.out[l]*numBeads(simdataarray2, array=i),0)
#            ind <- sample(1:numBeads(simdataarray2, array=i), nsamp)
        ## gross error model!
#            ind <- as.logical(rbinom(numBeads(simdataarray2, array=i), prob = per.out[l], size = 1))
#            ind <- (1:numBeads(simdataarray2, array=i))[ind]
#            simdataarray2@beadData[[arraynms[i]]]$G[ind] <- 2^16
        ## problem: there are bead types where more than 50% of the values are contaminated
        ## => there is no meaningful estimator which can handle this!
        ## Hence, we contiminate bead-type wise
            sel <- which(simdataarray2@beadData[[arraynms[i]]]$ProbeID != 0)
            pr <- simdataarray2@beadData[[arraynms[i]]]$ProbeID[sel]
            probes <- sort(unique(pr))
            indices <- NULL
            for(j in seq(along = probes)){
                ind <- pr == probes[j]
                ## < 50% of values contaminated
                repeat{
                    cont <- as.logical(rbinom(sum(ind), prob = per.out[l], size = 1))
                    if(sum(cont) < sum(ind)/2) break
                }
                indices <- c(indices, sel[ind][cont])
            }
            simdataarray2@beadData[[arraynms[i]]]$G[indices] <- 2^16
        }
    }
    tmp.ill <- createBeadSummaryData(simdataarray2, method="illumina", log=FALSE, n=3, imagesPerArray=2)
    tmp.ill2 <- createBeadSummaryData(simdataarray2, method="illumina", log=FALSE, n=2, imagesPerArray=2) 
    tmp.med <- createBeadSummaryData(simdataarray2, method="median", log=FALSE, imagesPerArray=2)
    tmp.mean <- createBeadSummaryData(simdataarray2, method="mean", log=FALSE, imagesPerArray=2)
    tmp.trim <- createBeadSummaryData(simdataarray2, method="trim", trim=0.1, log=FALSE, imagesPerArray=2)
    tmp.win <- createBeadSummaryData(simdataarray2, method="winsorize", trim=0.1, log=FALSE, imagesPerArray=2)
    tmp.rmx <- robloxbioc(simdataarray2, imagesPerArray=2, eps = max(per.out[l], 0.05))
    ill.summ.array2$NoBeads[[l]] = NoBeads(tmp.ill)
    ill.summ.array2$exprs[[l]] = exprs(tmp.ill)
    ill.summ.array2$se.exprs[[l]] = se.exprs(tmp.ill)
    ill.summ2.array2$NoBeads[[l]] = NoBeads(tmp.ill2)
    ill.summ2.array2$exprs[[l]] = exprs(tmp.ill2)
    ill.summ2.array2$se.exprs[[l]] = se.exprs(tmp.ill2)
    med.summ.array2$NoBeads[[l]] = NoBeads(tmp.med)
    med.summ.array2$exprs[[l]] = exprs(tmp.med)
    med.summ.array2$se.exprs[[l]] = se.exprs(tmp.med)
    mean.summ.array2$NoBeads[[l]] = NoBeads(tmp.mean)
    mean.summ.array2$exprs[[l]] = exprs(tmp.mean)
    mean.summ.array2$se.exprs[[l]] = se.exprs(tmp.mean)
    win.summ.array2$NoBeads[[l]] = NoBeads(tmp.win)
    win.summ.array2$exprs[[l]] = exprs(tmp.win)
    win.summ.array2$se.exprs[[l]] = se.exprs(tmp.win)
    trim.summ.array2$NoBeads[[l]] = NoBeads(tmp.trim)
    trim.summ.array2$exprs[[l]] = exprs(tmp.trim)
    trim.summ.array2$se.exprs[[l]] = se.exprs(tmp.trim)
    rmx.summ.array2$NoBeads[[l]] = NoBeads(tmp.rmx)
    rmx.summ.array2$exprs[[l]] = exprs(tmp.rmx)
    rmx.summ.array2$se.exprs[[l]] = se.exprs(tmp.rmx)
    rm(simdataarray2, tmp.ill, tmp.med, tmp.mean, tmp.trim, tmp.win, tmp.rmx)
}
save(ill.summ.array2, ill.summ2.array2, med.summ.array2, mean.summ.array2, 
     win.summ.array2, trim.summ.array2, rmx.summ.array2, file="sim.summary.bias.raw.rda")

load(file="sim.summary.bias.raw.rda")

## calculate bias - assume original, complete data set gave true values
bias.ill.array2 <- list(NoBeads=list(), exprs=list(), se.exprs=list())
bias.ill2.array2 <- bias.med.array2 <- bias.mean.array2 <- bias.trim.array2 <- bias.win.array2 <- bias.ill.array2
bias.rmx.array2 <- bias.ill.array2
for(l in 1:length(per.out)) {
    cat(l, "\n")
    bias.ill.array2$exprs[[l]] <- (ill.summ.array2$exprs[[l]]-mean.summ.array2$exprs[[1]])^2
    bias.ill.array2$se.exprs[[l]] <- ill.summ.array2$se.exprs[[l]]^2*ill.summ.array2$NoBeads[[l]]
    bias.ill.array2$NoBeads[[l]] <- mean.summ.array2$NoBeads[[1]]-ill.summ.array2$NoBeads[[l]]

    bias.ill2.array2$exprs[[l]] <- (ill.summ2.array2$exprs[[l]]-mean.summ.array2$exprs[[1]])^2
    bias.ill2.array2$se.exprs[[l]] <- ill.summ2.array2$se.exprs[[l]]^2*ill.summ2.array2$NoBeads[[l]]
    bias.ill2.array2$NoBeads[[l]] <- mean.summ.array2$NoBeads[[1]]-ill.summ2.array2$NoBeads[[l]]

    bias.med.array2$exprs[[l]] <- (med.summ.array2$exprs[[l]]-mean.summ.array2$exprs[[1]])^2
    bias.med.array2$se.exprs[[l]] <- med.summ.array2$se.exprs[[l]]^2*med.summ.array2$NoBeads[[l]]
    bias.med.array2$NoBeads[[l]] <- mean.summ.array2$NoBeads[[1]]-med.summ.array2$NoBeads[[l]]
  
    bias.mean.array2$exprs[[l]] <- (mean.summ.array2$exprs[[l]]-mean.summ.array2$exprs[[1]])^2
    bias.mean.array2$se.exprs[[l]] <- mean.summ.array2$se.exprs[[l]]^2*mean.summ.array2$NoBeads[[l]]
    bias.mean.array2$NoBeads[[l]] <- mean.summ.array2$NoBeads[[1]]-mean.summ.array2$NoBeads[[l]]

    bias.trim.array2$exprs[[l]] <- (trim.summ.array2$exprs[[l]]-mean.summ.array2$exprs[[1]])^2
    bias.trim.array2$se.exprs[[l]] <- trim.summ.array2$se.exprs[[l]]^2*trim.summ.array2$NoBeads[[l]]
    bias.trim.array2$NoBeads[[l]] <- mean.summ.array2$NoBeads[[1]]-trim.summ.array2$NoBeads[[l]]

    bias.win.array2$exprs[[l]] <- (win.summ.array2$exprs[[l]]-mean.summ.array2$exprs[[1]])^2
    bias.win.array2$se.exprs[[l]] <- win.summ.array2$se.exprs[[l]]^2*win.summ.array2$NoBeads[[l]]
    bias.win.array2$NoBeads[[l]] <- mean.summ.array2$NoBeads[[1]]-win.summ.array2$NoBeads[[l]]

    bias.rmx.array2$exprs[[l]] <- (rmx.summ.array2$exprs[[l]]-mean.summ.array2$exprs[[1]])^2
    bias.rmx.array2$se.exprs[[l]] <- rmx.summ.array2$se.exprs[[l]]^2*rmx.summ.array2$NoBeads[[l]]
    bias.rmx.array2$NoBeads[[l]] <- mean.summ.array2$NoBeads[[1]]-rmx.summ.array2$NoBeads[[l]]
}

## calculate average bias from `replicate arrays
avebias.ill.array2 <- avebias.ill2.array2 <- avebias.mean.array2 <- avebias.med.array2 <- avebias.trim.array2 <- avebias.win.array2 <- NULL
avebias.rmx.array2 <- NULL
for(l in 1:length(per.out)) {
    avebias.ill.array2[l] <- mean(as.vector(bias.ill.array2$exprs[[l]]), na.rm=TRUE)
    avebias.ill2.array2[l] <- mean(as.vector(bias.ill2.array2$exprs[[l]]), na.rm=TRUE)
    avebias.med.array2[l] <- mean(as.vector(bias.med.array2$exprs[[l]]), na.rm=TRUE)
    avebias.mean.array2[l] <- mean(as.vector(bias.mean.array2$exprs[[l]]), na.rm=TRUE)
    avebias.win.array2[l] <- mean(as.vector(bias.win.array2$exprs[[l]]), na.rm=TRUE)
    avebias.trim.array2[l] <- mean(as.vector(bias.trim.array2$exprs[[l]]), na.rm=TRUE)
    avebias.rmx.array2[l] <- mean(as.vector(bias.rmx.array2$exprs[[l]]), na.rm=TRUE)
}

## calculate variance
avevar.ill.array2 <- avevar.ill2.array2 <- avevar.mean.array2 <- avevar.med.array2 <- avevar.trim.array2 <- avevar.win.array2 <- NULL
avevar.rmx.array2 <- NULL
for(l in 1:length(per.out)) {
    avevar.ill.array2[l] <- mean(as.vector(bias.ill.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.ill2.array2[l] <- mean(as.vector(bias.ill2.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.med.array2[l] <- mean(as.vector(bias.med.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.mean.array2[l] <- mean(as.vector(bias.mean.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.win.array2[l] <- mean(as.vector(bias.win.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.trim.array2[l] <- mean(as.vector(bias.trim.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.rmx.array2[l] <- mean(as.vector(bias.rmx.array2$se.exprs[[l]]), na.rm=TRUE)
}

## calculate MSE
avemse.ill.array2 <- avemse.ill2.array2 <- avemse.mean.array2 <- avemse.med.array2 <- avemse.trim.array2 <- avemse.win.array2 <- NULL
avemse.rmx.array2 <- NULL
for(l in 1:length(per.out)) {
    avemse.ill.array2[l] <- mean(as.vector(bias.ill.array2$se.exprs[[l]])+as.vector(bias.ill.array2$exprs[[l]]), na.rm=TRUE)
    avemse.ill2.array2[l] <- mean(as.vector(bias.ill2.array2$se.exprs[[l]])+as.vector(bias.ill2.array2$exprs[[l]]), na.rm=TRUE)
    avemse.med.array2[l] <- mean(as.vector(bias.med.array2$se.exprs[[l]])+as.vector(bias.med.array2$exprs[[l]]), na.rm=TRUE)
    avemse.mean.array2[l] <- mean(as.vector(bias.mean.array2$se.exprs[[l]])+as.vector(bias.mean.array2$exprs[[l]]), na.rm=TRUE)
    avemse.win.array2[l] <- mean(as.vector(bias.win.array2$se.exprs[[l]])+as.vector(bias.win.array2$exprs[[l]]), na.rm=TRUE)
    avemse.trim.array2[l] <- mean(as.vector(bias.trim.array2$se.exprs[[l]])+as.vector(bias.trim.array2$exprs[[l]]), na.rm=TRUE)
    avemse.rmx.array2[l] <- mean(as.vector(bias.rmx.array2$se.exprs[[l]])+as.vector(bias.rmx.array2$exprs[[l]]), na.rm=TRUE)
}

lwd <- 2
sel <- 1:9
x <- per.out[sel]*100
myCol <- brewer.pal(5, "Set1")

pdf("Figure2.pdf", width=8, height=5)
par(mfrow=c(1,3))
par(mar=c(3,4,1.5,0.15), oma=c(1,0,0,0))
plot(x, log2(avebias.ill.array2[sel]), lwd=lwd, type="l", xlab="", ylab=expression(log[2]("average bias^2")), 
     main="Bias^2", ylim=c(0,30), col = myCol[1], 
     panel.first = abline(v = c(0, 10, 20, 30, 40), h = seq(0, 30, by = 5), lty = 2, col = "grey"))
points(x, log2(avebias.med.array2[sel]), type="l", col=myCol[2], lwd=lwd)
#points(x, avebias.trim.array2[sel], type="l", col=myCol[3], lwd=lwd)
points(x, log2(avebias.mean.array2[sel]), type="l", col=myCol[4], lwd=lwd)
points(x, log2(avebias.rmx.array2[sel]), type="l", col=myCol[5], lwd=lwd)
legend("bottomright",legend=c("Illumina", "median", "mean", "rmx"), col=myCol[c(1,2,4,5)], lwd=2)
#plot(x, numout.ill.array2[sel], lwd=lwd, type="l", xlab="", ylab="% observations removed by summary method", main="(b)", ylim=c(0,40))
#points(x, numout.med.array2[sel], type="l", col=2, lwd=lwd)
#points(x, numout.trim.array2[sel], type="l", col=3, lwd=lwd)
#points(x, numout.mean.array2[sel], type="l", col=4, lwd=lwd)
#abline(0,1,col="gray", lty=2)
plot(x, log2(avevar.ill.array2)[sel], lwd=lwd, type="l", xlab="", ylab=expression(log[2]("average variance")), 
     main="Variance", ylim=c(0,30),col = myCol[1], 
     panel.first = abline(v = c(0, 10, 20, 30, 40), h = seq(0, 30, by = 5), lty = 2, col = "grey"))
points(x, log2(avevar.med.array2)[sel], type="l", col=myCol[2], lwd=lwd)
#points(x, log2(avevar.trim.array2)[sel], type="l", col=myCol[3], lwd=lwd)
points(x, log2(avevar.mean.array2)[sel], type="l", col=myCol[4], lwd=lwd)
points(x, log2(avevar.rmx.array2)[sel], type="l", col=myCol[5], lwd=lwd)
plot(x, log2(avemse.ill.array2)[sel], lwd=lwd, type="l", xlab="", ylab=expression(log[2]("average MSE")), 
     main="MSE", ylim=c(0,30),col = myCol[1], 
     panel.first = abline(v = c(0, 10, 20, 30, 40), h = seq(0, 30, by = 5), lty = 2, col = "grey"))
points(x, log2(avemse.med.array2)[sel], type="l", col=myCol[2], lwd=lwd)
#points(x, log2(avemse.trim.array2)[sel], type="l", col=myCol[3], lwd=lwd)
points(x, log2(avemse.mean.array2)[sel], type="l", col=myCol[4], lwd=lwd)
points(x, log2(avemse.rmx.array2)[sel], type="l", col=myCol[5], lwd=lwd)
mtext("% outliers simulated", side=1, outer=TRUE)
dev.off()


## log-scale
## Saturated
ill.log.summ.array2 <- list(NoBeads=list(), exprs=list(), se.exprs=list()) 
ill.log.summ2.array2 <- med.log.summ.array2 <- mean.log.summ.array2 <- trim.log.summ.array2 <- win.log.summ.array2 <- ill.log.summ.array2
ill.log.rmx.array2 <- ill.log.summ.array2

set.seed(06092007)
for(l in 1:length(per.out)) {
    cat(l, "\n")
    simdataarray2 <- copyBeadLevelList(bld.sharpen.array2)
    if(per.out[l]!=0) {
        for(i in 1:nchips) {
            cat(i, "of 12\n")
        ## no gross error model!
#            nsamp <- round(per.out[l]*numBeads(simdataarray2, array=i),0)
#            ind <- sample(1:numBeads(simdataarray2, array=i), nsamp)
        ## gross error model!
#            ind <- as.logical(rbinom(numBeads(simdataarray2, array=i), prob = per.out[l], size = 1))
#            ind <- (1:numBeads(simdataarray2, array=i))[ind]
#            simdataarray2@beadData[[arraynms[i]]]$G[ind] <- 2^16
        ## problem: there are bead types where more than 50% of the values are contaminated
        ## => there is no meaningful estimator which can handle this!
            sel <- which(simdataarray2@beadData[[arraynms[i]]]$ProbeID != 0)
            pr <- simdataarray2@beadData[[arraynms[i]]]$ProbeID[sel]
            probes <- sort(unique(pr))
            indices <- NULL
            for(j in seq(along = probes)){
                ind <- pr == probes[j]
                ## < 50% of values contaminated
                repeat{
                    cont <- as.logical(rbinom(sum(ind), prob = per.out[l], size = 1))
                    if(sum(cont) < sum(ind)/2) break
                }
                indices <- c(indices, sel[ind][cont])
            }
            simdataarray2@beadData[[arraynms[i]]]$G[indices] <- 2^16
        }
    }
    tmp.ill.log <- createBeadSummaryData(simdataarray2, method="illumina", log=TRUE, n=3, imagesPerArray=2)
    tmp.ill2.log <- createBeadSummaryData(simdataarray2, method="illumina", log=TRUE, n=2, imagesPerArray=2) 
    tmp.med.log <- createBeadSummaryData(simdataarray2, method="median", log=TRUE, imagesPerArray=2)
    tmp.mean.log <- createBeadSummaryData(simdataarray2, method="mean", log=TRUE, imagesPerArray=2)
    tmp.trim.log <- createBeadSummaryData(simdataarray2, method="trim", trim=0.1, log=TRUE, imagesPerArray=2)
    tmp.win.log <- createBeadSummaryData(simdataarray2, method="winsorize", trim=0.1, log=TRUE, imagesPerArray=2)
    tmp.rmx.log <- robloxbioc(simdataarray2, log = TRUE, imagesPerArray=2, eps = max(per.out[l], 0.05))
    ill.log.summ.array2$NoBeads[[l]] <- NoBeads(tmp.ill.log)
    ill.log.summ.array2$exprs[[l]] <- exprs(tmp.ill.log)
    ill.log.summ.array2$se.exprs[[l]] <- se.exprs(tmp.ill.log)
    ill.log.summ2.array2$NoBeads[[l]] <- NoBeads(tmp.ill2.log)
    ill.log.summ2.array2$exprs[[l]] <- exprs(tmp.ill2.log)
    ill.log.summ2.array2$se.exprs[[l]] <- se.exprs(tmp.ill2.log)
    med.log.summ.array2$NoBeads[[l]] <- NoBeads(tmp.med.log)
    med.log.summ.array2$exprs[[l]] <- exprs(tmp.med.log)
    med.log.summ.array2$se.exprs[[l]] <- se.exprs(tmp.med.log)
    mean.log.summ.array2$NoBeads[[l]] <- NoBeads(tmp.mean.log)
    mean.log.summ.array2$exprs[[l]] <- exprs(tmp.mean.log)
    mean.log.summ.array2$se.exprs[[l]] <- se.exprs(tmp.mean.log)
    win.log.summ.array2$NoBeads[[l]] <- NoBeads(tmp.win.log)
    win.log.summ.array2$exprs[[l]] <- exprs(tmp.win.log)
    win.log.summ.array2$se.exprs[[l]] <- se.exprs(tmp.win.log)
    trim.log.summ.array2$NoBeads[[l]] <- NoBeads(tmp.trim.log)
    trim.log.summ.array2$exprs[[l]] <- exprs(tmp.trim.log)
    trim.log.summ.array2$se.exprs[[l]] <- se.exprs(tmp.trim.log)
    rmx.log.summ.array2$NoBeads[[l]] <- NoBeads(tmp.rmx.log)
    rmx.log.summ.array2$exprs[[l]] <- exprs(tmp.rmx.log)
    rmx.log.summ.array2$se.exprs[[l]] <- se.exprs(tmp.rmx.log)
    rm(simdataarray2, tmp.ill.log, tmp.med.log, tmp.mean.log, tmp.trim.log, tmp.win.log, tmp.rmx.log)
}
save(ill.log.summ.array2, ill.log.summ2.array2, med.log.summ.array2, mean.log.summ.array2, 
     win.log.summ.array2, trim.log.summ.array2, rmx.log.summ.array2, file="sim.summary.bias.log.rda")

bias.ill.log.array2 <- list(NoBeads=list(), exprs=list(), se.exprs=list())
bias.ill2.log.array2 <- bias.med.log.array2 <- bias.mean.log.array2 <- bias.trim.log.array2 <- bias.win.log.array2 <- bias.ill.log.array2
bias.rmx.log.array2 <- bias.ill.log.array2
for(l in 1:length(per.out)) {
    cat(l, "\n")
    bias.ill.log.array2$exprs[[l]] <- ill.log.summ.array2$exprs[[l]]-mean.log.summ.array2$exprs[[1]]
    bias.ill.log.array2$se.exprs[[l]] <- ill.log.summ.array2$se.exprs[[l]]^2*ill.log.summ.array2$NoBeads[[l]]
    bias.ill.log.array2$NoBeads[[l]] <- mean.log.summ.array2$NoBeads[[1]]-ill.log.summ.array2$NoBeads[[l]]

    bias.ill2.log.array2$exprs[[l]] <- ill.log.summ2.array2$exprs[[l]]-mean.log.summ.array2$exprs[[1]]
    bias.ill2.log.array2$se.exprs[[l]] <- ill.log.summ2.array2$se.exprs[[l]]^2*ill.log.summ2.array2$NoBeads[[l]]
    bias.ill2.log.array2$NoBeads[[l]] <- mean.log.summ.array2$NoBeads[[1]]-ill.log.summ2.array2$NoBeads[[l]]

    bias.med.log.array2$exprs[[l]] <- med.log.summ.array2$exprs[[l]]-mean.log.summ.array2$exprs[[1]]
    bias.med.log.array2$NoBeads[[l]] <- mean.log.summ.array2$NoBeads[[1]]-med.log.summ.array2$NoBeads[[l]]
    
    bias.mean.log.array2$exprs[[l]] <- mean.log.summ.array2$exprs[[l]]-mean.log.summ.array2$exprs[[1]]
    bias.mean.log.array2$se.exprs[[l]] <- mean.log.summ.array2$se.exprs[[l]]^2*mean.log.summ.array2$NoBeads[[l]]
    bias.mean.log.array2$NoBeads[[l]] <- mean.log.summ.array2$NoBeads[[1]]-mean.log.summ.array2$NoBeads[[l]]

    bias.trim.log.array2$exprs[[l]] <- trim.log.summ.array2$exprs[[l]]-mean.log.summ.array2$exprs[[1]]
    bias.trim.log.array2$se.exprs[[l]] <- trim.log.summ.array2$se.exprs[[l]]^2*trim.log.summ.array2$NoBeads[[l]]
    bias.trim.log.array2$NoBeads[[l]] <- mean.log.summ.array2$NoBeads[[1]]-trim.log.summ.array2$NoBeads[[l]]

    bias.win.log.array2$exprs[[l]] <- win.log.summ.array2$exprs[[l]]-mean.log.summ.array2$exprs[[1]]
    bias.win.log.array2$se.exprs[[l]] <- win.log.summ.array2$se.exprs[[l]]^2*win.log.summ.array2$NoBeads[[l]]
    bias.win.log.array2$NoBeads[[l]] <- mean.log.summ.array2$NoBeads[[1]]-win.log.summ.array2$NoBeads[[l]]

    bias.rmx.log.array2$exprs[[l]] <- rmx.log.summ.array2$exprs[[l]]-mean.log.summ.array2$exprs[[1]]
    bias.rmx.log.array2$se.exprs[[l]] <- rmx.log.summ.array2$se.exprs[[l]]^2*rmx.log.summ.array2$NoBeads[[l]]
    bias.rmx.log.array2$NoBeads[[l]] <- mean.log.summ.array2$NoBeads[[1]]-rmx.log.summ.array2$NoBeads[[l]]
}

avebias.ill.log.array2 <- avebias.ill2.log.array2 <- avebias.mean.log.array2 <- avebias.med.log.array2 <- avebias.trim.log.array2 <- avebias.win.log.array2 <- NULL
avebias.rmx.log.array2 <- NULL
for(l in 1:length(per.out)) {
    avebias.ill.log.array2[l] <- mean(as.vector(bias.ill.log.array2$exprs[[l]]), na.rm=TRUE)
    avebias.ill2.log.array2[l] <- mean(as.vector(bias.ill2.log.array2$exprs[[l]]), na.rm=TRUE)
    avebias.med.log.array2[l] <- mean(as.vector(bias.med.log.array2$exprs[[l]]), na.rm=TRUE)
    avebias.mean.log.array2[l] <- mean(as.vector(bias.mean.log.array2$exprs[[l]]), na.rm=TRUE)
    avebias.win.log.array2[l] <- mean(as.vector(bias.win.log.array2$exprs[[l]]), na.rm=TRUE)
    avebias.trim.log.array2[l] <- mean(as.vector(bias.trim.log.array2$exprs[[l]]), na.rm=TRUE)
    avebias.rmx.log.array2[l] <- mean(as.vector(bias.rmx.log.array2$exprs[[l]]), na.rm=TRUE)
}

numout.ill.log.array2 <- numout.ill2.log.array2 <- numout.mean.log.array2 <- numout.med.log.array2 <- numout.trim.log.array2 <- numout.win.log.array2 <- NULL
numout.rmx.log.array2 <- NULL
for(l in 1:length(per.out)) {
    numout.ill.log.array2[l] <- mean(as.vector(bias.ill.log.array2$NoBeads[[l]]/mean.log.summ.array2$NoBeads[[1]]), na.rm=TRUE)*100
    numout.ill2.log.array2[l] <- mean(as.vector(bias.ill2.log.array2$NoBeads[[l]]/mean.log.summ.array2$NoBeads[[1]]), na.rm=TRUE)*100
    numout.med.log.array2[l] <- mean(as.vector(bias.med.log.array2$NoBeads[[l]]/mean.log.summ.array2$NoBeads[[1]]), na.rm=TRUE)*100
    numout.mean.log.array2[l] <- mean(as.vector(bias.mean.log.array2$NoBeads[[l]]/mean.log.summ.array2$NoBeads[[1]]), na.rm=TRUE)*100
    numout.win.log.array2[l] <- mean(as.vector(bias.win.log.array2$NoBeads[[l]]/mean.log.summ.array2$NoBeads[[1]]), na.rm=TRUE)*100
    numout.trim.log.array2[l] <- mean(as.vector(bias.trim.log.array2$NoBeads[[l]]/mean.log.summ.array2$NoBeads[[1]]), na.rm=TRUE)*100
    numout.rmx.log.array2[l] <- mean(as.vector(bias.rmx.log.array2$NoBeads[[l]]/mean.log.summ.array2$NoBeads[[1]]), na.rm=TRUE)*100
}

avevar.ill.log.array2 <- avevar.ill2.log.array2 <- avevar.mean.log.array2 <- avevar.med.log.array2 <- avevar.trim.log.array2 <- avevar.win.log.array2 <- NULL
avevar.rmx.log.array2 <- NULL
for(l in 1:length(per.out)) {
    avevar.ill.log.array2[l] <- mean(as.vector(bias.ill.log.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.ill2.log.array2[l] <- mean(as.vector(bias.ill2.log.array2$se.exprs[[l]]), na.rm=TRUE)
#    avevar.med.log.array2[l] <- mean(as.vector(bias.med.log.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.mean.log.array2[l] <- mean(as.vector(bias.mean.log.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.win.log.array2[l] <- mean(as.vector(bias.win.log.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.trim.log.array2[l] <- mean(as.vector(bias.trim.log.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.rmx.log.array2[l] <- mean(as.vector(bias.rmx.log.array2$se.exprs[[l]]), na.rm=TRUE)
}

pdf("Figure2b.pdf", width=12, height=8)
par(mfrow=c(1,2))
par(mar=c(3,4,1.5,0.15), oma=c(1,0,0,0))
plot(x, avebias.ill.log.array2[sel], lwd=lwd, type="l", xlab="", ylab="average bias", main="A", ylim=c(0,4))
points(x, avebias.med.log.array2[sel], type="l", col=2, lwd=lwd)
points(x, avebias.trim.log.array2[sel], type="l", col=3, lwd=lwd)
points(x, avebias.mean.log.array2[sel], type="l", col=4, lwd=lwd)
points(x, avebias.rmx.log.array2[sel], type="l", col=5, lwd=lwd)
legend("topleft",legend=c("Illumina", "median", "trimmed mean", "mean", "rmx"), col=1:5, lwd=2)
#plot(x, numout.ill.log.array2[sel], lwd=lwd, type="l", xlab="", ylab="% observations removed by summary method", main="(b)", ylim=c(0,40))
#points(x, numout.med.log.array2[sel], type="l", col=2, lwd=lwd)
#points(x, numout.trim.log.array2[sel], type="l", col=3, lwd=lwd)
#points(x, numout.mean.log.array2[sel], type="l", col=4, lwd=lwd)
#abline(0,1,col="gray", lty=2)
plot(x, log2(avevar.ill.log.array2)[sel], lwd=lwd, type="l", xlab="", ylab=expression(log[2]("average variance")), main="B", ylim=c(-5,5))
points(x, log2(avevar.trim.log.array2)[sel], type="l", col=3, lwd=lwd)
points(x, log2(avevar.mean.log.array2)[sel], type="l", col=4, lwd=lwd)
points(x, log2(avevar.rmx.log.array2)[sel], type="l", col=5, lwd=lwd)
mtext("% outliers simulated", side=1, outer=TRUE)
dev.off()

###########################################################
## Figure 3 - plot intensities and variances for spike probes only
pal <- brewer.pal(6, "Set1")

## Read in spike intensities (file supplied by Illumina) to retrieve IDs
spikecsv <- read.csv("spikeins_profile.csv")
spikeIDs <- as.character(spikecsv$ProbeID)

nspikes <- length(spikeIDs)
nprobes <- nrow(exprs(bsd.sharpen.subtract))

spikeInd <- NULL
for(i in 1:33) {
  spikeInd <- c(spikeInd, seq(1:nprobes)[rownames(exprs(bsd.sharpen.subtract))==spikeIDs[i]])
 }
  

concs <- rep(c(1000,300,100,30,10,3),4)
concs <- c(concs, rep(c(1,0.3,0.1,0.03,0.01,0),4))
concs <- rep(concs,3)

o <- order(concs, decreasing=TRUE)

allExprs <- cbind(exprs(bsd.sharpen.nobgc), exprs(bsd.sharpen.subtract),exprs(bsd.sharpen.normexp))

#se.exprs(bsd.nosharp.nobgc) = assayData(bsd.nosharp.nobgc)$BeadStDev
#se.exprs(bsd.sharpen.subtract) = assayData(bsd.sharpen.subtract)$BeadStDev

allVar <- cbind(NoBeads(bsd.sharpen.nobgc)*(se.exprs(bsd.sharpen.nobgc))^2,       
                NoBeads(bsd.sharpen.subtract)*(se.exprs(bsd.sharpen.subtract))^2,
                NoBeads(bsd.sharpen.normexp)*(se.exprs(bsd.sharpen.normexp))^2)

allVar <- log2(allVar)


pdf("Figure3.pdf", width=12, height=16)

par(mfrow=c(2,1), mar=c(2,4,1.5,0.15))
boxplot(as.data.frame(allExprs[spikeInd,o]), col=c(rep(pal[1],4), rep(pal[2],4),rep(pal[3],4)), 
        outline=FALSE, xlab="", ylab="Non-normalised log2 intensity", names=NULL, xaxt="n", 
        ylim=c(4,16), main="A",font.main=2,medlwd=1)
abline(v=c(13,25,37,49,61,73,85,97,109,121,133)-0.5, lty=2)
legend("topright", c("No background","Subtract","Normexp"), fill=c(pal[1], pal[2], pal[3]), bg="white")
axis(side=1,at=c(13,25,37,49,61,73,85,97,109,121,133,145)-6, 
     labels=c("1000pM", "300pM", "100pM", "30pM", "10pM", "3pM", "1pM", "0.3pM", "0.1pM", "0.03pM", "0.01pM", "0pM"))
boxplot(as.data.frame(allVar[spikeInd,o]), col=c(rep(pal[1],4), rep(pal[2],4),rep(pal[3],4)), 
        outline=FALSE, xlab="", ylab="Non-normalised log2 variance",names=NULL, xaxt="n", 
        main="B", font.main=2,medlwd=1)
abline(v=c(13,25,37,49,61,73,85,97,109,121,133)-0.5, lty=2)
axis(side=1,at=c(13,25,37,49,61,73,85,97,109,121,133,145)-6, 
     labels=c("1000pM", "300pM", "100pM", "30pM", "10pM", "3pM", "1pM", "0.3pM", "0.1pM", "0.03pM", "0.01pM", "0pM"))
dev.off()


###########################################################
## Figure 4 - plot intensities of the spike probes at each concentration
## Average over replicate arrays, separate line for each probe
## Use sharpened, background corrected intensities
targets <- read.table("spike_targets.txt",header=TRUE)
## set up colours
cls <- rainbow(33)

arraynms <- as.character(targets$ArrayNo)
narrays <- length(arraynms)

## set up design matrix for linear model
design <- cbind(as.numeric(targets$SpikeConc==1000),
                as.numeric(targets$SpikeConc==300),
                as.numeric(targets$SpikeConc==100),
                as.numeric(targets$SpikeConc==30),
                as.numeric(targets$SpikeConc==10),
                as.numeric(targets$SpikeConc==3),
                as.numeric(targets$SpikeConc==1),
                as.numeric(targets$SpikeConc==0.3),
                as.numeric(targets$SpikeConc==0.1),
                as.numeric(targets$SpikeConc==0.03),
                as.numeric(targets$SpikeConc==0.01),
                as.numeric(targets$SpikeConc==0))[seq(1,narrays, by=2),]

colnames(design) <- paste("pm", unique(targets$SpikeConc), sep="")

## set up contrasts matrix for the most similar concentration differences
conts <- makeContrasts(pm1000-pm300, pm300-pm100, pm100-pm30, pm30-pm10, 
                       pm10-pm3, pm3-pm1, pm1-pm0.3, pm0.3-pm0.1, pm0.3-pm0.1, 
                       pm0.1-pm0.03,pm0.03-pm0.01, levels=design)

spikecsv <- read.csv("spikeins_profile.csv")
spikeIDs <- as.character(spikecsv$ProbeID)
targetnames <- as.character(spikecsv$TargetID)

nspikes <- length(spikeIDs)
nprobes <- nrow(exprs(bsd.sharpen.subtract))


spikeInd <- NULL
for(i in 1:nspikes) {
    spikeInd <- c(spikeInd, seq(1:nprobes)[rownames(exprs(bsd.sharpen.subtract))==spikeIDs[i]])
}

 
concs <- rep(c(1000,300,100,30,10,3),4)
concs <- c(concs, rep(c(1,0.3,0.1,0.03,0.01,0),4))
concs <- rep(c(1000,300,100,30,10,3),4)
concs <- c(concs, rep(c(1,0.3,0.1,0.03,0.01,0),4))
 
o <- order(concs, decreasing=TRUE)

l1  <- lmFit(exprs(bsd.sharpen.subtract), design=design)

pdf("Figure4.pdf", width=12, height=12)
par(mar=c(4,4,1.5,0.15))
plot(l1$coeff[spikeInd[1],], ylim=range(5.5,16), type="n", ylab="log2 intensity", 
     xlab="Concentration",xaxt="n", cex=0.6, main="Non-normalised intensities of the spike genes")
for(i in 1:33){
    text(labels=targetnames[i],x=1:12,y=as.numeric(l1$coeff[spikeInd[i],]),col=cls[i],pch=16,cex=0.8)
    plot.smooth.line(1:12, as.numeric(l1$coeff[spikeInd[i],]), col=cls[i])
}
axis(side=1, at=1:12, labels=c("1000pM", "300pM", "100pM", "30pM", "10pM", "3pM", "1pM", 
                                "0.3pM", "0.1pM", "0.03pM", "0.01pM", "0pM"))
dev.off()


###########################################################
## Figure 5 - MA plots of data
## without background adjustment, with background subtracted data
## and background corrected and normalised data
spikeInd <-  NULL
for(i in 1:nspikes) {
    spikeInd <- c(spikeInd, seq(1:nprobes)[rownames(exprs(bsd.sharpen.subtract))==spikeIDs[i]])
}

controlInfo <- read.table("SpikesAveragedControls.txt", sep="\t", header=T)

controlInfo <- controlInfo[,grep("Signal", colnames(controlInfo))]

# 7th row are theaverages of the negative control beads
negControls <- controlInfo[7,]
negControls <- as.numeric(negControls)
negControlsMatrix <- matrix(rep(negControls, each=nprobes), nprobes, narrays/2)

B <- exprs(bsd.sharpen.subtract.raw) - negControlsMatrix
bsd.bgnorm <- assayDataElementReplace(bsd.sharpen.subtract.raw, "exprs",log2(B))

MAplot <- function(bsd, a1, a2,...){
    M <- exprs(bsd)[,a1] - exprs(bsd)[,a2]
    A <- 0.5*(exprs(bsd)[,a1] + exprs(bsd)[,a2])

    smoothScatter(A, M, pch=16, cex=0.5,...)
    points(A[spikeInd], M[spikeInd], col="red", pch=16)

    abline(h=c(log2(3),0), lty=2)
}

xlim <- c(-1,16)
ylim <- c(-3,3)


pdf("Figure5.pdf", width=12, height=8)
par(mfrow=c(1,3))
par(oma=c(1.1,1.1,0,0))
par(mar=c(2.2,2.2,1.5,0.15))
  MAplot(bsd.sharpen.nobgc, 24, 25,xlim=xlim, ylim=ylim, ylab="", xlab="")
  mtext("A", side=3, font=2)

  MAplot(bsd.sharpen.subtract,24, 25, xlim=xlim, ylim=ylim, ylab="", xlab="")
  mtext("B", side=3, font=2)
 
  MAplot(bsd.bgnorm, 24, 25, xlim=xlim, ylim=ylim, ylab="", xlab="", nbin=512)
  mtext("C", side=3, font=2)

  mtext("log-ratios", side=2, outer=TRUE)
  mtext("average expression", side=1, outer=TRUE)
dev.off()


## percentage negative values after background normalisation
apply(B<0, 2, FUN="sum", na.rm=TRUE)/dim(B)[1]*100

round(quantile(apply(B<0, 2, FUN="sum", na.rm=TRUE)/dim(B)[1]*100, c(0,0.25,0.5,0.75,1)),2)


###########################################################
## Fitting the linear model to separate background methods

## Make design matrix with 12 columns, one for each concentration
targets <- read.table("spike_targets.txt",header=TRUE)

arraynms <- as.character(targets$ArrayNo)
narrays <- length(arraynms)

design <- cbind(as.numeric(targets$SpikeConc==1000),
                as.numeric(targets$SpikeConc==300),
                as.numeric(targets$SpikeConc==100),
                as.numeric(targets$SpikeConc==30),
                as.numeric(targets$SpikeConc==10),
                as.numeric(targets$SpikeConc==3),
                as.numeric(targets$SpikeConc==1),
                as.numeric(targets$SpikeConc==0.3),
                as.numeric(targets$SpikeConc==0.1),
                as.numeric(targets$SpikeConc==0.03),
                as.numeric(targets$SpikeConc==0.01),
                as.numeric(targets$SpikeConc==0))[seq(1,narrays, by=2),]

colnames(design) <- paste("pm", unique(targets$SpikeConc), sep="")

## Make pairwise contrasts between all concentrations
contr.1000 <- makeContrasts(pm1000-pm300, pm1000-pm100, pm1000-pm30, pm1000-pm10, pm1000-pm3, 
                            pm1000-pm1, pm1000-pm0.3, pm1000-pm0.1, pm1000-pm0.03, 
                            pm1000-pm0.01, pm1000-pm0, levels=design)
contr.300 <- makeContrasts(pm300-pm100, pm300-pm30, pm300-pm10, pm300-pm3, pm300-pm1, 
                           pm300-pm0.3, pm300-pm0.1, pm300-pm0.03, pm300-pm0.01, pm300-pm0, 
                           levels=design)
contr.100 <- makeContrasts(pm100-pm30, pm100-pm10, pm100-pm3, pm100-pm1, pm100-pm0.3, 
                           pm100-pm0.1, pm100-pm0.03, pm100-pm0.01, pm100-pm0, levels=design)
contr.30 <- makeContrasts(pm30-pm10, pm30-pm3, pm30-pm1, pm30-pm0.3, pm30-pm0.1, pm30-pm0.03, 
                          pm30-pm0.01, pm30-pm0, levels=design)
contr.10 <- makeContrasts(pm10-pm3, pm10-pm1, pm10-pm0.3, pm10-pm0.1, pm10-pm0.03,levels=design)
contr.3 <- makeContrasts(pm3-pm1, pm3-pm0.3, pm3-pm0.1, pm3-pm0.03, levels=design)
contr.1 <- makeContrasts(pm1-pm0.3, pm1-pm0.1, pm1-pm0.03, levels=design)
contr.0.3 <- makeContrasts(pm0.3-pm0.1, pm0.3-pm0.03, levels=design)
contr.0.1 <- makeContrasts(pm0.1-pm0.03, levels=design)

contr.all <- cbind(contr.1000,contr.300,contr.100,contr.30,contr.10, contr.3, contr.1, contr.0.3, contr.0.1)


## Normalise and fit the model to sharpened subtracted data
fit.sharpen.subtract <- lmFit(normalizeQuantiles(exprs(bsd.sharpen.subtract)), design=design)
contr.sharpen.subtract <- contrasts.fit(fit.sharpen.subtract, contr.all)

## Empirical bayes smoothing of variances
contr.sharpen.subtract <- eBayes(contr.sharpen.subtract)

## Calculate a weight for each bead type using the standard error
## and number of replicate observations
w.probes <- se.exprs(bsd.sharpen.subtract)^2*bsd.sharpen.subtract@assayData$NoBeads
w.probes <- 1/w.probes

##Re-fit the model using weights
fit.pwts.subtract <- lmFit(normalizeQuantiles(exprs(bsd.sharpen.subtract)),
                           design=design, weights=w.probes)
contr.pwts <- contrasts.fit(fit.pwts.subtract, contr.all)
contr.pwts.subtract <- eBayes(contr.pwts)

fit.sharpen.nobgc <- lmFit(normalizeQuantiles(exprs(bsd.sharpen.nobgc)), design=design)
contr.sharpen.nobgc <- contrasts.fit(fit.sharpen.nobgc, contr.all)
contr.sharpen.nobgc <- eBayes(contr.sharpen.nobgc)
w.probes <- se.exprs(bsd.sharpen.nobgc)^2*bsd.sharpen.nobgc@assayData$NoBeads

w.probes <- 1/w.probes
fit.pwts.nobgc  <- lmFit(normalizeQuantiles(exprs(bsd.sharpen.nobgc)), design=design, weights=w.probes)

contr.pwts <- contrasts.fit(fit.pwts.nobgc, contr.all)
contr.pwts.nobgc <- eBayes(contr.pwts)

fit.sharpen.normexp <- lmFit(normalizeQuantiles(exprs(bsd.sharpen.normexp)), design=design)
contr.sharpen.normexp <- contrasts.fit(fit.sharpen.normexp, contr.all)
contr.sharpen.normexp <- eBayes(contr.sharpen.normexp)
w.probes <- se.exprs(bsd.sharpen.normexp)^2*bsd.sharpen.normexp@assayData$NoBeads

fit.pwts.normexp <- lmFit(normalizeQuantiles(exprs(bsd.sharpen.normexp)), design=design, weights=w.probes)

contr.pwts <- contrasts.fit(fit.pwts.normexp, contr.all)
contr.pwts.normexp <- eBayes(contr.pwts)

spikeInd <- match(spikeIDs, rownames(contr.sharpen.subtract$t))
nonspikes <- setdiff(1:nrow(exprs(bsd.sharpen.subtract)), spikeInd)


###########################################################
## Figure 6 - log-odds and log-ratios for 
## contrast between 3pM and 1pM concentrations

methods <- c("no bg", "no bg + weights", "subtract", "subtract + weights", 
              "normexp", "normexp + weights")
col <- c(pal[1],pal[1],pal[2],pal[2],pal[3],pal[3])

transcol <- NULL
for(i in 1:6){
    temp <- col2rgb(col[i])
    transcol[i] <- rgb(temp[1]/256, temp[2]/256,temp[3]/256, alpha=0.3)
}


pdf("Figure6.pdf",width=12, height=8)
par(mar=c(8,4,1,0.15))
par(mfrow=c(1,2))
i <- 44
cont <- cbind(contr.sharpen.nobgc$lods[,i],contr.pwts.nobgc$lods[,i],contr.sharpen.subtract$lods[,i], contr.pwts.subtract$lods[,i],contr.sharpen.normexp$lods[,i], contr.pwts.normexp$lods[,i])
boxplot(as.data.frame(cont[spikeInd,]), main="A", ylim=c(-10,95), names=methods, las=2, 
        ylab="log-odds", col=col,cex.lab=0.8)
boxplot(as.data.frame(cont[nonspikes,]), main="",  names=rep("",6), las=2, col=transcol, 
        add=TRUE, pch="x")

cont <- cbind(contr.sharpen.nobgc$coeff[,i], contr.pwts.nobgc$coeff[,i], contr.sharpen.subtract$coeff[,i],    
              contr.pwts.subtract$coeff[,i], contr.sharpen.normexp$coeff[,i], contr.pwts.normexp$coeff[,i])
boxplot(as.data.frame(cont[spikeInd,]), main="B", # log-ratios for contrast 3pM - 1pM", 
        ylim=c(-0.2,2.5), names=methods, las=2, ylab="log-ratios", col=col,cex.lab=0.8)
boxplot(as.data.frame(cont[nonspikes,]), main="",  names=rep("",6), las=2, add=TRUE,pch="x",col=transcol)
abline(h=log2(3),lty=2)
dev.off()


###########################################################
## Figure 7 - plot intensities of negative control probes

## Read in control probe summary profile (output by BeadStudio)
controlInfo <- read.table("FullSpikeControls.txt", sep="\t",header=TRUE)

## select negatives
nControls <- which(controlInfo[,1] == "negative")
controlInfo <- controlInfo[nControls,]

sigCols <- grep("Signal", colnames(controlInfo))

set.seed(13092007)
s <- sample(1:nrow(controlInfo), 50)
subset <- controlInfo[s,sigCols]

require(Biobase)
med <- rowMedians(subset.1)
ordmed <- order(med)
subset <- subset[ordmed,]
s <- s[ordmed]

pdf("Figure7.pdf",width=12, height=8)
par(mar=c(5.2,4,1.5,0.15))
boxplot(as.data.frame(log2(t(subset))), ylab="log2 intensity",
        las=2, xlab="", names=controlInfo[s,2], col=pal[1], 
        main="Distribution of intensities for 50 negative controls across all arrays")
dev.off()


###########################################################
## Figure 8 - Plotting array intensities in terms of base composition
## quantile normalised, background corrected data is used

## Read annotation file and process sequences (stored in column 8)

## Plot for probes on strip 1 only
M1 <- read.table("NewAnnotationM1.txt", sep="\t", header=T,fill=TRUE,quote="")

M1probestring <- gsub(", ","",toString(as.character(M1[,8])))
M1probestring <- (strsplit(M1probestring,split=NULL))[[1]]

## Ade, Gdes, Cdes, Tdes are the A, G, C, T
## matrices referred to in Method
## ie 48,000 rows by 50 columns
M1Ades <- as.numeric(M1probestring=="A")
M1Ades <- matrix(M1Ades,ncol=50,byrow=T)
M1Gdes <- as.numeric(M1probestring=="G")
M1Gdes <- matrix(M1Gdes,ncol=50,byrow=T)
M1Cdes <- as.numeric(M1probestring=="C")
M1Cdes <- matrix(M1Cdes,ncol=50,byrow=T)
M1Tdes <- as.numeric(M1probestring=="T")
M1Tdes <- matrix(M1Tdes,ncol=50,byrow=T)
rownames(M1Ades) <- rownames(M1Gdes) <- rownames(M1Cdes) <- rownames(M1Tdes) <- M1[,3]

## Match the ids on the array to their position
## in the annotation file
m <- match(M1[,3],rownames(exprs(bsd.sharpen.subtract)))


M1Acount <- rowSums((M1Ades)
M1Gcount <- rowSums(M1Gdes)
M1Ccount <- rowSums(M1Cdes)
M1Tcount <- rowSums(M1Tdes)
M1GCcount <- M1Gcount + M1Ccount
M1GCdes <- M1Gdes + M1Cdes

##Get probes found on the first strip only
M1oddstrips <- which(nchar(M1[,3]) < 9)

temp <- normalizeQuantiles(exprs(bsd.sharpen.subtract))[m,]

var <- log2(NoBeads(bsd.sharpen.subtract)*se.exprs(bsd.sharpen.subtract)^2)[m,1]

pdf("Figure8.pdf", width=12, height=14)
par(mar=c(4,4,1,0.15))
par(mfrow=c(4,2))
boxplot(temp[M1oddstrips,1]~M1Acount[M1oddstrips], col=pal[1], outline=FALSE, 
        ylim=c(4,11), main="Number of A bases", ylab="Normalised log2 intensity",
        varwidth=TRUE)
abline(h=6.38)
boxplot(temp[M1oddstrips,1]~M1Tcount[M1oddstrips], col=pal[1], outline=FALSE, 
        ylim=c(4,11), main="Number of T bases", ylab="Normalised log2 intensity",
        varwidth=TRUE)
abline(h=6.38)
boxplot(temp[M1oddstrips,1]~M1Gcount[M1oddstrips], col=pal[1], outline=FALSE, 
        ylim=c(4,11), main="Number of G bases", ylab="Normalised log2 intensity",
        varwidth=TRUE)
abline(h=6.38)
boxplot(temp[M1oddstrips,1]~M1Ccount[M1oddstrips], col=pal[1], outline=FALSE, 
        ylim=c(4,11), main="Number of C bases", ylab="Normalised log2 intensity",
        varwidth=TRUE)
abline(h=6.38)
boxplot(temp[M1oddstrips,1]~M1GCcount[M1oddstrips], col=pal[1], outline=FALSE, 
        ylim=c(4,11), main="GC count", ylab="Normalised log2 intensity",
        varwidth=TRUE)
abline(h=6.38)
boxplot(var[M1oddstrips]~M1GCcount[M1oddstrips], col=pal[1], outline=FALSE, 
        main="GC Count", ylab="log2 variance", varwidth=TRUE)
abline(h=-2.17)

## Fit the model to estimate fhe relative effect of having
## an A, C, G at each position
l <- lm(temp[M1oddstrips,1]~M1Ades[M1oddstrips,]+M1Cdes[M1oddstrips,]+M1Gdes[M1oddstrips,])
plot(1:50,l$coefficients[2:51], type="n", xlab="Base Position", 
     main="Effect of base on log2 intensity relative to T", ylab="Effect Size", 
     ylim=c(-0.4,0.4))
text(1:50,l$coefficients[2:51], labels="A", col="red")
text(1:50,l$coefficients[52:101], labels="C", col="blue")
text(1:50,l$coefficients[102:151], labels="G", col="green")
abline(h=0)

l <- lm(var[M1oddstrips]~M1Ades[M1oddstrips,]+M1Cdes[M1oddstrips,]+M1Gdes[M1oddstrips,])
plot(1:50,l$coefficients[2:51], type="n", xlab="Base Position", 
     main="Effect of base on log2 variance relative to T", ylim=c(-0.4,0.4), ylab="Effect Size")
text(1:50,l$coefficients[2:51], labels="A", col="red")
text(1:50,l$coefficients[52:101], labels="C", col="blue")
text(1:50,l$coefficients[102:151], labels="G", col="green")
abline(h=0)

dev.off()


###########################################################
## Plot the discrepancy in ranking of log-odds due to GC count
pdf("Figure9.pdf", width=12, height=8)
par(mar=c(4,4,1.5,0.15))

## Choose the contrast 3pM - 1pM
colnames(contr.sharpen.subtract$lods)[44]
i <- 44
rk <- rank(contr.sharpen.subtract$lods[,i])[m]
boxplot(rk[M1oddstrips]~M1GCcount[M1oddstrips], col=pal[1], 
        ylab="Ranking based on log-odds", xlab="GC count", 
        main="3pM - 1pM", varwidth=TRUE, outline=TRUE, pch="x") 
dev.off()

###############################################################################
################## Annotation Details ##########################
###############################################################################
## Bead types on Strip 2 have 9 character probe IDs, everything else is on strip1

M1oddstrips <- which(nchar(M1[,3]) < 9)
M1evenstrips <- which(nchar(M1[,3]) == 9)


## Column 21 contains comments about the mapping of the probe
## -including intronic and intergenic matches
intronicGenes <- grep("Intronic", M1[,21])
intergenicGenes <- grep("Intergenic", M1[,21])

probeCats <- matrix(nrow=3, ncol=2)
colnames(probeCats) <- c("First Strip", "Second Strip")
rownames(probeCats) <- c("Total","Intronic", "Intergenic")

probeCats[1,1] <- length(M1oddstrips)
probeCats[1,2] <- length(M1evenstrips)

probeCats[2,1] <- length(intersect(intronicGenes, M1oddstrips))
probeCats[2,2] <- length(intersect(intronicGenes, M1evenstrips))

probeCats[3,1] <- length(intersect(intergenicGenes, M1oddstrips))
probeCats[3,2] <- length(intersect(intergenicGenes, M1evenstrips))

## Percentage of intronic probes on 1st strip
probeCats[2,1]/probeCats[1,1]*100

## Percentage of intergenic probes on 1st strip
probeCats[3,1]/probeCats[1,1]*100

## Percentage of intronic probes on 2nd strip
probeCats[2,2]/probeCats[1,2]*100

## Percentage of intergenic probes on 2nd strip
probeCats[3,2]/probeCats[1,2]*100

##Numbers of genes found on strip 1 and strip 2 and their intersection
strip1Transcripts <- M1[M1oddstrips,5]
strip2Transcripts <- M1[M1evenstrips,5]

strip1Genes <- unique(M1[M1oddstrips,4])
strip2Genes <- unique(M1[M1evenstrips,4])

length(intersect(strip1Genes, strip2Genes))


###############################################################################
######Probe Thermodynamics calculations
###############################################################################
## thermodynamic properties calculated using Citation: Kibbe WA. 
## 'OligoCalc: an online oligonucleotide properties calculator'. (2007) Nucleic Acids Res. 35(webserver issue): May 25.  

## Load pre-calculated thermodynamic properties
spikeTP <- read.csv("spikeThermoProperties.csv")

## calculate thermodynamic properties
deltaH <- c(-8,-5.6,-6.6,-8.2,-6.6,-8.8,-9.4,-11.8,-10.5,-10.9)
deltaS <- c(-21.9,-15.2,-18.4,-21,-16.4,-23.5,-25.5,-29,-26.4,-28.4)
deltaG <- c(-1.2,-0.9,-0.9,-1.7,-1.5,-1.5,-1.5,-2.8,-2.3,-2.1)

deltaH <- deltaH[c(1,2,3,5,4,7,6,8,10,9)]
deltaS <- deltaS[c(1,2,3,5,4,7,6,8,10,9)]
deltaG <- deltaG[c(1,2,3,5,4,7,6,8,10,9)]

listH <- rep(NA,17)
listS <- rep(NA,17)
listG <- rep(NA,17)

for(k in 1:17){
    strvec <- (strsplit(as.character(spikeTP$Seq[k]),"")[[1]])
    dh <- 0
    ds <- 0
    dg <- 0
    for(i in 1:(length(strvec)-1)){
        if(strvec[i]=="G"){
            if(strvec[i+1]=="G"){
                ind <- 9
            }
            if(strvec[i+1]=="C"){
                ind <- 10
            }
            if(strvec[i+1]=="A"){
                ind <- 7
            }
            if(strvec[i+1]=="T"){
                ind<-6
            }
        }
        if(strvec[i]=="C"){
            if(strvec[i+1]=="G"){
                ind <- 8
            }
            if(strvec[i+1]=="C"){
                ind <- 9
            }
            if(strvec[i+1]=="A"){
                ind <- 5
            }
            if(strvec[i+1]=="T"){
                ind <- 4
            }
        }
        if(strvec[i]=="A"){
            if(strvec[i+1]=="G"){
                ind <- 4
            }
            if(strvec[i+1]=="C"){
                ind <- 6
            }
            if(strvec[i+1]=="A"){
                ind <- 1
            }
            if(strvec[i+1]=="T"){
                ind <- 2
            }
        }
        if(strvec[i]=="T"){
            if(strvec[i+1]=="G"){
                ind <- 5
            }
            if(strvec[i+1]=="C"){
                ind <- 7
            }
            if(strvec[i+1]=="A"){
                ind <- 3
            }
            if(strvec[i+1]=="T"){
                ind <- 1
            }
        }
        dh <- dh+deltaH[ind]
        dg <- dg+deltaG[ind]
        ds <- ds+deltaS[ind]
    }
#dh<-dh+0.6
#dg<-dg+3.4
#ds<-ds-9

    listH[k] <- dh
    listS[k] <- ds
    listG[k] <- dg
}

listH <- -listH
listG <- -listG
listS <- -listS

##check both sets agree
listS - spiketherm[,8]
listH - spiketherm[,7]
listG - spiketherm[,6]


###########################################################
## Draw figure 10
pdf("Figure10.pdf", width=10, height=10)
par(mar=c(5.2,4,1.5,0.15))
plot(cor(cbind(spikeTP[,2:13], spikeTP[,15:19]))[1:12,14], axes=F, xlab="concentration", 
      ylab="correlation of intensity with ...", pch=16)
points(1:12,cor(cbind(spikeTP[,2:13],spikeTP[,15:19]))[1:12,15],pch=16,col="purple")
points(1:12,cor(cbind(spikeTP[,2:13],spikeTP[,15:19]))[1:12,16],pch=16,col="green")
points(1:12,cor(cbind(spikeTP[,2:13],spikeTP[,15:19]))[1:12,17],pch=16,col="blue")
lines(1:12,cor(cbind(spikeTP[,2:13],spikeTP[,15:19]))[1:12,14],lty=2)
lines(1:12,cor(cbind(spikeTP[,2:13],spikeTP[,15:19]))[1:12,15],lty=2,col="purple")
lines(1:12,cor(cbind(spikeTP[,2:13],spikeTP[,15:19]))[1:12,16],lty=2,col="green")
lines(1:12,cor(cbind(spikeTP[,2:13],spikeTP[,15:19]))[1:12,17],lty=2,col="blue")
axis(1,at=1:12,labels=c(1000,300,100,30,10,3,1,0.3,0.1,0.03,0.01,0))
axis(2)
box()
legend(9,.6,legend=c("melting temp","free energy","enthalpy","entropy"),fill=c("black","purple","green","blue"))
dev.off()

## Read in the negative control data
#negcon<-read.csv("negativecontrols.csv")
#colnames(negcon)<-c("species","empty","ID","seq","olicode","group","groupID")

negcon <- read.table("ControlSequences.txt", sep="\t",header=T)
colnames(negcon) <- c("species","empty","ID","seq","olicode","group","groupID")
negperf <- read.table("FullSpikeControls.txt", sep="\t",header=T)

#negperf<-read.table("ControlsExpress.txt",sep="\t",header=T)
negperf <- negperf[match(negcon$olicode,negperf$ProbeID),]


## calculate thermodynamic values for negative controls
negconT <- matrix(NA,ncol=3,nrow=1688)
for(k in 1:1688){
    strvec <- (strsplit(as.character(negcon$seq[k]),"")[[1]])
    dh <- 0
    ds <- 0
    dg <- 0

    for(i in 1:(length(strvec)-1)){
        if(strvec[i]=="G"){
            if(strvec[i+1]=="G"){
                ind <- 9
            }
            if(strvec[i+1]=="C"){
                ind <- 10
            }
            if(strvec[i+1]=="A"){
                ind <- 7
            }
            if(strvec[i+1]=="T"){
                ind <- 6
            }
        }
        if(strvec[i]=="C"){
            if(strvec[i+1]=="G"){
                ind <- 8
            }
            if(strvec[i+1]=="C"){
                ind <- 9
            }
            if(strvec[i+1]=="A"){
                ind <- 5
            }
            if(strvec[i+1]=="T"){
                ind <- 4
            }
        }
        if(strvec[i]=="A"){
            if(strvec[i+1]=="G"){
                ind <- 4
            }
            if(strvec[i+1]=="C"){
                ind <- 6
            }
            if(strvec[i+1]=="A"){
                ind <- 1
            }
            if(strvec[i+1]=="T"){
                ind <- 2
            }
        }
        if(strvec[i]=="T"){
            if(strvec[i+1]=="G"){
                ind <- 5
            }
            if(strvec[i+1]=="C"){
                ind <- 7
            }
            if(strvec[i+1]=="A"){
                ind <- 3
            }
            if(strvec[i+1]=="T"){
                ind <- 1
            }
        }
        dh <- dh+deltaH[ind]
        dg <- dg+deltaG[ind]
        ds <- ds+deltaS[ind]
    }
    dh <- dh+0.6
    dg <- dg+3.4
    ds <- ds-9
    
#if((strvec[length(strvec)]=="A")|(strvec[length(strvec)]=="T")){
#dh<-dh+3.72
#dg<-dg+0.45
#ds<-ds+10.5
#}

    negconT[k,1] <- dh
    negconT[k,2] <- ds
    negconT[k,3] <- dg
}
colnames(negconT) <- c("H","S","G")


## correlations
summary((cor(cbind(negconT[(negcon$olicode<100000000)&(negcon$groupID=="permuted_negative"),],
                   log2(negperf[(negcon$olicode<100000000)&(negcon$groupID=="permuted_negative"),seq(3,97,2)]))))[4:51,1])
summary((cor(cbind(negconT[(negcon$olicode<100000000)&(negcon$groupID=="permuted_negative"),],
                   log2(negperf[(negcon$olicode<100000000)&(negcon$groupID=="permuted_negative"),seq(3,97,2)]))))[4:51,2])
summary((cor(cbind(negconT[(negcon$olicode<100000000)&(negcon$groupID=="permuted_negative"),],
                   log2(negperf[(negcon$olicode<100000000)&(negcon$groupID=="permuted_negative"),seq(3,97,2)]))))[4:51,3])


###############################################################################
#######Supplementary Plot of false discovery curves####
###############################################################################

## Number of false discoveries is defined as number of non-spikes 
nprobes <- nrow(exprs(bsd.sharpen.subtract))
fp <- rep(1, nprobes)
fp[spikeInd] <- 0
sel <- 1:50
ylim <- c(0,50)
nmethods <- 6

methods <- c("no bg", "no bg + weights", "subtract", "subtract + weights", "normexp", "normexp + weights")
ncontr <- ncol(contr.all)
col <- c(pal[1],pal[1],pal[2],pal[2],pal[3],pal[3])

pdf("fdr.pdf")
for(i in 1:ncontr) {
    plot(sel, cumsum(fp[order(abs(contr.sharpen.nobgc$t)[,i], decreasing=TRUE)][sel]), 
         xlab="Number selected", ylab="Number of false discoveries", 
         main=paste(colnames(contr.all)[i], "fold-change"), type="l", col=col[1])
    points(sel, cumsum(fp[order(abs(contr.pwts.nobgc$t)[,i], decreasing=TRUE)][sel]), type="l", col=col[2],lty=2)
    points(sel, cumsum(fp[order(abs(contr.sharpen.subtract$t)[,i], decreasing=TRUE)][sel]), type="l", col=col[3])
    points(sel, cumsum(fp[order(abs(contr.pwts.subtract$t)[,i], decreasing=TRUE)][sel]), type="l", col=col[4],lty=2)
    points(sel, cumsum(fp[order(abs(contr.sharpen.normexp$t)[,i], decreasing=TRUE)][sel]), type="l", col=col[5])	
    points(sel, cumsum(fp[order(abs(contr.pwts.normexp$t)[,i], decreasing=TRUE)][sel]), type="l", col=col[6],lty=2)	
    legend("topleft", legend=methods, col=col, cex=0.5, lty=c(1,2))  
}
dev.off()


library(geneplotter)
library(beadarray)

## Diagnostic plots of bead level data
## Image plots of foreground
load("bld.sharpen.nobgc.rda")

arraynms <- arrayNames(bld.sharpen.nobgc)
narrays <- length(arraynms)

zlim <- c(8.5,10.5)
pdf("imageplots.G.pdf", width=4, height=10)
par(mfrow=c(12,1), mai=c(0.3,0.1,0.3,0.1), oma=c(0,0,0,0))
for(i in 1:narrays) {
    imageplot(bld.sharpen.nobgc, array=i, main=arraynms[i], log=TRUE, ncol=100, nrow=50, 
              high="blue", low="white", whatToPlot="G", zlim=zlim)
    cat(i, " ")
}
dev.off()

zlim <- c(9,9.5)
pdf("imageplots.Gb.pdf", width=4, height=10)
par(mfrow=c(12,1), mai=c(0.3,0.1,0.3,0.1), oma=c(0,0,0,0))
for(i in 1:narrays) {
    imageplot(bld.sharpen.nobgc, array=i, main=arraynms[i], log=TRUE, ncol=100, nrow=50, 
              high="blue", low="white", whatToPlot="Gb", zlim=zlim)
    cat(i, " ")
}
dev.off()


arraycols <- rep(rainbow(8), each=12)
ylim <- c(8, 11)
png("boxplots.G.png", width=1024, height=640)
boxplotBeads(bld.sharpen.nobgc, names=seq(1:narrays), las=2, outline=FALSE, 
             col=arraycols, main="Raw Foreground", what="G", ylim=ylim, 
             ylab=expression(log[2](intensity)), medlwd=1)
dev.off()

png("boxplots.Gb.png", width=1024, height=640)
boxplotBeads(bld.sharpen.nobgc, names=seq(1:narrays), las=2, outline=FALSE, 
             col=arraycols, main="Raw Background", what="Gb", ylim=ylim, 
             ylab=expression(log[2](intensity)), medlwd=1)
dev.off()

rm(bld.sharpen.nobgc)
gc()

## Diagnostic plots for bead summary data  
load("bsd.log.ill.rda")

narrays <- ncol(exprs(bsd.sharpen.subtract))
ngenes <- nrow(exprs(bsd.sharpen.subtract))
arraynms <- gsub("_1", "", colnames(exprs(bsd.sharpen.subtract)))
medsig <- rowMedians(exprs(bsd.sharpen.subtract), na.rm=TRUE)

## scatter plots of intensities from individual arrays versus median array
pdf("smoothscatter.summarydata.versus.median.pdf", width=10, height=6)
par(mfrow=c(2,3), mai=c(0.5,0.5,0.3,0.1), oma=c(1.1,1.1,0,0))
for(i in 1:narrays) {
    smoothScatter(medsig, exprs(bsd.sharpen.subtract)[,i], main=arraynms[i], ylab="", xlab="")
    abline(0,1,col=2, lty=2, lwd=2)
    if(i%%6==0) {
        mtext(expression(log[2]("median intensity")), side=1, outer=TRUE)
        mtext(expression(log[2](intensity)), side=2, outer=TRUE)
    }
}
dev.off()

arraycols <- rep(rainbow(8), each=6)

## boxplots of summarised data
pdf("boxplots.summarised.data.pdf", width=10, height=6)
boxplot(as.data.frame(exprs(bsd.sharpen.subtract)), names=seq(1:narrays), xlab="array", 
        ylab=expression(log[2](intensity)), main="Summarised expression data per array", 
        las=2, col=arraycols)
boxplot(as.data.frame(se.exprs(bsd.sharpen.subtract)), names=seq(1:narrays), xlab="array", 
        ylab="Standard error", main="Standard errors per array", las=2, ylim=c(0,0.5), 
        col=arraycols)
boxplot(as.data.frame(NoBeads(bsd.sharpen.subtract)), names=seq(1:narrays), xlab="array", 
        ylab="Number of beads per bead type", main="Number of beads per array", las=2, 
        col=arraycols)
dev.off()

png("boxplots.means.png", width=1024, height=640)
boxplot(as.data.frame(exprs(bsd.sharpen.subtract)), names=seq(1:narrays), xlab="array", 
        ylab=expression(log[2](intensity)), main="Summarised expression data per array", 
        las=2, col=arraycols)
dev.off()

png("boxplots.standarderrors.png", width=1024, height=640)
boxplot(as.data.frame(se.exprs(bsd.sharpen.subtract)), names=seq(1:narrays), xlab="array", 
        ylab="Standard error", main="Standard errors per array", las=2, ylim=c(0,0.5), 
        col=arraycols)
dev.off()

png("boxplots.numbeads.png", width=1024, height=640)
boxplot(as.data.frame(NoBeads(bsd.sharpen.subtract)), names=seq(1:narrays), xlab="array", 
        ylab="Number of beads per bead type", main="Number of beads per array", las=2, 
        col=arraycols)
dev.off()
