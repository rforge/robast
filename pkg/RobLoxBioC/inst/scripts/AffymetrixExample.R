###############################################################################
## Affymetrix example
###############################################################################

###############################################################################
## References:
## L.M. Cope1 , R.A. Irizarry, H.A. Jaffee, Z. Wu and T.P. Speed (2004):
## "A benchmark for Affymetrix GeneChip expression measures". 
## Bioinformatics 20(3): 323-331
##
## R.A. Irizarry1, Z. Wu and H.A. Jaffee (2006):
## "Comparison of Affymetrix GeneChip expression measures"
## Bioinformatics 22(7): 789-794
###############################################################################

###############################################################################
## Data
## Spike-in hgu95a data:
## http://www.biostat.jhsph.edu/~ririzarr/affycomp/spikein.tgz
##
## Spike-in hgu133a data:
## http://www.biostat.jhsph.edu/~ririzarr/affycomp/hgu133spikein.tgz
##
## Dilution data: 
## The links:
## http://qolotus02.genelogic.com/datasets.nsf/               (no response)
## and
## http://www.genelogic.com/media/studies/dilution.cfm        (not found)
## seem not to lead to the data any longer.
## An email to the support of genelogic is still unanswered ...
###############################################################################


library(affy)
library(affycomp) ## Version 1.19.4

###################
## replace with your path to hgu95a data!!!
PATH <- "./spikein"
###################

fn <- list.celfiles(path = PATH, full.names=TRUE)
data(spikein.phenodata)
(spikein.hgu95a <- read.affybatch(filenames = fn, 
                                  phenoData = spikein.phenodata))

###################
## replace with your path to hgu133a data!!!
PATH <- "./hgu133spikein"
###################

fn <- list.celfiles(path = PATH, full.names=TRUE)
fn <- fn[c(seq(1,40,3), seq(2, 41, 3), seq(3, 42, 3))]
fn <- fn[c(6:14, 1:5, 20:28, 15:19, 34:42, 29:33)]
data(hgu133a.spikein.phenodata)

## Attention:
## Order of filenames in fn has to be identical to 
## sampleNames(hgu133a.spikein.phenodata)!!!

(spikein.hgu133a <- read.affybatch(filenames = fn,
                                  phenoData = hgu133a.spikein.phenodata))


###########################################################
## assessments for MAS 5.0 and RMA including dilution data from package affycomp
###########################################################
data(mas5.assessment)
data(rma.assessment)
data(mas5.assessment.133)
data(rma.assessment.133)

res.rma <- rma(spikein.hgu133a)

library(RobLoxBioC)
###########################################################
## Example 1: Analogous to "classical" MAS 5.0
## computation takes about 38 sec on Intel P9500 (64bit Linux, 4 GByte RAM)
###########################################################
## hgu95a
system.time(eset.hgu95a <- robloxbioc(spikein.hgu95a, normalize = TRUE, add.constant = 0))
eset.hgu95a.log2 <- eset.hgu95a
exprs(eset.hgu95a.log2) <- log2(exprs(eset.hgu95a))
roblox.hgu95a <- assessSpikeIn(eset.hgu95a.log2, method.name = "roblox")

## hgu133a
system.time(eset.hgu133a <- robloxbioc(spikein.hgu133a, normalize = TRUE, add.constant = 0))
eset.hgu133a.log2 <- eset.hgu133a
exprs(eset.hgu133a.log2) <- log2(exprs(eset.hgu133a))
roblox.hgu133a <- assessSpikeIn(eset.hgu133a.log2, method.name = "roblox")


###########################################################
## Example 2: MAS 5.0 + 32
## computation takes about 38 sec on Intel P9500 (64bit Linux, 4 GByte RAM)
###########################################################
## hgu95a
system.time(eset.hgu95a32 <- robloxbioc(spikein.hgu95a, normalize = TRUE, add.constant = 32))
eset.hgu95a.log232 <- eset.hgu95a32
exprs(eset.hgu95a.log232) <- log2(exprs(eset.hgu95a32))
roblox.hgu95a32 <- assessSpikeIn(eset.hgu95a.log232, method.name = "roblox + 32")

## hgu133a
system.time(eset.hgu133a32 <- robloxbioc(spikein.hgu133a, normalize = TRUE, add.constant = 32))
eset.hgu133a.log232 <- eset.hgu133a32
exprs(eset.hgu133a.log232) <- log2(exprs(eset.hgu133a32))
roblox.hgu133a32 <- assessSpikeIn(eset.hgu133a.log232, method.name = "roblox + 32")


###########################################################
## Example 3: MAS 5.0 with PM only
## computation takes about 19 sec on Intel P9500 (64bit Linux, 4 GByte RAM)
###########################################################
## hgu95a
system.time(eset.hgu95a.pmonly <- robloxbioc(spikein.hgu95a, bg.correct = TRUE, 
                                             pmcorrect = FALSE, normalize = TRUE, 
                                             add.constant = 0))
eset.hgu95a.log2.pmonly <- eset.hgu95a.pmonly
exprs(eset.hgu95a.log2.pmonly) <- log2(exprs(eset.hgu95a.pmonly))
roblox.hgu95a.pmonly <- assessSpikeIn(eset.hgu95a.log2.pmonly, method.name = "roblox (PM)")

## hgu133a
system.time(eset.hgu133a.pmonly <- robloxbioc(spikein.hgu133a, bg.correct = TRUE, 
                                             pmcorrect = FALSE, normalize = TRUE, 
                                             add.constant = 0))
eset.hgu133a.log2.pmonly <- eset.hgu133a.pmonly
exprs(eset.hgu133a.log2.pmonly) <- log2(exprs(eset.hgu133a.pmonly))
roblox.hgu133a.pmonly <- assessSpikeIn(eset.hgu133a.log2.pmonly, method.name = "roblox (PM)")


###########################################################
## Example 4: MAS 5.0 + 32 with PM only
## computation takes about 19 sec on Intel P9500 (64bit Linux, 4 GByte RAM)
###########################################################
## hgu95a
system.time(eset.hgu95a.pmonly32 <- robloxbioc(spikein.hgu95a, bg.correct = TRUE, 
                                             pmcorrect = FALSE, normalize = TRUE, 
                                             add.constant = 32))
eset.hgu95a.log2.pmonly32 <- eset.hgu95a.pmonly32
exprs(eset.hgu95a.log2.pmonly32) <- log2(exprs(eset.hgu95a.pmonly32))
roblox.hgu95a.pmonly32 <- assessSpikeIn(eset.hgu95a.log2.pmonly32, method.name = "roblox + 32 (PM)")

## hgu133a
system.time(eset.hgu133a.pmonly32 <- robloxbioc(spikein.hgu133a, bg.correct = TRUE, 
                                             pmcorrect = FALSE, normalize = TRUE, 
                                             add.constant = 32))
eset.hgu133a.log2.pmonly32 <- eset.hgu133a.pmonly32
exprs(eset.hgu133a.log2.pmonly32) <- log2(exprs(eset.hgu133a.pmonly32))
roblox.hgu133a.pmonly32 <- assessSpikeIn(eset.hgu133a.log2.pmonly32, method.name = "roblox + 32 (PM)")


###############################################################################
## Figure 1: The MA plot shows log fold change as a function of mean log 
## expression level. A set of 14 arrays representing a single experiment from 
## the Affymetrix spike-in data are used for this plot. A total of 13 sets of 
## fold changes are generated by comparing the first array in the set to each 
## of the others. Genes are symbolized by numbers representing the nominal 
## log2 fold change for the gene. Non-differentially expressed genes with 
## observed fold changes larger than 2 are plotted in red. All other probesets 
## are represented with black dots.
###############################################################################
## hgu95a
par(mfrow = c(3, 2))
affycompPlot(roblox.hgu95a$MA)
affycompPlot(roblox.hgu95a32$MA)
affycompPlot(roblox.hgu95a.pmonly$MA)
affycompPlot(roblox.hgu95a.pmonly32$MA)
affycompPlot(mas5.assessment$MA)
affycompPlot(rma.assessment$MA)

## hgu133a
par(mfrow = c(3, 2))
affycompPlot(roblox.hgu133a$MA)
affycompPlot(roblox.hgu133a32$MA)
affycompPlot(roblox.hgu133a.pmonly$MA)
affycompPlot(roblox.hgu133a.pmonly32$MA)
affycompPlot(mas5.assessment.133$MA)
affycompPlot(rma.assessment.133$MA)


###############################################################################
## Figure 4a: Average observed log2 intensity plotted against nominal log2 
## concentration for each spiked-in gene for all arrays in Affymetrix spike-In 
## experiment
###############################################################################
## hgu95a
par(mfrow = c(3, 2))
affycomp.figure4a(roblox.hgu95a$Signal)
affycomp.figure4a(roblox.hgu95a32$Signal)
affycomp.figure4a(roblox.hgu95a.pmonly$Signal)
affycomp.figure4a(roblox.hgu95a.pmonly32$Signal)
affycomp.figure4a(mas5.assessment$Signal)
affycomp.figure4a(rma.assessment$Signal)

## hgu133a
par(mfrow = c(3, 2))
affycomp.figure4a(roblox.hgu133a$Signal)
affycomp.figure4a(roblox.hgu133a32$Signal)
affycomp.figure4a(roblox.hgu133a.pmonly$Signal)
affycomp.figure4a(roblox.hgu133a.pmonly32$Signal)
affycomp.figure4a(mas5.assessment.133$Signal)
affycomp.figure4a(rma.assessment.133$Signal)

## Comparison plot
## hgu95a
affycomp.compfig4a(list(roblox.hgu95a$Signal, 
                        roblox.hgu95a32$Signal,
                        roblox.hgu95a.pmonly$Signal,
                        roblox.hgu95a.pmonly32$Signal,
                        mas5.assessment$Signal, 
                        rma.assessment$Signal), 
                   method.names = c("roblox", "roblox + 32", "roblox (PM)", 
                                    "roblox + 32 (PM)", "MAS 5.0", "RMA"))

## hgu133a
affycomp.compfig4a(list(roblox.hgu133a$Signal, 
                        roblox.hgu133a32$Signal,
                        roblox.hgu133a.pmonly$Signal,
                        roblox.hgu133a.pmonly32$Signal,
                        mas5.assessment.133$Signal, 
                        rma.assessment.133$Signal), 
                   method.names = c("roblox", "roblox + 32", "roblox (PM)", 
                                    "roblox + 32 (PM)", "MAS 5.0", "RMA"))


###############################################################################
## Figure 5: A typical identification rule for differential expression filters 
## genes with fold change exceeding a given threshold. This figure shows 
## average ROC curves which offer a graphical representation of both 
## specificity and sensitivity for such a detection rule. 
## a) Average ROC curves based on comparisons with nominal fold changes 
## ranging from 2 to 4096. 
## b) As a) but with nominal fold changes equal to 2.
###############################################################################
## Figure 5a:
## hgu95a
par(mfrow = c(3, 2))
affycomp.figure5a(roblox.hgu95a$FC)
affycomp.figure5a(roblox.hgu95a32$FC)
affycomp.figure5a(roblox.hgu95a.pmonly$FC)
affycomp.figure5a(roblox.hgu95a.pmonly32$FC)
affycomp.figure5a(mas5.assessment$FC)
affycomp.figure5a(rma.assessment$FC)

## hgu133a
par(mfrow = c(3, 2))
affycomp.figure5a(roblox.hgu133a$FC)
affycomp.figure5a(roblox.hgu133a32$FC)
affycomp.figure5a(roblox.hgu133a.pmonly$FC)
affycomp.figure5a(roblox.hgu133a.pmonly32$FC)
affycomp.figure5a(mas5.assessment.133$FC)
affycomp.figure5a(rma.assessment.133$FC)

## Comparison plot
## hgu95a
affycomp.compfig5a(list(roblox.hgu95a$FC,
                        roblox.hgu95a32$FC,
                        roblox.hgu95a.pmonly$FC,
                        roblox.hgu95a.pmonly32$FC,
                        mas5.assessment$FC, 
                        rma.assessment$FC), 
                   method.names = c("roblox", "roblox + 32", "roblox (PM)", 
                                    "roblox + 32 (PM)", "MAS 5.0", "RMA"))

## hgu133a
affycomp.compfig5a(list(roblox.hgu133a$FC, 
                        roblox.hgu133a32$FC,
                        roblox.hgu133a.pmonly$FC,
                        roblox.hgu133a.pmonly32$FC,
                        mas5.assessment$FC, 
                        rma.assessment$FC), 
                   method.names = c("roblox", "roblox + 32", "roblox (PM)", 
                                    "roblox + 32 (PM)", "MAS 5.0", "RMA"))

## Figure 5b:
## hgu95a
par(mfrow = c(3, 2))
affycomp.figure5b(roblox.hgu95a$FC)
affycomp.figure5b(roblox.hgu95a32$FC)
affycomp.figure5b(roblox.hgu95a.pmonly$FC)
affycomp.figure5b(roblox.hgu95a.pmonly32$FC)
affycomp.figure5b(mas5.assessment$FC)
affycomp.figure5b(rma.assessment$FC)

## hgu133a
par(mfrow = c(3, 2))
affycomp.figure5b(roblox.hgu133a$FC)
affycomp.figure5b(roblox.hgu133a32$FC)
affycomp.figure5b(roblox.hgu133a.pmonly$FC)
affycomp.figure5b(roblox.hgu133a.pmonly32$FC)
affycomp.figure5b(mas5.assessment.133$FC)
affycomp.figure5b(rma.assessment.133$FC)

## Comparison plot
## hgu95a
affycomp.compfig5b(list(roblox.hgu95a$FC, 
                        roblox.hgu95a32$FC,
                        roblox.hgu95a.pmonly$FC,
                        roblox.hgu95a.pmonly32$FC,
                        mas5.assessment$FC, 
                        rma.assessment$FC), 
                   method.names = c("roblox", "roblox + 32", "roblox (PM)", 
                                    "roblox + 32 (PM)", "MAS 5.0", "RMA"))

## hgu133a
affycomp.compfig5b(list(roblox.hgu133a$FC, 
                        roblox.hgu133a32$FC,
                        roblox.hgu133a.pmonly$FC,
                        roblox.hgu133a.pmonly32$FC,
                        mas5.assessment$FC, 
                        rma.assessment$FC), 
                   method.names = c("roblox", "roblox + 32", "roblox (PM)", 
                                    "roblox + 32 (PM)", "MAS 5.0", "RMA"))


###############################################################################
## Figure 6:            
## a) Observed log fold changes plotted against nominal log fold changes. The
## dashed lines represent highest, 25th highest, 100th highest, 25th 
## percentile, 75th percentile, smallest 100th, smallest 25th, and smallest 
## log fold change for the genes that were not differentially expressed. 
## b) Like a) but the observed fold changes were calculated for spiked in 
## genes with nominal concentrations no higher than 2pM.
###############################################################################
## Figure 6a:
## hgu95a
par(mfrow = c(3, 2))
affycomp.figure6a(roblox.hgu95a$FC)
affycomp.figure6a(roblox.hgu95a32$FC)
affycomp.figure6a(roblox.hgu95a.pmonly$FC)
affycomp.figure6a(roblox.hgu95a.pmonly32$FC)
affycomp.figure6a(mas5.assessment$FC)
affycomp.figure6a(rma.assessment$FC)

## hgu133a
par(mfrow = c(3, 2))
affycomp.figure6a(roblox.hgu133a$FC)
affycomp.figure6a(roblox.hgu133a32$FC)
affycomp.figure6a(roblox.hgu133a.pmonly$FC)
affycomp.figure6a(roblox.hgu133a.pmonly32$FC)
affycomp.figure6a(mas5.assessment.133$FC)
affycomp.figure6a(rma.assessment.133$FC)

## Figure 6b:
## hgu95a
par(mfrow = c(3, 2))
affycomp.figure6b(roblox.hgu95a$FC)
affycomp.figure6b(roblox.hgu95a32$FC)
affycomp.figure6b(roblox.hgu95a.pmonly$FC)
affycomp.figure6b(roblox.hgu95a.pmonly32$FC)
affycomp.figure6b(mas5.assessment$FC)
affycomp.figure6b(rma.assessment$FC)

## hgu133a
par(mfrow = c(3, 2))
affycomp.figure6b(roblox.hgu133a$FC)
affycomp.figure6b(roblox.hgu133a32$FC)
affycomp.figure6b(roblox.hgu133a.pmonly$FC)
affycomp.figure6b(roblox.hgu133a.pmonly32$FC)
affycomp.figure6b(mas5.assessment.133$FC)
affycomp.figure6b(rma.assessment.133$FC)


###############################################################################
## Table
###############################################################################
## hgu95a
round(tableAll(roblox.hgu95a, roblox.hgu95a32, roblox.hgu95a.pmonly, 
               roblox.hgu95a.pmonly32, mas5.assessment, rma.assessment), 4)

## hgu133a
round(tableAll(roblox.hgu133a, roblox.hgu133a32, roblox.hgu133a.pmonly, 
               roblox.hgu133a.pmonly32, mas5.assessment.133, rma.assessment.133), 4)

## smaller table, more informative ...
## hgu95a
affycompTable(roblox.hgu95a, mas5.assessment, rma.assessment)

## hgu133a
affycompTable(roblox.hgu133a, mas5.assessment.133, rma.assessment.133)
