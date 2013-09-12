pkgname <- "RobAStBase"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('RobAStBase')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("0RobAStBase-package")
### * 0RobAStBase-package

flush(stderr()); flush(stdout())

### Name: RobAStBase-package
### Title: Robust Asymptotic Statistics
### Aliases: RobAStBase-package RobAStBase
### Keywords: package

### ** Examples

library(RobAStBase)

## some L2 differentiable parametric family from package distrMod, e.g.
B <- BinomFamily(size = 25, prob = 0.25) 

## classical optimal IC
IC0 <- optIC(model = B, risk = asCov())
plot(IC0) # plot IC
checkIC(IC0, B)



cleanEx()
nameEx("ALEstimate-class")
### * ALEstimate-class

flush(stderr()); flush(stdout())

### Name: ALEstimate-class
### Title: ALEstimate-class.
### Aliases: ALEstimate-class pIC pIC,ALEstimate-method asbias
###   asbias,ALEstimate-method show,ALEstimate-method
###   confint,ALEstimate,missing-method
###   confint,ALEstimate,symmetricBias-method
###   confint,ALEstimate,onesidedBias-method
###   confint,ALEstimate,asymmetricBias-method
### Keywords: classes

### ** Examples

## prototype
new("ALEstimate")



cleanEx()
nameEx("BdStWeight-class")
### * BdStWeight-class

flush(stderr()); flush(stdout())

### Name: BdStWeight-class
### Title: Robust Weight classes for bounded, standardized weights
### Aliases: BdStWeight-class stand,BdStWeight-method
###   stand<-,BdStWeight-method
### Keywords: classes

### ** Examples

## prototype
new("BdStWeight")



cleanEx()
nameEx("BoundedWeight-class")
### * BoundedWeight-class

flush(stderr()); flush(stdout())

### Name: BoundedWeight-class
### Title: Robust Weight classes for bounded weights
### Aliases: BoundedWeight-class clip,BoundedWeight-method
###   clip<-,BoundedWeight-method
### Keywords: classes

### ** Examples

## prototype
new("BoundedWeight")



cleanEx()
nameEx("ContIC-class")
### * ContIC-class

flush(stderr()); flush(stdout())

### Name: ContIC-class
### Title: Influence curve of contamination type
### Aliases: ContIC-class CallL2Fam<-,ContIC-method cent cent,ContIC-method
###   cent<- cent<-,ContIC-method clip,ContIC-method clip<-
###   clip<-,ContIC-method lowerCase<- lowerCase<-,ContIC-method stand<-
###   stand<-,ContIC-method neighbor,ContIC-method
###   generateIC,ContNeighborhood,L2ParamFamily-method show,ContIC-method
### Keywords: classes

### ** Examples

IC1 <- new("ContIC")
plot(IC1)



cleanEx()
nameEx("ContIC")
### * ContIC

flush(stderr()); flush(stdout())

### Name: ContIC
### Title: Generating function for ContIC-class
### Aliases: ContIC
### Keywords: robust

### ** Examples

IC1 <- ContIC()
plot(IC1)



cleanEx()
nameEx("ContNeighborhood-class")
### * ContNeighborhood-class

flush(stderr()); flush(stdout())

### Name: ContNeighborhood-class
### Title: Contamination Neighborhood
### Aliases: ContNeighborhood-class
### Keywords: classes models

### ** Examples

new("ContNeighborhood")



cleanEx()
nameEx("ContNeighborhood")
### * ContNeighborhood

flush(stderr()); flush(stdout())

### Name: ContNeighborhood
### Title: Generating function for ContNeighborhood-class
### Aliases: ContNeighborhood
### Keywords: models

### ** Examples

ContNeighborhood()

## The function is currently defined as
function(radius = 0){ 
    new("ContNeighborhood", radius = radius) 
}



cleanEx()
nameEx("FixRobModel-class")
### * FixRobModel-class

flush(stderr()); flush(stdout())

### Name: FixRobModel-class
### Title: Robust model with fixed (unconditional) neighborhood
### Aliases: FixRobModel-class neighbor<-,FixRobModel-method
###   show,FixRobModel-method
### Keywords: classes models

### ** Examples

new("FixRobModel")



cleanEx()
nameEx("FixRobModel")
### * FixRobModel

flush(stderr()); flush(stdout())

### Name: FixRobModel
### Title: Generating function for FixRobModel-class
### Aliases: FixRobModel
### Keywords: models

### ** Examples

(M1 <- FixRobModel())

## The function is currently defined as
function(center = ParamFamily(), neighbor = ContNeighborhood()){
    new("FixRobModel", center = center, neighbor = neighbor)
}



cleanEx()
nameEx("HampIC-class")
### * HampIC-class

flush(stderr()); flush(stdout())

### Name: HampIC-class
### Title: Influence curve of Hampel type
### Aliases: HampIC-class lowerCase lowerCase,HampIC-method neighborRadius
###   neighborRadius,HampIC-method neighborRadius<-
###   neighborRadius<-,HampIC-method stand stand,HampIC-method
###   weight,HampIC-method biastype,HampIC-method normtype,HampIC-method
### Keywords: classes

### ** Examples

IC1 <- new("HampIC")
plot(IC1)



cleanEx()
nameEx("HampelWeight-class")
### * HampelWeight-class

flush(stderr()); flush(stdout())

### Name: HampelWeight-class
### Title: Robust Weight classes for weights of Hampel type
### Aliases: HampelWeight-class cent,HampelWeight-method
###   cent<-,HampelWeight-method
### Keywords: classes

### ** Examples

## prototype
new("HampelWeight")



cleanEx()
nameEx("IC-class")
### * IC-class

flush(stderr()); flush(stdout())

### Name: IC-class
### Title: Influence curve
### Aliases: IC-class CallL2Fam CallL2Fam,IC-method CallL2Fam<-
###   CallL2Fam<-,IC-method modifyIC modifyIC,IC-method
###   checkIC,IC,missing-method checkIC,IC,L2ParamFamily-method
###   evalIC,IC,numeric-method evalIC,IC,matrix-method show,IC-method
### Keywords: classes robust

### ** Examples

IC1 <- new("IC")
plot(IC1)



cleanEx()
nameEx("IC")
### * IC

flush(stderr()); flush(stdout())

### Name: IC
### Title: Generating function for IC-class
### Aliases: IC
### Keywords: robust

### ** Examples

IC1 <- IC()
plot(IC1)



cleanEx()
nameEx("InfRobModel-class")
### * InfRobModel-class

flush(stderr()); flush(stdout())

### Name: InfRobModel-class
### Title: Robust model with infinitesimal (unconditional) neighborhood
### Aliases: InfRobModel-class neighbor<-,InfRobModel-method
###   show,InfRobModel-method
### Keywords: classes models

### ** Examples

new("InfRobModel")



cleanEx()
nameEx("InfRobModel")
### * InfRobModel

flush(stderr()); flush(stdout())

### Name: InfRobModel
### Title: Generating function for InfRobModel-class
### Aliases: InfRobModel
### Keywords: models

### ** Examples

(M1 <- InfRobModel())

## The function is currently defined as
function(center = L2ParamFamily(), neighbor = ContNeighborhood()){
    new("InfRobModel", center = center, neighbor = neighbor)
}



cleanEx()
nameEx("InfluenceCurve-class")
### * InfluenceCurve-class

flush(stderr()); flush(stdout())

### Name: InfluenceCurve-class
### Title: Influence curve
### Aliases: InfluenceCurve-class addInfo<- addInfo<-,InfluenceCurve-method
###   addRisk<- addRisk<-,InfluenceCurve-method Curve
###   Curve,InfluenceCurve-method Domain,InfluenceCurve-method Infos
###   Infos,InfluenceCurve-method Infos<- Infos<-,InfluenceCurve-method
###   Map,InfluenceCurve-method name,InfluenceCurve-method
###   name<-,InfluenceCurve-method Range,InfluenceCurve-method Risks
###   Risks,InfluenceCurve-method Risks<- Risks<-,InfluenceCurve-method
###   show,InfluenceCurve-method
### Keywords: classes robust

### ** Examples

new("InfluenceCurve")



cleanEx()
nameEx("InfluenceCurve")
### * InfluenceCurve

flush(stderr()); flush(stdout())

### Name: InfluenceCurve
### Title: Generating function for InfluenceCurve-class
### Aliases: InfluenceCurve
### Keywords: robust

### ** Examples

InfluenceCurve()

## The function is currently defined as
InfluenceCurve <- function(name, Curve = EuclRandVarList(EuclRandVariable(Domain = Reals())), 
                           Risks, Infos){
    if(missing(name))
        name <- "influence curve"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                     dimnames=list(character(0), c("method", "message")))
    
    return(new("InfluenceCurve", name = name, Curve = Curve, 
               Risks = Risks, Infos = Infos))
}



cleanEx()
nameEx("MEstimate-class")
### * MEstimate-class

flush(stderr()); flush(stdout())

### Name: MEstimate-class
### Title: MEstimate-class.
### Aliases: MEstimate-class Mroot Mroot,MEstimate-method
###   show,MEstimate-method
### Keywords: classes

### ** Examples

## prototype
new("MEstimate")



cleanEx()
nameEx("RobAStBaseMASK")
### * RobAStBaseMASK

flush(stderr()); flush(stdout())

### Name: RobAStBaseMASK
### Title: Masking of/by other functions in package "RobAStBase"
### Aliases: RobAStBaseMASK MASKING
### Keywords: programming distribution documentation

### ** Examples

RobAStBaseMASK()



cleanEx()
nameEx("RobAStBaseOptions")
### * RobAStBaseOptions

flush(stderr()); flush(stdout())

### Name: RobAStBaseOptions
### Title: Function to change the global variables of the package
###   'RobAStBase'
### Aliases: RobAStBaseOptions getRobAStBaseOption kStepUseLast
### Keywords: misc robust

### ** Examples

RobAStBaseOptions()
RobAStBaseOptions("kStepUseLast")
RobAStBaseOptions("kStepUseLast" = TRUE)
# or
RobAStBaseOptions(kStepUseLast = 1e-6)
getRobAStBaseOption("kStepUseLast")



cleanEx()
nameEx("RobWeight-class")
### * RobWeight-class

flush(stderr()); flush(stdout())

### Name: RobWeight-class
### Title: Robust Weight classes
### Aliases: RobWeight-class name,RobWeight-method name<-,RobWeight-method
###   weight weight,RobWeight-method weight<- weight<--methods
###   weight<-,RobWeight-method
### Keywords: classes

### ** Examples

## prototype
new("RobWeight")



cleanEx()
nameEx("TotalVarIC-class")
### * TotalVarIC-class

flush(stderr()); flush(stdout())

### Name: TotalVarIC-class
### Title: Influence curve of total variation type
### Aliases: TotalVarIC-class CallL2Fam<-,TotalVarIC-method clipLo
###   clip,TotalVarIC-method clipLo,TotalVarIC-method clipLo<-
###   clipLo<-,TotalVarIC-method clipUp clipUp,TotalVarIC-method clipUp<-
###   clipUp<-,TotalVarIC-method lowerCase<-,TotalVarIC-method
###   neighbor,TotalVarIC-method show,TotalVarIC-method
###   stand<-,TotalVarIC-method
###   generateIC,TotalVarNeighborhood,L2ParamFamily-method
### Keywords: classes robust

### ** Examples

IC1 <- new("TotalVarIC")
plot(IC1)



cleanEx()
nameEx("TotalVarIC")
### * TotalVarIC

flush(stderr()); flush(stdout())

### Name: TotalVarIC
### Title: Generating function for TotalVarIC-class
### Aliases: TotalVarIC
### Keywords: robust

### ** Examples

IC1 <- TotalVarIC()
plot(IC1)



cleanEx()
nameEx("TotalVarNeighborhood-class")
### * TotalVarNeighborhood-class

flush(stderr()); flush(stdout())

### Name: TotalVarNeighborhood-class
### Title: Total variation neighborhood
### Aliases: TotalVarNeighborhood-class
### Keywords: classes models

### ** Examples

new("TotalVarNeighborhood")



cleanEx()
nameEx("TotalVarNeighborhood")
### * TotalVarNeighborhood

flush(stderr()); flush(stdout())

### Name: TotalVarNeighborhood
### Title: Generating function for TotalVarNeighborhood-class
### Aliases: TotalVarNeighborhood
### Keywords: models

### ** Examples

TotalVarNeighborhood()

## The function is currently defined as
function(radius = 0){ 
    new("TotalVarNeighborhood", radius = radius) 
}



cleanEx()
nameEx("checkIC")
### * checkIC

flush(stderr()); flush(stdout())

### Name: checkIC
### Title: Generic Function for Checking ICs
### Aliases: checkIC
### Keywords: robust

### ** Examples

IC1 <- new("IC")
checkIC(IC1)



cleanEx()
nameEx("comparePlot")
### * comparePlot

flush(stderr()); flush(stdout())

### Name: comparePlot-methods
### Title: Compare - Plots
### Aliases: comparePlot comparePlot-methods comparePlot,IC,IC-method
### Keywords: robust

### ** Examples

if(require(ROptEst)){

N0 <- NormLocationScaleFamily(mean=0, sd=1) 
N0.Rob1 <- InfRobModel(center = N0, neighbor = ContNeighborhood(radius = 0.5))

IC1 <- optIC(model = N0, risk = asCov())
IC2 <- optIC(model = N0.Rob1, risk = asMSE())

comparePlot(IC1,IC2)

set.seed(12); data <- r(N0)(20)
comparePlot(IC1, IC2, data=data, with.lab = TRUE,
            which.lbs = c(1:4,15:20),
            which.Order = 1:6,
            return.Order = TRUE)

## don't run to reduce check time on CRAN
## Not run: 
##D ## selection of subpanels for plotting
##D par(mfrow=c(1,1))
##D comparePlot(IC1, IC2 ,mfColRow = FALSE, to.draw.arg=c("mean"),
##D             panel.first= grid(),ylim=c(-4,4),xlim=c(-6,6))
##D ## matrix-valued ylim
##D comparePlot(IC1, IC2, panel.first= grid(),ylim=c(-4,4,0,4),xlim=c(-6,6))
##D 
##D x <- c(data,-12,10)
##D comparePlot(IC1, IC2, data=x, which.Order=10,
##D             panel.first= grid(), ylim=c(-4,4,0,4), xlim=c(-6,6))
##D 
##D Y <- Chisq(df=1)* DiscreteDistribution(c(-1,1))
##D comparePlot(IC1, IC2, data=x, which.Order=10,
##D             scaleX = TRUE, scaleX.fct=pnorm, scaleX.inv=qnorm,
##D             scaleY = TRUE, scaleY.fct=p(Y), scaleY.inv=q(Y),
##D             panel.first= grid(), ylim=c(-4,4,0,4), xlim=c(-6,6))
##D comparePlot(IC1, IC2, data=x, which.Order=10,
##D             scaleX = TRUE, scaleX.fct=pnorm, scaleX.inv=qnorm,
##D             scaleY = TRUE, scaleY.fct=p(Y), scaleY.inv=q(Y),
##D             x.ticks = c(-Inf, -10, -1,0,1,10,Inf),
##D             y.ticks = c(-Inf, -5, -1,0,1,5,Inf),
##D             panel.first= grid(), ylim=c(-4,4,0,4), xlim=c(-6,6))
##D 
##D ## with use of trafo-matrix:
##D G <- GammaFamily(scale = 1, shape = 2)
##D ## explicitely transforming to
##D ## MASS parametrization:
##D mtrafo <- function(x){
##D      nms0 <- names(c(main(param(G)),nuisance(param(G))))
##D      nms <- c("shape","rate")
##D      fval0 <- c(x[2], 1/x[1])
##D      names(fval0) <- nms
##D      mat0 <- matrix( c(0, -1/x[1]^2, 1, 0), nrow = 2, ncol = 2,
##D                      dimnames = list(nms,nms0))                          
##D      list(fval = fval0, mat = mat0)}
##D G2 <- G
##D trafo(G2) <- mtrafo
##D G2
##D G2.Rob1 <- InfRobModel(center = G2, neighbor = ContNeighborhood(radius = 0.5))
##D system.time(IC1 <- optIC(model = G2, risk = asCov()))
##D system.time(IC2 <- optIC(model = G2.Rob1, risk = asMSE()))
##D system.time(IC2.i <- optIC(model = G2.Rob1, risk = asMSE(normtype=InfoNorm())))
##D system.time(IC2.s <- optIC(model = G2.Rob1, risk = asMSE(normtype=SelfNorm())))
##D 
##D comparePlot(IC1,IC2, IC2.i, IC2.s)
## End(Not run)

}



cleanEx()
nameEx("cutoff-class")
### * cutoff-class

flush(stderr()); flush(stdout())

### Name: cutoff-class
### Title: Cutoff class for distance-distance plots
### Aliases: cutoff-class cutoff.quantile<-,cutoff-method cutoff.quantile<-
###   cutoff.quantile,cutoff-method cutoff.quantile name,cutoff-method
###   fct,cutoff-method
### Keywords: classes

### ** Examples

cutoff()



cleanEx()
nameEx("cutoff")
### * cutoff

flush(stderr()); flush(stdout())

### Name: cutoff
### Title: Generating function(s) for class 'cutoff'
### Aliases: cutoff cutoff.sememp cutoff.chisq
### Keywords: hplot

### ** Examples

cutoff()
cutoff.sememp()
cutoff.chisq()



cleanEx()
nameEx("ddPlot-methods")
### * ddPlot-methods

flush(stderr()); flush(stdout())

### Name: ddPlot-methods
### Title: Methods for Function ddPlot in Package 'RobAStBase'
### Aliases: ddPlot ddPlot-methods ddPlot,matrix-method
###   ddPlot,numeric-method ddPlot,data.frame-method
### Keywords: methods hplot

### ** Examples

MX <- matrix(rnorm(1500),nrow=6)
QM <- matrix(rnorm(36),nrow=6); QM <- QM %*% t(QM)
ddPlot(data=MX, dist.y=QFNorm(QuadF=PosSemDefSymmMatrix(QM)))



cleanEx()
nameEx("getRiskFctBV-methods")
### * getRiskFctBV-methods

flush(stderr()); flush(stdout())

### Name: getRiskFctBV-methods
### Title: Methods for Function getRiskFctBV in Package 'RobAStBase'
### Aliases: getRiskFctBV getRiskFctBV-methods
###   getRiskFctBV,asGRisk,ANY-method getRiskFctBV,asMSE,ANY-method
###   getRiskFctBV,asSemivar,onesidedBias-method
###   getRiskFctBV,asSemivar,asymmetricBias-method
### Keywords: classes

### ** Examples

myrisk <- asMSE()
getRiskFctBV(myrisk)



cleanEx()
nameEx("infoPlot")
### * infoPlot

flush(stderr()); flush(stdout())

### Name: infoPlot
### Title: Plot absolute and relative information
### Aliases: infoPlot infoPlot-methods infoPlot,IC-method
### Keywords: robust

### ** Examples

N <- NormLocationScaleFamily(mean=0, sd=1) 
IC1 <- optIC(model = N, risk = asCov())
infoPlot(IC1)

## don't run to reduce check time on CRAN
## Not run: 
##D ## selection of subpanels for plotting
##D par(mfrow=c(1,2))
##D infoPlot(IC1, mfColRow = FALSE, to.draw.arg=c("Abs","sd"))
##D infoPlot(IC1, mfColRow = FALSE, to.draw.arg=c("Abs","sd"), log="y")
##D 
##D infoPlot(IC1, mfColRow = FALSE, to.draw.arg=c("Abs","mean"), 
##D               panel.first= grid(), ylim = c(0,4), xlim = c(-6,6))
##D infoPlot(IC1, mfColRow = FALSE, to.draw.arg=c("Abs","mean"), 
##D               panel.first= grid(), ylim = c(0,4,-3,3), xlim = c(-6,6))
##D 
##D par(mfrow=c(1,3))
##D infoPlot(IC1, mfColRow = FALSE, panel.first= grid(),
##D          ylim = c(0,4,0,.3,0,.8), xlim=c(-6,6))
##D par(mfrow=c(1,1))
##D 
##D data <- r(N)(20)
##D par(mfrow=c(1,3))
##D infoPlot(IC1, data=data, mfColRow = FALSE, panel.first= grid(),
##D          with.lab = TRUE, cex.pts=2,
##D          which.lbs = c(1:4,15:20), which.Order = 1:6,
##D          return.Order = TRUE)
##D infoPlot(IC1, data=data[1:10], mfColRow = FALSE, panel.first= grid(),
##D          with.lab = TRUE, cex.pts=0.7)
##D par(mfrow=c(1,1))
## End(Not run)




cleanEx()
nameEx("internals_ddPlot")
### * internals_ddPlot

flush(stderr()); flush(stdout())

### Name: internals_for_RobAStBase_ddPlot
### Title: Internal / Helper functions of package RobAStBase for ddPlot
### Aliases: internals_for_RobAStBase_ddPlot .ddPlot.MatNtNtCoCo
### Keywords: internal hplot

### ** Examples

MX <- matrix(rnorm(1500),nrow=6)
QM <- matrix(rnorm(36),nrow=6); QM <- QM %*% t(QM)
RobAStBase:::.ddPlot.MatNtNtCoCo(data=MX, 
        dist.y=QFNorm(QuadF=PosSemDefSymmMatrix(QM)),
        xlab="Norm.x",ylab="Norm.y", cex.idn = 1.3, offset=0,
        lwd=2, lwd.cutoff=4, lty=2, col.cutoff =2, col.idn="green",
        col = "blue", adj=0.4, pos=4,id.n = sample(1:200,size=100),
        lab.pts=letters,log="x", main="GA", sub="NO",cex.sub=0.2)



cleanEx()
nameEx("interpolRisk-class")
### * interpolRisk-class

flush(stderr()); flush(stdout())

### Name: interpolRisk-class
### Title: Interpolated Risks
### Aliases: interpolRisk-class OMSRRisk-class RMXRRisk-class MBRRisk-class
###   OMSRRisk RMXRRisk MBRRisk
### Keywords: classes

### ** Examples

new("OMSRRisk")
OMSRRisk()
RMXRRisk()
MBRRisk()
myrisk <- MBRRisk(samplesize=100)
samplesize(myrisk)
samplesize(myrisk) <- 20



cleanEx()
nameEx("kStepEstimator")
### * kStepEstimator

flush(stderr()); flush(stdout())

### Name: kStepEstimator
### Title: Function for the computation of k-step estimates
### Aliases: kStepEstimator
### Keywords: univar robust

### ** Examples

if(require(ROptEst)){
## 1. generate a contaminated sample
ind <- rbinom(100, size=1, prob=0.05)
x <- rnorm(100, mean=0, sd=(1-ind) + ind*9)

## 2. Kolmogorov(-Smirnov) minimum distance estimator
(est0 <- MDEstimator(x=x, NormLocationScaleFamily()))

## 3. k-step estimation: radius known
N1 <- NormLocationScaleFamily(mean=estimate(est0)["mean"], sd=estimate(est0)["sd"])
N1.Rob <- InfRobModel(center = N1, neighbor = ContNeighborhood(radius = 0.5))
IC1 <- optIC(model = N1.Rob, risk = asMSE())
(est1 <- kStepEstimator(x, IC1, est0, steps = 3, withPIC = TRUE))
estimate(est1)
ksteps(est1)
pICList(est1)
start(est1)

## don't run to reduce check time on CRAN
## Not run: 
##D ## a transformed model
##D tfct <- function(x){
##D     nms0 <- c("mean","sd")
##D     nms  <- "comb"
##D     fval0 <- x[1]+2*x[2]
##D     names(fval0) <- nms
##D     mat0 <- matrix(c(1,2), nrow = 1, dimnames = list(nms,nms0))
##D     return(list(fval = fval0, mat = mat0))
##D }
##D 
##D N1.traf <- N1; trafo(N1.traf) <- tfct
##D N1R.traf <- N1.Rob; trafo(N1R.traf) <- tfct
##D IC1.traf <- optIC(model = N1R.traf, risk = asMSE())
##D (est0.traf <- MDEstimator(x, N1.traf))
##D (est1.traf <- kStepEstimator(x, IC1.traf, est0, steps = 3,
##D                 withIC = TRUE, withPIC = TRUE, withUpdateInKer = FALSE))
##D (est1a.traf <- kStepEstimator(x, IC1.traf, est0, steps = 3,
##D                 withIC = TRUE, withPIC = TRUE, withUpdateInKer = TRUE))
##D estimate(est1.traf)
##D ksteps(est1.traf)
##D pICList(est1.traf)
##D startval(est1.traf)
##D 
##D untransformed.estimate(est1.traf)
##D uksteps(est1.traf)
##D ICList(est1.traf)
##D ustartval(est1.traf)
##D 
##D estimate(est1a.traf)
##D ksteps(est1a.traf)
##D pICList(est1a.traf)
##D startval(est1a.traf)
##D 
##D untransformed.estimate(est1a.traf)
##D uksteps(est1a.traf)
##D ICList(est1a.traf)
##D ustartval(est1a.traf)
## End(Not run)
}



cleanEx()
nameEx("makeIC-methods")
### * makeIC-methods

flush(stderr()); flush(stdout())

### Name: makeIC-methods
### Title: Generic Function for making ICs consistent at a possibly
###   different model
### Aliases: makeIC makeIC-methods makeIC,IC,missing-method
###   makeIC,IC,L2ParamFamily-method makeIC,list,L2ParamFamily-method
###   makeIC,function,L2ParamFamily-method
### Keywords: robust

### ** Examples

## default IC
IC1 <- new("IC")

## L2-differentiable parametric family
B <- BinomFamily(13, 0.3)

## check IC properties
checkIC(IC1, B)

## make IC
IC2 <- makeIC(IC1, B)

## check IC properties
checkIC(IC2)

## slot modifyIC is filled in case of IC2
IC3 <- modifyIC(IC2)(BinomFamily(13, 0.2), IC2)
checkIC(IC3)
## identical to
checkIC(IC3, BinomFamily(13, 0.2))

IC4 <- makeIC(sin, B)
checkIC(IC4)

(IC5 <- makeIC(list(function(x)x^3), B, name="a try"))
plot(IC5)
checkIC(IC5)

## don't run to reduce check time on CRAN
## Not run: 
##D N0 <- NormLocationScaleFamily()
##D IC6 <- makeIC(list(sin,cos),N0)
##D plot(IC6)
##D checkIC(IC6)
##D 
##D getRiskIC(IC6,risk=trAsCov())$trAsCov$value
##D getRiskIC(IC6,risk=asBias(),neighbor=ContNeighborhood())$asBias$value
## End(Not run)




cleanEx()
nameEx("optIC")
### * optIC

flush(stderr()); flush(stdout())

### Name: optIC
### Title: Generic function for the computation of optimally robust ICs
### Aliases: optIC optIC-methods optIC,L2ParamFamily,asCov-method
### Keywords: robust

### ** Examples

B <- BinomFamily(size = 25, prob = 0.25) 

## classical optimal IC
IC0 <- optIC(model = B, risk = asCov())
plot(IC0) # plot IC
checkIC(IC0, B)



cleanEx()
nameEx("outlyingPlotIC")
### * outlyingPlotIC

flush(stderr()); flush(stdout())

### Name: outlyingPlotIC
### Title: Function outlyingPlotIC in Package 'RobAStBase'
### Aliases: outlyingPlotIC
### Keywords: hplot

### ** Examples

if(require(ROptEst)){
## generates normal location and scale family with mean = -2 and sd = 3
N0 <- NormLocationScaleFamily()
N0.IC0 <- optIC(model = N0, risk = asCov())
N0.Rob1 <- InfRobModel(center = N0, neighbor = ContNeighborhood(radius = 0.5))
N0.IC1 <- optIC(model = N0.Rob1, risk = asMSE())
xn <- c(rnorm(100),rcauchy(20)+20)
outlyingPlotIC(xn, IC.x=N0.IC0)
outlyingPlotIC(xn, IC.x=N0.IC1)
}



cleanEx()
nameEx("plot-methods")
### * plot-methods

flush(stderr()); flush(stdout())

### Name: plot-methods
### Title: Methods for Function plot in Package 'RobAStBase'
### Aliases: plot plot-methods plot,IC,missing-method
###   plot,IC,numeric-method
### Keywords: methods distribution

### ** Examples

IC1 <- new("IC")
plot(IC1)
plot(IC1, main = TRUE, panel.first= grid(),
     col = "blue", cex.main = 2, cex.inner = 1)

### selection of subpanels for plotting
N <- NormLocationScaleFamily(mean=0, sd=1) 
IC2 <- optIC(model = N, risk = asCov())
par(mfrow=c(1,1))
plot(IC2, main = TRUE, panel.first= grid(),
     col = "blue", cex.main = 2, cex.inner = 0.6,
     mfColRow = FALSE, to.draw.arg=c("sd"))

## xlim and ylim arguments
plot(IC2, main = TRUE, panel.first= grid(), 
     ylim=c(-3,3), xlim=c(-2,3))
plot(IC2, main = TRUE, panel.first= grid(), 
     ylim=c(-3,3,-1,3), xlim=c(-2,3),
     with.legend = TRUE)

data <- r(N)(30)
plot(IC2, data, panel.first= grid(),
     ylim = c(-3,3,-1,3), xlim=c(-2,3),
     cex.pts = 3, pch.pts = 1:2, col.pts="green",
     with.lab = TRUE, which.lbs = c(1:4,15:20),
     which.Order = 1:6, return.Order = TRUE)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("qqplot")
### * qqplot

flush(stderr()); flush(stdout())

### Name: qqplot
### Title: Methods for Function qqplot in Package 'RobAStBase'
### Aliases: qqplot qqplot-methods qqplot,ANY,RobModel-method
###   qqplot,ANY,InfRobModel-method qqplot,ANY,kStepEstimate-method
### Keywords: hplot distribution

### ** Examples

qqplot(r(Norm(15,sqrt(30)))(40), Chisq(df=15))
RobM <- InfRobModel(center = NormLocationFamily(mean=13,sd=sqrt(28)),
                    neighbor = ContNeighborhood(radius = 0.4))
x <- r(Norm(15,sqrt(30)))(20)
qqplot(x, RobM)
qqplot(x, RobM, alpha.CI=0.9)
## further examples for ANY,kStepEstimator-method
## in example to roptest() in package ROptEst



cleanEx()
nameEx("samplesize-methods")
### * samplesize-methods

flush(stderr()); flush(stdout())

### Name: samplesize-methods
### Title: Methods for Function samplesize in Package 'RobAStBase'
### Aliases: samplesize,interpolRisk-method samplesize
###   samplesize<-,interpolRisk-method samplesize<-
### Keywords: classes

### ** Examples

myrisk <- MBRRisk(samplesize=100)
samplesize(myrisk)
samplesize(myrisk) <- 20



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
