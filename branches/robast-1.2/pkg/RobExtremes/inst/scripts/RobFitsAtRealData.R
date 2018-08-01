#######################################################################
# example of robust fits at real data not included in Rd due to timings
#######################################################################

require(RobExtremes)
require(ismev)
require(fitdistrplus) ## for dataset groundbeef


help(package="RobExtremes")
help("RobExtremes-package")

#----------------------------------------
## data sets
data(groundbeef)
data(portpirie)
data(rain)
detach(package:ismev)
detach(package:fitdistrplus)
#----------------------------------------
## create contaminated datasets
raini <- rain[rain>10]
rainc <- c(raini,1000,10000)
portpiriei <- portpirie[,2]
portpiriec <- c(portpiriei,100)
grbsi <- groundbeef$serving
grbsc <- c(grbsi,10000,1000)


options("newDevice"=TRUE) ## opens a new plot for each graphic
##--------------------------
## GEV
##--------------------------

## compare the fit from ismev:
ppfiti <- ismev::gev.fit(portpiriei)
gev.diag(ppfiti)
  ##
GEVFam <- GEVFamilyMuUnknown(withPos=FALSE)

mlEi <- MLEstimator(portpiriei, GEVFam)
## MLE as asy lin estimator with IC optimal as to asCov() - risk:
system.time(mlEiALE <- roptest(portpiriei, GEVFam,risk=asCov()))
## with 10 steps and with
system.time(mlEi10ALE <- roptest(portpiriei, GEVFam,risk=asCov(),steps=10))
system.time(MBRi <- MBREstimator(portpiriei, GEVFam))
## synonymous to
## system.time(MBRi0 <- roptest(portpiriei, GEVFam,risk=MBRRisk()))
system.time(RMXi <- RMXEstimator(portpiriei, GEVFam))
## synonymous to
## system.time(RMXi <- roptest(portpiriei, GEVFam,risk=RMXRRisk()))
## in fact the precision of the pIC is not too good, but the resp. estimate only differs
## little to the situation where we enforce IC conditions
checkIC(pIC(RMXi))
system.time(RMXiw <- RMXEstimator(portpiriei, GEVFam,withMakeIC=TRUE))
checkIC(pIC(RMXiw))
estimate(RMXi)
estimate(RMXiw)

## our output:
mlEi
mlEiALE
MBRi
### comparison of the estimates
ppfiti$mle
estimate(mlEi)
estimate(mlEiALE) ### somewhat different after 1st step
estimate(mlEi10ALE) ### similar after 10 steps
estimate(MBRi)
estimate(RMXi)
estimate(RMXiw)
### where do the robust estimators spend their time?
attr(MBRi, "timings")

## our return values can be plugged into ismev-diagnostics:
devNew()
gev.diag(mlEi)
devNew()
gev.diag(MBRi)
devNew()
gev.prof(mlEi, m = 10, 4.1, 5)
devNew()
gev.profxi(MBRi, -0.3, 0.3)

## this is how the MLEi IC looks like
devNew()
plot(pIC(mlEi10ALE),portpiriei)
## with additional arguments
devNew()
plot(pIC(mlEi10ALE),portpiriei,
     cex.lbs=1.8,with.lab=TRUE, which.Order=1:3, col.pts="red",
     adj.lbs=c(1.3,1.3), pch.pts=20,pch.npts=21,
     cex.pts=2,cex.npts=.6,lab.pts=c("H","G","I"))
## this is how the MBRi IC looks like
devNew()
plot(pIC(MBRi),portpiriei,
     cex.lbs=1.8,with.lab=TRUE, which.Order=1:3, col.pts="red",
     adj.lbs=c(1.3,1.3), pch.pts=20,pch.npts=21,
     cex.pts=2,cex.npts=.6,lab.pts=c("H","G","I"))
## compare the 4 ICs
devNew()
comparePlot(obj1=pIC(mlEi10ALE),obj2=pIC(MBRi),obj3=pIC(RMXi),obj4=pIC(mlEiALE),
     portpiriei, with.legend=TRUE,# lwd=1, lty=1:4,
     cex.lbs=1.8,with.lab=TRUE, which.Order=1:3, col.pts="red",
     adj.lbs=c(1.3,1.3), pch.pts=20,pch.npts=21,
     cex.pts=2,cex.npts=.6,lab.pts=c("H","G","I"))

devNew()
comparePlot(obj1=pIC(mlEi10ALE),obj2=pIC(MBRi),obj3=pIC(RMXi),
     with.legend=TRUE, ylim=matrix(c(-1.5,1.5,-1,3,-5,5),2,3))

devNew()
infoPlot(pIC(MBRi), portpiriei, with.legend=TRUE, lty=1:4, log="y",
     cex.lbs=1.8,with.lab=TRUE, which.Order=1:3, col.pts="red",
     adj.lbs=c(1.3,1.3), pch.pts=20,pch.npts=21,
     cex.pts=2,cex.npts=.6,lab.pts=c("H","G","I"))

devNew()
infoPlot(pIC(RMXi), portpiriei, with.legend=TRUE, lty=1:4, log="y",
     cex.lbs=1.8,with.lab=TRUE, which.Order=1:3, col.pts="red",
     adj.lbs=c(1.3,1.3), pch.pts=20,pch.npts=21,
     cex.pts=2,cex.npts=.6,lab.pts=c("H","G","I"))

## rescale y axes according to different probit scales
devNew()
plot(pIC(mlEi10ALE),portpiriei, with.lab=TRUE,
     which.lbs=c(1:10), attr.pre=FALSE, col.pts="blue", col.npts="green",
     scaleY = TRUE, scaleY.fct=list(function(x)pnorm(x,sd=4),
     function(x)pnorm(x,sd=10),function(x)pnorm(x,sd=30)),
     scaleY.inv=list(function(x)qnorm(x,sd=4),
     function(x)qnorm(x,sd=10),function(x)qnorm(x,sd=30)),scaleN=20)
## contaminated:
ppfitc <- ismev::gev.fit(portpiriec)
devNew()
gev.diag(ppfitc)
  ##
mlEc <- MLEstimator(portpiriec, GEVFamilyMuUnknown(withPos=FALSE))
system.time(MBRc <- MBREstimator(portpiriec, GEVFamilyMuUnknown(withPos=FALSE)))
system.time(RMXc <- RMXEstimator(portpiriec, GEVFamilyMuUnknown(withPos=FALSE)))
## our output:
mlEc
MBRc

## compare the estimates
# ideal situation -- only minimal differences
ppfiti$mle
estimate(mlEi)
estimate(MBRi)
estimate(RMXi)
# contaminated situation -- the mle reacts sharply, the robust est's don't
ppfitc$mle
estimate(mlEc)
estimate(MBRc)
estimate(RMXc)

## diagnostics from ismev applied to our estimators:
devNew()
gev.diag(mlEc)
devNew()
gev.diag(MBRc)
devNew()
gev.prof(mlEc, m = 10, 4.1, 5)
devNew()
gev.profxi(mlEc, -0.3, 0.3)

## diagnostics from pkg 'distrMod'/'RobAStBase'
devNew()
qqplot(portpiriec,MBRc)
devNew()
qqplot(portpiriec,MBRc,ylim=c(3.5,5))
devNew()
returnlevelplot(portpiriec,MBRc)
devNew()
returnlevelplot(portpiriec,MBRc,ylim=c(3.5,5))

## here the MBR-IC looks as follows
devNew()
plot(pIC(MBRc),portpiriec,log="x",xlim=c(1,1.5*max(portpiriec)))

##--------------------------
## GPD
##--------------------------

rnfiti <- ismev::gpd.fit(rain,10)
gpd.diag(rnfiti)
  ##
mlE2i <- MLEstimator(raini, GParetoFamily(loc=10))
gpd.diag(mlE2i)
system.time(MBR2i <- MBREstimator(raini, GParetoFamily(loc=10)))
system.time(RMX2i <- RMXEstimator(raini, GParetoFamily(loc=10)))
mlE2i
MBR2i
estimate(mlE2i)
estimate(MBR2i)
devNew()
plot(MBR2i@pIC)
devNew()
gpd.diag(mlE2i)
devNew()
gpd.diag(MBR2i)
devNew()
gpd.prof(mlE2i, m = 10, 55, 77)
devNew()
gpd.profxi(MBR2i, -0.02, 0.02)
GP <- GParetoFamily(scale=rnfiti$mle[1],shape=rnfiti$mle[2],loc=10)
returnlevelplot(rain, GP, MaxOrPOT="POT", xlim=c(1e-1,1e3))

## contaminated:

rnfitc <- ismev::gpd.fit(c(rain,1000,10000),10)
devNew()
gpd.diag(rnfitc)
  ##
mlE2c <- MLEstimator(rainc, GParetoFamily(loc=10))
devNew()
gpd.diag(mlE2c)
system.time(MBR2c <- MBREstimator(rainc, GParetoFamily(loc=10)))
system.time(RMX2c <- RMXEstimator(rainc, GParetoFamily(loc=10)))

## again a comparison, and again MLE is shuttered, the robust ones keep cool
estimate(mlE2i)
estimate(MBR2i)
estimate(RMX2i)
estimate(mlE2c)
estimate(MBR2c)
estimate(RMX2c)
devNew()
gpd.diag(mlE2c)
devNew()
gpd.diag(MBR2c)
devNew()
gpd.diag(RMX2c)
devNew()
gpd.prof(mlE2c, m = 10, 55, 77)
devNew()
gpd.profxi(mlE2c, -0.02, 0.02)
devNew()
plot(pIC(MBR2c))

devNew()
qqplot(rainc,MBR2c)
devNew()
qqplot(rainc,MBR2c,ylim=c(5,100))
devNew()
qqplot(rainc,MBR2c,xlim=c(5,100),ylim=c(5,100),log="xy")
devNew()
qqplot(rainc,MBR2c,xlim=c(5,100),ylim=c(5,100),log="xy",
       cex.pts=2,col.pts="blue",with.lab=TRUE,cex.lbs=.9,which.Order=1:3)

devNew()
returnlevelplot(raini,MBR2i,MaxOrPot="POT",threshold=0)
devNew()
returnlevelplot(raini,MBR2i,MaxOrPot="POT",threshold=0, withLab=TRUE, cex.lbl=0.8)
devNew()
returnlevelplot(rainc,MBR2c,MaxOrPot="POT",threshold=0)
devNew()
returnlevelplot(rainc,MBR2c,ylim=c(10,100),MaxOrPot="POT",threshold=0)
#
L2F <- eval(MBR2c@pIC@CallL2Fam)
dI2c <- L2F@distribution
devNew()
qqplot(rainc,dI2c)
rainc.10 <- rainc-10
devNew()
qqplot(rainc.10,dI2c-10)
devNew()
returnlevelplot(rainc.10,dI2c-10,MaxOrPot="POT",threshold=0)

## wrong data set
dI2i <- distribution(eval(MBR2i@pIC@CallL2Fam))
loc(dI2i) <- 0
devNew()
qqplot(portpiriei-10,dI2i)
devNew()
qqplot(portpiriec,MBR2c)
### all points are red

## right data set
devNew()
qqplot(raini-10,dI2i)
devNew()
qqplot(rainc,MBR2c)


#######################################################
# synthetic Pareto data
#######################################################
set.seed(20180723)
x <- actuar::rpareto1(100, min=2, shape=0.4)
xc <- c(x, 1e15,1e12,1e69, 2.001,2.00000001)
PM <- ParetoFamily(Min=2)
mlE3i <- MLEstimator(x,PM)
mlE3c <- MLEstimator(xc,PM)
devNew()
qqplot(x, mlE3i, log="xy")
devNew()
qqplot(xc, mlE3c, log="xy")
devNew()
returnlevelplot(x, mlE3i, MaxOrPOT="POT",ylim=c(1,1e5),log="y")

system.time(MBR3i <- MBREstimator(x, PM))
system.time(RMX3i <- RMXEstimator(x, PM))
system.time(MBR3c <- MBREstimator(xc, PM))
system.time(RMX3c <- RMXEstimator(xc, PM))
estimate(mlE3i)
estimate(MBR3i)
estimate(RMX3i)
estimate(mlE3c)
estimate(MBR3c)
estimate(RMX3c)
devNew()
plot(pIC(MBR3i))
devNew()
plot(pIC(RMX3i))

#######################################################
# Weibull data
#######################################################
WF <- WeibullFamily()
system.time(mlE4i <- MLEstimator(grbsi, WF))
system.time(MBR4i <- MBREstimator(grbsi, WF))
system.time(OMS4i <- OMSEstimator(grbsi, WF))
## synonymous to
## system.time(OMS4i <- roptest(grbsi, WF, risk= OMSRRisk()))
system.time(RMX4i <- RMXEstimator(grbsi, WF))
system.time(mlE4c <- MLEstimator(grbsc, WF))
system.time(MBR4c <- MBREstimator(grbsc, WF))
system.time(OMS4c <- OMSEstimator(grbsc, WF))
system.time(RMX4c <- RMXEstimator(grbsc, WF))
estimate(mlE4i)
estimate(MBR4i)
estimate(RMX4i)
estimate(OMS4i)
estimate(mlE4c)
estimate(MBR4c)
estimate(OMS4c)
estimate(RMX4c)
devNew()
plot(pIC(MBR4i))
devNew()
plot(pIC(RMX4i))
devNew()
qqplot(grbsi, RMX4i)
devNew()
qqplot(grbsc, RMX4c, log="xy")

#######################################################
# Gamma data
#######################################################

GF <- GammaFamily()
system.time(mlE5i <- MLEstimator(grbsi, GF))
system.time(OMS5i <- MBREstimator(grbsi, GF))
system.time(RMX5i <- OMSEstimator(grbsi, GF))
system.time(MBR5i <- RMXEstimator(grbsi, GF))
system.time(mlE5c <- MLEstimator(grbsc, GF))
system.time(OMS5c <- MBREstimator(grbsc, GF))
system.time(RMX5c <- OMSEstimator(grbsc, GF))
system.time(MBR5c <- RMXEstimator(grbsc, GF))
estimate(mlE5i)
estimate(RMX5i)
estimate(OMS5i)
estimate(MBR5i)
estimate(mlE5c)
estimate(OMS5c)
estimate(RMX5c)
estimate(MBR5c)
devNew()
plot(pIC(OMS5i))
devNew()
plot(pIC(RMX5i))
devNew()
plot(pIC(MBR5i))
devNew()
qqplot(grbsi, RMX5i)
devNew()
qqplot(grbsc, RMX5c, log="xy")
