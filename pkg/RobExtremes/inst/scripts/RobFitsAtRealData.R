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
system.time(MBRi <- roptest(portpiriei, GEVFam,risk=MBRRisk()))
system.time(RMXi <- roptest(portpiriei, GEVFam,risk=RMXRRisk()))
## in fact the precision of the pIC is not too good, but the resp. estimate only differs
## little to the situation where we enforce IC conditions
checkIC(pIC(RMXi))
system.time(RMXiw <- roptest(portpiriei, GEVFam,risk=RMXRRisk(),withMakeIC=TRUE))
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
gev.diag(mlEi)
gev.diag(MBRi)
gev.prof(mlEi, m = 10, 4.1, 5)
gev.profxi(MBRi, -0.3, 0.3)

## this is how the MLEi IC looks like
plot(pIC(mlEi10ALE),portpiriei)
## with additional arguments
plot(pIC(mlEi10ALE),portpiriei,
     cex.lbs=1.8,with.lab=TRUE, which.Order=1:3, c
     ol.pts="red", adj.lbs=c(1.3,1.3), pch.pts=20,pch.npts=21,
     cex.pts=2,cex.npts=.6,lab.pts=c("H","G","I"))
## this is how the MBRi IC looks like
plot(pIC(MBRi))
## rescale y axes according to different probit scales
plot(pIC(mlEi10ALE),portpiriei, with.lab=TRUE,
     which.lbs=c(1:10), attr.pre=FALSE, col.pts="blue", col.npts="green",
     scaleY = TRUE, scaleY.fct=list(function(x)pnorm(x,sd=4),
     function(x)pnorm(x,sd=10),function(x)pnorm(x,sd=30)),
     scaleY.inv=list(function(x)qnorm(x,sd=4),
     function(x)qnorm(x,sd=10),function(x)qnorm(x,sd=30)),scaleN=20)
## contaminated:
ppfitc <- ismev::gev.fit(portpiriec)
gev.diag(ppfitc)
  ##
mlEc <- MLEstimator(portpiriec, GEVFamilyMuUnknown(withPos=FALSE))
system.time(MBRc <- roptest(portpiriec, GEVFamilyMuUnknown(withPos=FALSE),risk=MBRRisk()))
system.time(RMXc <- roptest(portpiriec, GEVFamilyMuUnknown(withPos=FALSE),risk=RMXRRisk()))
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
gev.diag(mlEc)
gev.diag(MBRc)
gev.prof(mlEc, m = 10, 4.1, 5)
gev.profxi(mlEc, -0.3, 0.3)

## diagnostics from pkg 'distrMod'/'RobAStBase'
qqplot(portpiriec,MBRc)
qqplot(portpiriec,MBRc,ylim=c(3.5,5))
returnlevelplot(portpiriec,MBRc)
returnlevelplot(portpiriec,MBRc,ylim=c(3.5,5))

## here the MBR-IC looks as follows
plot(pIC(MBRc),portpiriec,log="x",xlim=c(1,1.5*max(portpiriec)))

##--------------------------
## GPD
##--------------------------

rnfiti <- ismev::gpd.fit(rain,10)
gpd.diag(rnfiti)
  ##
mlE2i <- MLEstimator(raini, GParetoFamily(loc=10))
gpd.diag(mlE2i)
system.time(MBR2i <- roptest(raini, GParetoFamily(loc=10),risk=MBRRisk()))
system.time(RMX2i <- roptest(raini, GParetoFamily(loc=10),risk=RMXRRisk()))
mlE2i
MBR2i
estimate(mlE2i)
estimate(MBR2i)
plot(MBR2i@pIC)
gpd.diag(mlE2i)
gpd.diag(MBR2i)
gpd.prof(mlE2i, m = 10, 55, 77)
gpd.profxi(MBR2i, -0.02, 0.02)
GP <- GParetoFamily(scale=rnfiti$mle[1],shape=rnfiti$mle[2],loc=10)
returnlevelplot(rain, GP, MaxOrPOT="POT", xlim=c(1e-1,1e3))

## contaminated:

rnfitc <- ismev::gpd.fit(c(rain,1000,10000),10)
gpd.diag(rnfitc)
  ##
mlE2c <- MLEstimator(rainc, GParetoFamily(loc=10))
gpd.diag(mlE2c)
system.time(MBR2c <- roptest(rainc, GParetoFamily(loc=10),risk=MBRRisk()))
system.time(RMX2c <- roptest(rainc, GParetoFamily(loc=10),risk=RMXRRisk()))

## again a comparison, and again MLE is shuttered, the robust ones keep cool
estimate(mlE2i)
estimate(MBR2i)
estimate(RMX2i)
estimate(mlE2c)
estimate(MBR2c)
estimate(RMX2c)
gpd.diag(mlE2c)
gpd.diag(MBR2c)
gpd.diag(RMX2c)
gpd.prof(mlE2c, m = 10, 55, 77)
gpd.profxi(mlE2c, -0.02, 0.02)
plot(pIC(MBR2c))

qqplot(rainc,MBR2c)
qqplot(rainc,MBR2c,ylim=c(5,100))
qqplot(rainc,MBR2c,xlim=c(5,100),ylim=c(5,100),log="xy")
qqplot(rainc,di,xlim=c(5,100),ylim=c(5,100),log="xy",
       cex.pts=2,col.pts="blue",with.lab=TRUE,cex.lbs=.9,which.Order=1:3)

returnlevelplot(raini,MBR2i,MaxOrPot="POT",threshold=0)
returnlevelplot(raini,MBR2i,MaxOrPot="POT",threshold=0, withLab=TRUE, cex.lbl=0.8)
returnlevelplot(rainc,MBR2c,MaxOrPot="POT",threshold=0)
returnlevelplot(rainc,MBR2c,ylim=c(10,100),MaxOrPot="POT",threshold=0)
#
L2F <- eval(MBR2c@pIC@CallL2Fam)
dI2c <- L2F@distribution
#loc(dI2c) <- 0
qqplot(rainc,dI2c)
rainc.10 <- rainc-10
qqplot(rainc.10,dI2c-10)
## to be fixed
returnlevelplot(rainc.10,dI2c-10,MaxOrPot="POT",threshold=0)
## wrong data set
dI2i <- distribution(eval(MBR2i@pIC@CallL2Fam))
loc(dI2i) <- 0
qqplot(portpiriei-10,dI2i)
qqplot(portpiriec,MBR2c)
## right data set
qqplot(raini-10,dI2i)
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
qqplot(x, mlE3i, log="xy")
qqplot(xc, mlE3c, log="xy")
returnlevelplot(x, mlE3i, MaxOrPOT="POT",ylim=c(1,1e5),log="y")

system.time(MBR3i <- roptest(x, PM,risk=MBRRisk()))
system.time(RMX3i <- roptest(x, PM,risk=RMXRRisk()))
system.time(MBR3c <- roptest(xc, PM,risk=MBRRisk()))
system.time(RMX3c <- roptest(xc, PM,risk=RMXRRisk()))
estimate(mlE3i)
estimate(MBR3i)
estimate(RMX3i)
estimate(mlE3c)
estimate(MBR3c)
estimate(RMX3c)
plot(pIC(MBR3i))
plot(pIC(RMX3i))

#######################################################
# Weibull data
#######################################################
WF <- WeibullFamily()
system.time(mlE4i <- MLEstimator(grbsi, WF))
system.time(MBR4i <- roptest(grbsi, WF,risk=MBRRisk()))
system.time(OMS4i <- roptest(grbsi, WF,risk=OMSRRisk()))
system.time(RMX4i <- roptest(grbsi, WF,risk=RMXRRisk()))
system.time(mlE4c <- MLEstimator(grbsc, WF))
system.time(MBR4c <- roptest(grbsc, WF,risk=MBRRisk()))
system.time(OMS4c <- roptest(grbsc, WF,risk=OMSRRisk()))
system.time(RMX4c <- roptest(grbsc, WF,risk=RMXRRisk()))
estimate(mlE4i)
estimate(MBR4i)
estimate(RMX4i)
estimate(OMS4i)
estimate(mlE4c)
estimate(MBR4c)
estimate(OMS4c)
estimate(RMX4c)
plot(pIC(MBR4i))
plot(pIC(RMX4i))
qqplot(grbsi, RMX4i)
qqplot(grbsc, RMX4c, log="xy")

#######################################################
# Gamma data
#######################################################

GF <- GammaFamily()
system.time(mlE5i <- MLEstimator(grbsi, GF))
system.time(OMS5i <- roptest(grbsi, GF,risk=OMSRRisk()))
system.time(RMX5i <- roptest(grbsi, GF,risk=RMXRRisk()))
system.time(mlE5c <- MLEstimator(grbsc, GF))
system.time(OMS5c <- roptest(grbsc, GF,risk=OMSRRisk()))
system.time(RMX5c <- roptest(grbsc, GF,risk=RMXRRisk()))
estimate(mlE5i)
estimate(RMX5i)
estimate(OMS5i)
estimate(mlE5c)
estimate(OMS5c)
estimate(RMX5c)
plot(pIC(OMS5i))
plot(pIC(RMX5i))
qqplot(grbsi, RMX5i)
qqplot(grbsc, RMX5c, log="xy")
