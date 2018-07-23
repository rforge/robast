#######################################################################
# example of robust fits at real data not included in Rd due to timings
#######################################################################

require(RobExtremes)
require(ismev)
data(portpirie)
data(rain)
detach(package:ismev)
raini <- rain[rain>10]
rainc <- c(raini,1000,10000)
portpiriei <- portpirie[,2]
portpiriec <- c(portpiriei,100)

##--------------------------
## GEV
##--------------------------

ppfiti <- ismev::gev.fit(portpiriei)
gev.diag(ppfiti)
  ##
mlEi <- MLEstimator(portpiriei, GEVFamilyMuUnknown(withPos=FALSE))
system.time(MBRi <- roptest(portpiriei, GEVFamilyMuUnknown(withPos=FALSE),risk=MBRRisk()))
system.time(RMXi <- roptest(portpiriei, GEVFamilyMuUnknown(withPos=FALSE),risk=RMXRRisk()))
## in fact the precision of the pIC is not too good, but the resp. estimate only differs
## little to the situation where we enforce IC conditions
checkIC(pIC(RMXi))
system.time(RMXiw <- roptest(portpiriei, GEVFamilyMuUnknown(withPos=FALSE),risk=RMXRRisk(),withMakeIC=TRUE))
checkIC(pIC(RMXiw))
estimate(RMXi)
estimate(RMXiw)

mlEi
MBRi
estimate(mlEi)
estimate(MBRi)
estimate(RMXi)
estimate(RMXiw)
attr(MBRi, "timings")
gev.diag(mlEi)
gev.diag(MBRi)
gev.prof(mlEi, m = 10, 4.1, 5)
gev.profxi(MBRi, -0.3, 0.3)
plot(pIC(MBRi))

## contaminated:
ppfitc <- ismev::gev.fit(portpiriec)
gev.diag(ppfitc)
  ##
mlEc <- MLEstimator(portpiriec, GEVFamilyMuUnknown(withPos=FALSE))
system.time(MBRc <- roptest(portpiriec, GEVFamilyMuUnknown(withPos=FALSE),risk=MBRRisk()))
system.time(RMXc <- roptest(portpiriec, GEVFamilyMuUnknown(withPos=FALSE),risk=RMXRRisk()))
mlEc
MBRc
estimate(mlEi)
estimate(MBRi)
estimate(RMXi)
estimate(mlEc)
estimate(MBRc)
estimate(RMXc)
gev.diag(mlEc)
gev.diag(MBRc)
gev.prof(mlEc, m = 10, 4.1, 5)
gev.profxi(mlEc, -0.3, 0.3)
qqplot(portpiriec,MBRc)
qqplot(portpiriec,MBRc,ylim=c(3.5,5))
returnlevelplot(portpiriec,MBRc)
returnlevelplot(portpiriec,MBRc,ylim=c(3.5,5))

plot(pIC(MBRc))

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
mlE2c
MBR2c
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
qqplot(rainc,MBR2c,xlim=c(5,100),ylim=c(5,100),log="xy",cex.pch=8,exp.fadcol.pch = .55,pch=19)

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
dI2i <- distribution(eval(MBR2i@pIC@CallL2Fam))
loc(dI2i) <- 0
## wrong data set
qqplot(portpiriei-10,dI2i)
qqplot(portpiriec,MBR2c)


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
