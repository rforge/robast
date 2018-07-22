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
mlEi
MBRi
estimate(mlEi)
estimate(MBRi)
attr(MBRi, "timings")
gev.diag(mlEi)
gev.diag(MBRi)
gev.prof(mlEi, m = 10, 4.1, 5)
gev.profxi(MBRi, -0.3, 0.3)
plot(MBRi@pIC)

## contaminated:
ppfitc <- ismev::gev.fit(portpiriec)
gev.diag(ppfitc)
  ##
mlEc <- MLEstimator(portpiriec, GEVFamilyMuUnknown(withPos=FALSE))
system.time(MBRc <- roptest(portpiriec, GEVFamilyMuUnknown(withPos=FALSE),risk=MBRRisk()))
mlEc
MBRc
estimate(mlEc)
estimate(MBRc)
estimate(mlEi)
estimate(MBRi)
gev.diag(mlEc)
gev.diag(MBRc)
gev.prof(mlEc, m = 10, 4.1, 5)
gev.profxi(mlEc, -0.3, 0.3)
qqplot(portpiriec,MBRc)
qqplot(portpiriec,MBRc,ylim=c(3.5,5))
returnlevelplot(portpiriec,MBRc)
returnlevelplot(portpiriec,MBRc,ylim=c(3.5,5))

plot(MBRc@pIC)

##--------------------------
## GPD
##--------------------------

rnfiti <- ismev::gpd.fit(rain,10)
gpd.diag(rnfiti)
  ##
mlE2i <- MLEstimator(raini, GParetoFamily(loc=10))
gpd.diag(mlE2i)
system.time(MBR2i <- roptest(raini, GParetoFamily(loc=10),risk=MBRRisk()))
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
mlE2c
MBR2c
estimate(mlE2c)
estimate(MBR2c)
estimate(mlE2i)
estimate(MBR2i)
gpd.diag(mlE2c)
gpd.diag(MBR2c)
gpd.prof(mlE2c, m = 10, 55, 77)
gpd.profxi(mlE2c, -0.02, 0.02)
plot(MBR2c@pIC)

qqplot(rainc,MBR2c)
qqplot(rainc,MBR2c,ylim=c(5,100))
qqplot(rainc,MBR2c,xlim=c(5,100),ylim=c(5,100),log="xy")

returnlevelplot(raini,MBR2i,MaxOrPot="POT",threshold=0)
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
