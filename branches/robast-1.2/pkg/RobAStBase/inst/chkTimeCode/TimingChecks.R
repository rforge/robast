###
# to check timings carefully do:
#
##  remove all comment tags ##-t-## from kStepEstimator.R
##   -> this costs about > 1sec per step in kStep
##      but reveals where time is spent (better than too detailed R-prof
#
###
###
## this file is used for testing only
require(ROptEst)
ind <- rbinom(100, size=1, prob=0.05)
x <- rnorm(100, mean=0, sd=(1-ind) + ind*9)
## 2. Kolmogorov(-Smirnov) minimum distance estimator
(est0 <- MDEstimator(x=x, NormLocationScaleFamily()))
## 3. k-step estimation: radius known
N1 <- NormLocationScaleFamily(mean=estimate(est0)["mean"], sd=estimate(est0)["sd"])
N1.Rob <- InfRobModel(center = N1, neighbor = ContNeighborhood(radius = 0.5))
IC1 <- optIC(model = N1.Rob, risk = asMSE())
syest1 <- system.time(est1 <- kStepEstimator(x, IC1, est0, steps = 3, withPIC = TRUE))
attr(est1,"timings")
## a transformed model
tfct <- function(x){
    nms0 <- c("mean","sd")
    nms  <- "comb"
    fval0 <- x[1]+2*x[2]
    names(fval0) <- nms
    mat0 <- matrix(c(1,2), nrow = 1, dimnames = list(nms,nms0))
    return(list(fval = fval0, mat = mat0))
}
N1.traf <- N1; trafo(N1.traf) <- tfct
N1R.traf <- N1.Rob; trafo(N1R.traf) <- tfct
IC1.traf <- optIC(model = N1R.traf, risk = asMSE())
syest0.traf <- system.time(est0.traf <- MDEstimator(x, N1.traf))
syest1.traf <- system.time(est1.traf <- kStepEstimator(x, IC1.traf, est0, steps = 3,
                withIC = TRUE, withPIC = TRUE, withUpdateInKer = FALSE))
syest1a.traf <- system.time(est1a.traf <- kStepEstimator(x, IC1.traf, est0, steps = 3,
                withIC = TRUE, withPIC = TRUE, withUpdateInKer = TRUE))
syest0.traf
syest1.traf
colSums(attr(est1.traf,"timings"))
syest1a.traf
colSums(attr(est1a.traf,"timings"))
syest1.0 <- system.time(est1.0 <- kStepEstimator(x, IC1, est0, steps = 3, withPIC = TRUE))
syest1
colSums(attr(est1,"timings"))
syest1.0
colSums(attr(est1.0,"timings"))

syest1
syest1.0
syest0.traf
syest1.traf
syest1a.traf
attr(est1,"timings")
attr(est1.0,"timings")
attr(est1.traf,"timings")
attr(est1a.traf,"timings")
set.seed(123)
ind <- rbinom(100, size=1, prob=0.05)
x <- rbinom(100, size=25, prob=(1-ind)*0.25 + ind*0.9)
## ML-estimate
sy.MLE <- system.time(MLE.bin <- MLEstimator(x, BinomFamily(size = 25)))
## compute optimally robust estimators
sy.OMSE <- system.time(OMSE.bin <- OMSEstimator(x, BinomFamily(size = 25), steps = 3))
sy.MLE
sy.OMSE
attr(OMSE.bin,"kStepTimings")
colSums(attr(OMSE.bin,"kStepTimings"))

require(RobExtremes)
require(ismev)
require(fitdistrplus) ## for dataset groundbeef
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

GEVFam <- GEVFamilyMuUnknown(withPos=FALSE)
mlEi <- MLEstimator(portpiriei, GEVFam)
## MLE as asy lin estimator with IC optimal as to asCov() - risk:
syMLEiALE <- system.time(mlEiALE <- roptest(portpiriei, GEVFam,risk=asCov()))
## with 10 steps and with
syMLEi10ALE <- system.time(mlEi10ALE <- roptest(portpiriei, GEVFam,risk=asCov(),steps=10))
syMBRi <- system.time(MBRi <- MBREstimator(portpiriei, GEVFam))
## synonymous to
## system.time(MBRi0 <- roptest(portpiriei, GEVFam,risk=MBRRisk()))
syRMXi <- system.time(RMXi <- RMXEstimator(portpiriei, GEVFam))
## synonymous to
## system.time(RMXi <- roptest(portpiriei, GEVFam,risk=RMXRRisk()))
## in fact the precision of the pIC is not too good, but the resp. estimate only differs
## little to the situation where we enforce IC conditions
syRMXiw <- system.time(RMXiw <- RMXEstimator(portpiriei, GEVFam,withMakeIC=TRUE))

syMLEiALE
attr(mlEiALE,"kStepTimings")
syMLEi10ALE
attr(mlEi10ALE,"kStepTimings")
syMBRi
attr(MBRi,"kStepTimings")
syRMXi
attr(RMXi,"kStepTimings")
syRMXiw
attr(RMXiw,"kStepTimings")
###
 IC <- RMXi@pIC@Curve
 ICw <- RMXiw@pIC@Curve
 L2F <- eval(CallL2Fam(RMXi@pIC))
 IC1 <- as(diag(3) %*% IC, "EuclRandVariable")
 IC1w <- as(diag(3) %*% ICw, "EuclRandVariable")
 system.time(E(L2F, IC1 %*% t(IC1)))
 system.time(E(L2F, IC1w %*% t(IC1w)))

require(ROptEst)
set.seed(123)
ind <- rbinom(100, size=1, prob=0.05)
x0 <- rpois(100,1)
x <- c(rpois(100, lambda=(1-ind)*1 + ind*100),-(1:4))
## ML-estimate
sy.MLE <- system.time(MLE.bin <- MLEstimator(x, PoisFamily()))
## compute optimally robust estimators
sy.OMSE <- system.time(OMSE.bin <- OMSEstimator(x, PoisFamily(), steps = 3))
sy.OMSE0 <- system.time(OMSE0.bin <- OMSEstimator(x0, PoisFamily(), steps = 3))

sapply(c(x,1e9),pIC(OMSE0.bin)@Curve[[1]]@Map[[1]])
sapply(c(x,1e9),pIC(OMSE.bin)@Curve[[1]]@Map[[1]])
system.time(print(CVM<-CvMMDEstimator(x,PoisFamily())))
system.time(print(CVM<-CvMMDEstimator(x,PoisFamily(),muDatOrMod="Dat")))
