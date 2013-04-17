####################################################
### Tests fuer InterpolRisiken
####################################################
PFam <- NULL
mytest <- function(PF = GParetoFamily, xi = 0.5, seed=130313, beta=1){
### arguments
## PF: generating function of the family
## xi: shape parameter
## seed: seed fuer den Zufallszahlengenerator
   PFam <- PF(shape=xi,scale=beta)
   cat("\n\n\n---------------------------------\n")
   cat("  ", name(PFam),"  ")
   cat("\n---------------------------------\n")
   set.seed(seed)
   dat0 <- r(PFam)(100)
   cat("\n\n\n---------------------------------\n")
   cat("RMXE")
   cat("\n---------------------------------\n")
   try({
   print(system.time(
   {re1<-roptest(dat0,PFam,risk=RMXRRisk())}
   ))
   print(re1)
   print(checkIC(pIC(re1)))
   },silent=TRUE)
   cat("\n\n\n---------------------------------\n")
   cat("OMSE")
   cat("\n---------------------------------\n")
   try({
   system.time(re2<-roptest(dat0,PFam,risk=OMSRRisk()))
   print(re2)
   print(checkIC(pIC(re2)))
   },silent=TRUE)
   cat("\n\n\n---------------------------------\n")
   cat("MBRE")
   cat("\n---------------------------------\n")
   try({
   system.time(re3<-roptest(dat0,PFam,risk=MBRRisk()))
   print(re3)
   print(checkIC(pIC(re3)))
   },silent=TRUE)
}
mytest(PF=GParetoFamily)
mytest(PF=GEVFamily)
mytest(PF=GammaFamily)
mytest(PF=WeibullFamily)
mytest(PF=GParetoFamily,xi=1)
mytest(PF=GEVFamily,xi=1)
mytest(PF=GammaFamily,xi=1)
mytest(PF=WeibullFamily,xi=1)
mytest(PF=GParetoFamily,xi=0.1)
mytest(PF=GEVFamily,xi=0.1)
mytest(PF=GammaFamily,xi=0.1)
mytest(PF=WeibullFamily,xi=0.1)
mytest(PF=GParetoFamily,xi=10)
mytest(PF=GEVFamily,xi=10)
mytest(PF=GammaFamily,xi=10)
mytest(PF=WeibullFamily,xi=10)
