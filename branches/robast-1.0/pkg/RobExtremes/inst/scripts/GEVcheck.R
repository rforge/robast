require(RobExtremes)
set.seed(123)

x <- rgev(100,shape=2,scale=30,loc=4) ### x still causes problems....
x1 <- rgev(100,shape=-.2,scale=30,loc=4)
gev0 <- GEVFamilyMuUnknown(withPos=FALSE)
gev1 <- GEVFamily(loc=quantile(x,exp(-1)),withPos=FALSE)
MLEstimator(x,gev0)
MLEstimator(x,gev1)
MLEstimator(x1,gev0)
MLEstimator(x1,gev1)
fucheck <- function(x,mu,sigma,xi){
   1+min((x-mu)/sigma*xi)
}
(est1 <- RobExtremes:::.getMuBetaXiGEV(x,withPos=FALSE))
fucheck(x,mu=est1[1],sigma=est1[2],xi=est1[3])
(est2 <- RobExtremes:::.getBetaXiGEV(x,mu=quantile(x,exp(-1)),withPos=FALSE))
fucheck(x,mu=quantile(x,exp(-1)),sigma=est2[1],xi=est2[2])

(est3 <- RobExtremes:::.getMuBetaXiGEV(x1,withPos=FALSE))
fucheck(x1,mu=est3[1],sigma=est3[2],xi=est3[3])
(est4 <- RobExtremes:::.getBetaXiGEV(x1,mu=quantile(x,exp(-1)),withPos=FALSE))
fucheck(x1,mu=quantile(x,exp(-1)),sigma=est4[1],xi=est4[2])
