### getLMs for Pareto
### produced Jul 20, 2018 P.Ruckdeschel

## note:: the only param to be estimated in the smooth model framework
##        is the shape k
##        now, independent of fixed parameter Min
#         L2deriv in this case is distributed as (1-Exp(1))/k
##        so shape k acts as a scaling parameter, and we can compute
##        LMs under reference value k = 1, Min = 1,
##        the respective for general k OptICs then are simply  k \times OptIC(k=1)

.checkIC <- function(IC,L2deriv,Distr){
   IC. <- as(diag(dimension(IC@Curve)) %*% IC@Curve, "EuclRandVariable")
   IC.. <- function(x) sapply(x,function(y) c(IC.@Map[[1]](q.l(Distr)(y))))
   IC2. <- IC. %*% t(L2deriv)
   IC2.. <- function(x) sapply(x,function(y) c(IC2.@Map[[1]](q.l(Distr)(y))))
   IV <- integrate(IC..,0,1)$value
   IV2 <- integrate(IC2..,0,1)$value - 1
   return(c(cent=IV,consist=IV2))
}
..getLM <- function(IC){
   return(c(b=clip(IC),a=cent(IC),aw=cent(weight(IC)),A=stand(IC),Aw=stand(weight(IC))))
}

PF <- ParetoFamily(shape=1)
## necessary Jul 20, 2018 due to bug:
PF@L2derivDistr[[1]] <- 1-Exp(1)

L2deriv. <- as(diag(1) %*% PF@L2deriv, "EuclRandVariable")
D1. <- PF@distribution

ICRMX <- radiusMinimaxIC(L2Fam=PF, neighbor= ContNeighborhood(), risk = asMSE(), verbose = TRUE, loRad = 0, upRad = Inf, z.start = 0, A.start = 2, upper = 1e7, lower = 1e-7, OptOrIter = "iterate", maxiter = 150, tol = .Machine$double.eps^0.7, loRad0 = 1e-3)
.checkIC(ICRMX,L2deriv.,D1.)
(RMXw <- ..getLM(ICRMX))

RobPF <- InfRobModel(center = PF, neighbor = ContNeighborhood(radius = 0.5))
ICMBR <- optIC(model = RobPF, risk = asBias(), verbose = TRUE, z.start = 0, A.start = 2, upper = 1e7, lower = 1e-7, OptOrIter = "iterate", maxiter = 150, tol = .Machine$double.eps^0.7)
.checkIC(ICMBR,L2deriv.,D1.)
(MBRw <- ..getLM(ICMBR))

asM <- asMSE()
ICOMS <- optIC(model = RobPF, risk = asM, verbose = TRUE, z.start = 0, A.start = 2, upper = 1e7, lower = 1e-7, OptOrIter = "iterate", maxiter = 150, tol = .Machine$double.eps^0.7)
.checkIC(ICRMX,L2deriv.,D1.)
(OMSw <- ..getLM(ICOMS))

.ParetoLM <- list()
.ParetoLM$RMX <- RMXw
.ParetoLM$OMS <- OMSw
.ParetoLM$MBR <- MBRw

baseP <- "C:/rtest/RobASt/branches/robast-1.1/pkg"
bufferP <- "RobExtremesBuffer"
RobExP <- "RobExtremes/R"
unzF <- file.path(baseP,bufferP,"Paretosysdata.rda")
zF <-   file.path(baseP,bufferP,"ParetoZipsysdata.rda")
rdaF <- file.path(baseP,RobExP,"sysdata.rda")

save(.ParetoLM, file=unzF)
save(.ParetoLM, file=zF, compress="xz")
save(.ParetoLM, file=rdaF, compress="xz")
