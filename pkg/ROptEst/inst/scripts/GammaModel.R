###############################################################################
## Example: Gamma Family
###############################################################################
require(ROptEst)

## generates Gamma Family with 
## scale = 1 and shape = 1
G <- GammaFamily(scale = 1, shape = 2)
G       # show G
plot(G) # plot of Gammad(scale = 1, shape = 2) and L_2 derivative
checkL2deriv(G)

## classical optimal IC
IC0 <- optIC(model = G, risk = asCov())
IC0       # show IC
system.time(checkIC(IC0), gcFirst = TRUE)
Risks(IC0)
plot(IC0) # plot IC

## L_2 family + infinitesimal neighborhood
RobG1 <- InfRobModel(center = G, neighbor = ContNeighborhood(radius = 0.5))
RobG1     # show RobB1

## MSE solution
system.time(IC1 <- optIC(model=RobG1, risk=asMSE()), gcFirst = TRUE)
IC1
checkIC(IC1)
Risks(IC1)
plot(IC1)
x11()
infoPlot(IC1)

## lower case solutions
system.time(IC2 <- optIC(model=RobG1, risk=asBias(), tol = 1e-10), gcFirst = TRUE)
IC2
checkIC(IC2)
Risks(IC2)
plot(IC2)
x11()
infoPlot(IC2)

## Hampel solution
system.time(IC3 <- optIC(model=RobG1, risk=asHampel(bound=clip(IC1))), gcFirst = TRUE)
IC3
checkIC(IC3)
Risks(IC3)
plot(IC3)
x11()
infoPlot(IC3)

## radius minimax IC
## takes quite some time - about 30 min.
system.time(IC4 <- radiusMinimaxIC(L2Fam=G, neighbor=ContNeighborhood(), 
                risk=asMSE(), loRad=0, upRad=Inf), gcFirst = TRUE)

## least favorable radius
## takes quite some time - several hours!
#system.time(r.rho1 <- leastFavorableRadius(L2Fam=G, neighbor=ContNeighborhood(),
#                    risk=asMSE(), rho=0.5))

## one-step estimation
## 1. generate a contaminated sample
ind <- rbinom(100, size=1, prob=0.05) 
x <- (1-ind)*rgamma(100, scale = 1, shape = 2) + ind*10

## 2. Kolmogorov(-Smirnov) minimum distance estimator
(est0 <- ksEstimator(x=x, Gammad()))

## 3. one-step estimation: radius known
RobG3 <- InfRobModel(center=GammaFamily(scale = est0$scale, shape = est0$shape), 
                neighbor=ContNeighborhood(radius=0.5))
IC9 <- optIC(model=RobG3, risk=asMSE())
(est1 <- oneStepEstimator(x, IC=IC9, start=est0))
