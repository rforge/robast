###############################################################################
## Example: Gamma Family
###############################################################################
require(ROptEst)

## generates Gamma Family with 
## scale = 2 and shape = 0.5
G <- GammaFamily(scale = 2, shape = 0.5)
G       # show G
plot(G) # plot of Gammad(scale = 2, shape = 0.5) and L_2 derivative
distrExOptions(ErelativeTolerance = 1e-8) # increase precision for E
checkL2deriv(G)

## classical optimal IC
IC0 <- optIC(model = G, risk = asCov())
IC0       # show IC
checkIC(IC0) # seems to be a numerical problem!?
Risks(IC0)
plot(IC0) # plot IC

## L_2 family + infinitesimal neighborhood
RobG1 <- InfRobModel(center = G, neighbor = ContNeighborhood(radius = 0.5))
RobG1     # show RobB1

## MSE solution
system.time(IC1 <- optIC(model=RobG1, risk=asMSE()))
IC1
checkIC(IC1)
Risks(IC1)
plot(IC1)
devNew()
infoPlot(IC1)

## lower case solutions
system.time(IC2 <- optIC(model=RobG1, risk=asBias(), tol = 1e-10))
IC2
checkIC(IC2)
Risks(IC2)
plot(IC2)
devNew()
infoPlot(IC2)

## Hampel solution
system.time(IC3 <- optIC(model=RobG1, risk=asHampel(bound=clip(IC1))))
IC3
checkIC(IC3)
Risks(IC3)
plot(IC3)
x11()
infoPlot(IC3)

## radius minimax IC
## takes quite some time - about 30 min.
system.time(IC4 <- radiusMinimaxIC(L2Fam=G, neighbor=ContNeighborhood(), 
            risk=asMSE(), loRad=0, upRad=Inf))

## least favorable radius
## takes really long time - several hours!
#system.time(r.rho1 <- leastFavorableRadius(L2Fam=G, neighbor=ContNeighborhood(),
#                    risk=asMSE(), rho=0.5))

## one-step estimation
## 1. generate a contaminated sample
ind <- rbinom(100, size=1, prob=0.05) 
x <- (1-ind)*rgamma(100, scale = 1, shape = 2) + ind*rgamma(100, scale = 3, shape = 5)

## 2. Kolmogorov(-Smirnov) minimum distance estimator
(est0 <- MDEstimator(x=x, GammaFamily()))

## non-robust ML estimator
MLEstimator(x=x, GammaFamily())

## 3. one-step estimation: radius known
RobG3 <- InfRobModel(center=GammaFamily(scale = estimate(est0)["scale"], 
                     shape = estimate(est0)["shape"]), 
                neighbor=ContNeighborhood(radius=0.5))
IC9 <- optIC(model=RobG3, risk=asMSE())
(est1 <- oneStepEstimator(x, IC=IC9, start=estimate(est0)))
