###############################################################################
## Example: Normal location and scale
###############################################################################
require(ROptEst)
options("newDevice"=TRUE)

## generates normal location and scale family with mean = -2 and sd = 3
N0 <- NormLocationScaleFamily(mean=-2, sd=3) 
N0        # show G0
plot(N0)  # plot of Norm(mean = -2, sd = 3) and L_2 derivative
checkL2deriv(N0)

## classical optimal IC
N0.IC0 <- optIC(model = N0, risk = asCov())
N0.IC0       # show IC
checkIC(N0.IC0)
Risks(N0.IC0)
plot(N0.IC0) # plot IC

## L_2 family + infinitesimal neighborhood
N0.Rob1 <- InfRobModel(center = N0, neighbor = ContNeighborhood(radius = 0.5))
N0.Rob1     # show N0.Rob1

## MSE solution
system.time(N0.IC1 <- optIC(model = N0.Rob1, risk = asMSE()))
checkIC(N0.IC1)
Risks(N0.IC1)
plot(N0.IC1)
infoPlot(N0.IC1)

system.time(N0.IC1.i <- optIC(model = N0.Rob1, risk = asMSE(normtype=InfoNorm())))
checkIC(N0.IC1.i)
Risks(N0.IC1.i)
plot(N0.IC1.i)
infoPlot(N0.IC1.i)

system.time(N0.IC1.s <- optIC(model = N0.Rob1, risk = asMSE(normtype=SelfNorm())))
checkIC(N0.IC1.s)
Risks(N0.IC1.s)
plot(N0.IC1.s)
infoPlot(N0.IC1.s)

comparePlot(N0.IC1,N0.IC1.i,N0.IC1.s)


## lower case solutions
(N0.IC2 <- optIC(model = N0.Rob1, risk = asBias(), tol = 1e-10))
checkIC(N0.IC2)
Risks(N0.IC2)
plot(N0.IC2)
infoPlot(N0.IC2)

(N0.IC2.i <- optIC(model = N0.Rob1, risk = asBias(normtype=InfoNorm()), tol = 1e-10))
checkIC(N0.IC2.i)
Risks(N0.IC2.i)
plot(N0.IC2.i)
infoPlot(N0.IC2.i)

## Hampel solution
(N0.IC3 <- optIC(model = N0.Rob1, risk = asHampel(bound = clip(N0.IC1))))
checkIC(N0.IC3)
Risks(N0.IC3)
plot(N0.IC3) 
infoPlot(N0.IC3)

## radius minimax IC
## (may take quite some time!)
system.time(N0.IC4 <- radiusMinimaxIC(L2Fam=N0, neighbor=ContNeighborhood(), 
                risk=asMSE(), loRad=0, upRad=Inf))
checkIC(N0.IC4)
Risks(N0.IC4)
plot(N0.IC4) 
infoPlot(N0.IC4)

(N0.IC4.i <- radiusMinimaxIC(L2Fam=N0, neighbor=ContNeighborhood(), 
                risk=asMSE(normtype=InfoNorm()), loRad=0, upRad=Inf))
checkIC(N0.IC4.i)
Risks(N0.IC4.i)
plot(N0.IC4.i) 
infoPlot(N0.IC4.i)

## takes extremely long time:
(N0.IC4.s <- radiusMinimaxIC(L2Fam=N0, neighbor=ContNeighborhood(), 
                risk=asMSE(normtype=SelfNorm()), loRad=0, upRad=Inf))
checkIC(N0.IC4.s)
Risks(N0.IC4.s)
plot(N0.IC4.s) 
infoPlot(N0.IC4.s)


## least favorable radius
## (may take quite some time!)
N0.r.rho1 <- leastFavorableRadius(L2Fam=N0, neighbor=ContNeighborhood(),
                    risk=asMSE(), rho=0.5)

## one-step estimation
## 1. generate a contaminated sample
ind <- rbinom(100, size=1, prob=0.05) 
x <- rnorm(100, mean=0, sd=(1-ind) + ind*9)

## 2. Kolmogorov(-Smirnov) minimum distance estimator
(est0 <- MDEstimator(x=x, NormLocationScaleFamily()))

## 3. k-step estimation: radius known
N1 <- NormLocationScaleFamily(mean=estimate(est0)["mean"], sd=estimate(est0)["sd"])
N1.Rob <- InfRobModel(center = N1, neighbor = ContNeighborhood(radius = 0.5))
IC1 <- optIC(model = N1.Rob, risk = asMSE())
(est1 <- kStepEstimator(x, IC1, est0, steps = 3))

## 4. k-step estimation: radius unknown
## rough estimate: 1-10% contamination
## => r\in[0.1,1.0]

## takes some time
IC2 <- radiusMinimaxIC(L2Fam=N1, neighbor=ContNeighborhood(),risk=asMSE(), 
                       loRad=0.1, upRad=1.0) 
(est2 <- oneStepEstimator(x, IC2, est0))


## a simple example
library(MASS)
data(chem)
initial.est <- c(median(chem), mad(chem))
system.time(ROest1 <- roptest(chem, NormLocationScaleFamily(), eps.upper = 0.1, steps = 3L, 
                              initial.est = initial.est, useLast = TRUE))

## optimized for speed
library(RobLox)
system.time(ROest2 <- roblox(chem, eps.upper = 0.1, k = 3, returnIC = TRUE))

## comparison
estimate(ROest1)
estimate(ROest2)

## confidence intervals
confint(ROest1, symmetricBias())
confint(ROest2, symmetricBias())
