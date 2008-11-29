###############################################################################
## Example: Binomial Family
###############################################################################
require(ROptEst)

## generates Binomial Family with 
## m = 25 and probability of success theta = 0.25
B <- BinomFamily(size = 25, prob = 0.25) 
B       # show B 
plot(B) # plot of Binom(size = 25, prob = 0.25) and L_2 derivative
checkL2deriv(B)

## classical optimal IC
IC0 <- optIC(model = B, risk = asCov())
IC0       # show IC
plot(IC0) # plot IC
checkIC(IC0)
Risks(IC0)

## lower case radius
lowerCaseRadius(L2Fam = B, neighbor = ContNeighborhood(), risk = asMSE())
lowerCaseRadius(L2Fam = B, neighbor = TotalVarNeighborhood(), risk = asMSE())

## L_2 family + infinitesimal neighborhood
RobB1 <- InfRobModel(center = B, neighbor = ContNeighborhood(radius = 0.5))
RobB1     # show RobB1
(RobB2 <- InfRobModel(center = B, neighbor = TotalVarNeighborhood(radius = 0.5)))

## MSE solution
system.time(IC1 <- optIC(model=RobB1, risk=asMSE()), gcFirst = TRUE)
IC1
checkIC(IC1)
Risks(IC1)
getRiskIC(IC1, asBias(), ContNeighborhood()) # standardized bias
getRiskIC(IC1, asMSE(), ContNeighborhood(radius = 0.5))

(Cov1 <- getRiskIC(IC1, asCov()))
(mse1 <- getRiskIC(IC1, asMSE(), TotalVarNeighborhood(radius = 0.5)))
(bias1 <- getRiskIC(IC1, asBias(), TotalVarNeighborhood()))
## only suboptimal -> ToDo-List
addRisk(IC1) <- list(Cov1, mse1, bias1)
Risks(IC1)
plot(IC1)

system.time(IC2 <- optIC(model=RobB2, risk=asMSE()), gcFirst = TRUE)
IC2
checkIC(IC2)
Risks(IC2)
getRiskIC(IC2, asMSE(), TotalVarNeighborhood(radius = 0.5))
getRiskIC(IC2, asBias(), TotalVarNeighborhood())
getRiskIC(IC2, asBias(), ContNeighborhood())
Cov2 <- getRiskIC(IC2, asCov())
addRisk(IC2) <- Cov2
Risks(IC2)
plot(IC2)

## lower case solutions
(IC3 <- optIC(model=RobB1, risk=asBias()))
checkIC(IC3)
Risks(IC3)
plot(IC3)

(IC4 <- optIC(model=RobB2, risk=asBias()))
checkIC(IC4)
Risks(IC4)
plot(IC4)


## Hampel solution
(IC5 <- optIC(model=RobB1, risk=asHampel(bound=clip(IC1))))
checkIC(IC5)
Risks(IC5)
plot(IC5)

(IC6 <- optIC(model=RobB2, risk=asHampel(bound=Risks(IC2)$asBias), maxiter = 200))
checkIC(IC6)
Risks(IC6)
plot(IC6)


## radius minimax IC
system.time(IC7 <- radiusMinimaxIC(L2Fam=B, neighbor=ContNeighborhood(), 
                        risk=asMSE(), loRad=0, upRad=1), gcFirst = TRUE)
IC7
checkIC(IC7)
Risks(IC7)
plot(IC7)

system.time(IC8 <- radiusMinimaxIC(L2Fam=B, neighbor=TotalVarNeighborhood(), 
                        risk=asMSE(), loRad=0, upRad=1), gcFirst = TRUE)
IC8
checkIC(IC8)
Risks(IC8)
plot(IC8)


## least favorable radius
## (may take quite some time!)
system.time(r.rho1 <- leastFavorableRadius(L2Fam=B, neighbor=ContNeighborhood(),
                    risk=asMSE(), rho=0.5), gcFirst = TRUE)
r.rho1
system.time(r.rho2 <- leastFavorableRadius(L2Fam=B, neighbor=TotalVarNeighborhood(),
                    risk=asMSE(), rho=0.5), gcFirst = TRUE)
r.rho2

## one-step estimation
## 1. generate a contaminated sample
ind <- rbinom(100, size=1, prob=0.05) 
x <- rbinom(100, size=25, prob=(1-ind)*0.25 + ind*0.75)

## 2. Kolmogorov(-Smirnov) minimum distance estimator
(est0 <- ksEstimator(x=x, Binom(size=25), param = "prob"))

## 3. one-step estimation: radius known
RobB3 <- InfRobModel(center=BinomFamily(size=25, prob=est0$prob), 
                neighbor=ContNeighborhood(radius=0.5))
IC9 <- optIC(model=RobB3, risk=asMSE())
(est1 <- oneStepEstimator(x, IC=IC9, start=est0$prob))

RobB4 <- InfRobModel(center=BinomFamily(size=25, prob=est0$prob), 
                neighbor=TotalVarNeighborhood(radius=0.25))
IC10 <- optIC(model=RobB4, risk=asMSE())
(est1 <- oneStepEstimator(x, IC=IC10, start=est0$prob))

## 4. one-step estimation: radius interval
IC11 <- radiusMinimaxIC(L2Fam=BinomFamily(size=25, prob=est0$prob),
                neighbor=ContNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)
(est2 <- oneStepEstimator(x, IC=IC11, start=est0$prob))

IC12 <- radiusMinimaxIC(L2Fam=BinomFamily(size=25, prob=est0$prob),
                neighbor=TotalVarNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)
(est2 <- oneStepEstimator(x, IC=IC12, start=est0$prob))
