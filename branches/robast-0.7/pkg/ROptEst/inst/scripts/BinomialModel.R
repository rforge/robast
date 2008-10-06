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
system.time(IC1 <- optIC(model=RobB1, risk=asMSE()))
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

system.time(IC2 <- optIC(model=RobB2, risk=asMSE()))
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

(IC6 <- optIC(model=RobB2, risk=asHampel(bound=Risks(IC2)$asBias$value), maxiter = 200))
checkIC(IC6)
Risks(IC6)
plot(IC6)


## radius minimax IC
system.time(IC7 <- radiusMinimaxIC(L2Fam=B, neighbor=ContNeighborhood(), 
                        risk=asMSE(), loRad=0, upRad=1))
IC7
checkIC(IC7)
Risks(IC7)
plot(IC7)

system.time(IC8 <- radiusMinimaxIC(L2Fam=B, neighbor=TotalVarNeighborhood(), 
                        risk=asMSE(), loRad=0, upRad=1))
IC8
checkIC(IC8)
Risks(IC8)
plot(IC8)


## least favorable radius
## (may take quite some time!)
system.time(r.rho1 <- leastFavorableRadius(L2Fam=B, neighbor=ContNeighborhood(),
                    risk=asMSE(), rho=0.5))
r.rho1
system.time(r.rho2 <- leastFavorableRadius(L2Fam=B, neighbor=TotalVarNeighborhood(),
                    risk=asMSE(), rho=0.5))
r.rho2


###############################################################################
## k-step (k >= 1) estimation
###############################################################################

## one-step estimation
## 1. generate a contaminated sample
ind <- rbinom(100, size=1, prob=0.05) 
x <- rbinom(100, size=25, prob=(1-ind)*0.25 + ind*0.75)

## 2. Kolmogorov(-Smirnov) minimum distance estimator
(est0 <- MDEstimator(x=x, BinomFamily(size=25)))

## 3.1. one-step estimation: radius known
## ac) Define infinitesimal robust model
RobB3 <- InfRobModel(center=BinomFamily(size=25, prob=estimate(est0)),
                     neighbor=ContNeighborhood(radius=0.5))
## bc) Compute optimally robust IC
IC9 <- optIC(model=RobB3, risk=asMSE())
checkIC(IC9)
## cc) Determine 1-step estimate
(est1c <- oneStepEstimator(x, IC=IC9, start=est0))

## instead of ac)-cc) you can also use function roptest
est1c1 <- roptest(x, BinomFamily(size = 25), eps = 0.05, initial.est = est0)
checkIC(pIC(est1c1))
## you can also omit step 2
est1c2 <- roptest(x, BinomFamily(size = 25), eps = 0.05, distance = KolmogorovDist)
checkIC(pIC(est1c2))

## Using Cramer-von-Mises MD estimator (default)
est1c3 <- roptest(x, BinomFamily(size = 25), eps = 0.025)
checkIC(pIC(est1c3))

## comparison of estimates
estimate(est1c)
estimate(est1c1)
estimate(est1c2)
estimate(est1c3)


## av) Define infinitesimal robust model
RobB4 <- InfRobModel(center=BinomFamily(size=25, prob=estimate(est0)),
                neighbor=TotalVarNeighborhood(radius=0.25))
## bv) Compute optimally robust IC
IC10 <- optIC(model=RobB4, risk=asMSE())
checkIC(IC10)
## cv) Determine 1-step estimate
(est1v <- oneStepEstimator(x, IC=IC10, start=est0))

## instead of av)-cv) you can also use function roptest
est1v1 <- roptest(x, BinomFamily(size = 25), eps = 0.025, initial.est = est0,
                  neighbor = TotalVarNeighborhood())
checkIC(pIC(est1v1))
## you can also omit step 2
est1v2 <- roptest(x, BinomFamily(size = 25), eps = 0.025,
                  neighbor = TotalVarNeighborhood(), distance = KolmogorovDist)
checkIC(pIC(est1v2))

## Using Cramer-von-Mises MD estimator (default)
est1v3 <- roptest(x, BinomFamily(size = 25), eps = 0.025, neighbor = TotalVarNeighborhood())
checkIC(pIC(est1v3))

## comparison of estimates
estimate(est1v)
estimate(est1v1)
estimate(est1v2)
estimate(est1v3)

## 3.2. k-step estimation: radius known
IC9 <- optIC(model=RobB3, risk=asMSE())
(est2c <- kStepEstimator(x, IC=IC9, start=est0, steps = 3L))

est2c1 <- roptest(x, BinomFamily(size = 25), eps = 0.05, initial.est = est0, steps = 3L)
checkIC(pIC(est2c1))

est2c2 <- roptest(x, BinomFamily(size = 25), eps = 0.05, steps = 3L,
                  distance = KolmogorovDist)
checkIC(pIC(est2c2))

## Using Cramer-von-Mises MD estimator
est2c3 <- roptest(x, BinomFamily(size = 25), eps = 0.05, steps = 3L)
checkIC(pIC(est2c3))

## comparison of estimates
estimate(est2c)
estimate(est2c1)
estimate(est2c2)
estimate(est2c3)


IC10 <- optIC(model=RobB4, risk=asMSE())
(est2v <- kStepEstimator(x, IC=IC10, start=est0, steps = 3L))
checkIC(pIC(est2v))

est2v1 <- roptest(x, BinomFamily(size = 25), eps = 0.025, initial.est = est0, 
                  steps = 3L, neighbor = TotalVarNeighborhood())
checkIC(pIC(est2v1))

est2v2 <- roptest(x, BinomFamily(size = 25), eps = 0.025, steps = 3L,
                  distance = KolmogorovDist, neighbor = TotalVarNeighborhood())
checkIC(pIC(est2v2))

## Using Cramer-von-Mises MD estimator
est2v3 <- roptest(x, BinomFamily(size = 25), eps = 0.025, steps = 3L, 
                  neighbor = TotalVarNeighborhood())
checkIC(pIC(est2v3))

## comparison of estimates
estimate(est2v)
estimate(est2v1)
estimate(est2v2)
estimate(est2v3)


## 4.1. one-step estimation: radius interval
IC11 <- radiusMinimaxIC(L2Fam=BinomFamily(size=25, prob=estimate(est0)),
                neighbor=ContNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)
(est3c <- oneStepEstimator(x, IC=IC11, start=est0))
checkIC(pIC(est3c))

IC12 <- radiusMinimaxIC(L2Fam=BinomFamily(size=25, prob=estimate(est0)),
                neighbor=TotalVarNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)
(est3v <- oneStepEstimator(x, IC=IC12, start=est0))
checkIC(pIC(est3v))

## maximum radius for given sample size n: sqrt(n)*0.5
(est3c1 <- roptest(x, BinomFamily(size = 25), eps.upper = 0.5))
checkIC(pIC(est3c1))
(est3v1 <- roptest(x, BinomFamily(size = 25), eps.upper = 0.5, neighbor = TotalVarNeighborhood()))
checkIC(pIC(est3v1))

## 4.2. k-step estimation: radius interval
IC11 <- radiusMinimaxIC(L2Fam=BinomFamily(size=25, prob=estimate(est0)),
                neighbor=ContNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)
(est4c <- kStepEstimator(x, IC=IC11, start=est0, steps = 3L))
checkIC(pIC(est4c))

IC12 <- radiusMinimaxIC(L2Fam=BinomFamily(size=25, prob=estimate(est0)),
                neighbor=TotalVarNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)
(est4v <- kStepEstimator(x, IC=IC12, start=est0, steps = 3L))
checkIC(pIC(est4v))

## maximum radius for given sample size n: sqrt(n)*0.5
(est4c1 <- roptest(x, BinomFamily(size = 25), eps.upper = 0.5, steps = 3L))
checkIC(pIC(est4c1))
(est4v1 <- roptest(x, BinomFamily(size = 25), eps.upper = 0.5, neighbor = TotalVarNeighborhood(),
        steps = 3L))
checkIC(pIC(est4v1))
