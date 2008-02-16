###############################################################################
## Example: Exponential Scale Family
###############################################################################
require(ROptEst)

## generates Exponential Scale Family with rate = 1
E1 <- ExpScaleFamily(rate = 1) 
E1        # show E1
plot(E1)  # plot of Exp(rate = 1) and L_2 derivative
checkL2deriv(E1)

# classical optimal IC
E1.IC0 <- optIC(model = E1, risk = asCov())
E1.IC0       # show IC
checkIC(E1.IC0)
Risks(E1.IC0)
plot(E1.IC0) # plot IC

# L_2 family + infinitesimal neighborhood
E1.Rob1 <- InfRobModel(center = E1, neighbor = ContNeighborhood(radius = 0.5))
E1.Rob1     # show E1.Rob1
E1.Rob2 <- InfRobModel(center = E1, neighbor = TotalVarNeighborhood(radius = 0.5))

# MSE solution
(E1.IC1 <- optIC(model=E1.Rob1, risk=asMSE()))
checkIC(E1.IC1)
Risks(E1.IC1)
plot(E1.IC1)
(E1.IC2 <- optIC(model=E1.Rob2, risk=asMSE()))
checkIC(E1.IC2)
Risks(E1.IC2)
plot(E1.IC2)

# lower case solutions
(E1.IC3 <- optIC(model=E1.Rob1, risk=asBias()))
checkIC(E1.IC3)
Risks(E1.IC3)
plot(E1.IC3)
(E1.IC4 <- optIC(model=E1.Rob2, risk=asBias()))
checkIC(E1.IC4)
Risks(E1.IC4)
plot(E1.IC4)

# Hampel solution
(E1.IC5 <- optIC(model=E1.Rob1, risk=asHampel(bound=clip(E1.IC1))))
checkIC(E1.IC5)
Risks(E1.IC5)
plot(E1.IC5)
(E1.IC6 <- optIC(model=E1.Rob2, risk=asHampel(bound=Risks(E1.IC2)$asBias), maxiter = 200))
checkIC(E1.IC6)
Risks(E1.IC6)
plot(E1.IC6)

# radius minimax IC
(E1.IC7 <- radiusMinimaxIC(L2Fam=E1, neighbor=ContNeighborhood(), 
                risk=asMSE(), loRad=0, upRad=0.5))
checkIC(E1.IC7)
Risks(E1.IC7)
(E1.IC8 <- radiusMinimaxIC(L2Fam=E1, neighbor=TotalVarNeighborhood(), 
                risk=asMSE(), loRad=0, upRad=0.5))
checkIC(E1.IC8)
Risks(E1.IC8)

# least favorable radius
# (may take quite some time!)
(E1.r.rho1 <- leastFavorableRadius(L2Fam=E1, neighbor=ContNeighborhood(),
                    risk=asMSE(), rho=0.5))
(E1.r.rho2 <- leastFavorableRadius(L2Fam=E1, neighbor=TotalVarNeighborhood(),
                    risk=asMSE(), rho=1/3))

## one-step estimation
## 1. generate a contaminated sample
ind <- rbinom(1e2, size=1, prob=0.05) 
E1.x <- rexp(1e2, rate=(1-ind)*2+ind*10)

## 2. Kolmogorov(-Smirnov) minimum distance estimator
(E1.est0 <- ksEstimator(x=E1.x, Exp()))

## 3. one-step estimation: radius known
E1.Rob3 <- InfRobModel(center=ExpScaleFamily(rate=E1.est0$rate), 
                neighbor=ContNeighborhood(radius=0.5))
E1.IC9 <- optIC(model=E1.Rob3, risk=asMSE())
(E1.est1 <- oneStepEstimator(E1.x, IC=E1.IC9, start=E1.est0$rate))

## 4. one-step estimation: radius interval
E1.IC10 <- radiusMinimaxIC(L2Fam=ExpScaleFamily(rate=E1.est0$rate),
                neighbor=ContNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)
(E1.est2 <- oneStepEstimator(E1.x, IC=E1.IC10, start=E1.est0$rate))
