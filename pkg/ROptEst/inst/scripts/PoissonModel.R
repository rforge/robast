###############################################################################
## Example: Poisson Family
###############################################################################
require(ROptEst)

distroptions("TruncQuantile", 1e-10) # increases numerical support of Pois; 
                                     # i.e., increases precision of the 
                                     # computations
## generates Poisson Family with theta = 4.5
P <- PoisFamily(lambda = 4.5) 
P       # show P
plot(P) # plot of Pois(lambda = 4.5) and L_2 derivative
checkL2deriv(P)

# classical optimal IC
IC0 <- optIC(model = P, risk = asCov())
IC0       # show IC
checkIC(IC0)
Risks(IC0)
plot(IC0) # plot IC

# L_2 family + infinitesimal neighborhood
RobP1 <- InfRobModel(center = P, neighbor = ContNeighborhood(radius = 0.5))
RobP1     # show RobP1
(RobP2 <- InfRobModel(center = P, neighbor = TotalVarNeighborhood(radius = 0.5)))

## lower case radius
lowerCaseRadius(L2Fam = P, ContNeighborhood(radius = 0.5), risk = asMSE())
lowerCaseRadius(L2Fam = P, TotalVarNeighborhood(radius = 0.5), risk = asMSE())

# MSE solution
(IC1 <- optIC(model=RobP1, risk=asMSE()))
checkIC(IC1)
Risks(IC1)
plot(IC1)

(IC2 <- optIC(model=RobP2, risk=asMSE()))
checkIC(IC2)
Risks(IC2)
plot(IC2)


# lower case solutions
(IC3 <- optIC(model=RobP1, risk=asBias()))
checkIC(IC3)
Risks(IC3)
plot(IC3)

(IC4 <- optIC(model=RobP2, risk=asBias()))
checkIC(IC4)
Risks(IC4)
plot(IC4)

# Hampel solution
(IC5 <- optIC(model=RobP1, risk=asHampel(bound=clip(IC1))))
checkIC(IC5)
Risks(IC5)
plot(IC5)

(IC6 <- optIC(model=RobP2, risk=asHampel(bound=Risks(IC2)$asBias), maxiter = 200))
checkIC(IC6)
Risks(IC6)
plot(IC6)


# radius minimax IC
(IC7 <- radiusMinimaxIC(L2Fam=P, neighbor=ContNeighborhood(), 
                risk=asMSE(), loRad=0, upRad=0.5))
checkIC(IC7)
Risks(IC7)
plot(IC7)

(IC8 <- radiusMinimaxIC(L2Fam=P, neighbor=TotalVarNeighborhood(), 
                risk=asMSE(), loRad=0, upRad=Inf))
checkIC(IC8)
Risks(IC8)
plot(IC8)

# least favorable radius
# (may take quite some time!)
(r.rho1 <- leastFavorableRadius(L2Fam=P, neighbor=ContNeighborhood(),
                    risk=asMSE(), rho=0.5))
(r.rho2 <- leastFavorableRadius(L2Fam=P, neighbor=TotalVarNeighborhood(),
                    risk=asMSE(), rho=1/3))

## one-step estimation
## Example: Rutherford-Geiger (1910)
## cf. Feller~(1968), Section VI.7 (a)
x <- c(rep(0, 57), rep(1, 203), rep(2, 383), rep(3, 525), rep(4, 532), 
       rep(5, 408), rep(6, 273), rep(7, 139), rep(8, 45), rep(9, 27), 
       rep(10, 10), rep(11, 4), rep(12, 0), rep(13, 1), rep(14, 1))
       
## 0. mean (classical optimal)
(est0 <- mean(x))

## 1. Kolmogorov(-Smirnov) minimum distance estimator
(est1 <- ksEstimator(x=x, Pois()))

## 2. one-step estimation: radius interval
## 2.1 small amount of contamination < 2%
IC9 <- radiusMinimaxIC(L2Fam=PoisFamily(lambda=est1$lambda),
                neighbor=ContNeighborhood(), risk=asMSE(), loRad=0, upRad=1)
(est21 <- oneStepEstimator(x, IC=IC9, start=est1$lambda))
## 2.2 amount of contamination unknown
IC10 <- radiusMinimaxIC(L2Fam=PoisFamily(lambda=est1$lambda),
                neighbor=ContNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)
(est22 <- oneStepEstimator(x, IC=IC10, start=est1$lambda))

distroptions("TruncQuantile", 1e-5) # default
