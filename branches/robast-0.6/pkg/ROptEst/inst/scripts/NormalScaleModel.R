###############################################################################
## Example: Normal Scale
###############################################################################
require(ROptEst)

## generates Normal Scale Family with scale = 3
N0 <- NormScaleFamily(mean=2, sd=3) 
N0        # show G0
plot(N0)  # plot of Norm(mean = 2, sd = 3) and L_2 derivative
checkL2deriv(N0)

## classical optimal IC
N0.IC0 <- optIC(model = N0, risk = asCov())
N0.IC0       # show IC
plot(N0.IC0) # plot IC
checkIC(N0.IC0)
Risks(N0.IC0)

## L_2 family + infinitesimal neighborhood
N0.Rob1 <- InfRobModel(center = N0, neighbor = ContNeighborhood(radius = 0.5))
N0.Rob1     # show N0.Rob1
N0.Rob2 <- InfRobModel(center = N0, neighbor = TotalVarNeighborhood(radius = 0.5))

## MSE solution
(N0.IC1 <- optIC(model=N0.Rob1, risk=asMSE()))
checkIC(N0.IC1)
Risks(N0.IC1)
plot(N0.IC1)

(N0.IC2 <- optIC(model=N0.Rob2, risk=asMSE()))
checkIC(N0.IC2)
Risks(N0.IC2)
plot(N0.IC2)

## lower case solutions
(N0.IC3 <- optIC(model=N0.Rob1, risk=asBias()))
checkIC(N0.IC3)
Risks(N0.IC3)
plot(N0.IC3)
(N0.IC4 <- optIC(model=N0.Rob2, risk=asBias()))
checkIC(N0.IC4)
Risks(N0.IC4)
plot(N0.IC4)

## Hampel solution
(N0.IC5 <- optIC(model=N0.Rob1, risk=asHampel(bound=clip(N0.IC1))))
checkIC(N0.IC5)
Risks(N0.IC5)
plot(N0.IC5)
(N0.IC6 <- optIC(model=N0.Rob2, risk=asHampel(bound=Risks(N0.IC2)$asBias$value), maxiter = 200))
checkIC(N0.IC6)
Risks(N0.IC6)
plot(N0.IC6)

## radius minimax IC
(N0.IC7 <- radiusMinimaxIC(L2Fam=N0, neighbor=ContNeighborhood(), 
                risk=asMSE(), loRad=0, upRad=Inf))
checkIC(N0.IC7)
Risks(N0.IC7)
plot(N0.IC7)
(N0.IC8 <- radiusMinimaxIC(L2Fam=N0, neighbor=TotalVarNeighborhood(), 
                risk=asMSE(), loRad=0, upRad=Inf))
checkIC(N0.IC8)
Risks(N0.IC8)
plot(N0.IC8)

## least favorable radius
## (may take quite some time!)
(N0.r.rho1 <- leastFavorableRadius(L2Fam=N0, neighbor=ContNeighborhood(),
                    risk=asMSE(), rho=0.5))
(N0.r.rho2 <- leastFavorableRadius(L2Fam=N0, neighbor=TotalVarNeighborhood(),
                    risk=asMSE(), rho=1/3))
