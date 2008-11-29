library(ROptEst)

###############################################################################
## start of tests
###############################################################################

## positive-definite, symmetric matrices
new("PosDefSymmMatrix", diag(2))
PosDefSymmMatrix(1)
PosDefSymmMatrix(diag(3))


## Distribution symmetry
S1 <- new("NoSymmetry")
type(S1)
NoSymmetry()
S2 <- new("EllipticalSymmetry", SymmCenter = 2)
type(S2)
EllipticalSymmetry(SymmCenter = 1)
S3 <- new("SphericalSymmetry", SymmCenter = -2)
type(S3)
SphericalSymmetry(SymmCenter = c(0,0))
new("DistrSymmList", list(S1, S2, S3))
DistrSymmList(S1, S2, S3)

## Distribution symmetry
S4 <- new("NonSymmetric")
type(S4)
NonSymmetric()
S5 <- new("EvenSymmetric", SymmCenter = -1)
type(S5)
EvenSymmetric(SymmCenter = 0)
S6 <- new("OddSymmetric", SymmCenter = 3)
type(S6)
OddSymmetric(SymmCenter = c(1,1))
new("FunSymmList", list(S4, S5, S6))
FunSymmList(S4, S5, S6)


## parametric family
(PF <- new("ParamFamily"))
plot(PF)
ParamFamily()


## L2-differentiable parametric family
(L2PF <- new("L2ParamFamily"))
plot(L2PF)
L2ParamFamily()


## simple L2-differentiable parametric families
BinomFamily()
BinomFamily(size = 10)
BinomFamily(prob = 0.4)
BinomFamily(size = 100, prob = 0.3)
BinomFamily(size = 50, prob = 0.8, trafo = matrix(-1))

PoisFamily()
PoisFamily(lambda = 10)
PoisFamily(lambda = 20, trafo = matrix(3))

NormLocationFamily()
NormLocationFamily(mean = 2)
NormLocationFamily(sd = 0.1)
NormLocationFamily(mean = -3, sd = 2)
NormLocationFamily(mean = 1, sd = 0.5, trafo = matrix(-1))

GumbelLocationFamily()
GumbelLocationFamily(loc = -2)
GumbelLocationFamily(scale = 2)
GumbelLocationFamily(loc = 1, scale = 0.5)
GumbelLocationFamily(loc = 10, scale = 10, trafo = matrix(0.5))

NormScaleFamily()
NormScaleFamily(sd = 3)
NormScaleFamily(mean = 5)
NormScaleFamily(sd = 0.1, mean = -3)
NormScaleFamily(sd = 2.5, mean = 1, trafo = matrix(0.1))

ExpScaleFamily()
ExpScaleFamily(rate = 0.5)
ExpScaleFamily(rate = 2, trafo = matrix(2))

LnormScaleFamily()
LnormScaleFamily(meanlog = 0.5)
LnormScaleFamily(sdlog = 0.1)
LnormScaleFamily(meanlog = -0.3, sdlog = 2)
LnormScaleFamily(meanlog = 2, sdlog = 1.2, trafo = matrix(2.5))

G1 <- GammaFamily()
name(G1)
name(G1) <- "standard Gamma family"
name(G1)
distribution(G1)
(old <- props(G1))
addProp(G1) <- "test"
props(G1)
props(G1) <- old
props(G1)
param(G1)
main(G1)
nuisance(G1)
trafo(G1)
L2deriv(G1)
L2derivDistr(G1)
L2derivSymm(G1)
FisherInfo(G1)
GammaFamily(scale = 2)
GammaFamily(shape = 0.75)
GammaFamily(scale = 1.5, shape = 2)
GammaFamily(scale = 3, shape = 1.5, trafo = matrix(c(3, 0, 0, 1), ncol = 2))

NormLocationScaleFamily()
NormLocationScaleFamily(mean = 1)
NormLocationScaleFamily(sd = 0.5)
NormLocationScaleFamily(mean = -3, sd = 2)
N1 <- NormLocationScaleFamily(mean = 2, sd = 0.1, trafo = matrix(c(1, 0), ncol = 2))
plot(N1)

## robust models
new("FixRobModel")
(RM1 <- FixRobModel(center = NormLocationFamily()))
FixRobModel(center = PoisFamily(), neighbor = TotalVarNeighborhood(radius = 0.5))
new("InfRobModel")
(RM2 <- InfRobModel(center = NormLocationScaleFamily()))
InfRobModel(center = BinomFamily(size=10), neighbor = TotalVarNeighborhood(radius = 0.2))


###############################################################################
## end of tests
###############################################################################

q("no")
