if(FALSE){
require(ROptEst)
options("newDevice"=TRUE)

## generates normal location and scale family with mean = -2 and sd = 3
 ### checks for lower case in various standardizations
N0 <- NormLocationScaleFamily(mean=-2, sd=3)
N0.Rob1<- InfRobModel(center = N0, neighbor = ContNeighborhood(radius = 15));
N0.IC2.MBRE <- optIC(model = N0.Rob1, risk = asBias(), tol = 1e-10);print(stand(N0.IC2));print(cent(N0.IC2));print(stand(N0.IC2)/max(stand(N0.IC2)));print(cent(N0.IC2)/max(stand(N0.IC2)));print(clip(N0.IC2))
plot(N0.IC2.MBRE)
N0.IC2.OMSE <- optIC(model = N0.Rob1, risk = asMSE(), tol = 1e-10);print(stand(N0.IC2.i)/max(stand(N0.IC2.i)));print(cent(N0.IC2.i)/max(stand(N0.IC2.i)));print(clip(N0.IC2.i));print(stand(N0.IC2)/max(stand(N0.IC2)));print(cent(N0.IC2)/max(stand(N0.IC2)));print(clip(N0.IC2))
plot(N0.IC2.OMSE)
N0.IC2.MBRE.i <- optIC(model = N0.Rob1, risk = asBias(normtype=InfoNorm()), tol = 1e-10);print(stand(N0.IC2.i));print(cent(N0.IC2.i));print(stand(N0.IC2.i)/max(stand(N0.IC2.i)));print(cent(N0.IC2.i)/max(stand(N0.IC2.i)));print(clip(N0.IC2.i));print(stand(N0.IC2)/max(stand(N0.IC2)));print(cent(N0.IC2)/max(stand(N0.IC2)));print(clip(N0.IC2))
plot(N0.IC2.MBRE.i)
N0.IC2.OMSE.i <- optIC(model = N0.Rob1, risk = asMSE(normtype=InfoNorm()), tol = 1e-10);print(stand(N0.IC2.i)/max(stand(N0.IC2.i)));print(cent(N0.IC2.i)/max(stand(N0.IC2.i)));print(clip(N0.IC2.i));print(stand(N0.IC2)/max(stand(N0.IC2)));print(cent(N0.IC2)/max(stand(N0.IC2)));print(clip(N0.IC2))
plot(N0.IC2.OMSE.i)
N0.IC2.MBRE.s <- optIC(model = N0.Rob1, risk = asBias(normtype=SelfNorm()), tol = 1e-10);print(stand(N0.IC2.i));print(cent(N0.IC2.i));print(stand(N0.IC2.i)/max(stand(N0.IC2.i)));print(cent(N0.IC2.i)/max(stand(N0.IC2.i)));print(clip(N0.IC2.i));print(stand(N0.IC2)/max(stand(N0.IC2)));print(cent(N0.IC2)/max(stand(N0.IC2)));print(clip(N0.IC2))
plot(N0.IC2.MBRE.s)
N0.IC2.OMSE.s <- optIC(model = N0.Rob1, risk = asMSE(normtype=SelfNorm()), tol = 1e-10);print(stand(N0.IC2.i)/max(stand(N0.IC2.i)));print(cent(N0.IC2.i)/max(stand(N0.IC2.i)));print(clip(N0.IC2.i));print(stand(N0.IC2)/max(stand(N0.IC2)));print(cent(N0.IC2)/max(stand(N0.IC2)));print(clip(N0.IC2))
plot(N0.IC2.OMSE.s)
}
