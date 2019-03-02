if(FALSE){
require(ROptEst)
options("newDevice"=TRUE)

## generates normal location and scale family with mean = -2 and sd = 3
 ### checks for lower case in various standardizations
N0 <- NormLocationScaleFamily(mean=-2, sd=3)
N0.Rob1<- InfRobModel(center = N0, neighbor = ContNeighborhood(radius = 15));
N0.IC2.MBRE <- optIC(model = N0.Rob1, risk = asBias(), tol = 1e-10);print(stand(N0.IC2.MBRE));print(cent(N0.IC2.MBRE));print(stand(N0.IC2.MBRE)/max(stand(N0.IC2.MBRE)));print(cent(N0.IC2.MBRE)/max(stand(N0.IC2.MBRE)));print(clip(N0.IC2.MBRE))
plot(N0.IC2.MBRE)
N0.IC2.OMSE <- optIC(model = N0.Rob1, risk = asMSE(), tol = 1e-10);print(stand(N0.IC2.OMSE)/max(stand(N0.IC2.OMSE)));print(cent(N0.IC2.OMSE)/max(stand(N0.IC2.OMSE)));print(clip(N0.IC2.OMSE));print(stand(N0.IC2.OMSE)/max(stand(N0.IC2.OMSE)));print(cent(N0.IC2.OMSE)/max(stand(N0.IC2.OMSE)));print(clip(N0.IC2.OMSE))
plot(N0.IC2.OMSE)
N0.IC2.MBRE.i <- optIC(model = N0.Rob1, risk = asBias(normtype=InfoNorm()), tol = 1e-10);print(stand(N0.IC2.MBRE.i));print(cent(N0.IC2.MBRE.i));print(stand(N0.IC2.MBRE.i)/max(stand(N0.IC2.MBRE.i)));print(cent(N0.IC2.MBRE.i)/max(stand(N0.IC2.MBRE.i)));print(clip(N0.IC2.MBRE.i));print(stand(N0.IC2.MBRE.i)/max(stand(N0.IC2.MBRE.i)));print(cent(N0.IC2.MBRE.i)/max(stand(N0.IC2.MBRE.i)));print(clip(N0.IC2.MBRE.i))
plot(N0.IC2.MBRE.i)
N0.IC2.OMSE.i <- optIC(model = N0.Rob1, risk = asMSE(normtype=InfoNorm()), tol = 1e-10);print(stand(N0.IC2.OMSE.i)/max(stand(N0.IC2.OMSE.i)));print(cent(N0.IC2.OMSE.i)/max(stand(N0.IC2.OMSE.i)));print(clip(N0.IC2.OMSE.i));print(stand(N0.IC2.OMSE.i)/max(stand(N0.IC2.OMSE.i)));print(cent(N0.IC2.OMSE.i)/max(stand(N0.IC2.OMSE.i)));print(clip(N0.IC2.OMSE.i))
plot(N0.IC2.OMSE.i)
N0.IC2.MBRE.s <- optIC(model = N0.Rob1, risk = asBias(normtype=SelfNorm()), tol = 1e-10);print(stand(N0.IC2.MBRE.s));print(cent(N0.IC2.MBRE.s));print(stand(N0.IC2.MBRE.s)/max(stand(N0.IC2.MBRE.s)));print(cent(N0.IC2.MBRE.s)/max(stand(N0.IC2.MBRE.s)));print(clip(N0.IC2.MBRE.s));print(stand(N0.IC2.MBRE.s)/max(stand(N0.IC2.MBRE.s)));print(cent(N0.IC2.MBRE.s)/max(stand(N0.IC2.MBRE.s)));print(clip(N0.IC2.MBRE.s))
plot(N0.IC2.MBRE.s)
N0.IC2.OMSE.s <- optIC(model = N0.Rob1, risk = asMSE(normtype=SelfNorm()), tol = 1e-10);print(stand(N0.IC2.OMSE.s)/max(stand(N0.IC2.OMSE.s)));print(cent(N0.IC2.OMSE.s)/max(stand(N0.IC2.OMSE.s)));print(clip(N0.IC2.OMSE.s));print(stand(N0.IC2.OMSE.s)/max(stand(N0.IC2.OMSE.s)));print(cent(N0.IC2.OMSE.s)/max(stand(N0.IC2.OMSE.s)));print(clip(N0.IC2.OMSE.s))
plot(N0.IC2.OMSE.s)
}
