require(RobExtremes)

checkSn <- function(xi, PdG=GPareto){
  dist <- PdG(shape=xi, scale= 1)
  s1 <- system.time(Sn1 <- Sn(dist))
  dist2 <- as(dist,"AbscontDistribution")
  s2 <- system.time(Sn2 <- Sn(dist2))
  c(Sn1,Sn2,Sn2-Sn1,(Sn2/Sn1-1)*100,s1[3],s2[3])
}

xiv <- c(0.001,0.2,0.7, 1.0, 2.0,5.0, 10.0,15)
system.time(sapply(xiv, function(xi){
       res <- rbind(GPD=checkSn(xi,GPareto),
                    GEV=checkSn(xi,GEV),
                    Gamma=checkSn(xi,Gammad),
                    Weibull= checkSn(xi,Weibull))
       colnames(res) <- c("Sn-int.", "Sn-ex.", "DeltaSn", "relDel",
                          "Time-int", "Time-ex.")
       cat("\nxi=",xi,":\n")
       print(round(res,4))
       return(invisible(NULL))
       }))

gSn <- function(xi, PF) Sn(PF(shape=xi,scale=1))
gSna <- function(xi, PF) Sn(as(PF(shape=xi,scale=1),"AbscontDistribution"))
xig <- seq(0.01,10,by=0.05)
system.time({ S1g <- sapply(xig, gSn, GPareto)
  S2g <- sapply(xig, gSn, GEV)
  S3g <- sapply(xig, gSn, Gammad)
  S4g <- sapply(xig, gSn, Weibull)
})
system.time({S1ga <- sapply(xig, gSna, GPareto)
S2ga <- sapply(xig, gSna, GEV)
S3ga <- sapply(xig, gSna, Gammad)
S4ga <- sapply(xig, gSna, Weibull)})
par(mfrow=c(2,2))
plot(xig, S1g, type="l")
lines(xig, S1ga, col="red")
plot(xig, S2g, type="l")
lines(xig, S2ga, col="red")
plot(xig, S3g, type="l")
lines(xig, S3ga, col="red")
plot(xig, S4g, type="l")
lines(xig, S4ga, col="red")
par(mfrow=c(1,1))
