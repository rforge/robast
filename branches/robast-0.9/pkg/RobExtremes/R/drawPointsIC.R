##########################################
##                                      ## 
##       Function for drawing IC        ##
##       on observations                ##
##                                      ##
##########################################

##@x - observations
##@IC - 1-dim influence function (object)
##@col - color of points
##@bg - background for semitranparence
##@cex - scaling of points


drawPointsIC = function(x,IC,col=NULL,bg=NULL,cex=NULL,xlab=NULL,ylab=NULL){

##evaluation of IC on observations
IC0 = evalIC(IC,as.matrix(x)) 

if (is.null(col)) col = rep(rgb(192,192,192,maxColorValue=255,alpha=0.8), length= dim(IC0)[1])
if (is.null(bg)) bg = rep(c("#70707050"), length = dim(IC0)[1])
if (is.null(xlab)) xlab = "x"
if (is.null(ylab)) {ylab = numeric()
for (i in 1:dim(IC0)[1]) ylab[i] = paste("IC of",names(L2fam@param@main[i]),"parameter")}


par(mfrow=c(dim(IC0)[1],1),mar=c(2,4,3,2))

for (i in 1:dim(IC0)[1]){

##plotting IC 
plot(x,IC0[i,],ylab=ylab[i],xlab=xlab,'l',main = ifelse(i==1,paste(IC@name)," "),col="gray")

##plotting observations onto IC
points(x,IC0[i,]
,pch = 21
,cex = (colSums(IC0^2))^(1/2)*sqrt(2)
,col = col[i]
,bg = bg[i]
)

grid()
 }
}

##Example

##GPD family
L2fam.GPD = GParetoFamily(shape = 0.7)
IC.GPD.MLE = optIC(L2fam.GPD,risk=asCov())
#IC.GPD.OMSE = optIC(L2fam.GPD,risk=asMSE())

xgpd = q(L2fam.GPD@distribution)(seq(0,1,length=50))

# par(mfcol=c(1,2))

pdf("IC_GPD.pdf")
drawPointsIC(xgpd,IC.GPD.MLE,col=c(rgb(7,154,121,maxColorValue=255,alpha=0.8),
rgb(7,154,121,maxColorValue=255,alpha=0.8)),bg=c(c("#00640040"),c("#00640040")))

# drawPointsIC(xgpd,IC.GPD.OMSE,col=c(rgb(7,154,121,maxColorValue=255,alpha=0.8),
# rgb(7,154,121,maxColorValue=255,alpha=0.8)),bg=c(c("#ff100040"),c("#ff100040")))
dev.off()

##GEV family
L2fam.GEV = GEVFamily(shape = 0.7)
IC.GEV.MLE = optIC(L2fam.GEV,risk=asCov())
#IC.GEV.OMSE = optIC(L2fam.GEV,risk=asMSE())

xgev = q(L2fam.GEV@distribution)(seq(0,1,length=50))

pdf("IC_GEV.pdf")
# par(mfcol=c(1,2))
drawPointsIC(xgev,IC.GEV.MLE,col=c(rgb(7,154,121,maxColorValue=255,alpha=0.8),
rgb(7,154,121,maxColorValue=255,alpha=0.8)),bg=c(c("#00640040"),c("#00640040")))

# drawPointsIC(xgpd,IC.GEV.OMSE,col=c(rgb(7,154,121,maxColorValue=255,alpha=0.8),
# rgb(7,154,121,maxColorValue=255,alpha=0.8)),bg=c(c("#ff100040"),c("#ff100040")))
dev.off()

