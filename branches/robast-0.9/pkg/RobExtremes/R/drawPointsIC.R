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

if (is.null(col)) col = rgb(192,192,192,maxColorValue=255,alpha=0.8) 
if (is.null(bg)) bg = c("#70707050")
if (is.null(bg)) xlab = "x"
if (is.null(ylab)) {ylab = numeric()
for (i in 1:dim(IC0)[1]) ylab[i] = paste("IC of",names(L2fam@param@main[i]),"parameter")}


par(mfrow=c(dim(IC0)[1],1),mar=c(2,4,3,2))

for (i in 1:dim(IC0)[1]){

##plotting IC 
plot(x,IC0[i,],ylab=ylab[i],xlab=xlab,'l',main = ifelse(i==1,paste(IC@name)," "),col="gray")

##plotting observations onto IC
points(x,IC0[i,]
,pch = 21
,cex = abs(IC0[i,]*sqrt(2))
,col = col
,bg = bg
)

grid()
 }
}

##Example
L2fam = GParetoFamily(shape = 0.7)
IC = optIC(L2fam,risk=asCov())
x = q(L2fam@distribution)(seq(0,1,length=50))

##with probit transformation
drawPointsIC(x,IC)
