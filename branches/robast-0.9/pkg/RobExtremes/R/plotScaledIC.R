##########################################
##                                      ## 
##       IC Scaling of axes             ##
##                                      ##
##########################################


##@IC - influence function
##@L2fam - L2 parametric family
##@probit - if probit transformation should be done, otherwise is a log-transformation of x-axis
##@xlim, ylim - plot limits
##@xticks, yticks - ticks for plot



plotScaledIC = function(IC,L2fam,probit=TRUE,xlim=NULL,ylim=NULL,xticks=NULL,yticks=NULL){

## get distribution from L2-diff family
X = L2fam@distribution

## get L2 derivative from L2-diff family
L2derivative = L2fam@L2deriv.fct(L2fam@param)

## dimension of derivative
dim = length(L2derivative)

xgrid = seq(1e-9,1-1e-9,length=1000)

## quantile transformation of x-scale
x = q(X)(xgrid)

## evaluation of IC on quantiles
ICeval = evalIC(IC,as.matrix(x))

ICtr = ICeval

## set plot limits in (0,1)-cub for probit-scaled IC
if(probit){
xlim = c(0,1)
ylim = c(0,1)}

minIC = min(ICtr)
imin = which(ICtr == minIC)
 
maxIC = max(ICtr)
imax = which(ICtr == maxIC)

## calculate reasonable plot limits for a non-transformed IC
if ((is.null(ylim))&&(is.null(xlim))){

if (!probit){

if (dim == 1){
  xmin = x[imin]
  xmax = x[imax]
}else{
  dist = numeric(length(xgrid))
for (i in 1:length(xgrid)){ 
  sum = 0
for (j in 1:dim){sum = ICtr[j,i]^2+sum}
  dist[i] = sqrt(sum)
 }
}

idist.max = which(dist == max(dist))
idist.min = which(dist == min(dist))

xmax.IC = x[idist.max]
xmin.IC = x[idist.min]

cat("\n IC takes its maximal values at",xmax.IC, "and minimal values at",xmin.IC,"\n")

xmin = min(xmin.IC,xmax.IC)
xmax = max(xmin.IC,xmax.IC)

spanx = (xmax-xmin)/length(x)
spany = (maxIC-minIC)/length(x) 

xlim = c(min(x),max(x))
ylim = c(minIC,maxIC)
 }
}

if (probit){
## optimal min bias IC
 robModel <- InfRobModel(center = L2fam, neighbor = ContNeighborhood(radius = 0.5))
 ICmbr = optIC(model = robModel, risk = asBias())
 sd = 2*ICmbr@clip 

## probit transformation of IC
 for (i in 1:dim) ICtr[i,] = pnorm(ICeval[i,],sd=sd)
}

if(is.null(xticks)&&is.null(yticks)){
## ticks for x,y-axis
if (is.null(xticks)) xticks = c(-Inf,-1000,-100,-10,-5,-1,-0.5,-0.1,0.0,0.1,0.5,1,5,10,100,1000,Inf)
if (is.null(yticks)) yticks = c(-Inf,-1000,-100,-10,-5,-1,-0.5,-0.1,0.0,0.1,0.5,1,5,10,100,1000,Inf)
}

if (probit){
## transformation of ticks for x-axis
 xax = round(p(X)(xticks),1)
## probit transformation of ticks for y-axis
 yax = round(pnorm(yticks,sd=sd),1)
}else{
 xax = round(q(X)(p(X)(xticks)),1)
 yax = round(yticks,1)
}

## do ticks lie to near to each other?
indy = !duplicated(yax)
yax = yax[indy]
yticks = yticks[indy]
indx = !duplicated(xax)
xax = xax[indx]
xticks = xticks[indx]

## function creating a grid on plot
fun = function(){for(i in 1:length(yax)){
 abline(v=xax[i],col="grey",lty="dashed",lwd=0.1)
 abline(h=yax[i],col="grey",lty="dashed",lwd=0.1)
 }
}

## plotting 
if (!probit){
 plot(x, ICtr[1,], 'l', ylim=ylim, xlim=xlim, main=" ", xlab="x", ylab="IC", axes=FALSE, col=1,log="x")}
else{
 plot(xgrid, ICtr[1,], 'l', ylim=ylim, xlim=xlim, main=" ", xlab="x", ylab="probit(IC)", axes=FALSE, col=1)}

##grid generation 
fun()

## plotting all of the IC components
if (dim > 1){
 if (!probit){
  for (i in 2:dim) lines(x, ICtr[i,], col = i)}
   else{
  for (i in 2:dim) lines(xgrid, ICtr[i,], col = i)
 }
}

## drawing the ticks and labels on the plot
axis(side=1, at= xax, labels = c(expression(paste(-infinity)),xticks[2:(length(xticks)-1)],expression(paste(infinity))), pos=0)
axis(side=2, at= yax, labels = c(expression(paste(-infinity)),yticks[2:(length(yticks)-1)],expression(paste(infinity))))
}

##Example

L2fam = GParetoFamily(shape = 0.7)
IC = optIC(L2fam,risk=asCov())

##with probit transformation
plotScaledIC(IC,L2fam,probit = TRUE)

##without probit transformation
plotScaledIC(IC,L2fam,probit = FALSE,xlim=c(0.01,100),ylim=c(-10,10),yticks=c(-Inf,-10,-5,-1,0,1,5,10,Inf),xticks=c(0.01,1,5,10,100,Inf))