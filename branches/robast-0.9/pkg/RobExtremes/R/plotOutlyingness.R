##########################################
##                                      ## 
##    Wrapper for outlyingnessPlot.R    ##
##                                      ##
##                                      ##
##########################################

if(FALSE){
##projection distance
  qfun = function(x){p0 = p(X)(x); q0 = q(X)(p0)}
  QProj <- function(){new("NormType", name="Quantiles", fct=qfun)}

##@x - dataset
##@X - random variable
##@fam - parameter family
##@alpha - confidence level for quantile
#
plotOutlyingness = function(x,alpha=0.99,X=GPareto(),fam=GParetoFamily()){

  ##logarithmic representation (for distributions with positive support)
  fam@distribution = log(fam@distribution)

  ##classical IC
  ICmle <- optIC(model=fam,risk=asCov())

  ##parameter for plotting
  par(cex=1,bty="n")

##call of routine from RobAStBase  
outlyingPlotIC(x
 ,IC.x = ICmle
 ,IC.y = ICmle
 ,dist.x = QProj()
 #NormType() - Euclidean norm, default - Mahalanobis norm
 #,dist.y = NormType()
 ,adj = 0.1
 ,pch = 16
 ,col = rgb(152,152,152,maxColorValue=255)
 ,col.idn = rgb(102,102,102,maxColorValue=255)
 ,cex.idn = 1.7
 ,col.cutoff = rgb(202,202,202,maxColorValue=255) 
 ,offset = 0
 ,cutoff.quantile.y = 0.99
 ,cutoff.quantile.x = 0.99
 ,cutoff.x = cutoff()
 ,cutoff.y = cutoff.sememp()
 ,robCov.x = TRUE
 ,robCov.y = TRUE
 ,tf.x = function(x)log(x)
 ,cex.main = 1.5
 ,cex.lab = 1.5
 ,cex = 1.5
 #,col.lab=FhGred
 ,lwd.cutoff = 3
 #,jitt.fac = 300
 ,col.abline = rgb(102,102,102,maxColorValue=255)
 ,cex.abline = 1.5
 ,main = ""#"Outlyingness Plot"
 ,xlab="Theoretical log-quantiles"
 ,ylab="Mahalanobis distance"
)
}

##Example
X= GPareto()
fam = GParetoFamily()
x = r(X)(1000)
plotOutlyingness(x,alpha=0.9,X=X,fam=fam)
}