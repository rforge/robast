##########################################
##                                      ## 
##    Wrapper for outlyingnessPlot.R    ##
##                                      ##
##                                      ##
##########################################

##projection distance
qfun = function(x){p0 = p(X)(x); q0 = q(X)(p0)}
QProj <- function(){new("NormType", name="Quantiles", fct=qfun)}

##@x - dataset
##@X - random variable
##@fam - parameter family
##@alpha - confidence level for quantile
#
plotOutlyingness = function(x,alpha=0.99,fam=GParetoFamily(),...,alpha.trsp = 100, with.legend = TRUE){
  
  mc <- as.list(match.call(expand.dots = FALSE))[-1]
  dots <- mc$"..."
  if(is.null(mc$alpha.trsp)) mc$alpha.trsp <- 100
  if(is.null(mc$with.legend)) mc$with.legend <- TRUE
  if(missing(x)) stop("Argument 'x' must be given as argument to 'plotOutlyingness'")
  if(missing(alpha)) stop("Argument 'alpha' must be given as argument to 'plotOutlyingness'")
  if(missing(fam)) stop("Argument 'fam' must be given as argument to 'plotOutlyingness'")

  ##logarithmic representation (for distributions with positive support)
  fam@distribution = log(fam@distribution - q(fam@distribution)(0))

  ##classical IC
  ICmle <- optIC(model=fam,risk=asCov())

  ##parameter for plotting
  if(mc$with.legend)
  {par(cex=1,bty="n", col = addAlphTrsp2col(rgb(0,0,0,maxColorValue=255), mc$alpha.trsp),
       col.main = "black", col.lab = "black")}
  else
  {par(cex=1,bty="n", col = addAlphTrsp2col(rgb(0,0,0,maxColorValue=255), mc$alpha.trsp), 
       col.main = "white", col.lab = "white")}
  
  cutoff.quantile.x = alpha
  cutoff.quantile.y = alpha
  
##call of routine from RobAStBase  
outlyingPlotIC(x
 ,IC.x = ICmle
 ,IC.y = ICmle
 ,dist.x = QProj()
 #NormType() - Euclidean norm, default - Mahalanobis norm
 #,dist.y = NormType()
 ,adj = 0.1
 ,pch = 21
 ,col.idn = rgb(102,102,102,maxColorValue=255)
 ,cex.idn = 1.7
 ,col.cutoff = rgb(202,202,202,maxColorValue=255) 
 ,offset = 0
 ,cutoff.quantile.x = cutoff.quantile.x
 ,cutoff.quantile.y = cutoff.quantile.y
 ,cutoff.x = cutoff()
 ,cutoff.y = cutoff.sememp()
 ,robCov.x = TRUE
 ,robCov.y = TRUE
 ,tf.x = function(x)log(x - q(fam@distribution)(0))
 ,cex.main = 1.5
 ,cex.lab = 1.5
 ,cex = 1.5
 #,col.lab=FhGred
 ,lwd.cutoff = 3
 #,jitt.fac = 300
 ,col.abline = rgb(52,52,52,maxColorValue=255)
 ,cex.abline = 1.2
 ,adj.abline = c(0.8, 0.2)
 ,main = ""#"Outlyingness Plot"
 ,xlab="Theoretical log-quantiles"
 ,ylab="Mahalanobis distance"
)
}

##Examples
require(RobExtremes)
require(distr)

# # GPD
# fam = GParetoFamily()
# X=distribution(fam)
# x = r(X)(1000)
# plotOutlyingness(x,alpha=0.99,fam=fam,alpha.trsp=50, with.legend = FALSE)

# GEV
fam = GEVFamily()
X=distribution(fam)
x = r(X)(1000)
plotOutlyingness(x,alpha=0.95,fam=fam,alpha.trsp=50, with.legend = TRUE)

# # Gamma
# fam = GammaFamily()
# X=distribution(fam)
# x = r(X)(1000)
# plotOutlyingness(x,alpha=0.95,fam=fam,alpha.trsp=50, with.legend = TRUE)

# # Weibull
# fam = WeibullFamily()
# X=distribution(fam)
# x = r(X)(1000)
# plotOutlyingness(x,alpha=0.95,fam=fam,alpha.trsp=50, with.legend = TRUE)

