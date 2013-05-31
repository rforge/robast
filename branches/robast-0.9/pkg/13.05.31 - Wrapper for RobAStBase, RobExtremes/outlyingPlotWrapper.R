##########################################
##                                      ## 
##    Wrapper for outlyingPlot.R        ##
##                                      ##
##                                      ##
##########################################

## projection distance
qfun = function(x){p0 = p(X)(x); q0 = q(X)(p0)}
QProj <- function(){new("NormType", name="Quantiles", fct=qfun)}

##@x - dataset
##@X - random variable
##@fam - parameter family
##@alpha - confidence level for quantile
## alpha.trsp - optional transparency of the plot
## with.legend - optional legend indicator
## withCall - optional indicator of the function call
#
outlyingPlotWrapper = function(x,alpha,fam,...,alpha.trsp = 100, with.legend = TRUE, withCall = TRUE){
  ###
  ### 1. grab the dots (and manipulate it within the wrapper function)
  ###
  ###
  ### do something to fix the good default arguments
  ###
  mc <- as.list(match.call(expand.dots = FALSE))[-1]
  dots <- mc$"..."
  if(is.null(mc$alpha.trsp)) alpha.trsp <- 100
  if(is.null(mc$with.legend)) mc$with.legend <- TRUE
  if(is.null(mc$withCall)) mc$withCall <- TRUE
  if(missing(x)) stop("Argument 'x' must be given as argument to 'outlyingPlotWrapper'")
  if(missing(alpha)) stop("Argument 'alpha' must be given as argument to 'outlyingPlotWrapper'")
  if(missing(fam)) stop("Argument 'fam' must be given as argument to 'outlyingPlotWrapper'")
  ###
  ### 2. build up the argument list for the (powerful/fullfledged)
  ### graphics/diagnostics function;
  ##

  argsList <- list(data = substitute(x)
                   ,IC.x = substitute(optIC(model = fam, risk = asCov()))
                   ,IC.y = substitute(optIC(model = fam, risk = asCov()))
                   ,dist.x = substitute(QProj())
                   #NormType() - Euclidean norm, default - Mahalanobis norm
                   ,dist.y = substitute(NormType())
                   ,adj = 0.1
                   ,pch = 21
                   ,col.idn = substitute(rgb(102,102,102,maxColorValue=255))
                   ,cex.idn = substitute(1.7)
                   ,col.cutoff = substitute(rgb(202,202,202,maxColorValue=255))
                   ,offset = substitute(0)
                   ,cutoff.quantile.x = substitute(alpha)
                   ,cutoff.quantile.y = substitute(alpha)
                   ,cutoff.x = substitute(cutoff())
                   ,cutoff.y = substitute(cutoff.sememp())
                   ,robCov.x = substitute(TRUE)
                   ,robCov.y = substitute(TRUE)
                   ,tf.x = substitute(function(x)log(x - q(fam@distribution)(0)))
                   ,cex.main = substitute(1.5)
                   ,cex.lab = substitute(1.5)
                   ,cex = substitute(1.5)
                   ,lwd.cutoff = substitute(3)
                   ,col.abline = substitute(rgb(52,52,52,maxColorValue=255))
                   ,cex.abline = substitute(1.2)
                   ,adj.abline = substitute(c(0.8, 0.2))
                   ,main = ""#"Outlyingness Plot"
                   ,xlab=substitute("Theoretical log-quantiles")
                   ,ylab=substitute("Mahalanobis distance")
                   ,bty = substitute("n")
                   ,col = substitute(addAlphTrsp2col(rgb(0,0,0,maxColorValue=255), substitute(alpha.trsp)))
                   )
  
  ##parameter for plotting
  if(mc$with.legend)
  {
    argsList$col.main <- "black"
    argsList$col.lab <- "black"
  }
  else
  {
    argsList$col.main <- "white"
    argsList$col.lab <- "white"
  }
  
  args <- c(argsList, dots)
  ###
  ### 3. build up the call but grab it and write it into an object
  ###
  cl <- substitute(do.call(outlyingPlotIC,args0), list(args0=args))
  ### manipulate it so that the wrapper do.call is ommitted
  cl0 <- as.list(cl)[-1]
  mycall <- c(cl0[1],unlist(cl0[-1]))
  mycall <- as.call(mycall)
  ###
  ### 4. evaluate the call (i.e., produce the graphic)
  ###
  eval(mycall)
  ###
  ### 5. return the call (if withCall==TRUE)
  ###
  if(mc$withCall) print(mycall)
  
}

##Examples
require(RobExtremes)
require(distr)

# GPD
dev.new()
fam = GParetoFamily()
X=distribution(fam)
x = r(X)(1000)
outlyingPlotWrapper(x,alpha=0.99,fam=fam, alpha.trsp=50, with.legend = FALSE)

# GEV
dev.new()
fam = GEVFamily()
X=distribution(fam)
x = r(X)(1000)
outlyingPlotWrapper(x,alpha=0.95,fam=fam, with.legend = TRUE, withCall = TRUE)

# Gamma
dev.new()
fam = GammaFamily()
X=distribution(fam)
x = r(X)(1000)
outlyingPlotWrapper(x,alpha=0.95,fam=fam, alpha.trsp=70)

# Weibull
dev.new()
fam = WeibullFamily()
X=distribution(fam)
x = r(X)(1000)
outlyingPlotWrapper(x,alpha=0.95,fam=fam, alpha.trsp=50, with.legend = TRUE, withCall = FALSE)

