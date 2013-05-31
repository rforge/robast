##########################################
##                                      ## 
##    Wrapper for infoPlot.R            ##
##    (infoPlot method for IC)          ##
##                                      ##
##########################################

##IC - influence curve
##data - dataset
## with.legend - optional legend indicator
## withCall - optional indicator of the function call
#
infoPlotWrapper = function(IC, data,...,alpha.trsp = 100,with.legend = TRUE, withCall = TRUE){
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
  if(missing(IC)) stop("Argument 'IC' must be given as argument to 'ICAllPlotWrapper'")
  if(missing(data)) data <- NULL
  ###
  ### 2. build up the argument list for the (powerful/fullfledged)
  ### graphics/diagnostics function;
  ##  
  
    argsList <- list(object = substitute(IC)
                     ,data = substitute(data)
                     ,withSweave = substitute(getdistrOption("withSweave"))
                     ,lwd = substitute(par("lwd"))
                     ,lty = substitute("solid")
                     ,colI = substitute(grey(0.5))
                     ,lwdI = substitute(0.7*par("lwd"))
                     ,ltyI = substitute("dotted")
                     ,main = substitute(FALSE)
                     ,inner = substitute(TRUE)
                     ,sub = substitute(FALSE)
                     ,col.inner = substitute(par("col.main"))
                     ,cex.inner = substitute(0.8)
                     ,bmar = substitute(par("mar")[1])
                     ,tmar = substitute(par("mar")[3])
                     ,with.legend = substitute(TRUE)
                     ,legend = substitute(NULL)
                     ,legend.bg = substitute("white")
                     ,legend.location = substitute("bottomright")
                     ,legend.cex = substitute(0.8)
                     ,scaleX = substitute(FALSE)
                     ,scaleX.fct = substitute(p(eval(object@CallL2Fam)))
                     ,scaleX.inv = substitute(q(eval(object@CallL2Fam)))
                     ,scaleY = substitute(FALSE)
                     ,scaleY.fct = substitute(pnorm)
                     ,scaleY.inv=substitute(qnorm)
                     ,scaleN = substitute(9)
                     ,x.ticks = substitute(NULL)
                     ,y.ticks = substitute(NULL)
                     ,mfColRow = substitute(TRUE)
                     ,to.draw.arg = substitute(NULL)
                     ,cex.pts = substitute(1)
                     ,col.pts = substitute(addAlphTrsp2col(rgb(0,255,0,maxColorValue=255), substitute(alpha.trsp)))
                     ,pch.pts = substitute(1)
                     ,jitter.fac = substitute(1)
                     ,with.lab = substitute(FALSE)
                     ,lab.pts = substitute(NULL)
                     ,lab.font = substitute(NULL)
                     ,alpha.trsp = substitute(NA)
                     ,which.lbs = substitute(NULL)
                     ,which.Order  = substitute(NULL)
                     ,return.Order = substitute(FALSE)
                     ,ylab.abs = substitute("absolute information")
                     ,ylab.rel= substitute("relative information")
                     ,adj = substitute(0.1)
                     ,cex.main = substitute(1.5)
                     ,cex.lab = substitute(1.5)
                     ,cex = substitute(1.5)
                     ,bty = substitute("n")
                     ,panel.first= substitute(grid())
                     ,col = substitute("blue")
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
  cl <- substitute(do.call(infoPlot,args0), list(args0=args))
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
fam = GParetoFamily()
IC <- optIC(model = fam, risk = asCov())
Y=distribution(fam)
data = r(Y)(1000)
dev.new()
infoPlotWrapper(IC, alpha.trsp=50, with.legend = FALSE)
dev.new()
infoPlotWrapper(IC, data, alpha.trsp=50, with.legend = FALSE)


# GEV
fam = GEVFamily()
IC <- optIC(model = fam, risk = asCov())
Y=distribution(fam)
data = r(Y)(1000)
dev.new()
infoPlotWrapper(IC, with.legend = TRUE, withCall = TRUE)
dev.new()
infoPlotWrapper(IC, data, with.legend = TRUE, withCall = TRUE)

# Gamma
fam = GammaFamily()
IC <- optIC(model = fam, risk = asCov())
Y=distribution(fam)
data = r(Y)(1000)
dev.new()
infoPlotWrapper(IC, alpha.trsp=70)
dev.new()
infoPlotWrapper(IC, data, alpha.trsp=70)

# Weibull
fam = WeibullFamily()
IC <- optIC(model = fam, risk = asCov())
Y=distribution(fam)
data = r(Y)(1000)
dev.new()
infoPlotWrapper(IC, alpha.trsp=50, with.legend = TRUE, withCall = FALSE)
dev.new()
infoPlotWrapper(IC, data, alpha.trsp=50, with.legend = TRUE, withCall = FALSE)


