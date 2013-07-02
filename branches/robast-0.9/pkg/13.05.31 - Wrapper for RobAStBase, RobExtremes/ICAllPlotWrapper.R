##########################################
##                                      ## 
##    Wrapper for AllPlot.R             ##
##    (plot method for IC)              ##
##                                      ##
##########################################

##IC - influence curve
##y - dataset
## with.legend - optional legend indicator
## withCall - optional indicator of the function call
#
ICAllPlotWrapper = function(IC, y,...,alpha.trsp = 100, with.legend = TRUE, rescale = FALSE ,withCall = TRUE){
  ###
  ### 1. grab the dots (and manipulate it within the wrapper function)
  ###
  ###
  ### do something to fix the good default arguments
  ###
  if(missing(IC)) stop("Argument 'IC' must be given as argument to 'ICAllPlotWrapper'")
  mc <- as.list(match.call(expand.dots = FALSE))[-1]
  dots <- mc$"..."
  if(missing(y)){
    alpha.trsp <- 100
  } else {
    if(is.null(mc$alpha.trsp)){
      alpha.trsp <- 30
      if(length(y) < 1000){
        alpha.trsp <- 50
      }
      if(length(y) < 100){
        alpha.trsp <- 100
      }
    }
  }
  if(is.null(mc$with.legend)) mc$with.legend <- TRUE
  if(is.null(mc$rescale)) mc$rescale <- FALSE
  if(is.null(mc$withCall)) mc$withCall <- TRUE
  ###
  ### 2. build up the argument list for the (powerful/fullfledged)
  ### graphics/diagnostics function;
  ##
  
  ## For GEVFamily we have to scale the axes
  if((mc$rescale) & (as.list(IC@CallL2Fam)[[1]] == "GEVFamily")){
    if(missing(y)){
      scaleList <- list(scaleX = substitute(TRUE)
                      ,scaleX.fct = substitute(p(eval(IC@CallL2Fam)))
                      ,scaleX.inv = substitute(q(eval(IC@CallL2Fam)))
                      ,scaleY = substitute(TRUE)
                      ,scaleY.fct = substitute(pnorm)
                      ,scaleY.inv = substitute(qnorm)
                      ,x.ticks = substitute(NULL)
                      ,y.ticks = substitute(NULL)
      )
    } else {
      scaleList <- list(scaleX = substitute(TRUE)
                        ,scaleY = substitute(TRUE)
      )
    }
  }else{
    if(missing(y)){
      scaleList <- list(scaleX = substitute(FALSE)
                      ,scaleX.fct = substitute(p(eval(IC@CallL2Fam)))
                      ,scaleX.inv = substitute(q(eval(IC@CallL2Fam)))
                      ,scaleY = substitute(FALSE)
                      ,scaleY.fct = substitute(pnorm)
                      ,scaleY.inv=substitute(qnorm)
                      ,x.ticks = substitute(NULL)
                      ,y.ticks = substitute(NULL)
      )
    } else {
      scaleList <- list(scaleX = substitute(FALSE)
                        ,scaleY = substitute(FALSE)
      )
    }
  }

  
  if(missing(y)){
    argsList <- c(list(x = substitute(IC)
                     ,withSweave = substitute(getdistrOption("withSweave"))
                     ,main = substitute(FALSE)
                     ,inner = substitute(TRUE)
                     ,sub = substitute(FALSE)
                     ,col.inner = substitute(par("col.main"))
                     ,cex.inner = substitute(0.8)
                     ,bmar = substitute(par("mar")[1])
                     ,tmar = substitute(par("mar")[3])
                     ,with.legend = substitute(FALSE)
                     ,legend = substitute(NULL)
                     ,legend.bg = substitute("white")
                     ,legend.location = substitute("bottomright")
                     ,legend.cex = substitute(0.8)
                     ,withMBR = substitute(FALSE)
                     ,MBRB = substitute(NA)
                     ,MBR.fac = substitute(2)
                     ,col.MBR = substitute(par("col"))
                     ,lty.MBR = substitute("dashed")
                     ,lwd.MBR = substitute(0.8)
                     ,scaleN = substitute(9)
                     ,mfColRow = substitute(TRUE)
                     ,to.draw.arg = substitute(NULL)
                     ,adj = substitute(0.5)
                     ,cex.main = substitute(1.5)
                     ,cex.lab = substitute(1.5)
                     ,cex = substitute(1.5)
                     ,bty = substitute("o")
                     ,panel.first= substitute(grid())
                     ,col = substitute("blue")
    ), scaleList)
  }else{
    argsList <- c(list(x = substitute(IC)
                     ,y = substitute(y)
                     ,cex.pts = substitute(0.3)
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
                     ,scaleN = substitute(9)
                     ,adj = substitute(0.5)
                     ,cex.main = substitute(1.5)
                     ,cex.lab = substitute(1.5)
                     ,cex = substitute(1.5)
                     ,bty = substitute("o")
                     ,panel.first= substitute(grid())
                     ,col = substitute("blue")
    ), scaleList)
  }


  
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
  cl <- substitute(do.call(plot,args0), list(args0=args))
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
y = r(Y)(1000)
dev.new()
ICAllPlotWrapper(IC, with.legend = FALSE)
dev.new()
ICAllPlotWrapper(IC, y, with.legend = FALSE)

# GEV
fam = GEVFamily()
IC <- optIC(model = fam, risk = asCov())
Y=distribution(fam)
y = r(Y)(1000)
dev.new()
ICAllPlotWrapper(IC, alpha.trsp=100, with.legend = TRUE, rescale = TRUE, withCall = TRUE)
dev.new()
ICAllPlotWrapper(IC, y, alpha.trsp=100, with.legend = TRUE, rescale = TRUE, withCall = TRUE)

# Gamma
fam = GammaFamily()
IC <- optIC(model = fam, risk = asCov())
Y=distribution(fam)
y = r(Y)(1000)
dev.new()
ICAllPlotWrapper(IC)
dev.new()
ICAllPlotWrapper(IC, y)

# Weibull
fam = WeibullFamily()
IC <- optIC(model = fam, risk = asCov())
Y=distribution(fam)
y = r(Y)(1000)
dev.new()
ICAllPlotWrapper(IC, alpha.trsp=30, with.legend = TRUE, withCall = FALSE)
dev.new()
ICAllPlotWrapper(IC, y, alpha.trsp=30, with.legend = TRUE, withCall = FALSE)

