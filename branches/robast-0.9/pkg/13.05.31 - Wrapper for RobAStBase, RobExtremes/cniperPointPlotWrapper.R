##########################################
##                                      ## 
##    Wrapper for cniperPointPlot.R     ##
##                                      ##
##                                      ##
##########################################

### aditional function
merge.lists <- function(a, b){
  a.names <- names(a)
  b.names <- names(b)
  m.names <- sort(unique(c(a.names, b.names), fromLast = TRUE))
  sapply(m.names, function(i) {
    if (is.list(a[[i]]) & is.list(b[[i]])) merge.lists(a[[i]], b[[i]])
    else if (i %in% b.names) b[[i]]
    else a[[i]]
  }, simplify = FALSE)
}

##@fam - parameter family
## lower - left point of the x-axis
## upper - right point of the x-axis
## alpha.trsp - optional transparency of the plot
## with.legend - optional legend indicator
## withCall - optional indicator of the function call
#
cniperPointPlotWrapper = function(fam,...
                                  ,lower = getdistrOption("DistrResolution"), upper=1-getdistrOption("DistrResolution")
                                  ,with.legend = TRUE, withCall = TRUE){
  ###
  ### 1. grab the dots (and manipulate it within the wrapper function)
  ###
  ###
  ### do something to fix the good default arguments
  ###
  mc <- as.list(match.call(expand.dots = FALSE))[-1]
  dots <- mc$"..."
  if(is.null(mc$lower)) lower <- getdistrOption("DistrResolution")
  if(is.null(mc$upper)) upper <- 1-getdistrOption("DistrResolution")
  if(is.null(mc$with.legend)) mc$with.legend <- TRUE
  if(is.null(mc$withCall)) mc$withCall <- TRUE
  if(missing(fam)) stop("Argument 'fam' must be given as argument to 'cniperPointPlotWrapper'")
  ###
  ### 2. build up the argument list for the (powerful/fullfledged)
  ### graphics/diagnostics function;
  ##  

  argsList <- list(L2Fam = substitute(fam)
                   ,data = substitute(NULL)
                   ,neighbor = substitute(ContNeighborhood(radius = 0.5))
                   ,risk = substitute(asMSE())
                   ,lower = substitute(lower)
                   ,upper = substitute(upper)
                   ,n = substitute(101)
                   ,withMaxRisk = substitute(TRUE)
                   ,scaleX = substitute(FALSE)
                   ,scaleX.fct = substitute(p(fam))
                   ,scaleX.inv = substitute(q(fam))
                   ,scaleY = substitute(FALSE)
                   ,scaleY.fct = substitute(pnorm)
                   ,scaleY.inv = substitute(qnorm)
                   ,scaleN = substitute(9)
                   ,x.ticks = substitute(NULL)
                   ,y.ticks = substitute(NULL)
                   ,cex.pts = substitute(1)
                   ,col.pts = substitute(par("col"))
                   ,pch.pts = substitute(1)
                   ,jitter.fac = substitute(1)
                   ,with.lab = substitute(FALSE)
                   ,lab.pts = substitute(NULL)
                   ,lab.font = substitute(NULL)
                   ,alpha.trsp = substitute(NA)
                   ,which.lbs = substitute(NULL)
                   ,which.Order  = substitute(NULL)
                   ,return.Order = substitute(FALSE)
                   ,adj = 0.5
                   ,cex.main = substitute(1.5)
                   ,cex.lab = substitute(1.5)
                   ,main = ""#"Outlyingness Plot"
                   ,xlab=substitute("Dirac point")
                   ,ylab=substitute("Asymptotic Risk difference (classic - robust)")
                   ,bty = substitute("o")
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
  
  args <- merge.lists(argsList, dots)
  ###
  ### 3. build up the call but grab it and write it into an object
  ###
  cl <- substitute(do.call(cniperPointPlot,args0), list(args0=args))
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

# WRite the correct path to the modified file cniperCont.R from the ROptEst package
source("D:/Dropbox/My Mathematics/Researches Misha/Current Research/11.06 - KL PhD/PhD Thesis/Reports for Project/13.07.16 - Wrapper for RobAStBase, RobExtremes/cniperCont.R")

# GPD
dev.new()
fam = GParetoFamily()
cniperPointPlotWrapper(fam=fam, main = "GPD", lower = 0, upper = 10, withCall = FALSE)

# GEV
dev.new()
fam = GEVFamily()
cniperPointPlotWrapper(fam=fam, main = "GEV", lower = 0, upper = 5, withCall = FALSE)

# Gamma
dev.new()
fam = GammaFamily()
cniperPointPlotWrapper(fam=fam, main = "Gamma", lower = 0, upper = 5, withCall = FALSE)

# Weibull
dev.new()
fam = WeibullFamily()
cniperPointPlotWrapper(fam=fam, main = "Weibull", withCall = FALSE)



