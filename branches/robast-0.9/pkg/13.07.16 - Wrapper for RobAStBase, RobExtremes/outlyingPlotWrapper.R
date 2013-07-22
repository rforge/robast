##########################################
##                                      ## 
##    Wrapper for outlyingPlot.R        ##
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

## projection distance
qfun = function(x){p0 = p(X)(x); q0 = q(X)(p0)}
QProj <- function(){new("NormType", name="Quantiles", fct=qfun)}


##############################################################
#' Wrapper function for plot method for IC
#'  
#' The Wrapper takes most of arguments to the plot method
#' by default and gives a user possibility to run the 
#' function with low number of arguments
#' 
#' @param x data coercable to matrix. The data at which to produce the ddPlot
#' 
#' @param alpha confidence level for quantile (0 to 1)   
#' 
#' @param fam object of class L2ParamFamily
#' 
#' @param ... additional parameters (in particular to be passed on to \code{plot})
#'
#' @param alpha.trsp the transparency argument (0 to 100) for ploting the data
#' 
#' @param with.legend the flag for showing the legend of the plot
#'
#' @param withCall the flag for the call output 
#'
#' @usage outlyingPlotWrapper(x,alpha,fam,...,alpha.trsp = 100, with.legend = TRUE, withCall = TRUE)
#' 
#' @return produces an outlyingness plot based on distances applied to ICs. If withCall = TRUE, the call of the function is returned
#' 
#' @export
#' @docType function
#' @rdname outlyingPlotWrapper
#'
#' @import phylobase
#' @import vegan
#' @import igraph
#' @importFrom multtest mt.maxT
#' @importFrom multtest mt.minP
#'
#' @examples
#' # GPD
#' dev.new()
#' fam = GParetoFamily()
#' X=distribution(fam)
#' x = r(X)(1000)
#' outlyingPlotWrapper(x,alpha=0.99,fam=fam, main = "GPD", withCall = FALSE)
#' 
#' # GEV
#' dev.new()
#' fam = GEVFamily()
#' X=distribution(fam)
#' x = r(X)(1000)
#' outlyingPlotWrapper(x,alpha=0.95,fam=fam, main = "GEV", withCall = FALSE)
#' 
#' # Gamma
#' dev.new()
#' fam = GammaFamily()
#' X=distribution(fam)
#' x = r(X)(1000)
#' outlyingPlotWrapper(x,alpha=0.95,fam=fam, main = "Gamma", withCall = FALSE)
#' 
#' # Weibull
#' dev.new()
#' fam = WeibullFamily()
#' X=distribution(fam)
#' x = r(X)(1000)
#' outlyingPlotWrapper(x,alpha=0.95,fam=fam, main = "Weibull", withCall = FALSE)
##############################################################

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
  if(missing(x)) stop("Argument 'x' must be given as argument to 'outlyingPlotWrapper'")
  if(missing(alpha)) stop("Argument 'alpha' must be given as argument to 'outlyingPlotWrapper'")
  if(missing(fam)) stop("Argument 'fam' must be given as argument to 'outlyingPlotWrapper'")
  mc <- as.list(match.call(expand.dots = FALSE))[-1]
  dots <- mc$"..."
  if(missing(x)){
    alpha.trsp <- 100
  } else {
    if(is.null(mc$alpha.trsp)){
      alpha.trsp <- 30
      if(length(x) < 1000){
        alpha.trsp <- 50
      }
      if(length(x) < 100){
        alpha.trsp <- 100
      }
    }
  }
  if(is.null(mc$with.legend)) mc$with.legend <- TRUE
  if(is.null(mc$withCall)) mc$withCall <- TRUE

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
                   ,adj = substitute(0.5)
                   ,pch = substitute(21)
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
                   ,bty = substitute("o")
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
  
  args <- merge.lists(argsList, dots)
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
outlyingPlotWrapper(x,alpha=0.99,fam=fam, main = "GPD", withCall = FALSE)

# GEV
dev.new()
fam = GEVFamily()
X=distribution(fam)
x = r(X)(1000)
outlyingPlotWrapper(x,alpha=0.95,fam=fam, main = "GEV", withCall = FALSE)

# Gamma
dev.new()
fam = GammaFamily()
X=distribution(fam)
x = r(X)(1000)
outlyingPlotWrapper(x,alpha=0.95,fam=fam, main = "Gamma", withCall = FALSE)

# Weibull
dev.new()
fam = WeibullFamily()
X=distribution(fam)
x = r(X)(1000)
outlyingPlotWrapper(x,alpha=0.95,fam=fam, main = "Weibull", withCall = FALSE)
