### * <HEADER>
###
attach(NULL, name = "CheckExEnv")
assign("nameEx", 
       local({
	   s <- "__{must remake R-ex/*.R}__"
           function(new) {
               if(!missing(new)) s <<- new else s
           }
       }),
       pos = "CheckExEnv")
## Add some hooks to label plot pages for base and grid graphics
assign("base_plot_hook",
       function() {
           pp <- par(c("mfg","mfcol","oma","mar"))
           if(all(pp$mfg[1:2] == c(1, pp$mfcol[2]))) {
               outer <- (oma4 <- pp$oma[4]) > 0; mar4 <- pp$mar[4]
               mtext(sprintf("help(\"%s\")", nameEx()), side = 4,
                     line = if(outer)max(1, oma4 - 1) else min(1, mar4 - 1),
              outer = outer, adj = 1, cex = .8, col = "orchid", las=3)
           }
       },
       pos = "CheckExEnv")
assign("grid_plot_hook",
       function() {
           pushViewport(viewport(width=unit(1, "npc") - unit(1, "lines"),
                                 x=0, just="left"))
           grid.text(sprintf("help(\"%s\")", nameEx()),
                     x=unit(1, "npc") + unit(0.5, "lines"),
                     y=unit(0.8, "npc"), rot=90,
                     gp=gpar(col="orchid"))
       },
       pos = "CheckExEnv")
setHook("plot.new",     get("base_plot_hook", pos = "CheckExEnv"))
setHook("persp",        get("base_plot_hook", pos = "CheckExEnv"))
setHook("grid.newpage", get("grid_plot_hook", pos = "CheckExEnv"))
assign("cleanEx",
       function(env = .GlobalEnv) {
	   rm(list = ls(envir = env, all.names = TRUE), envir = env)
           RNGkind("default", "default")
	   set.seed(1)
   	   options(warn = 1)
	   .CheckExEnv <- as.environment("CheckExEnv")
	   delayedAssign("T", stop("T used instead of TRUE"),
		  assign.env = .CheckExEnv)
	   delayedAssign("F", stop("F used instead of FALSE"),
		  assign.env = .CheckExEnv)
	   sch <- search()
	   newitems <- sch[! sch %in% .oldSearch]
	   for(item in rev(newitems))
               eval(substitute(detach(item), list(item=item)))
	   missitems <- .oldSearch[! .oldSearch %in% sch]
	   if(length(missitems))
	       warning("items ", paste(missitems, collapse=", "),
		       " have been removed from the search path")
       },
       pos = "CheckExEnv")
assign("ptime", proc.time(), pos = "CheckExEnv")
grDevices::postscript("ROptEst-Ex.ps")
assign("par.postscript", graphics::par(no.readonly = TRUE), pos = "CheckExEnv")
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"), pager="console")
options(warn = 1)    
library('ROptEst')

assign(".oldSearch", search(), pos = 'CheckExEnv')
assign(".oldNS", loadedNamespaces(), pos = 'CheckExEnv')
cleanEx(); nameEx("getL1normL2deriv")
### * getL1normL2deriv

flush(stderr()); flush(stdout())

### Name: getL1normL2deriv
### Title: Calculation of L1 norm of L2derivative
### Aliases: getL1normL2deriv getL1normL2deriv-methods
###   getL1normL2deriv,UnivariateDistribution-method
###   getL1normL2deriv,RealRandVariable-method


### ** Examples

##



cleanEx(); nameEx("getL2normL2deriv")
### * getL2normL2deriv

flush(stderr()); flush(stdout())

### Name: getL2normL2deriv
### Title: Calculation of L2 norm of L2derivative
### Aliases: getL2normL2deriv


### ** Examples

##



cleanEx(); nameEx("leastFavorableRadius")
### * leastFavorableRadius

flush(stderr()); flush(stdout())

### Name: leastFavorableRadius
### Title: Generic Function for the Computation of Least Favorable Radii
### Aliases: leastFavorableRadius leastFavorableRadius-methods
###   leastFavorableRadius,L2ParamFamily,UncondNeighborhood,asGRisk-method


### ** Examples

N <- NormLocationFamily(mean=0, sd=1) 
leastFavorableRadius(L2Fam=N, neighbor=ContNeighborhood(),
                     risk=asMSE(), rho=0.5)



cleanEx(); nameEx("lowerCaseRadius")
### * lowerCaseRadius

flush(stderr()); flush(stdout())

### Name: lowerCaseRadius
### Title: Computation of the lower case radius
### Aliases: lowerCaseRadius lowerCaseRadius-methods
###   lowerCaseRadius,L2ParamFamily,ContNeighborhood,asMSE,BiasType-method
###   lowerCaseRadius,L2ParamFamily,TotalVarNeighborhood,asMSE,BiasType-met
###   hod


### ** Examples

lowerCaseRadius(BinomFamily(size = 10), ContNeighborhood(), asMSE())
lowerCaseRadius(BinomFamily(size = 10), TotalVarNeighborhood(), asMSE())



cleanEx(); nameEx("optIC")
### * optIC

flush(stderr()); flush(stdout())

### Name: optIC
### Title: Generic function for the computation of optimally robust ICs
### Aliases: optIC optIC-methods optIC,L2ParamFamily,asCov-method
###   optIC,InfRobModel,asRisk-method optIC,InfRobModel,asUnOvShoot-method
###   optIC,FixRobModel,fiUnOvShoot-method


### ** Examples

B <- BinomFamily(size = 25, prob = 0.25) 

## classical optimal IC
IC0 <- optIC(model = B, risk = asCov())
plot(IC0) # plot IC
checkIC(IC0, B)



cleanEx(); nameEx("optRisk")
### * optRisk

flush(stderr()); flush(stdout())

### Name: optRisk
### Title: Generic function for the computation of the minimal risk
### Aliases: optRisk optRisk-methods optRisk,L2ParamFamily,asCov-method
###   optRisk,InfRobModel,asRisk-method
###   optRisk,FixRobModel,fiUnOvShoot-method


### ** Examples

optRisk(model = NormLocationScaleFamily(), risk = asCov())



cleanEx(); nameEx("radiusMinimaxIC")
### * radiusMinimaxIC

flush(stderr()); flush(stdout())

### Name: radiusMinimaxIC
### Title: Generic function for the computation of the radius minimax IC
### Aliases: radiusMinimaxIC radiusMinimaxIC-methods
###   radiusMinimaxIC,L2ParamFamily,UncondNeighborhood,asGRisk-method


### ** Examples

N <- NormLocationFamily(mean=0, sd=1) 
radiusMinimaxIC(L2Fam=N, neighbor=ContNeighborhood(), 
                risk=asMSE(), loRad=0.1, upRad=0.5)



cleanEx(); nameEx("trAsCov-class")
### * trAsCov-class

flush(stderr()); flush(stdout())

### Name: trAsCov-class
### Title: Trace of asymptotic covariance
### Aliases: trAsCov-class
### Keywords: classes

### ** Examples

new("trAsCov")



cleanEx(); nameEx("trFiCov-class")
### * trFiCov-class

flush(stderr()); flush(stdout())

### Name: trFiCov-class
### Title: Trace of finite-sample covariance
### Aliases: trFiCov-class
### Keywords: classes

### ** Examples

new("trFiCov")



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
