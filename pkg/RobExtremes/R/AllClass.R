.onLoad <- function(lib, pkg){
#    require("methods", character = TRUE, quietly = TRUE)
}


.onAttach <- function(library, pkg)
{
#unlockBinding(".RobExtremesOptions", asNamespace("RobExtremes"))
#msga <- gettext("Note: Packages \"e1071\", \"moments\", \"fBasics\" should be attached ")
#msgb <- gettext("/before/ package \"distrEx\". See distrExMASK().")
buildStartupMessage(pkg = "RobExtremes", msga="", msgb="",
                    library = library, packageHelp = TRUE
#                    MANUAL="http://www.uni-bayreuth.de/departments/math/org/mathe7/DISTR/distr.pdf",
#VIGNETTE = gettext("Package \"distrDoc\" provides a vignette to this package as well as to several related packages; try vignette(\"distr\").")
)
  invisible()
}


.onUnload <- function(libpath)
{
    library.dynam.unload("RobExtremes", libpath)
}


#setClassUnion("ParamWithLocAndScaleAndShapeFamParameterUnion",
#               c("ParamWithScaleFamParameter",
#                 "ParamWithShapeFamParameter")
#         )

setClass("ParamWithLocAndScaleAndShapeFamParameter",
            contains = c("ParamWithScaleAndShapeFamParameter")
)


# parameter of Gumbel distribution
setClass("GumbelParameter", representation(loc = "numeric", 
                                           scale = "numeric"), 
            prototype(name = gettext("parameter of a Gumbel distribution"),
                      loc = 0, scale = 1),
            contains = "Parameter",
            validity = function(object){
                if(length(object@scale) != 1)
                    stop("length of 'scale' is not equal to 1")
                if(length(object@loc) != 1)
                    stop("length of 'loc' is not equal to 1")
                if(object@scale <= 0)
                    stop("'scale' has to be positive")
                else return(TRUE)
            })

# Gumbel distribution
setClass("Gumbel", 
            prototype = prototype(r = function(n){ evd::rgumbel(n, loc = 0, scale = 1) },
                                  d = function(x, log){ evd::dgumbel(x, loc = 0, scale = 1, log = FALSE) },
                                  p = function(q, lower.tail = TRUE, log.p = FALSE){ 
                                         p0 <- evd::pgumbel(q, loc = 0, scale = 1, lower.tail = lower.tail)
                                         if(log.p) return(log(p0)) else return(p0) 
                                  },
                                  q = function(p, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE){
                                      ## P.R.: changed to vectorized form 
                                      p1 <- if(log.p) exp(p) else p
                                                                                      
                                      in01 <- (p1>1 | p1<0)
                                      i01 <- .isEqual01(p1) 
                                      i0 <- (i01 & p1<1)   
                                      i1 <- (i01 & p1>0)
                                      ii01 <- .isEqual01(p1) | in01
                                                    
                                      p0 <- p
                                      p0[ii01] <- if(log.p) log(0.5) else 0.5
                                                    
                                      q1 <- evd::qgumbel(p0, loc = 0, scale = 1,
                                                    lower.tail = lower.tail) 
                                      q1[i0] <- if(lower.tail) -Inf else Inf
                                      q1[i1] <- if(!lower.tail) -Inf else Inf
                                      q1[in01] <- NaN
                                      
                                      return(q1)  
                                      },
                                  img = new("Reals"),
                                  param = new("GumbelParameter"),
                                  .logExact = FALSE,
                                  .lowerExact = TRUE),
            contains = "AbscontDistribution")


###### Pareto distribution by Nataliya Horbenko, ITWM, 18-03-09
## Class: ParetoParameter
setClass("ParetoParameter", 
          representation = representation(shape = "numeric",
                                          Min = "numeric"
                                          ), 
          prototype = prototype(shape = 1, Min = 1, name = 
                      gettext("Parameter of a Pareto distribution")
                      ), 
          contains = "Parameter"
          )

## Class: Pareto distribution
setClass("Pareto",  
          prototype = prototype(
                      r = function(n){ rpareto1(n, shape = 1, min = 1) },
                      d = function(x, log = FALSE){ 
                              dpareto1(x, shape = 1, min = 1, log = log) 
                                          },
                      p = function(q, lower.tail = TRUE, log.p = FALSE ){ 
                              ppareto1(q, shape = 1, min = 1, 
                                     lower.tail = lower.tail, log.p = log.p) 
                                          },
                      q = function(p, lower.tail = TRUE, log.p = FALSE ){ 
                        ## P.R.: changed to vectorized form 
                               p1 <- if(log.p) exp(p) else p
                                                                               
                               in01 <- (p1>1 | p1<0)
                               i01 <- .isEqual01(p1) 
                               i0 <- (i01 & p1<1)   
                               i1 <- (i01 & p1>0)
                               ii01 <- .isEqual01(p1) | in01
                                             
                               p0 <- p
                               p0[ii01] <- if(log.p) log(0.5) else 0.5
                                             
                               q1 <- qpareto1(p0, shape = 1,  min =  1, 
                                           lower.tail = lower.tail, log.p = log.p) 
                               q1[i0] <- if(lower.tail) -Inf else Inf
                               q1[i1] <- if(!lower.tail) -Inf else Inf
                               q1[in01] <- NaN
                               
                               return(q1)  
                            },
                      param = new("ParetoParameter"),
                      img = new("Reals"),
                      .logExact = TRUE,
                      .lowerExact = TRUE),
          contains = "AbscontDistribution"
          )

## Class: GParetoParameter
setClass("GParetoParameter", 
          representation = representation(loc = "numeric", scale = "numeric", shape = "numeric"
                                          ), 
          prototype = prototype(loc = 0, scale = 1, shape = 0, name = 
                      gettext("Parameter of a generalized Pareto distribution")
                      ), 
          contains = "Parameter"
          )
## Class: Generalized Pareto distribution
setClass("GPareto",  
          prototype = prototype(
                      r = function(n){ rgpd(n,loc = 0, scale = 1, shape = 1) },
                      d = function(x, log = FALSE){ 
                              dgpd(x, loc = 0, scale = 1, shape = 1, log = log) 
                                          },
                      p = function(q, lower.tail = TRUE, log.p = FALSE ){ 
                              p0 <- pgpd(q, loc = 0, scale = 1, shape = 1)
                              if(!lower.tail ) p0 <- 1-p0
                              if(log.p) p0 <- log(p0)
                              return(p0)},
                      q = function(p, lower.tail = TRUE, log.p = FALSE ){ 
                        ## P.R.: changed to vectorized form 
                               p1 <- if(log.p) exp(p) else p
                               if(!lower.tail) p1 <- 1-p1
                                                                               
                               in01 <- (p1>1 | p1<0)
                               i01 <- .isEqual01(p1) 
                               i0 <- (i01 & p1<1)   
                               i1 <- (i01 & p1>0)
                               ii01 <- .isEqual01(p1) | in01
                                             
                               p0 <- p
                               p0[ii01] <- if(log.p) log(0.5) else 0.5
                                             
                               q1 <- qgpd(p0,loc=0, scale = 1, shape = 1) 
                               q1[i0] <- if(lower.tail) -Inf else Inf
                               q1[i1] <- if(!lower.tail) -Inf else Inf
                               q1[in01] <- NaN
                               
                               return(q1)  
                            },
                      param = new("GParetoParameter"),
                      img = new("Reals"),
                      .withArith = FALSE,
                      .withSim = FALSE,
                      .logExact = TRUE,
                      .lowerExact = TRUE),
          contains = "AbscontDistribution"
          )


## Class: GEVParameter
setClass("GEVParameter", 
          representation = representation(loc = "numeric", scale = "numeric", shape = "numeric"
                                          ), 
          prototype = prototype(loc = 0, scale = 1, shape = 0.5, name = 
                      gettext("Parameter of a generalized extreme value distribution")
                      ), 
          contains = "Parameter"
          )
## Class: Generalized extreme value distribution
setClass("GEV",  
          prototype = prototype(
                      r = function(n){ rgev(n,loc = 0, scale = 1, shape = 0.5) },
                      d = function(x, log = FALSE){ 
                              dgev(x, loc = 0, scale = 1, shape = 0.5, log = log)
                                          },
                      p = function(q, lower.tail = TRUE, log.p = FALSE ){ 
                              p0 <- pgev(q, loc = 0, scale = 1, shape = 0.5)
                              if(!lower.tail ) p0 <- 1-p0
                              if(log.p) p0 <- log(p0)
                              return(p0)},
                      q = function(p, lower.tail = TRUE, log.p = FALSE ){ 
                        ## analogous to GPD
                               p1 <- if(log.p) exp(p) else p
                               if(!lower.tail) p1 <- 1-p1
                                                                               
                               in01 <- (p1>1 | p1<0)
                               i01 <- .isEqual01(p1) 
                               i0 <- (i01 & p1<1)   
                               i1 <- (i01 & p1>0)
                               ii01 <- .isEqual01(p1) | in01
                                             
                               p0 <- p
                               p0[ii01] <- if(log.p) log(0.5) else 0.5
                                             
                               q1 <- qgev(p0,loc=0, scale = 1, shape = 0.5) 
                               q1[i0] <- if(lower.tail) -Inf else Inf
                               q1[i1] <- if(!lower.tail) -Inf else Inf
                               q1[in01] <- NaN
                               
                               return(q1)  
                            },
                      param = new("GEVParameter"),
                      img = new("Reals"),
                      .withArith = FALSE,
                      .withSim = FALSE,
                      .logExact = TRUE,
                      .lowerExact = TRUE),
          contains = "AbscontDistribution"
          )
## Gumbel location family
setClass("GumbelLocationFamily",
          contains = "L2LocationFamily")

### for integration:
setClassUnion("DistributionsIntegratingByQuantiles",
               c("Weibull", "GEV", "GPareto", "Pareto", "Gammad"))


## models:
setClass("ParetoFamily", contains="L2ParamFamily")
setClass("GParetoFamily", contains="L2ScaleShapeUnion")
setClass("GEVFamily", contains="L2ScaleShapeUnion")
setClass("WeibullFamily", contains="L2ScaleShapeUnion")

## virtual in-between class for common parts in modifyModel - method
setClass("L2LocScaleShapeUnion", representation(locscaleshapename = "character"),
         contains = c("L2GroupParamFamily","VIRTUAL")
        )

setClass("GEVFamilyMuUnknown", contains="L2LocScaleShapeUnion")


setClass("LDEstimate",
         representation(location = "numeric",
                        dispersion = "numeric"
                        ),
         prototype(name = "LD estimate",
                   estimate = numeric(0),
                   samplesize = numeric(0),
                   completecases = logical(0),
                   asvar = NULL,
                   estimate.call = call("{}"),
                   location = 0,
                   dispersion = 1,
                   Infos = matrix(c(character(0),character(0)), ncol=2,
                                  dimnames=list(character(0), c("method", "message"))),
                   nuis.idx = NULL,
                   trafo = list(fct = function(x){
                                      list(fval = x, mat = matrix(1))},
                                mat = matrix(1))
                   ),
         contains = "Estimate")

setOldClass("gev.fit")
setOldClass("gpd.fit")
setClass("GPDEstimate", contains="Estimate")
setClass("GPDMCEstimate", contains=c("MCEstimate", "GPDEstimate"))
setClass("GPDLDEstimate", contains=c("LDEstimate", "GPDEstimate"))
setClass("GPDkStepEstimate", contains=c("kStepEstimate", "GPDEstimate"))
setClass("GPDORobEstimate", contains=c("ORobEstimate", "GPDkStepEstimate"))
setClass("GEVEstimate", contains="Estimate")
setClass("GEVLDEstimate", contains=c("LDEstimate", "GEVEstimate"))
setClass("GEVkStepEstimate", contains=c("kStepEstimate", "GEVEstimate"))
setClass("GEVORobEstimate", contains=c("ORobEstimate", "GEVkStepEstimate"))
setClass("GEVMCEstimate", contains=c("MCEstimate", "GEVEstimate"))
