.onLoad <- function(lib, pkg){
#    require("methods", character = TRUE, quietly = TRUE)
#    require("distr", character = TRUE, quietly = TRUE)
#    require("distrEx", character = TRUE, quietly = TRUE)
#    require("RandVar", character = TRUE, quietly = TRUE)
#    require("distrMod", character = TRUE, quietly = TRUE)
#    require("RobAStBase", character = TRUE, quietly = TRUE)
}


## asymptotic Anscombe risk
setClass("asAnscombe", representation(eff = "numeric"),
            prototype = prototype(eff = .95,
                             type = "optimal bias robust IC for given ARE in the ideal model"),
            contains = "asRiskwithBias",
            validity = function(object){
                if(any(object@eff <= 0|object@eff > 1))
                    stop("'eff' has to be in (0,1]")
                else TRUE
            })


## asymptotic L4 error
setClass("asL4", contains = "asGRisk",
            prototype = prototype(type = "asymptotic mean power 4 error"))
## asymptotic L1 error
setClass("asL1", contains = "asGRisk",
            prototype = prototype(type = "asymptotic mean absolute error"))

setClass("ORobEstimate",
         representation(roptestCall = "OptionalCall"),
         prototype(name = "Optimally robust asymptotically linear estimate",
                   estimate = numeric(0),
                   samplesize = numeric(0),
                   completecases = logical(0),
                   estimate.call = call("{}"),
                   steps = integer(0),
                   asvar = NULL,
                   asbias = NULL,
                   pIC = NULL,
                   pICList = NULL,
                   ICList = NULL,
                   ksteps = NULL,
                   uksteps = NULL,
                   start = matrix(0),
                   startval = matrix(0),
                   ustartval = matrix(0),
                   nuis.idx = NULL,
                   trafo = list(fct = function(x){
                                      list(fval = x, mat = matrix(1))},
                                mat = matrix(1)), ### necessary for comparison with unit matrix
                   Infos = matrix(c(character(0),character(0)), ncol=2,
                                  dimnames=list(character(0), c("method", "message"))),
                   untransformed.estimate = NULL,
                   untransformed.asvar = NULL,
                   robestCall = NULL,
                   roptestCall = NULL),
         contains = "kStepEstimate")
