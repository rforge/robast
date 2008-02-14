.onLoad <- function(lib, pkg){
    require("methods", character = TRUE, quietly = TRUE)
    require("distr", character = TRUE, quietly = TRUE)
    require("distrEx", character = TRUE, quietly = TRUE)
    require("distrMod", character = TRUE, quietly = TRUE)
    require("RandVar", character = TRUE, quietly = TRUE)
}

# risks (e.g., risk of estimator)
setClass("RiskType", representation(type = "character"), contains = "VIRTUAL")
# asymptotic risk
setClass("asRisk", contains = c("RiskType", "VIRTUAL"))
# asymptotic covariance
setClass("asCov", contains = "asRisk", 
            prototype = prototype(type = "asymptotic covariance"))
# trace of asymptotic covariance
setClass("trAsCov", contains = "asRisk", 
            prototype = prototype(type = "trace of asymptotic covariance"))
# asymptotic Hampel risk
setClass("asHampel", representation(bound = "numeric"), 
            prototype = prototype(bound = Inf, 
                             type = "trace of asymptotic covariance for given bias bound"),
            contains = "asRisk", 
            validity = function(object){
                if(any(object@bound <= 0))
                    stop("'bound' has to be positive")
                else TRUE
            })
# asymptotic bias
setClass("asBias", contains = "asRisk", 
            prototype = prototype(type = "asymptotic bias"))

# convex asymptotic risk
setClass("asGRisk", contains = c("asRisk", "VIRTUAL")) 
# asymptotic mean square error
setClass("asMSE", contains = "asGRisk", 
            prototype = prototype(type = "asymptotic mean square error"))
# asymptotic under-/overshoot probability
setClass("asUnOvShoot", representation(width = "numeric"), 
            prototype = prototype(type = "asymptotic under-/overshoot probability"),
            contains = "asGRisk",
            validity = function(object){
                if(length(object@width) != 1)
                    stop("length of 'width' has to be 1")
                if(any(object@width <= 0))
                    stop("'width' has to be positive")
                else TRUE
            })
# finite-sample risk
setClass("fiRisk", contains = c("RiskType", "VIRTUAL"))
# finite-sample covariance
setClass("fiCov", contains = "fiRisk", 
            prototype = prototype(type = "finite-sample covariance"))
# trace of finite-sample covariance
setClass("trFiCov", contains = "fiRisk", 
            prototype = prototype(type = "trace of finite-sample covariance"))
# finite-sample Hampel risk
setClass("fiHampel", representation(bound = "numeric"),
            prototype = prototype(bound = Inf, 
                             type = "finite-sample variance for given bias bound"),
            contains = "fiRisk", 
            validity = function(object){
                if(any(object@bound <= 0))
                    stop("'bound' has to be positive")
                else TRUE
            })
# finite-sample mean square error
setClass("fiMSE", contains = "fiRisk", 
            prototype = prototype(type = "finite-sample mean square error"))
# finite-sample bias
setClass("fiBias", contains = "fiRisk", 
            prototype = prototype(type = "finite-sample bias"))
# finite-sample under-/overshoot probability
setClass("fiUnOvShoot", representation(width = "numeric"), 
            prototype = prototype(type = "finite-sample under-/overshoot probability"),
            contains = "fiRisk",
            validity = function(object){
                if(length(object@width) != 1)
                    stop("length of 'width' has to be 1")
                if(any(object@width <= 0))
                    stop("'width' has to be positive")
                else TRUE
            })

################################################################################
# Bias Classes
################################################################################

### from session 10-01-08 : class Bias type
setClass("BiasType", representation(name = "character"),
          contains = "VIRTUAL")

setClass("symmetricBiasType", prototype = prototype(name = "symmetric bias"),
          contains = "BiasType")

setClass("onesidedBiasType", representation(sign = "numeric"),
          prototype = prototype(name = "positive bias", sign = 1),
          contains = "BiasType")

setClass("asymmetricBiasType", 
          representation(nu = "numeric"), ### weights acc. to paper
          prototype = prototype(name = "asymmetric bias", nu = c(1,1)),
          contains = "BiasType")
