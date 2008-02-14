.onLoad <- function(lib, pkg){
    require("methods", character = TRUE, quietly = TRUE)
    require("distr", character = TRUE, quietly = TRUE)
    require("distrEx", character = TRUE, quietly = TRUE)
    require("distrMod", character = TRUE, quietly = TRUE)
    require("RandVar", character = TRUE, quietly = TRUE)
}

# neighborhood
setClass("Neighborhood",
            representation(type = "character",
                           radius = "numeric"), 
            contains = "VIRTUAL")
# unconditional (errors-in-variables) neighborhood
setClass("UncondNeighborhood", contains = c("Neighborhood", "VIRTUAL"))
# unconditional convex contamination neighborhood
setClass("ContNeighborhood", contains = "UncondNeighborhood",
            prototype = prototype(type = "(uncond.) convex contamination neighborhood",
                                  radius = 0))
# unconditional total variation neighborhood
setClass("TotalVarNeighborhood", contains = "UncondNeighborhood",
            prototype = prototype(type = "(uncond.) total variation neighborhood",
                                  radius = 0))
# robust model
setClass("RobModel",
            representation(center = "ProbFamily",
                           neighbor = "Neighborhood"),
            contains = "VIRTUAL")
# robust model with fixed (unconditional) neighborhood
setClass("FixRobModel",
            prototype = prototype(center = new("ParamFamily"),
                                  neighbor = new("ContNeighborhood")),
            contains = "RobModel",
            validity = function(object){
                if(!is(object@neighbor, "UncondNeighborhood"))
                    stop("'neighbor' is no unconditional neighborhood")
                if(any(object@neighbor@radius < 0 || object@neighbor@radius > 1))
                    stop("neighborhood radius has to be in [0, 1]")
                else return(TRUE)
            })
# robust model with infinitesimal (unconditional) neighborhood
setClass("InfRobModel",
            prototype = prototype(center = new("L2ParamFamily"),
                                  neighbor = new("ContNeighborhood")),
            contains = "RobModel",
            validity = function(object){
                if(!is(object@neighbor, "UncondNeighborhood"))
                    stop("'neighbor' is no unconditional neighborhood")
                if(any(object@neighbor@radius < 0))
                    stop("'radius' has to be in [0, Inf]")
                else return(TRUE)
            })
# Influence curve/function with domain: EuclideanSpace
setClass("InfluenceCurve", 
            representation(name = "character", 
                           Curve = "EuclRandVarList", 
                           Risks = "list",
                           Infos = "matrix"),
            validity = function(object){
                if(!is(Domain(object@Curve[[1]]), "EuclideanSpace"))
                    stop("The domain of 'Curve' has to be a Euclidean space")
                if(!is.character(object@Infos))
                    stop("'Infos' contains no matrix of characters")
                for(char in names(object@Risks))
                    if(!extends(char, "RiskType"))
                        stop(paste(char, "is no valid 'RiskType'"))
                if(ncol(object@Infos)!=2)
                    stop("'Infos' must have two columns")
                else TRUE
            })
# partial incluence curve
setClass("IC", representation(CallL2Fam = "call"),
            prototype(name = "square integrable (partial) influence curve",
                      Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x}), 
                                               Domain = Reals())),
                      Risks = list(),
                      Infos = matrix(c(character(0),character(0)), ncol=2,
                                dimnames=list(character(0), c("method", "message"))),
                      CallL2Fam = call("L2ParamFamily")),
            contains = "InfluenceCurve",
            validity = function(object){
                L2Fam <- eval(object@CallL2Fam)
                trafo <- L2Fam@param@trafo
                if(nrow(trafo) != dimension(object@Curve))
                    stop("wrong dimension of 'Curve'")
                if(dimension(Domain(L2Fam@L2deriv[[1]])) != dimension(Domain(object@Curve[[1]])))
                    stop("dimension of 'Domain' of 'L2deriv' != dimension of 'Domain' of 'Curve'")

                return(TRUE)
            })
# (partial) influence curve of contamination type
setClass("ContIC", 
            representation(clip = "numeric",
                           cent = "numeric",
                           stand = "matrix",
                           lowerCase = "OptionalNumeric",
                           neighborRadius = "numeric"), 
            prototype(name = "IC of contamination type",
                      Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x}), 
                                                    Domain = Reals())),
                      Risks = list(),
                      Infos = matrix(c(character(0),character(0)), ncol=2,
                                dimnames=list(character(0), c("method", "message"))),
                      CallL2Fam = call("L2ParamFamily"),
                      clip = Inf, cent = 0, stand = as.matrix(1),
                      lowerCase = NULL,
                      neighborRadius = 0),
            contains = "IC",
            validity = function(object){
                if(any(object@neighborRadius < 0)) # radius vector?!
                    stop("'neighborRadius' has to be in [0, Inf]")
                if(length(object@cent) != nrow(object@stand))
                    stop("length of centering constant != nrow of standardizing matrix")
                if((length(object@clip) != 1) && (length(object@clip) != length(object@Curve)))
                    stop("length of clipping bound != 1 and != length of 'Curve'")
                if(!is.null(object@lowerCase))
                    if(length(object@lowerCase) != nrow(object@stand))
                        stop("length of 'lowerCase' != nrow of standardizing matrix")
                L2Fam <- eval(object@CallL2Fam)
                if(!identical(dim(L2Fam@param@trafo), dim(object@stand)))
                    stop(paste("dimension of 'trafo' of 'param' != dimension of 'stand'"))
                return(TRUE)
            })
# (partial) influence curve of total variation type
setClass("TotalVarIC",
            representation(clipLo = "numeric",
                           clipUp = "numeric",
                           stand = "matrix",
                           lowerCase = "OptionalNumeric",
                           neighborRadius = "numeric"),
            prototype(name = "IC of total variation type",
                      Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x}),
                                                               Domain = Reals())),
                      Risks = list(),
                      Infos = matrix(c(character(0),character(0)), ncol=2,
                                dimnames=list(character(0), c("method", "message"))),
                      CallL2Fam = call("L2ParamFamily"),
                      clipLo = -Inf, clipUp = Inf, stand = as.matrix(1),
                      lowerCase = NULL,
                      neighborRadius = 0),
            contains = "IC",
            validity = function(object){
                if(any(object@neighborRadius < 0)) # radius vector?!
                    stop("'neighborRadius' has to be in [0, Inf]")
                if((length(object@clipLo) != 1) && (length(object@clipLo) != length(object@Curve)))
                    stop("length of lower clipping bound != 1 and != length of 'Curve'")
                if((length(object@clipLo) != 1) && (length(object@clipLo) != length(object@Curve)))
                    stop("length of upper clipping bound != 1 and != length of 'Curve'")
                L2Fam <- eval(object@CallL2Fam)
                if(!identical(dim(L2Fam@param@trafo), dim(object@stand)))
                    stop(paste("dimension of 'trafo' of 'param' != dimension of 'stand'"))
                return(TRUE)
            })
