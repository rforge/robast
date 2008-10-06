###############################################################################
## Functions and methods for "ALEstimate" classes and subclasses
###############################################################################

setMethod("pIC", "ALEstimate", function(object) object@pIC)
setMethod("asbias", "ALEstimate", function(object) object@asbias)
setMethod("steps", "kStepEstimate", function(object) object@steps)
setMethod("Mroot", "MEstimate", function(object) object@Mroot)
