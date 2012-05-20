
################################
##
## Class: ParetoParameter
##                             
################################

## Access Methods
setMethod("shape", "ParetoParameter", function(object) object@shape)
setMethod("Min", "ParetoParameter", function(object) object@Min)

## Replace Methods
setReplaceMethod("shape", "ParetoParameter", 
                  function(object, value){ object@shape <- value; object})
setReplaceMethod("Min", "ParetoParameter", 
                  function(object, value){ object@Min <- value; object})

setValidity("ParetoParameter", function(object){
  if(length(shape(object)) != 1)
    stop("shape has to be a numeric of length 1")    
  if(shape(object) <= 0)
    stop("shape has to be positive")
  if(length(Min(object)) != 1)
    stop("Min has to be a numeric of length 1")    
  if(Min(object) <= 0)
    stop("Min has to be positive")
  else return(TRUE)
})

################################
##            .Object@img <- new("Naturals")

## Class: Pareto distribution
##
################################

Pareto <- function(shape = 1, Min = 1) 
               new("Pareto", shape = shape, Min = Min)

## wrapped access methods
setMethod("shape", "Pareto", function(object) shape(param(object)))
setMethod("Min", "Pareto", function(object) Min(param(object)))

## wrapped replace methods
setMethod("shape<-", "Pareto", function(object, value) 
           new("Pareto", shape = value, Min = Min(object)))
setMethod("Min<-", "Pareto", function(object, value) 
           new("Pareto", shape = shape(object), Min = value))

