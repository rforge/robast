
################################
##
## Class: ParetoParameter
##                             
################################

## Access Methods
setMethod("shape", "ParetoParameter", function(object) object@shape)
setMethod("Min", "ParetoParameter", function(object) object@Min)
setMethod("scale", "ParetoParameter",
           function(x, center = TRUE, scale = TRUE) x@Min)
### odd arg-list due to existing function in base package

## Replace Methods
setReplaceMethod("shape", "ParetoParameter", 
                  function(object, value){ 
                    if(length(value) != 1 || value <= 0)
                        stop("'value' has to be a single positive number")
                    object@shape <- value; object})
setReplaceMethod("Min", "ParetoParameter", 
                  function(object, value){ 
                    if(length(value) != 1 || value <= 0)
                        stop("'value' has to be a single positive number")
                    object@Min <- value; object})

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
setMethod("scale", "Pareto",
           function(x, center = TRUE, scale = TRUE) Min(param(x)))

## wrapped replace methods
setMethod("shape<-", "Pareto", function(object, value)
           new("Pareto", shape = value, Min = Min(object)))
setMethod("Min<-", "Pareto", function(object, value) 
           new("Pareto", shape = shape(object), Min = value))

setMethod("*", c("Pareto","numeric"),
          function(e1, e2){
            if (length(e2)>1) stop("length of operator must be 1")
            if (isTRUE(all.equal(e2,0)))
                return(new("Dirac", location = 0, .withArith = TRUE))
            Pareto <- new("Pareto", Min=Min(e1)* abs(e2),
                                 shape=shape(e1))
            if(e2<0) Pareto <- (-1)*as(Pareto,"AbscontDistribution")
            return(Pareto)
          })
