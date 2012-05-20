###########################################
## Class: Generalized Pareto Distribution  
##  
## @param: location, scale, shape 
###########################################

## access methods
setMethod("location", "GParetoParameter", function(object) object@loc)
setMethod("loc", "GParetoParameter", function(object) object@loc)
setMethod("scale", "GParetoParameter", 
           function(x, center = TRUE, scale = TRUE) x@scale)
           ### odd arg-list due to existing function in base package 
setMethod("shape", "GParetoParameter", function(object) object@shape)

## replace Methods
setReplaceMethod("loc", "GParetoParameter", 
    function(object, value){ object@loc <- value; object })
setReplaceMethod("location", "GParetoParameter", 
    function(object, value){ object@loc <- value; object })
setReplaceMethod("scale", "GParetoParameter", 
    function(object, value){
        if(length(value) != 1 || value <= 0)
            stop("'value' has to be a single positive real number!")
        object@scale <- value; object})
setReplaceMethod("shape", "GParetoParameter", 
    function(object, value){ object@shape <- value; object})

## wrapped access methods
setMethod("location", "GPareto", function(object) loc(object@param))
setMethod("loc", "GPareto", function(object) loc(object@param))
setMethod("scale", "GPareto", 
           function(x, center = TRUE, scale = TRUE)
           scale(x@param))
setMethod("shape", "GPareto", function(object) shape(object@param))


## wrapped replace methods
setMethod("loc<-", "GPareto", function(object, value) 
           new("GPareto", loc = value, scale = scale(object), shape = shape(object)))
setMethod("location<-", "GPareto", function(object, value) 
           new("GPareto", loc = value, scale = scale(object), shape = shape(object)))
setMethod("scale<-", "GPareto", function(object, value) 
           new("GPareto", loc = loc(object), scale = value, shape = shape(object)))
setMethod("shape<-", "GPareto", function(object, value) 
           new("GPareto", loc = loc(object), scale = scale(object), shape = value))

setValidity("GParetoParameter", function(object){
  if(length(loc(object)) != 1)
    stop("location has to be a numeric of length 1") 
  if(length(scale(object)) != 1)
    stop("scale has to be a numeric of length 1")    
  if(scale(object) <= 0)
    stop("scale has to be positive")
  if(length(shape(object)) != 1)
    stop("shape has to be a numeric of length 1")    
#  if(shape(object) < 0)
#    stop("shape has to be non-negative")
  else return(TRUE)
})

## generating function
GPareto <- function(loc = 0, scale = 1, shape = 0, location = loc){ 
           if(!missing(loc)&&!missing(location)) 
              if(!isTRUE(all.equal(loc,location)))
                 stop("Only one of arguments 'loc' and 'location' may be used.")
           if(!missing(location)) loc <- location
           if(abs(shape) < .Machine$double.eps) return(loc+Exp(rate=1/scale))
           new("GPareto", loc = loc, scale = scale, shape = shape) }

## extra methods for GPareto distribution
setMethod("+", c("GPareto","numeric"),
          function(e1, e2){
            if (length(e2)>1) stop("length of operator must be 1")
            new("GPareto", loc = loc(e1) + e2, scale = scale(e1), shape=shape(e1)) 
          })

setMethod("*", c("GPareto","numeric"),
          function(e1, e2){
            if (length(e2)>1) stop("length of operator must be 1")
            if (isTRUE(all.equal(e2,0))) 
                return(new("Dirac", location = 0, .withArith = TRUE))
            GP <- new("GPareto", loc = loc(e1) * abs(e2), 
                                 scale = scale(e1)*abs(e2), shape=shape(e1))
            if(e2<0) GP <- (-1)*as(GP,"AbscontDistribution")
            return(GP)
          })


