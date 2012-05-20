###########################################
## Class: Generalized Extreme Value Distribution  
##  
## @param: location, scale, shape 
###########################################

## access methods
setMethod("location", "GEVParameter", function(object) object@loc)
setMethod("loc", "GEVParameter", function(object) object@loc)
setMethod("scale", "GEVParameter", 
           function(x, center = TRUE, scale = TRUE) x@scale)
           ### odd arg-list due to existing function in base package 
setMethod("shape", "GEVParameter", function(object) object@shape)

## replace Methods
setReplaceMethod("loc", "GEVParameter", 
    function(object, value){ object@loc <- value; object })
setReplaceMethod("location", "GEVParameter", 
    function(object, value){ object@loc <- value; object })
setReplaceMethod("scale", "GEVParameter", 
    function(object, value){
        if(length(value) != 1 || value <= 0)
            stop("'value' has to be a single positive real number!")
        object@scale <- value; object})
setReplaceMethod("shape", "GEVParameter", 
    function(object, value){ object@shape <- value; object})

## wrapped access methods
setMethod("location", "GEV", function(object) loc(object@param))
setMethod("loc", "GEV", function(object) loc(object@param))
setMethod("scale", "GEV", 
           function(x, center = TRUE, scale = TRUE)
           scale(x@param))
setMethod("shape", "GEV", function(object) shape(object@param))


## wrapped replace methods
setMethod("loc<-", "GEV", function(object, value) 
           new("GEV", loc = value, scale = scale(object), shape = shape(object)))
setMethod("location<-", "GEV", function(object, value) 
           new("GEV", loc = value, scale = scale(object), shape = shape(object)))
setMethod("scale<-", "GEV", function(object, value) 
           new("GEV", loc = loc(object), scale = value, shape = shape(object)))
setMethod("shape<-", "GEV", function(object, value) 
           new("GEV", loc = loc(object), scale = scale(object), shape = value))

setValidity("GEVParameter", function(object){
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
GEV <- function(loc = 0, scale = 1, shape = 0, location = loc){ 
           if(!missing(loc)&&!missing(location)) 
              if(!isTRUE(all.equal(loc,location)))
                 stop("Only one of arguments 'loc' and 'location' may be used.")
           if(!missing(location)) loc <- location
           #if(abs(shape) < .Machine$double.eps) return(Gumbel(loc=loc,scale=scale))
           new("GEV", loc = loc, scale = scale, shape = shape) }

## extra methods for GEV distribution
setMethod("+", c("GEV","numeric"),
          function(e1, e2){
            if (length(e2)>1) stop("length of operator must be 1")
            new("GEV", loc = loc(e1) + e2, scale = scale(e1), shape=shape(e1)) 
          })

setMethod("*", c("GEV","numeric"),
          function(e1, e2){
            if (length(e2)>1) stop("length of operator must be 1")
            if (isTRUE(all.equal(e2,0))) 
                return(new("Dirac", location = 0, .withArith = TRUE))
            GEV <- new("GEV", loc = loc(e1) * abs(e2), 
                                 scale = scale(e1)*abs(e2), shape=shape(e1))
            if(e2<0) GEV <- (-1)*as(GEV,"AbscontDistribution")
            return(GEV)
          })


