
if(!isGeneric("kMAD")){
   setGeneric("kMAD", function(x, k, ...) standardGeneric("kMAD"))
}


if(!isGeneric("loc")){
   setGeneric("loc", function(object) standardGeneric("loc"))
}

if(!isGeneric("loc<-")){
   setGeneric("loc<-", function(object, value) standardGeneric("loc<-"))
}


if(!isGeneric("Qn")){
   setGeneric("Qn", function(x, ...) standardGeneric("Qn"))
}

if(!isGeneric("Sn")){
   setGeneric("Sn", function(x, ...) standardGeneric("Sn"))
}

if(!isGeneric("dispersion")){
   setGeneric("dispersion", function(object, ...) standardGeneric("dispersion"))
}

if(!isGeneric(".loc")){
   setGeneric(".loc", function(L2Fam, ...) standardGeneric(".loc"))
}
