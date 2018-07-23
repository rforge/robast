
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
if(!isGeneric("gev.diag")){
   setGeneric("gev.diag", function(z) standardGeneric("gev.diag"))
}
if(!isGeneric("gev.prof")){
   setGeneric("gev.prof", function(z, ...) standardGeneric("gev.prof"))
}
if(!isGeneric("gev.profxi")){
   setGeneric("gev.profxi", function(z, ...) standardGeneric("gev.profxi"))
}
if(!isGeneric("gpd.diag")){
   setGeneric("gpd.diag", function(z,...) standardGeneric("gpd.diag"))
}
if(!isGeneric("gpd.prof")){
   setGeneric("gpd.prof", function(z, ...) standardGeneric("gpd.prof"))
}
if(!isGeneric("gpd.profxi")){
   setGeneric("gpd.profxi", function(z, ...) standardGeneric("gpd.profxi"))
}
if(!isGeneric("locscaleshapename")){
   setGeneric("locscaleshapename", function(object) standardGeneric("locscaleshapename"))
}
if(!isGeneric("locscaleshapename<-")){
   setGeneric("locscaleshapename<-", function(object,value) standardGeneric("locscaleshapename<-"))
}
if(!isGeneric("shapename")){
   setGeneric("shapename", function(object) standardGeneric("shapename"))
}
if(!isGeneric("locationname")){
   setGeneric("locationname", function(object) standardGeneric("locationname"))
}

