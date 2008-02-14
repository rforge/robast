if(!isGeneric("optIC")){
    setGeneric("optIC", function(model, risk, ...) standardGeneric("optIC"))
}
if(!isGeneric("getInfRobIC")){
    setGeneric("getInfRobIC", 
        function(L2deriv, risk, neighbor, ...) standardGeneric("getInfRobIC"))
}
if(!isGeneric("getFixRobIC")){
    setGeneric("getFixRobIC", 
        function(Distr, risk, neighbor, ...) standardGeneric("getFixRobIC"))
}
if(!isGeneric("getAsRisk")){
    setGeneric("getAsRisk", 
        function(risk, L2deriv, neighbor, ...) standardGeneric("getAsRisk"))
}
if(!isGeneric("getFiRisk")){
    setGeneric("getFiRisk", 
        function(risk, Distr, neighbor, ...) standardGeneric("getFiRisk"))
}
if(!isGeneric("getInfClip")){
    setGeneric("getInfClip", 
        function(clip, L2deriv, risk, neighbor, ...) standardGeneric("getInfClip"))
}
if(!isGeneric("getFixClip")){
    setGeneric("getFixClip", 
        function(clip, Distr, risk, neighbor, ...) standardGeneric("getFixClip"))
}
if(!isGeneric("getInfGamma")){
    setGeneric("getInfGamma", 
        function(L2deriv, risk, neighbor, ...) standardGeneric("getInfGamma"))
}
if(!isGeneric("getInfCent")){
    setGeneric("getInfCent", 
        function(L2deriv, neighbor, ...) standardGeneric("getInfCent"))
}
if(!isGeneric("getInfStand")){
    setGeneric("getInfStand", 
        function(L2deriv, neighbor, ...) standardGeneric("getInfStand"))
}
if(!isGeneric("getRiskIC")){
    setGeneric("getRiskIC", 
        function(IC, risk, neighbor, L2Fam, ...) standardGeneric("getRiskIC"))
}
if(!isGeneric("optRisk")){
    setGeneric("optRisk", function(model, risk, ...) standardGeneric("optRisk"))
}
if(!isGeneric("radiusMinimaxIC")){
    setGeneric("radiusMinimaxIC", function(L2Fam, neighbor, risk, ...) 
            standardGeneric("radiusMinimaxIC"))
}
if(!isGeneric("getIneffDiff")){
    setGeneric("getIneffDiff", function(radius, L2Fam, neighbor, risk, ...) 
            standardGeneric("getIneffDiff"))
}
if(!isGeneric("leastFavorableRadius")){
    setGeneric("leastFavorableRadius", function(L2Fam, neighbor, risk, ...) 
            standardGeneric("leastFavorableRadius"))
}
if(!isGeneric("lowerCaseRadius")){
    setGeneric("lowerCaseRadius", function(L2Fam, neighbor, risk, ...) standardGeneric("lowerCaseRadius"))
}
