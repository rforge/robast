if(!isGeneric("param<-")){ 
    setGeneric("param<-", function(object, value) standardGeneric("param<-"))
}
if(!isGeneric("radius")){ 
    setGeneric("radius", function(object) standardGeneric("radius"))
}
if(!isGeneric("center")){ 
    setGeneric("center", function(object) standardGeneric("center"))
}
if(!isGeneric("center<-")){
    setGeneric("center<-", function(object, value) standardGeneric("center<-"))
}
if(!isGeneric("neighbor")){
    setGeneric("neighbor", function(object) standardGeneric("neighbor"))
}
if(!isGeneric("neighbor<-")){
    setGeneric("neighbor<-", function(object, value) standardGeneric("neighbor<-"))
}
if(!isGeneric("bound")){ 
    setGeneric("bound", function(object) standardGeneric("bound"))
}
if(!isGeneric("width")){ 
    setGeneric("width", function(object) standardGeneric("width"))
}
if(!isGeneric("Curve")){
    setGeneric("Curve", function(object) standardGeneric("Curve"))
}
if(!isGeneric("Risks")){
    setGeneric("Risks", function(object) standardGeneric("Risks"))
}
if(!isGeneric("Risks<-")){
    setGeneric("Risks<-", function(object, value) standardGeneric("Risks<-"))
}
if(!isGeneric("addRisk<-")){
    setGeneric("addRisk<-", function(object, value) standardGeneric("addRisk<-"))
}
if(!isGeneric("Infos")){
    setGeneric("Infos", function(object) standardGeneric("Infos"))
}
if(!isGeneric("Infos<-")){
    setGeneric("Infos<-", function(object, value) standardGeneric("Infos<-"))
}
if(!isGeneric("addInfo<-")){
    setGeneric("addInfo<-", function(object, value) standardGeneric("addInfo<-"))
}
if(!isGeneric("CallL2Fam")){ 
    setGeneric("CallL2Fam", function(object) standardGeneric("CallL2Fam"))
}
if(!isGeneric("CallL2Fam<-")){ 
    setGeneric("CallL2Fam<-", function(object, value) standardGeneric("CallL2Fam<-"))
}
if(!isGeneric("generateIC")){
    setGeneric("generateIC", function(neighbor, L2Fam, ...) standardGeneric("generateIC"))
}
if(!isGeneric("checkIC")){
    setGeneric("checkIC", function(IC, L2Fam, ...) standardGeneric("checkIC"))
}
if(!isGeneric("evalIC")){
    setGeneric("evalIC", function(IC, x) standardGeneric("evalIC"))
}
if(!isGeneric("clip")){
    setGeneric("clip", function(object) standardGeneric("clip"))
}
if(!isGeneric("clip<-")){
    setGeneric("clip<-", function(object, value) standardGeneric("clip<-"))
}
if(!isGeneric("cent")){
    setGeneric("cent", function(object) standardGeneric("cent"))
}
if(!isGeneric("cent<-")){
    setGeneric("cent<-", function(object, value) standardGeneric("cent<-"))
}
if(!isGeneric("stand")){
    setGeneric("stand", function(object) standardGeneric("stand"))
}
if(!isGeneric("stand<-")){
    setGeneric("stand<-", function(object, value) standardGeneric("stand<-"))
}
if(!isGeneric("lowerCase")){
    setGeneric("lowerCase", function(object) standardGeneric("lowerCase"))
}
if(!isGeneric("lowerCase<-")){
    setGeneric("lowerCase<-", function(object, value) standardGeneric("lowerCase<-"))
}
if(!isGeneric("neighborRadius")){
    setGeneric("neighborRadius", function(object) standardGeneric("neighborRadius"))
}
if(!isGeneric("neighborRadius<-")){
    setGeneric("neighborRadius<-", function(object, value) standardGeneric("neighborRadius<-"))
}
if(!isGeneric("clipLo")){
    setGeneric("clipLo", function(object) standardGeneric("clipLo"))
}
if(!isGeneric("clipLo<-")){
    setGeneric("clipLo<-", function(object, value) standardGeneric("clipLo<-"))
}
if(!isGeneric("clipUp")){
    setGeneric("clipUp", function(object) standardGeneric("clipUp"))
}
if(!isGeneric("clipUp<-")){
    setGeneric("clipUp<-", function(object, value) standardGeneric("clipUp<-"))
}
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
if(!isGeneric("ksEstimator")){
    setGeneric("ksEstimator", 
        function(x, distribution, ...) standardGeneric("ksEstimator"))
}
#if(!isGeneric("cvmEstimator")){
#    setGeneric("cvmEstimator", 
#        function(x, distribution, weight, ...) standardGeneric("cvmEstimator"))
#}
if(!isGeneric("oneStepEstimator")){
    setGeneric("oneStepEstimator", 
        function(x, IC, start) standardGeneric("oneStepEstimator"))
}
if(!isGeneric("locMEstimator")){
    setGeneric("locMEstimator", function(x, IC, ...) standardGeneric("locMEstimator"))
}
if(!isGeneric("infoPlot")){
    setGeneric("infoPlot", function(object) standardGeneric("infoPlot"))
}
if(!isGeneric("sign")){
    setGeneric("sign", function(x) standardGeneric("sign"))
}
if(!isGeneric("sign<-")){
    setGeneric("sign<-", function(object,value) standardGeneric("sign<-"))
}
if(!isGeneric("nu")){
    setGeneric("nu", function(object) standardGeneric("nu"))
}
if(!isGeneric("nu<-")){
    setGeneric("nu<-", function(object,value) standardGeneric("nu<-"))
}
