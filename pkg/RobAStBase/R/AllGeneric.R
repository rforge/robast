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
if(!isGeneric("optIC")){
    setGeneric("optIC", function(model, risk, ...) standardGeneric("optIC"))
}
