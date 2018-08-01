setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GParetoFamily",estimator="Estimate"),
              function(PFam, estimator) as(estimator,"GPDEstimate"))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GParetoFamily",estimator="LDEstimate"),
              function(PFam, estimator) as(estimator,"GPDLDEstimate"))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GParetoFamily",estimator="kStepEstimate"),
              function(PFam, estimator) as(estimator,"GPDkStepEstimate"))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GParetoFamily",estimator="ORobEstimate"),
              function(PFam, estimator) as(estimator,"GPDORobEstimate"))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GParetoFamily",estimator="MCEstimate"),
              function(PFam, estimator) as(estimator,"GPDMCEstimate"))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamily",estimator="Estimate"),
              function(PFam, estimator) as(estimator,"GEVEstimate"))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamily",estimator="LDEstimate"),
              function(PFam, estimator) as(estimator,"GEVLDEstimate"))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamily",estimator="kStepEstimate"),
              function(PFam, estimator) as(estimator,"GEVkStepEstimate"))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamily",estimator="ORobEstimate"),
              function(PFam, estimator) as(estimator,"GEVORobEstimate"))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamily",estimator="MCEstimate"),
              function(PFam, estimator) as(estimator,"GEVMCEstimate"))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamilyMuUnknown",estimator="Estimate"),
              function(PFam, estimator) as(estimator,"GEVEstimate"))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamilyMuUnknown",estimator="LDEstimate"),
              function(PFam, estimator) as(estimator,"GEVLDEstimate"))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamilyMuUnknown",estimator="kStepEstimate"),
              function(PFam, estimator) as(estimator,"GEVkStepEstimate"))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamilyMuUnknown",estimator="ORobEstimate"),
              function(PFam, estimator) as(estimator,"GEVORobEstimate"))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamilyMuUnknown",estimator="MCEstimate"),
              function(PFam, estimator) as(estimator,"GEVMCEstimate"))
