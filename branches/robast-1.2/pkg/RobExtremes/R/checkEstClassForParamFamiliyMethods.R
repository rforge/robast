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
              function(PFam, estimator){# ret0 <- as(estimator,"GPDMCEstimate")
                 fromSlotNames <- slotNames(class(estimator))
                 to <- new("GPDMCALEstimate")
                 for(item in fromSlotNames) slot(to, item) <- slot(estimator,item)
                 to@pIC <- substitute(getPIC(estimator0), list(estimator0=estimator))
                 return(to)
              })
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GParetoFamily",estimator="MLEstimate"),
              getMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GParetoFamily",estimator="MCEstimate")))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GParetoFamily",estimator="MDEstimate"),
              getMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GParetoFamily",estimator="MCEstimate")))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GParetoFamily",estimator="CvMMDEstimate"),
              getMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GParetoFamily",estimator="MCEstimate")))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GParetoFamily",estimator="MCALEstimate"),
              function(PFam, estimator) as(estimator,"GPDMCALEstimate"))


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
              function(PFam, estimator){ #ret0 <- as(estimator,"GEVMCEstimate")
                 fromSlotNames <- slotNames(class(estimator))
                 to <- new("GEVMCALEstimate")
                 for(item in fromSlotNames) slot(to, item) <- slot(estimator,item)
                 to@pIC <- substitute(getPIC(estimator0), list(estimator0=estimator))
                 return(to)
              })
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamily",estimator="MLEstimate"),
              getMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamily",estimator="MCEstimate")))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamily",estimator="MDEstimate"),
              getMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamily",estimator="MCEstimate")))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamily",estimator="CvMMDEstimate"),
              getMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamily",estimator="MCEstimate")))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamily",estimator="MCALEstimate"),
              function(PFam, estimator) as(estimator,"GEVMCALEstimate"))


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
              function(PFam, estimator){ #ret0 <- as(estimator,"GEVMCEstimate")
                 fromSlotNames <- slotNames(class(estimator))
                 to <- new("GEVMCALEstimate")
                 for(item in fromSlotNames) slot(to, item) <- slot(estimator,item)
                 to@pIC <- substitute(getPIC(estimator0), list(estimator0=estimator))
                 return(to)
              })
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamilyMuUnknown",estimator="MLEstimate"),
              getMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamilyMuUnknown",estimator="MCEstimate")))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamilyMuUnknown",estimator="MDEstimate"),
              getMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamilyMuUnknown",estimator="MCEstimate")))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamilyMuUnknown",estimator="CvMMDEstimate"),
              getMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamilyMuUnknown",estimator="MCEstimate")))
setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="GEVFamilyMuUnknown",estimator="MCALEstimate"),
              function(PFam, estimator) as(estimator,"GEVMCALEstimate"))
