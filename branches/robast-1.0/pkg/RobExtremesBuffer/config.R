# ROBAST_BASE_DIR = if ("massini" == Sys.info()['login']){
# 	paste0(
# 			 if(.Platform$OS.type=="windows") "P:" else "/p/fm",
# 			 "/EugenMassini/robast"
# 	)
# }else{
#   "C:/rtest/RobASt"
# }
# 
# 
# # INTERPOLATION_DIR <- "branches/robast-1.0/pkg/RobExtremes/inst/AddMaterial/interpolation"
# # INTERPOLATION_FILE <- "plotInterpol.R"
# 
# INTERPOLATION_DIR <- "branches/robast-1.0/pkg/RobExtremesBuffer"
# INTERPOLATION_FILE <- "plotInterpolSimple.R"


TEST.save.grid <- TRUE

REQUIRED_PACKAGES <- c("shiny", "shinyjs", "ROptEst")