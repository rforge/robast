### rescale function

## famName - text argument stating the name of the family for which
## the rescaling is done
## dataFlag - flag, whether the the data is plotted or not
## rescaleFlag

# function returns the list of rescaling arguments to be passed on the 
# corresponding diagnostic function 

rescaleFunction <- function(famName, dataFlag, rescaleFlag){
  ### No rescaling ########################################
  
  if(!rescaleFlag){
    if(dataFlag){
      scaleList <- list(scaleX = substitute(FALSE)
                        ,scaleY = substitute(FALSE)
      )
    } else {
      scaleList <- list(scaleX = substitute(FALSE)
                        ,scaleX.fct = substitute(p(eval(IC@CallL2Fam)))
                        ,scaleX.inv = substitute(q(eval(IC@CallL2Fam)))
                        ,scaleY = substitute(FALSE)
                        ,scaleY.fct = substitute(pnorm)
                        ,scaleY.inv=substitute(qnorm)
                        ,x.ticks = substitute(NULL)
                        ,y.ticks = substitute(NULL)
      )
    }
  }
    
  ### GEV #################################################
  
  if((rescaleFlag) & (famName == "GEVFamily")){
    if(dataFlag){
      scaleList <- list(scaleX = substitute(TRUE)
                        ,scaleY = substitute(TRUE)
      )
    } else {
      scaleList <- list(scaleX = substitute(TRUE)
                        ,scaleX.fct = substitute(p(eval(IC@CallL2Fam)))
                        ,scaleX.inv = substitute(q(eval(IC@CallL2Fam)))
                        ,scaleY = substitute(TRUE)
                        ,scaleY.fct = substitute(pnorm)
                        ,scaleY.inv = substitute(qnorm)
                        ,x.ticks = substitute(NULL)
                        ,y.ticks = substitute(NULL)
      )
    }
  }
  
  ### GPD ################################################
  
  if((rescaleFlag) & (famName == "GParetoFamily")){
    if(dataFlag){
      scaleList <- list(scaleX = substitute(TRUE)
                        ,scaleY = substitute(TRUE)
      )
    } else {
      scaleList <- list(scaleX = substitute(TRUE)
                        ,scaleX.fct = substitute(p(eval(IC@CallL2Fam)))
                        ,scaleX.inv = substitute(q(eval(IC@CallL2Fam)))
                        ,scaleY = substitute(TRUE)
                        ,scaleY.fct = substitute(pnorm)
                        ,scaleY.inv = substitute(qnorm)
                        ,x.ticks = substitute(NULL)
                        ,y.ticks = substitute(NULL)
      )
    }
  }
  
  return(scaleList)
}

rescaleFunction("GEVFamily", TRUE, FALSE)