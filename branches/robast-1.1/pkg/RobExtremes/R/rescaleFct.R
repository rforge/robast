### rescale function

## dataFlag - flag, whether the the data is plotted or not
## rescaleFlag

# function returns the list of rescaling arguments to be passed on the
# corresponding diagnostic function

setMethod("rescaleFunction", signature(L2Fam="GParetoFamily"),
   function(L2Fam, dataFlag, rescaleFlag){
   if(!rescaleFlag)
       return(getMethod("rescaleFunction", "ANY")(L2Fam))
    if(dataFlag){
      scaleList <- list(scaleX = substitute(TRUE)
                        ,scaleY = substitute(TRUE)
      )
    } else {
      distr <- distribution(L2Fam)
      scaleList <- list(scaleX = substitute(TRUE)
                        ,scaleX.fct = substitute(p(distr))
                        ,scaleX.inv = substitute(q.l(distr))
                        ,scaleY = substitute(TRUE)
                        ,scaleY.fct = substitute(pnorm)
                        ,scaleY.inv = substitute(qnorm)
                        ,x.ticks = substitute(NULL)
                        ,y.ticks = substitute(NULL)
      )
    }
    return(scaleList)
})

setMethod("rescaleFunction", signature(L2Fam="GEVFamily"),
   getMethod("rescaleFunction", signature(L2Fam="GParetoFamily")))

setMethod("rescaleFunction", signature(L2Fam="GEVFamilyMuUnknown"),
   getMethod("rescaleFunction", signature(L2Fam="GParetoFamily")))

