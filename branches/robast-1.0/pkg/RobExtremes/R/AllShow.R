
setMethod("show", "LDEstimate",
    function(object){
       digits <- getOption("digits")
       show(as(object,"Estimate"))
       if(getdistrModOption("show.details")!="minimal"){
        cat("Location:", object@location, "\n")
        cat("Dispersion:", object@dispersion, "\n")
       }
    })


