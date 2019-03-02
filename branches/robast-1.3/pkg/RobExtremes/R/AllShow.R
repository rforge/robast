
setMethod("show", "LDEstimate",
    function(object){
       digits <- getOption("digits")
       getMethod("show","Estimate")(object)
       if(getdistrModOption("show.details")!="minimal"){
        cat("Location:", object@location, "\n")
        cat("Dispersion:", object@dispersion, "\n")
       }
    })


