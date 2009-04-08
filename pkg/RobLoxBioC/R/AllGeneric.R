############# preparations ################
.onLoad <- function(lib, pkg) {
    require("methods", character = TRUE, quietly = TRUE)
}

if(!isGeneric("robloxbioc")){
    setGeneric("robloxbioc", 
        function(x, ...) standardGeneric("robloxbioc"))
}

if(!isGeneric("KolmogorovMinDist")){
    setGeneric("KolmogorovMinDist", 
        function(x, D, ...) standardGeneric("KolmogorovMinDist"))
}
