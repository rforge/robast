############# preparations ################
.onLoad <- function(lib, pkg) {
    require("methods", character = TRUE, quietly = TRUE)
}

if(!isGeneric("robloxbioc")){
    setGeneric("robloxbioc", 
        function(x, ...) standardGeneric("robloxbioc"))
}
