# generating function for class 'PosDefSymmMatrix'
PosDefSymmMatrix <- function(mat){ 
    if(!is.matrix(mat)) mat <- as.matrix(mat)
    new("PosDefSymmMatrix", mat)
}
