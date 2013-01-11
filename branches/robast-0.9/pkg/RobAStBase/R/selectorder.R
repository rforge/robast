.SelectOrderData <- function(data, fct, which.lbs, which.Order){

     n <- if(is.null(dim(data))) nrow(data) else length(data)
     
     ind <- 1:n
     
     ### selection
     if(is.null(which.lbs)) which.lbs <- 1:n
     which.lbs0 <- (1:n) %in% which.lbs
     n <- sum(which.lbs0)
     which.lbx <- rep(which.lbs0, length.out=length(data))
     data <- data[which.lbx]
     ind <- ind[which.lbs0]

     ### function evaluation
     y <- sapply(data,fct)

     ## ordering
     oN <- order(y)
     ind1 <- rev(ind[oN])
     
     ## selection of ordered
     if(is.null(which.Order))
          which.Order <- 1:n
     oN <-  oN[(n+1)-which.Order]
     data <- if(is.null(dim(data))) data[oN,] else data[oN]
     y <- y[oN]
     ind <- ind[oN]

     return(list(data=data, y=y, ind=ind, ind1=ind1))
}


