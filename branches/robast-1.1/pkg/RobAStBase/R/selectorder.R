.SelectOrderData <- function(data, fct, which.lbs, which.Order, which.nonlbs = NULL){
   ## for data to be plotted in; performs two selections:
   ## on unordered (original) data (acc. to which.lbs)
   ## on data ordered acc. to fct a selection acc. to which.Order is done
   ## return value: list with elements
   #      data, the selected/thinned out data,
   #      y = fct(data)
   #      ind the indices of the selected data in the original data
   #      ind1 the indices of the data selected by which.lbs in the original data
     dimL <- !is.null(dim(data))
     d1  <- if(dimL) dim(data) else 1
     n   <- if(dimL) nrow(data) else length(data)
     ind <- 1:n

     ### function evaluation
     y <- if(dimL) apply(data, 1, fct) else sapply(data,fct)

#------------------------------------------------------------------------------
     ## selected data : data.t
#------------------------------------------------------------------------------

     ### first selection
     if(is.null(which.lbs)) which.lbs <- 1:n
     ## which.lbs0 is a logical of length the original data set, selecting
     ##     the remaining obs after first selection
     which.lbs0 <- ind %in% which.lbs
     # the remaining nb of obs after first selection
     n.s <- sum(which.lbs0) 
     ## produce index for shown data after first selection
     ind.s <- ind[which.lbs0]
     ## function values after first selection
     y.s <- y[ind.s]

     ### ordering 
     oN.s <- order(y.s)
     ## indices remaining after first selection ordered 
     ##         from largest function value to smallest
     ind1.s <- rev(ind[oN.s])

     ### second selection
     ## selection of ordered
     if(is.null(which.Order))
          which.Order <- 1:n.s ## if no 2nd selection performed use all remaining obs.

     ## from ranks in remaining selection pick out those in which.order
     in.t <- (n.s+1)-which.Order
     in.t <- in.t[in.t>0]
     oN.t <-  oN.s[in.t] ## use largest ones in this order
     oN.t <- oN.t[!is.na(oN.t)]

     ## remaining number of observations after 2nd selection
     n.t <- length(oN.t)
     ## observations indices after 2nd selection
     ind.t <- ind.s[oN.t]  
     ind.t <- ind.t[!is.na(ind.t)]
     ## function values after 2nd selection
     y.t <- y.s[oN.t]
     y.t <- y.t[!is.na(y.t)]

     ## produce data set of same dimensions as original one after 2nd selection
     which.lbt <- ind %in% ind.t
     ## if needed recycle arg which.lbs0 (e.g. if several columns)
     which.lbx <- rep(which.lbt, length.out=length(data))
     ## data after both selections
     data.t <- data[which.lbx]
     # if needed recast it to matrix/array
     if(dimL) dim(data.t) <- c(n.t,d1[-1])

#------------------------------------------------------------------------------
     ## data not labelled: data.ns
#------------------------------------------------------------------------------
     if(is.null(which.nonlbs)) which.nonlbs <- 1:n
     #### non selected obs' indices after 1st selection
     ind.ns0 <- ind[!which.lbs0]
     #### non selected obs' indices in 2nd selection
     ind.nt <- if(length(oN.t)) ind.s[-oN.t] else numeric(0)
     #### non selected obs' in total is the union of both non-selected ones
     ind.ns1 <- unique(sort(c(ind.ns0, ind.nt)))
     ind.ns <- ind.ns1[ind.ns1 %in% which.nonlbs]
     ## number of non-selected obs'
     n.ns <- length(ind.ns)

     which.lbns0 <-ind %in% ind.ns
     which.lbnx <- rep(which.lbns0, length.out=length(data))

     ## non selected data
     data.ns <- data[which.lbnx]
     # if needed recast it to matrix
     if(dimL) dim(data.ns) <- c(n.ns,d1[-1])

     y.ns <- y[ind.ns]

     return(list(data=data.t, y=y.t, ind=ind.t, ind1=ind1.s, data.ns=data.ns, y.ns=y.ns, ind.ns = ind.ns))
}

if(FALSE){
.SelectOrderData <- function(data, fct, which.lbs, which.Order){
   ## for data to be plot in performs two selections:
   ## on unordered (original) data (acc. to which.lbs)
   ## on data ordered acc. to fct a selection acc. to which.Order is done
   ## return value: list with elements
   #      data, the selected/thinned out data,
   #      y = fct(data)
   #      ind the indices of the selected data in the original data
   #      ind1 the indices of the data selected by which.lbs in the original data
     dimL <- !is.null(dim(data))
     d1  <- if(dimL) dim(data) else 1
     n   <- if(dimL) nrow(data) else length(data)
     ind <- 1:n
     
     ### selection
     if(is.null(which.lbs)) which.lbs <- 1:n
     ## which.lbs0 is a logical of length the original data set, selecting 
     ##     the remaining obs after first selection
     which.lbs0 <- (1:n) %in% which.lbs
     n0 <- n # n0 is the original nb of obs
     
     n <- sum(which.lbs0) # n now is the remaining nb of obs after first selection
     ## if needed recycle arg which.lbs0 (e.g. if several columns)
     which.lbx <- rep(which.lbs0, length.out=length(data))

     ## data not shown: data.ns 
     data.ns <- data[!which.lbx] ## select data not shown
     # if needed recast it to matrix
     if(dimL) dim(data.ns) <- c(n0-n,d1[-1])
     
     ## remaining data becomes data
     data <- data[which.lbx]
     # if needed recast it to matrix
     if(dimL) dim(data) <- c(n,d1[-1])

     ## produce indices for shown and non-shown data
     ind.ns <- ind[!which.lbs0]
     ind <- ind[which.lbs0]

     ### function evaluation
     y <- if(dimL) apply(data, 1, fct) else sapply(data,fct)
     y.ns <- if(dimL) apply(data.ns, 1, fct) else sapply(data.ns,fct)
     
     ## ordering
     oN <- order(y)
     ind1 <- rev(ind[oN])
     
     ## selection of ordered
     if(is.null(which.Order))
          which.Order <- 1:n ## if no 2nd selection performed use all remaining obs.
     
     oN <-  oN[(n+1)-which.Order] ## use largest ones in this order
     data <- if(dimL) data[oN,] else data[oN]
     y <- y[oN]
     ind <- ind[oN]

     return(list(data=data, y=y, ind=ind, ind1=ind1, data.ns=data.ns, y.ns=y.ns, ind.ns = ind.ns))
}
}

if(FALSE){
x <- rnorm(1000)
.SelectOrderData(x, function(x1)x1^2, 1:100, 1:5)
}
