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


#------------------------------------------------------------------------------
     ## firt selection: selected data in first : data.s
#------------------------------------------------------------------------------

     ### first selection
     if(is.null(which.lbs)) which.lbs <- 1:n
     ## which.lbs0 is a logical of length the original data set, selecting
     ##     the remaining obs after first selection
     which.lbs0 <- ind %in% which.lbs
     # the remaining nb of obs after first selection
     n.s <- sum(which.lbs0) 
     i.s <- 1:n.s
     ## produce index for shown data after first selection
     ind.s <- ind1.s <- ind[which.lbs0]
     ## first selection
     data.s <- .SelectIndex(data,1,ind.s)

#------------------------------------------------------------------------------
     ### second selection
#------------------------------------------------------------------------------

     ind2 <- ind.s

     ## function values only after first selection
     ### function evaluation
     y.s <- if(dimL) apply(data.s, 1, fct) else sapply(data.s,fct)
     ## simpler with ranks, see distrMod:::.labelprep
     rky.s <- n.s+1-rank(y.s)
     y2.s <- y.s
     sel2 <- i.s
     data.t <- data.s

     ## selection of ordered
     if(!is.null(which.Order)){
        sel2 <- i.s[rky.s %in% which.Order]
        ind2 <- ind2[sel2]
        y2.s <- y2.s[sel2]
        data.t <- .SelectIndex(data.s,1,sel2)
     }

     ord2 <- order(y2.s, decreasing = TRUE)
     ind2.s <- ind2[ord2]
     sel2 <- sel2[ord2]
     data.t <- .SelectIndex(data.t,1,ord2)
     y.t <- y2.s[ord2]
#------------------------------------------------------------------------------
     ## data not labelled: data.ns
#------------------------------------------------------------------------------
     ind.ns <- ind[-ind2]
     if(length(ind.ns) && !is.null(which.nonlbs))
         ind.ns <- ind.ns[ind.ns%in%which.nonlbs]
     ## non selected data
     data.ns <- .SelectIndex(data,1,ind.ns)
     y.ns <- if(dimL) apply(data.ns, 1, fct) else sapply(data.ns,fct)

     return(list(data=data.t, y=y.t, ind=ind2.s, ind1=ind1.s, data.ns=data.ns, y.ns=y.ns, ind.ns = ind.ns))
}

.SelectIndex <- function(data,index,selection){
  dims <- dim(data)
  if(is.null(dims)) return(data[selection])
  datav <- data
  dimv <- dims
  if(index!=1){
     n <- length(dims)
     dims1 <- dims[-index]
     ind0 <- 1:n
     ind1 <- if(index<n) c((1:(index-1))+1,1,((index+1):n)) else c((1:(index-1))+1,1)
     ind2 <- c(index,ind0[-index])
     datav <- aperm(data,ind2)
     dimv <- dims[ind2]
  }
  len0 <- dimv[1]
  len1 <- prod(dimv[-1])
  lens <- length(selection)
  sel <- numeric(lens*len1)
  dimss <- dimv
  dimss[1] <- lens
  for(j in 1:len1)
     sel[1:lens+(j-1)*lens] <- selection+(j-1)*len0
  datas <- datav[sel]
  dim(datas) <- dimss
  if(index!=1){
     datas <- aperm(datas,ind1)
  }
  return(datas)
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

     ## produce indices for shown and non-shown data
     ind.ns <- ind[!which.lbs0]
     ind <- ind[which.lbs0]

     ## data not shown: data.ns
#     data.ns <- if(dimL) data[!which.lbs,] else data[!which.lbs] ## select data not shown
     data.ns <- .SelectIndex(data,1,ind.ns)

     # if needed recast it to matrix
     # if(dimL) dim(data) <- c(n,d1[-1])


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
#     data <- if(dimL) data[oN,] else data[oN]
     data <- .SelectIndex(data,1,oN)
     y <- y[oN]
     ind <- ind[oN]

     return(list(data=data, y=y, ind=ind, ind1=ind1, data.ns=data.ns, y.ns=y.ns, ind.ns = ind.ns))
}
}

if(FALSE){
x <- rnorm(1000)
.SelectOrderData(x, function(x1)x1^2, 1:100, 1:5)
}
