 .cexscale <- function(y, y1=y, maxcex=4,mincex=0.05,cex, fun=NULL){
         if(length(y)==0) return(NA)
         if(is.list(y) && is.null(y[[1]])) return(NA)
         if(is.null(fun)) fun <- function(x) log(1+abs(x))
         ly <- fun(y)
         ly1 <- fun(unique(c(y,y1)))
         my <- min(ly1,na.rm=TRUE)
         My <- max(ly1,na.rm=TRUE)
         ly0 <- (ly-my)/(My-my)
         ly1 <- ly0*(maxcex-mincex)+mincex
         return(cex*ly1)
 }
