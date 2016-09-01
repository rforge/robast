#------------------------------------
#### utilities copied from package distr v.2.6  svn-rev 943
#------------------------------------

.isEqual <- function(p0, p1, tol = min( getdistrOption("TruncQuantile")/2,
                                          .Machine$double.eps^.7
                                          ))
                abs(p0-p1)< tol
.isEqual01<- function(x) .isEqual(x,0)|.isEqual(x,1)
