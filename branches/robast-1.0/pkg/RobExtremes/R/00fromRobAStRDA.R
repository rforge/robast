#------------------------------------
#### utilities copied from package RobAStRDA v.1.0  svn-rev 767
#------------------------------------

.versionSuff <- function(name){
    paste(sep="", name, if(getRversion()<"2.16") ".O" else ".N")
}
.MakeSmoothGridList <- function(thGrid, Y, df = NULL,
                            gridRestrForSmooth = NULL){
   if(length(dim(Y))==3)
      LMGrid <- Y[,1,,drop=TRUE]
   else LMGrid <- Y[,drop=FALSE]

  if(!is.null(df)){
    df0 <- vector("list",ncol(LMGrid))
    if(is.numeric(df)){
    df <- rep(df,length.out=ncol(LMGrid))
    for(i in 1:ncol(LMGrid)) df0[[i]] <- df[i]
    df <- df0
    }
  }else{
    df0 <- vector("list",ncol(LMGrid)+1)
    df0[[ncol(LMGrid)+1]] <- NULL
    df <- df0
  }

   iNA <- apply(LMGrid,1, function(u) any(is.na(u)))
   LMGrid <- LMGrid[!iNA,,drop=FALSE]
   thGrid <- thGrid[!iNA]
   oG <- order(thGrid)
   thGrid <- thGrid[oG]
   LMGrid <- LMGrid[oG,,drop=FALSE]

   if(is.null(gridRestrForSmooth))
      gridRestrForSmooth <- as.data.frame(matrix(TRUE,nrow(LMGrid),ncol(LMGrid)))
   if((is.vector(gridRestrForSmooth)&&!is.list(gridRestrForSmooth))||
       is.matrix(gridRestrForSmooth))
      gridRestrForSmooth <- as.data.frame(gridRestrForSmooth)

   if(is.list(gridRestrForSmooth)){
      gm <- vector("list",ncol(LMGrid))
      idx <- rep(1:length(gridRestrForSmooth), length.out=ncol(LMGrid))
      for (i in 1:ncol(LMGrid)){
           if(!is.null(gridRestrForSmooth[[idx[i]]])){
               gm[[i]] <- gridRestrForSmooth[[idx[i]]]
           }else{
               gm[[i]] <- rep(TRUE,nrow(LMGrid))
           }
      }
      gridRestrForSmooth <- gm
   }

   for(i in 1:ncol(LMGrid)){
       gmi <- gridRestrForSmooth[[i]]
       if(is.null(df[[i]])){
            SmoothSpline <- smooth.spline(thGrid[gmi], LMGrid[gmi, i])
            LMGrid[, i] <- predict(SmoothSpline, thGrid)$y
       } else {
            SmoothSpline <- smooth.spline(thGrid[gmi], LMGrid[gmi, i],
                                          df = df[[i]])
            LMGrid[, i] <- predict(SmoothSpline, thGrid)$y
       }
   }
   return(cbind(xi=thGrid,LM=LMGrid))
}

