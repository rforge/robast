source("config.R")

ranges.to.indicies <- function(ranges, values){
  if (is.null(ranges))
    return(NULL)
  
  closest <- function(y)which.min(abs(y-values))
  indices <- sapply(ranges, closest)
  indices <- unique(indices)
  return(indices)
}

ranges.to.grid <- function(ranges, values){
  # converts a list of ranges to grid vector
  idx <- sapply(ranges, function(x)ranges.to.indicies(x, values), simplify=FALSE, USE.NAMES=FALSE)
  
  seq.generator <- function(x){ # generate the sequence for ranges. It may happen that the same grid value is
          # used for an range. I.e. there is only one range. So we do not need the sequence, but just the value itself.
    if(length(x)==2)return(seq(x[1], x[2], by=1))
    else return(x)
  }
  
  res <- sapply(idx, seq.generator, simplify=FALSE, USE.NAMES=FALSE)
  res <- unlist(res, use.names=FALSE)
  res <- as.vector(res)
  return(-res) # NOTE: take care about the Minus sign!!!
}

list.to.matrix <- function(lst){
  ## converts a list of vectors of equal length to a matrix
  ##
  ## list entries are mapped to rows
  return(do.call(rbind, lst))
}

matrix.to.list <- function(mtx){
  ## converts a matrix to a list of vectors
  ##
  ## splits by rows
  result <- split(mtx, c(row(mtx))) 
  names(result) <- NULL
  return(result)
}

create.list.of.empty.lists <- function(length){
  ## Create an list of size length which contains empty lists.
  res <- vector("list", length)
  for(i in 1:length){
    res[[i]] <- list()
  }
  return(res)
}

updateRanges <- function(ranges, new.entry, round.digits=0){
  # add rounded if they differs
  # new.entry <- round(new.entry, digits=round.digits)  # < we do not apply rounding here 
  if (diff(new.entry) == 0)
    return(ranges)
  
  # We always add sorted ranges (x1 < x2)
  ranges[[length(ranges)+1]] <- sort(new.entry)

  # sort all ranges wrt their left limit
  v <- list.to.matrix(ranges) #do.call(rbind, ranges)
  v <- v[order(v[,1]),]
  v <- matrix(v, ncol=2)
  
  # combine overlapping
  y <- 1
  max <- if(length(v)==2) 1 else dim(v)[1]
  while(y < max){
    if (v[y, 2] > v[y+1, 1]){
      v[y,2] <- v[y+1, 2]
      v <- v[-(y+1),]
      max <- max - 1
    }else{
      y = y+ 1
    }
  }
  
  v <- matrix(v, ncol=2)
#   
#   res <- split(v, c(row(v)))
#   names(res) <- NULL
  
  # return(res)
  res <- matrix.to.list(v)
  return(res)
}

can.push <- function(dest, src){
  ## Used by zoom.history to check if a value is valid for push.
  ## 
  ## The reason to use this function is that for each change of axis there are to events. Even they belong 
  ## to the same brush box.
  last.idx <- length(dest)
  if (last.idx == 0)
    return(TRUE)
  
  last.entry <- dest[[last.idx]]
  
  if (any(sapply(last.entry, function(x)is.null(x))))
    return(TRUE)
  
  return(any(sapply(names(dest), function(e)any(dest[[e]] != src[[e]]))))
}

print.named.list <- function(lst, prefix=NULL, sep=' '){
  ## Print a named list of values.
  ## 
  ## Example: 
  ## x <- list(a=1, b=2, c=c(3,4)); print.named.list(x)
  ## will print:
  ## x$a = 1
  ## x$b = 2
  ## x$c = 3, 4
  ## 
  ## One can control the space before and after equality by sep.
  ## 
  ## If you want a different mapping than the variable name, then use prefix argument.
  
  if(is.null(prefix))
    prefix=deparse(substitute(lst))
  
  list.names <- names(lst)
  for(i in 1:length(lst)){
    e.name <- paste0(prefix, "$", list.names[i])
    e.value <- paste(lst[[i]], collapse=", ")
  }
}

parent.path <- function(){
  ## Returns the path which contains the caller's source file. 
  ##
  ## The returned path is without trailig slash.
  ##
  ## e.g. if it is called in /path/to/my_script.R 
  ## then the returned value is "/path/to"
  frame_files <- lapply(sys.frames(), function(x)x$ofile)
  frame_files <- Filter(Negate(is.null), frame_files)
  
  return(dirname(frame_files[[length(frame_files)]]))
}


loadRDataFileToEnv <- function(filePath, destEnvironment=new.env(), 
                         onError=function(filePath)stop("Fehler mit checkout"))
{
  ## Loads an R data file into new environment
  ##
  ## Returns the environment with the data
  if(!file.exists(filePath)) 
    return(onError(filePath))
  
  # destEnvironment <- new.env()
  load(filePath, envir=destEnvironment)
  
  return(destEnvironment)
}


get.partial.matched <- function(entry, list, errormsg=""){
  ## Checked partial matching.
  ##
  ## If the entry is in from, then the full value is returned.
  result <- list[pmatch(entry, list)]
  if (is.na(result))
    stop(errormsg)
  
  return(result)
}


load.grids <- function(gridName, familyName, baseDir){
  # dataEnviron <- loadRDataFileToEnv(file.path(baseDir, "branches/robast-1.0/pkg/RobExtremesBuffer/sysdata.rda"))
  dataEnviron <- loadRDataFileToEnv("sysdata.rda")
  # dataEnviron <- loadRDataFileToEnv(file.path(baseDir, "branches/robast-1.0/pkg/RobAStRDA/R/sysdata.rda"))
  
  return(loadGridsIntoEnv(dataEnviron, gridName, familyName))
}

expr.str <- function(expr){
  q.expr <- deparse(substitute(expr))
}

getGridLookupName <- function(gridName){
  return(paste0(".", gridName))
}

getFamilyLookupName <- function(familyName){
  return(gsub(" ", "", familyName))
}

getLMName <- function(lm) {
  NAMES_LM <- c("b", "a.a[sc]", "a.a[sh]", "z.i[sc]", "z.i[sh]",
                "A.a[sc,sc]", "A.a[sc,sh]", "A.a[sh,sc]", "A.a[sh,sh]",
                "A.i[sc,sc]", "A.i[sc,sh]", "A.i[sh,sc]", "A.i[sh,sh]")
  
  name <- NAMES_LM[lm]
  name <- paste("LM", name)
  
  return(name) 
}

loadGridsIntoEnv <- function(env, gridName, familyName){
  grid.lookup <- getGridLookupName(gridName)
  family.lookup <- getFamilyLookupName(familyName)
  
  # according to Sec. 5, WriteUp-Interpolators.txt 
  # grid - original interpolation grids
  # gridS - smoothed grids
  if (exists(grid.lookup, envir=env)) {
    lookup <- get(grid.lookup, envir=env)
  } else {
    warning(paste0("Grid '", grid.lookup, "' does not exists."))
  }
  lookup.family <- lookup[[family.lookup]]
  if (is.null(lookup.family))
    warning(paste0("Family '", familyName, "' does not exists for grid '", gridName, "'"))
  
  grid.orig <- lookup.family[["grid"]]
  if (is.null(grid.orig))
    warning(paste0("original grid does not exists for family '", familyName, "' and grid '", gridName, "'"))
  
  grid.smoothed <- lookup.family[["gridS"]]  
  if (is.null(grid.smoothed))
      warning(paste0("original grid does not exists for family '", familyName, "' and grid '", gridName, "'"))
  
  res <- list(orig=grid.orig, smoothed=grid.smoothed)
  return(res)
}

matlines.internal <- function(grid, grid.orig, grid.smoothed, idxCol, restriction, lwd, lty, col, xlab, ylab, main, withLegend, ...){
  if(is.null(restriction))
    restriction <- 1:nrow(grid)
  
  y.plot <- cbind(grid.orig[restriction, idxCol], grid[restriction, idxCol])
  y.lines <- cbind(grid.orig[restriction, idxCol], grid.smoothed[restriction, idxCol], grid[restriction,idxCol])
  x <- grid[restriction, 1]
  # par(xpd=TRUE)
  # type == 'n': no plotting
  # x <- c(xxxxxxx)
  # y <- matrix(y1; y2; y3), where y_i are vectors
  # Plotting is one x value (is the xi) and couple of y values.
  matplot(x, y.plot, type="n", xlab=xlab, ylab=ylab, main=main,  ...)
  matlines(x, y.lines, lwd=lwd, lty=lty, col=col)
  if(withLegend)
    legend("top", c("Original", "Smoothed (gespeicherte)", "in Bearbeitung"), lty=lty, col=col, lwd=lwd)
}


log.value <- function(x){
  x.name <- deparse(substitute(x))
  if(is.vector(x)){
    x.value <- paste(x, collapse = ", ")
  }else
    x.value <- x
  print(paste(x.name, "=", x.value))
}


# src for function: P:\EugenMassini\robast\branches\robast-1.0\pkg\RobExtremes\R\00fromRobAStRDA.R

# Call: .MakeSmoothGridList(thGrid=result[,1], Y=result[,-1], df = input$df, gridRestrForSmooth = gridRestrForSmooth)
#  result <- grids()[["smoothed"]] OR grids()[["orig"]]
.MakeSmoothGridList <- function(thGrid, Y, df = NULL,
                                gridRestrForSmooth = NULL){
  
  ############################################
  ### create internal lm-grid: lmgrid
  if(length(dim(Y))==3)
    lmgrid <- Y[,1,,drop=TRUE]
  else 
    lmgrid <- Y[,drop=FALSE]
  
  if (is.vector(lmgrid) && !is.list(lmgrid) && !is.matrix(lmgrid))
    lmgrid <- as.matrix(lmgrid, ncol=1)
  
  iNA <- any(is.na(lmgrid))
  lmgrid <- lmgrid[!iNA,,drop=FALSE]
  thGrid <- thGrid[!iNA]
  oG <- order(thGrid)
  thGrid <- thGrid[oG]
  lmgrid <- lmgrid[oG,,drop=FALSE]
  
  ############################################
  ### Handling of df.
  ### Set of each Lagrange multiplier the same df.
  if(!is.null(df)){
    df0 <- vector("list",ncol(lmgrid))
    if(is.numeric(df)){
      df <- as.list(rep(df, length.out=ncol(lmgrid)))
    }
  }else{ # handling for NULL (create a list of NULL)
    df0 <- vector("list",ncol(lmgrid)+1)
    df0[[ncol(lmgrid)+1]] <- NULL
    df <- df0
  }
  
  
  ############################################
  ### Handling of gridRestrForSmooth
  if(is.null(gridRestrForSmooth)){
    gridRestrForSmooth <- as.data.frame(matrix(TRUE, nrow(lmgrid), ncol(lmgrid)))
  }
  if ( ( is.vector(gridRestrForSmooth) && !is.list(gridRestrForSmooth) )
     || is.matrix(gridRestrForSmooth)){
    gridRestrForSmooth <- as.data.frame(gridRestrForSmooth)
  }
  if(is.list(gridRestrForSmooth)){
    gm <- vector("list",ncol(lmgrid))
    idx <- rep(1:length(gridRestrForSmooth), length.out=ncol(lmgrid))
    for (i in 1:ncol(lmgrid)){
      if(!is.null(gridRestrForSmooth[[idx[i]]])){
        gm[[i]] <- gridRestrForSmooth[[idx[i]]]
      }else{
        gm[[i]] <- rep(TRUE,nrow(lmgrid))
      }
    }
    gridRestrForSmooth <- gm
  }
  
  
  ############################################
  ### Create plots
  for(i in 1:ncol(lmgrid)){
    gmi <- gridRestrForSmooth[[i]]
    
    if(is.null(df[[i]])){
      SmoothSpline <- smooth.spline(thGrid[gmi], lmgrid[gmi, i])
    } else {
      SmoothSpline <- smooth.spline(thGrid[gmi], lmgrid[gmi, i], df = df[[i]])
    }
    
    lmgrid[, i] <- predict(SmoothSpline, thGrid)$y   
  }
  return(cbind(xi=thGrid,LM=lmgrid))
}


createStorableGrid <- function(familyName, gridName, editingGrid, origSmoothedGrid, useExisting, dfs, ranges){
  numLMs <- getNumMultiplicators(gridName, familyName)
  xiValues <- editingGrid[,1]
  
  
  storeGridComponent <- function(comp){
    if(useExisting[[comp]]) {
      return(origSmoothedGrid[,comp+1])
    } else {
      restrictions <- get.restrictions.for.smooth(comp, ranges, xiValues)
      smoothed <- applySmoothing(editingGrid, dfs[[comp]], restrictions)
      return(smoothed[,comp+1])
    }
  }
  
  totalGrid <- sapply(1:numLMs, storeGridComponent)
  totalGrid <- cbind(xiValues, totalGrid)
  
  if(TEST.save.grid)
    smoothed.totalgrid <<- totalGrid
  
  return(totalGrid)
}

saveGridToCsv <- function(familyName, gridName, editingGrid, origSmoothedGrid, useExisting, dfs, ranges){
  require(ROptEst)
  totalGrid <- createStorableGrid(familyName, gridName, editingGrid, origSmoothedGrid, useExisting, dfs, ranges)
  
  destFileName <- paste0("interpol", familyName, gridName, ".csv")
  .saveGridToCSV(totalGrid, destFileName, gridName, paste0(".", familyName))
}

getMainTitle <- function(gridName, familyName){
  familyName1 <- gsub(" [F,f]amily","", familyName)
  
  return(paste(gridName, familyName1, sep="-"))
}

getNumMultiplicators <- function(gridName, family){
  if(gridName=="Sn")
    return(1)
  return(if(family == "GEVU Family") 25 else 13)
}


isSnGridWithInvalidLM <- function(gridName, lm) {
  return((gridName == "Sn") && lm != 1)
}


## Create smooth grid
applySmoothing <- function(grid, df, gridRestrictions){
  # grid[,1] - The grid positions
  # grid[,2:end] - the Lagrange multiplier values
  result <- .MakeSmoothGridList(grid[,1], grid[,-1], df=df, gridRestrForSmooth=gridRestrictions)
  return(result)
}


## Uses the list of all ranges and deletes the required one
## 
## For deletion is the string value of input$ranges used,
## i.e. one assumes the index of the entry is in the string which.to.delete. 
## The format of the string should be "Index : somevalues"
##
## If the index is invalid or null the full range is returned
update.ranges.after.delete <- function(allRanges, which.to.delete){
  result <- allRanges
  if (!is.null(which.to.delete)){
    idx.to.delete <- as.numeric(gsub(" : .*$", "", which.to.delete))
    if (!is.na(idx.to.delete)){
      result[[idx.to.delete]] <- NULL
    }
  }
  return(result)
}

## Create the list of range as strings for the list output
update.ranges.output <- function(ranges){
  if(is.null(ranges))
    return(NULL)
  
  rangeNums <- sapply(ranges, function(x)paste(round(x, digits=3), collapse=", "))
  result <- sapply(seq_along(rangeNums), function(i)paste(i, ':', rangeNums[[i]]))
  return(result)
}

zoomIn <- function(brush, zoomList, zoomHistory){
  # Unfortunately this method is called twice (probably for each coordinate)
  # Hence we need to do the if check for not adding an element twice
  # It uses the heuristics that two distinguishing boxes will differ in both, xlim and ylim values
  if (is.null(brush))
    return(NULL)
  
  # store to history if the last differs from current
  if (can.push(zoomHistory, zoomList)){
    last.idx <- length(zoomHistory)
    zoomHistory[[last.idx + 1]] <<- zoomList
    
    # set new values
    res <- list (xlim=c(brush$xmin, brush$xmax), ylim=c(brush$ymin, brush$ymax))
    return(res)
  }
  
  return(NULL)
}

zoom.out <- function(){
  idx.last <- length(zoomHistory)
  if(idx.last > 0){
    last <- zoomHistory[[idx.last]]
    zoomHistory <<- zoomHistory[1:(idx.last-1)]
    
    res <-list(xlim=c(last$xlim[1], last$xlim[2]), ylim=c(last$ylim[1], last$ylim[2]))
    return(res)
  }
  return(NULL)
}


delete.ranges <- function(whichLM, state.ranges, input.ranges){
  if(!is.null(input.ranges) && (prev.deleted != input.ranges)){
    res <- update.ranges.after.delete(allRanges=state.ranges[[whichLM]], which.to.delete=input.ranges)
    prev.deleted <<- input.ranges
    return(res)
  }else{
    prev.deleted <<- ""
    return(NULL)
  }
}

get.restrictions.for.smooth <- function(which, from, grid.param){
  if (length(from)==0)
    return(NULL)
  ranges <- from[[which]]
  if (!is.null(ranges) && length(ranges) > 0){
    return(ranges.to.grid(ranges, grid.param))
  }
  
  return(NULL)
}

###########################################################################
# history local grid save intermediate will be saved in history Rdata 
# (analogue as in sysdata.rda, see Sec. 5 WriteUp-Interpolators.txt)
# [OptCrit],
# > [model1],
# >> [timestamp]
# >>> [ranges]
# >>> [dfs]
###########################################################################
addToHistory <- function(familyName, gridName, dfs, ranges, useExisting){
  commitsEnv <- loadRDataFileToEnv(HISTORY_COMMITS_FILE, onError=function(x)new.env())
  
  # Get entry
  gridLookupName <- getGridLookupName(gridName)
  familyLookup <- getFamilyLookupName(familyName)
  
  if(exists(gridLookupName, envir=commitsEnv)){
    models <- get(gridLookupName, envir=commitsEnv)
  }else{
    models <- list()
  }
  
  # append
  if(is.null(models[[familyLookup]])){
    models[[familyLookup]] <- list()
  }
  
  timestamp = format(Sys.time())
  models[[familyLookup]][[timestamp]] <- list(dfs=dfs, ranges=ranges, useExisting=useExisting)
  
  assign(gridLookupName, value=models, envir=commitsEnv)
  
  names <- ls(commitsEnv, all.names=TRUE)
  save(list=names, file=HISTORY_COMMITS_FILE, envir=commitsEnv)
}



loadDataFromHistory <- function(familyName, gridName) {
  commitsEnv <- loadRDataFileToEnv(HISTORY_COMMITS_FILE, onError=function(x)new.env())
  
  # Get entry
  gridLookupName <- getGridLookupName(gridName)
  familyLookup <- getFamilyLookupName(familyName)
  
  if(exists(gridLookupName, envir=commitsEnv)){
    models <- get(gridLookupName, envir=commitsEnv)
    return(models[[familyLookup]])
  }else{
    return(NULL)
  }
  
  # append
  if(is.null(models[[familyLookup]])){
    return(NULL)
  }

}


loadFromHistory <- function(familyName, gridName) {
  data <- loadDataFromHistory(familyName, gridName)
  
  timestamps <- names(data)
  timestamps <- as.POSIXlt(timestamps)
  
  latestTimestamp <- max(timestamps)
  
  # Need again a string to be able to access the data  
  latestTimestamp <- as.character(latestTimestamp)
  result <- data[[latestTimestamp]]
  return(result)
}


hasHistory <- function(familyName, gridName) {
  loadedHistory <- loadDataFromHistory(familyName, gridName)
  return(!is.null(loadedHistory))
}


checkRequiredPackages <- function(packages=REQUIRED_PACKAGES) {
  inQuotes <- function(x) paste("\"", x, "\"", sep="")
  
  notInstalled <- ! (packages %in% installed.packages())
  notInstalledPackages <- packages[notInstalled]
    
  if(length(notInstalledPackages) > 0) {
    packagesToInstall <- paste(inQuotes(notInstalledPackages), collapse=", ")
    cat("------------------------------------------------------------------\n")
    cat(paste("Please run> install.packages(", packagesToInstall, ")\n", sep=""))
    cat("------------------------------------------------------------------\n")
    
    stopApp()
  }
}

toogleAvailabilityComponents <- function(components, enable) {
  if(enable){
    sapply(components, function(x)shinyjs::enable(x))
  } else {
    sapply(components, function(x)shinyjs::disable(x))
  }
  
}
