

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


load.file.to <- function(filePath, destEnvironment=new.env(), 
                         on.not.exist=function(filePath)stop("Fehler mit checkout"))
{
  ## Loads an R data file into new environment
  ##
  ## Returns the environment with the data
  if(!file.exists(filePath)) 
    return(on.not.exist(filePath))
  
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
  # dataEnviron <- load.file.to(file.path(baseDir, "branches/robast-1.0/pkg/RobExtremesBuffer/sysdata.rda"))
  dataEnviron <- load.file.to("sysdata.rda")
  # dataEnviron <- load.file.to(file.path(baseDir, "branches/robast-1.0/pkg/RobAStRDA/R/sysdata.rda"))
  
  return(load.grids.env(dataEnviron, gridName, familyName))
}

expr.str <- function(expr){
  q.expr <- deparse(substitute(expr))
}

get.grid.lookup <- function(gridName){
  return(paste0(".", gridName))
}

get.family.lookup <- function(familyName){
  return(gsub(" ", "", familyName))
}

load.grids.env <- function(env, gridName, familyName){
  grid.lookup <- get.grid.lookup(gridName)
  family.lookup <- get.family.lookup(familyName)
  
  # according to Sec. 5, WriteUp-Interpolators.txt 
  # grid - original interpolation grids
  # gridS - smoothed grids
  lookup <- get0(grid.lookup, envir=env, ifnotfound = quote(warning(paste0("Grid '", grid.lookup, "' does not exists."))))
  
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