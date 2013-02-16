.versionSuff <- ROptEst:::.versionSuff
getShapeGrid <- function(gridsize=1000,centralvalue=0.7,
                         withPos=TRUE, cutoff.at.0=1e-4, fac = 2){

 gridsize <- gridsize/2
 xi.a <- (1:round(gridsize))/(gridsize+1)*10
 p.u <- .25/gridsize
 q0 <- qnorm(seq(from=0.5,to=1-p.u,length=gridsize%/%2))
 q <- sort(c(-q0,q0))
 xi.g <- centralvalue + q*fac
 if(withPos){
    xi.g <- unique(pmax(xi.g,0))
    l.xi.g <- length(xi.g)
    if(l.xi.g<gridsize){
       n.d <- (gridsize-l.xi.g)%/%8
       n.u <- 7*n.d + (gridsize-l.xi.g)%%8
       xi.n <- (1:n.d)/(n.d+1)*min(xi.g[xi.g>1e-9])
       xi.g <- sort(unique(c(0,xi.n,xi.g[xi.g>1e-9])))
       if(n.u>0){
          p.u2 <- 1-p.u + p.u*(1:n.u)/(n.u+1)
          xi.u <- centralvalue + fac * qnorm(p.u2)
          xi.g <- c(xi.g,xi.u)
       }
    }
 }
 xi.g <- sort(xi.g[(abs(xi.g)>cutoff.at.0)])
 l.xi.g <- length(xi.g)
 n.d <- (gridsize-l.xi.g)%/%8
 n.u <- 7*n.d + (gridsize-l.xi.g)%%8

 if(n.u>0){
    p.u <- pnorm((xi.g[l.xi.g]-centralvalue)/fac, lower.tail=FALSE)
    p.u2 <- 1-p.u + p.u*(1:n.u)/(n.u+1)
    xi.u <- centralvalue + fac * qnorm(p.u2)
    xi.g <- c(xi.g,xi.u)
 }

 if(n.d > 0 && !withPos){
    m.g.p <- min(pmax(xi.g,0))
    m.g.m <- min(pmax(-xi.g,0))
    n2.m <- n.d%/%2
    n2.p <- n2.m + n.d%%2
    xi.n.p <- cutoff.at.0+ (0:(n2.p-1))/n2.p*(m.g.p-cutoff.at.0)
    xi.g <- c(xi.g, xi.n.p)
    if(n2.m>0){
       xi.n.m <- cutoff.at.0- (1:n2.m)/(n2.m+1)*(m.g.m-cutoff.at.0)
       xi.g <- c(xi.g, xi.n.m)
    }
 }
 if(n.d > 0  && withPos){
    m.g <- min(xi.g)
    n2 <- n.d
    xi.n <- cutoff.at.0+ (0:(n2-1))/n2*(m.g-cutoff.at.0)
    xi.g <- c(xi.n, xi.g)
 }
 xi.g <- sort(unique(c(xi.g,xi.a)))
 return(xi.g)
}

getSnGrid <- function(xiGrid = getShapeGrid(), PFam=GParetoFamily(), low=0,
                      upp=1.01, accuracy = 10000, GridFileName="SnGrid.Rdata",
                      withSmooth = TRUE, withPrint = FALSE, withCall = FALSE){
   call <- match.call()
   xiGrid <- unique(sort(xiGrid))
   itSn <- 0
   getSn <- function(xi){
               itSn <<- itSn + 1
               if(withPrint) cat("Evaluation Nr.", itSn," at xi = ",xi,"\n")
               distr <- PFam@modifyParam(theta=c("scale"=1,"shape"=xi))
               return(Sn(x=as(distr,"AbscontDistribution"),
                         accuracy = accuracy, low=low, upp = upp))
               }
   SnGrid <- sapply(xiGrid,getSn)
   if(GridFileName!="") save(SnGrid, file=GridFileName)
   rm(PFam)
   iNA <- is.na(SnGrid)
   SnGrid <- SnGrid[!iNA]
   xiGrid <- xiGrid[!iNA]
   if(withSmooth)
      SnGrid <- smooth.spline(xiGrid,SnGrid)$y
   fct0 <- splinefun(x=xiGrid,y=SnGrid)
   xm <- xiGrid[1]
   ym <- SnGrid[1]
   dym <- (SnGrid[2]-SnGrid[1])/(xiGrid[2]-xiGrid[1])

   xM <- (rev(xiGrid))[1]
   yM <- (rev(SnGrid))[1]
   dyM <- ((rev(SnGrid))[2]-(rev(SnGrid))[1])/((rev(xiGrid))[2]-(rev(xiGrid))[1])
   fct <- function(x){
          y0 <- fct0(x)
          y1 <- y0
          y1[x<xm] <- ym+dym*(x[x<xm]-xm)
          y1[x>xM] <- yM+dyM*(x[x>xM]-xM)
          if(any(is.na(y0)))
             warning("There have been xi-values out of range of the interpolation grid.")
          return(y1)
   }
   environment(fct) <- nE <- new.env()
   assign("fct0",fct0, envir=nE)
       assign("yM",yM, envir=nE)
       assign("ym",ym, envir=nE)
       assign("dyM",dyM, envir=nE)
       assign("dym",dym, envir=nE)
   rm(itSn,getSn,iNA,fct0,ym,yM,dym,dyM)
   if(withCall) rm(call)
   return(list(grid = cbind(xi=xiGrid,Sn=SnGrid),
               fct = fct, call = if(withCall) call else NULL))
}

.saveInterpGrid <- function(xiGrid = getShapeGrid(), PFam = GParetoFamily(),
                        sysRdaFolder, sysdataWriteFile = "sysdata.rda",
                        getFun = getSnGrid, ..., nameInSysdata = ".SnGrids",
                        GridFileName, withSmooth = TRUE, withPrint = TRUE,
                        withCall = FALSE,
                        Y = NULL, elseFun = NULL){
   ### changed defaults and argnames (for historical reasons):
   ROptEst:::.saveInterpGrid(thGrid = xiGrid, PFam = PFam,
            sysRdaFolder = sysRdaFolder, sysdataWriteFile = sysdataWriteFile,
            getFun = getFun, ..., modifyfct = NULL,
            nameInSysdata = nameInSysdata, GridFileName = GridFileName,
            withSmooth = withSmooth, withPrint = withPrint, withCall = withCall,
            Y = Y, elseFun = elseFun)}