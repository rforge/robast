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
                      withPrint = FALSE){
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
   return(cbind(xi=xiGrid,Sn=SnGrid))
}

.generateInterpGridSn <- function(xiGrid = getShapeGrid(500, cutoff.at.0=0.005),
                                  PFam = GParetoFamily(), withPrint = TRUE){
   to <- gsub("XXXX",gsub(" ","",name(PFam)), "interpolSnXXXX.csv")
   Grid <- getSnGrid(xiGrid, PFam = PFam, withPrint = withPrint,
                     GridFileName = "Sn.Rdata")
  .saveGridToCSV(Grid,to,name(PFam),".Sn")
  return(invisible(NULL))
}
