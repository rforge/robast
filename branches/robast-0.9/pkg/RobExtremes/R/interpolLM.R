.saveGridToCSV <- ROptEst:::.saveGridToCSV
.readGridFromCSV <- ROptEst:::.readGridFromCSV
.MakeSmoothGridList <- RobAStRDA:::.MakeSmoothGridList

.modify.xi.PFam.call <- function(xi, PFam){
      Param <- param(PFam)
      param <- main(Param)
      param["shape"] <- xi
      main(Param) <- param
      nModel <- modifyModel(PFam, Param)
      return(nModel)
}

.RMXE.th <- ROptEst:::.RMXE.th
.MBRE.th <- ROptEst:::.MBRE.th
.OMSE.th <- ROptEst:::.OMSE.th

.RMXE.xi <- function(xi, PFam) .RMXE.th(xi, PFam, .modify.xi.PFam.call)
.MBRE.xi <- function(xi, PFam) .MBRE.th(xi, PFam, .modify.xi.PFam.call)
.OMSE.xi <- function(xi, PFam) .OMSE.th(xi, PFam, .modify.xi.PFam.call)



.getLMGrid <- function(xiGrid = getShapeGrid(),
                      PFam = GParetoFamily(scale=1,shape=2),
                      optFct = .RMXE.xi, GridFileName="LMGrid.Rdata",
                      withPrint = FALSE){
   ### changed defaults and argnames (for historical reasons):
   ROptEst:::.getLMGrid(thGrid = xiGrid, PFam = PFam, optFct = optFct,
           modifyfct = NULL, GridFileName = GridFileName,
           withPrint = withPrint)}


.svInt <- function(optF = .RMXE.th, xiGrid = getShapeGrid(700, cutoff.at.0=0.005),
#.svInt <- function(optF = .RMXE.th, xiGrid = getShapeGrid(5, cutoff.at.0=0.005),
                   PFam = GParetoFamily(shape=1,scale=2), radius = 0.5,
                   upper = 1e4, lower = 1e-4, OptOrIter = "iterate",
                   maxiter = 150, tol = .Machine$double.eps^0.5,
                   loRad = 0, upRad = Inf, loRad0 = 1e-3,
                   withStartLM = TRUE){
             namF <- gsub("\\.th$","",paste(deparse(substitute(optF))))
             namF <- gsub("^\\.(.+)","\\1",namF)
             to <- gsub("XXXX",gsub(" ","",name(PFam)),
                    gsub("YYYY", namF, "interpolYYYYXXXX.csv"))
             print(to)
             ROptEst:::.generateInterpGrid(thGrid = xiGrid,
                  PFam = PFam, toFileCSV = to,
                  getFun =  ROptEst:::.getLMGrid,
                  modifyfct = .modify.xi.PFam.call, optFct = optF,
                  nameInSysdata = namF, withPrint = TRUE, radius = radius,
                  upper = upper, lower = lower, OptOrIter = OptOrIter,
                  maxiter = maxiter, tol = tol, loRad = loRad, upRad = upRad,
                  loRad0 = loRad0, withStartLM = withStartLM)
}


