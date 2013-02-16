.modify.xi.PFam.call <- function(xi, PFam){
      Param <- param(PFam)
      param <- main(Param)
      param["shape"] <- xi
      main(Param) <- param
      nModel <- modifyModel(PFam, Param)
}

.RMXE.xi <- function(xi, PFam) ROptEst:::.RMXE.th(xi, PFam, .modify.xi.PFam.call)
.MBRE.xi <- function(xi, PFam) ROptEst:::.MBRE.th(xi, PFam, .modify.xi.PFam.call)
.OMSE.xi <- function(xi, PFam) ROptEst:::.OMSE.th(xi, PFam, .modify.xi.PFam.call)

.getLMGrid <- function(xiGrid = getShapeGrid(),
                      PFam = GParetoFamily(scale=1,shape=2),
                      optFct = .RMXE.xi, GridFileName="LMGrid.Rdata",
                      withSmooth = TRUE, withPrint = FALSE, withCall = FALSE){
   ### changed defaults and argnames (for historical reasons):
   ROptEst:::.getLMGrid(thGrid = xiGrid, PFam = PFam, optFct = optFct,
           modifyfct = NULL, GridFileName = GridFileName,
           withSmooth = withSmooth, withPrint = withPrint, withCall = withCall)}

.MakeGridList <- ROptEst:::.MakeGridList

.svInt <- function(optF = .RMXE.xi, xiGrid = getShapeGrid(500, cutoff.at.0=0.005),
                   sysRdafolder, PFam = GParetoFamily(shape=1,scale=2),
                   nam = "GPD"){
             namF <- gsub("\\.xi$","",paste(deparse(substitute(optF))))
             if(missing(sysRdafolder))
                sysRdafolder <- "C:/rtest/RobASt/branches/robast-0.9/pkg/RobExtremesBuffer"
             to <- gsub("XXXX",nam,
                    gsub("YYYY", gsub("^\\.(.+)","\\1",namF), "sysdataYYYYXXXX.rda"))
             ROptEst:::.saveInterpGrid(thGrid = xiGrid,
                  PFam = PFam, sysRdaFolder=sysRdafolder, sysdataWriteFile = to,
                  getFun =  .getLMGrid, modifyfct = NULL, optFct = optF,
                  nameInSysdata = nam, withSmooth = TRUE, withPrint = TRUE)
}


