.RMXE.th <- function(th, PFam, modifyfct){
      PFam <- modifyfct(th,PFam)
      IC <- radiusMinimaxIC(L2Fam=PFam, neighbor= ContNeighborhood(),
                            risk = asMSE(), verbose = FALSE)
      return(c(b=clip(IC), a=cent(IC), a.w = cent(weight(IC)),
                           A=stand(IC),  A.w = stand(weight(IC))))
}

.MBRE.th <- function(th, PFam, modifyfct){
      PFam <- modifyfct(th,PFam)
      RobM <- InfRobModel(center = PFam, neighbor = ContNeighborhood(radius = 15))
      IC <- optIC(model = RobM, risk = asBias(), verbose = FALSE)
      mA <- max(stand(IC))
      mAw <- max(stand(weight(IC)))
      return(c(b=clip(IC), a=cent(IC), aw=cent(weight(IC)),
               A=stand(IC)/mA, Aw=stand(weight(IC))/mAw))
}

.OMSE.th <- function(th, PFam, modifyfct){
      PFam <- modifyfct(th,PFam)
      RobM <- InfRobModel(center = PFam, neighbor = ContNeighborhood(radius = .5))
      IC <- optIC(model = RobM, risk = asMSE(), verbose = FALSE)
      res=c(b=clip(IC), a=cent(IC), a.w = cent(weight(IC)),
                A=stand(IC), A.w = stand(weight(IC)))
      return(res)
}


.getLMGrid <- function(thGrid, PFam, optFct = .RMXE.th, modifyfct,
                       GridFileName="LMGrid.Rdata", withPrint = FALSE){
   wprint <- function(...){ if (withPrint) print(...)}
   thGrid <- unique(sort(thGrid))
   itLM <- 0
   getLM <- function(th){
               itLM <<- itLM + 1
               if(withPrint) cat("Evaluation Nr.", itLM," at th = ",th,"\n")
               a <- try(
               optFct(th=th,PFam=PFam,modifyfct=modifyfct) , silent=TRUE)
               if(is(a,"try-error")) a <- rep(NA,13)
               return(a)
               }

   distroptions.old <- distroptions()
   distrExOptions.old <- distrExOptions()
   distroptions("withgaps"=FALSE)
   distrExOptions( MCIterations=1e6,
                   GLIntegrateTruncQuantile=.Machine$double.eps,
                   GLIntegrateOrder=1000,
                   ElowerTruncQuantile=1e-7,
                   EupperTruncQuantile=1e-7,
                   ErelativeTolerance = .Machine$double.eps^0.4,
                   m1dfRelativeTolerance = .Machine$double.eps^0.4,
                   m2dfRelativeTolerance = .Machine$double.eps^0.4,
                   nDiscretize = 300, IQR.fac = 20)
   on.exit({do.call(distrExOptions,args=distrExOptions.old)
            do.call(distroptions,args=distroptions.old)
            })
   LMGrid <- t(sapply(thGrid,getLM))

   iNA <- apply(LMGrid,1, function(u) any(is.na(u)))
   LMGrid <- LMGrid[!iNA,,drop=FALSE]
   thGrid <- thGrid[!iNA]
   oG <- order(thGrid)
   thGrid <- thGrid[oG]
   LMGrid <- LMGrid[oG,,drop=FALSE]
   Grid <- cbind(xi=thGrid,LM=LMGrid)

   if(GridFileName!="") save(Grid, file=GridFileName)
   wprint(Grid)
   return(Grid)
}


.saveGridToCSV <- function(Grid, toFileCSV, namPFam, nameInSysdata){
   write.table(format(Grid,digits=21),
               file=toFileCSV, row.names=FALSE, col.names=FALSE)
   toFileTXT <- gsub("(.+\\.)csv$","\\1txt",toFileCSV)
   cat(file=toFileTXT,gsub(" ","",namPFam),"\n",nameInSysdata)
   return(invisible(NULL))
}

.readGridFromCSV <- function(fromFileCSV){
  rg <- read.table(CSVFiles[1], colClasses=rep("character",2), sep=" ", header=FALSE)
  nrg <- nrow(rg)
  Grid <- matrix(as.numeric(as.matrix(rg)),nrow=nrg)

  as.matrix(read.csv(fromFileCSV)); dimnames(Grid) <- NULL
  fromFileTXT <- gsub("(.+\\.)csv$","\\1txt",fromFileCSV)
  res2 <- scan(file=fromFileTXT, what=c("character","character"))
  return(list(Grid=Grid, namPFam=res2[1], namInSysdata=res2[2]))
}

.generateInterpGrid <- function(thGrid, PFam, toFileCSV = "temp.csv",
            getFun = .getLMGrid, ..., modifyfct, nameInSysdata,
            GridFileName, withPrint = TRUE){
  if(missing(GridFileName))
     GridFileName <- paste(sub("^\\.(.+)","\\1",nameInSysdata),".Rdata",sep="")
  Grid <- getFun(thGrid = thGrid, PFam = PFam, ..., modifyfct = modifyfct,
                 withPrint = withPrint, GridFileName = GridFileName)
  .saveGridToCSV(Grid,toFileCSV,name(PFam),nameInSysdata)
  return(invisible(NULL))
}

