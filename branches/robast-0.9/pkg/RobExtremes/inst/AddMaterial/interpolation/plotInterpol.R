plotLM <- function(Gridnam,Famnam,whichLM, baseDir="C:/rtest/robast",
               withSmooth=FALSE, gridRestriction = NULL, prehook={}, posthook={}, ...){
   ## Gridnam in (Sn,OMSE,RMXE,MBRE)
   ## Famnam in "Generalized Pareto Family",
   ##           "GEV Family",
   ##           "Gamma family",
   ##           "Weibull Family"
   ## whichLM  ignoriert für Gridnam == Sn
   #           in 1:13 (clip=b, cent.a=a1.a,a2.a, cent.i=a1.i,a2.i,
   ##                  stand.a=A.a=matrix(c(A11.a,(A12.a+A21.a)/2,
   #                                       (A12.a+A21.a)/2,A.22.a), 2, 2),
   ##                  stand.i=A.i=matrix(c(A11.i,(A12.i+A21.i)/2,
   #                                       (A12.i+A21.i)/2,A.22.i), 2, 2),
   ##                 und optIC = Y.a min(1,b/norm(Y.i)), Y.* = A.* Lambda - a.*
   ## basedir: Oberverzeichnis des r-forge svn checkouts
   file <- file.path(baseDir, "branches/robast-0.9/pkg/RobAStRDA/R/sysdata.rda")
   if(!file.exists(file)) stop("Fehler mit Checkout")
   nE <- new.env()
   load(file, envir=nE)
   Gnams <- c("Sn","OMSE","RMXE","MBRE")
   Fnams <- c("Generalized Pareto Family",
              "GEV Family",
              "Gamma family",
              "Weibull Family")
   Gridnam <- Gnams[pmatch(Gridnam, Gnams)]
   Famnam <- Fnams[pmatch(Famnam, Fnams)]
   if(! Gridnam %in% Gnams) stop("Falscher Gittername")
   if(! Famnam %in% Fnams) stop("Falscher Familienname")
   isSn <- (Gridnam == "Sn")
   Famnam0 <- gsub(" ","",Famnam)
   GN0 <- Gridnam; if(isSn) GN0 <- "SnGrids"
   GN <- paste(".",GN0,sep="")
   funN <- paste("fun",".",if(getRversion()<"2.16") "O" else "N",sep="")
   gN <- if(withSmooth) "gridS" else "grid"
   gr <- get(GN,envir=nE)[[Famnam0]][[gN]]
   if(is.null(gridRestriction)) gridRestriction <- rep(TRUE, nrow(gr))
   if(!isSn) if(whichLM!="all") if(whichLM<1 | whichLM>13) stop("Falsche Koordinate")
   if(!isSn) if(whichLM=="all"){
      eval(prehook)
      par(mfrow=c(4,4))
      for(i in 2:14)
          plot(gr[gridRestriction,1], gr[gridRestriction,i], ...)
      par(mfrow=c(1,1))
      eval(posthook)
   return(invisible(NULL))
   }
   if(isSn) whichLM <- 1
   wM <- whichLM + 1
   eval(prehook)
   plot(gr[gridRestriction,1], gr[gridRestriction,wM], ...)
   eval(posthook)
   return(invisible(NULL))
}

if(FALSE){
## Examples
plotLM("OMSE","Gamma","all", type="l", gridR=-(1:20))
plotLM("OMSE","Pareto","all", type="l", gridR=-(1:20))
plotLM("OMSE","Gener","all", type="l", gridR=-(1:20))
plotLM("OMSE","GEV","all", type="l", gridR=-(1:20))
plotLM("OMSE","Wei","all", type="l", gridR=-(1:20))
plotLM("MBRE","Wei","all", type="l", gridR=-(1:20))
plotLM("MBRE","GE","all", type="l", gridR=-(1:20))
plotLM("MBRE","Gene","all", type="l", gridR=-(1:20))
plotLM("MBRE","Gam","all", type="l", gridR=-(1:20))
plotLM("RMXE","Gam","all", type="l", gridR=-(1:20))
plotLM("RMXE","Wei","all", type="l", gridR=-(1:20))
plotLM("RMXE","Gene","all", type="l", gridR=-(1:20))
plotLM("RMXE","GE","all", type="l", gridR=-(1:20))
}